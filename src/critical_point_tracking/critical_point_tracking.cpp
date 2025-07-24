#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include <chrono>


#include <tbb/tbb.h>
#include <atomic>

#include "opts.h"

#include "block.hpp"

#include "critical_point/span_filter.hpp"
#include "save_control_data.hpp"

// #include "find_isocontour.h"
// #include"span_filter.h"
// #include "find_root.h"
// #include "../critical_point/find_all_root.h"
// #include "ridge_valley_graph.h"
// #include "find_root_h.h"

// #include "connect_rv_graph.h"

// #include "transfer_data.h"
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

#include "xy_critical_point_finding.h"
#include "xy_critical_point_tracking.h"

#include "tracking_utility.h"
// #include "trace.h"

// #include "../morse_smale/find_isocontour.h"
// #include"../morse_smale/span_filter.h"
// #include "../morse_smale/find_root.h"
// #include "../critical_point/find_all_root.h"
// #include "../morse_smale/ridge_valley_graph.h"
// #include "../morse_smale/find_root_h.h"

// #include "../morse_smale/connect_rv_graph.h"


using namespace std;


namespace {
    tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
}

int main(int argc, char** argv)
{

    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    // string inControlPoint = "derivative_control_point.dat";

    string inControlPoint = "derivative_control_point.dat";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;
    int max_step = 5000; //max point number in one isocontour

    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string cp_tracing_file = "cp_tracing.dat";

    string edge_type = "";


    std::vector<double> same_root_epsilon; // same_root_epsilon
    double time_step = 1e-3;
    
    double grad_threshold = 1e-5;

    int correction_max_itr = 30;




    real_t initial_point_finding_hessian_threshold = 1e-20;
    real_t root_finding_grad_epsilon = 1e-8;
    real_t hessian_threshold_for_cpt_tracking = 1e-12;

    string input_shrink_ratio = "0-1-0-1-0-1";
    real_t dxy_dt_gradient_epsilon = 1e-10;

    string singular_point_file = "singular_point.dat";

    int max_itr=40;

    real_t point_itr_threshold = 0.5;

    double  spatial_step_size = 1.0;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "cp_tracing_file", cp_tracing_file, " file name of cp_tracing");
    ops >> opts::Option('z', "time_step",    time_step,       " time step size");
    ops >> opts::Option('g', "spatial_step_size",    spatial_step_size,       " spatial step size");
    // ops >> opts::Option('g', "grad_threshold",    grad_threshold,       " gradient is smaller enough to change to permutate the point a little bit");
    ops >> opts::Option('q', "edge_type",    edge_type,       "edge type file, (pseudo) ridge/valley");

    ops >> opts::Option('a', "inControlPoint",  inControlPoint,  " diy input derivative control point file name");

    ops >> opts::Option('x', "root_finding_grad_epsilon",    root_finding_grad_epsilon,       "first root finding epsilon");

    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
    ops >> opts::Option('m', "max_itr", max_itr, " max_itr");
    ops >> opts::Option('s', "singular_point_file", singular_point_file, " singular point file name");

    ops >> opts::Option('p', "point_itr_threshold", point_itr_threshold, " stop iteration when point update is less than point_itr_threshold * step size");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }


    std::istringstream iss(input_shrink_ratio);
    std::vector<double> shrink_ratio;
    double number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    }  




    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks to be moved out of core
    diy::Master      master(world,
            -1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy,
            &storage,
            &Block<real_t>::save,
            &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()

     // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";


    std::vector<std::vector<std::vector<std::vector<double>>>> geo_control_point; //[deriv][vars][dom][...]
    std::vector<std::vector<MatrixX<double>>> sci_deriv_control_points;//[vars][partial_deriv][...]
    save_control_points::load_control_points(inControlPoint.c_str(),geo_control_point,sci_deriv_control_points);

    std::vector<std::vector<VectorX<double>>> domain_root(master.size()); //[blocks,]

    std::vector<VectorX<double>> degenerate_points;
    degenerate_case_tracing::read_degenerate_point(singular_point_file,degenerate_points);

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;

        std::vector<VectorX<double>> root_ori;

        double min_ = local_domain_range.minCoeff();
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        std::vector<double> step_size(span_num.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);

       
        step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time
   

        same_root_epsilon = step_size; // same_root_epsilon

        std::cout<<Span_size.transpose()<<std::endl;
        std::cout<<"spatial step size "<<step_size[0]<<" "<<"time step " <<step_size.back()<<std::endl;

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);


        // std::vector<double> shrink_ratio_in_unit;
        // for(int i=0;i<shrink_ratio.size()/2;i++)
        // {
        //     shrink_ratio_in_unit.emplace_back((shrink_ratio[2*i]-b->core_mins[i])/ local_domain_range[i]);
        //     shrink_ratio_in_unit.emplace_back((shrink_ratio[2*i+1]-b->core_mins[i])/ local_domain_range[i]);
        // }

        std::vector<VectorX<real_t>> root; //the inner vector store the root in a span

        std::vector<VectorX<real_t>> root_record;

        std::vector<std::vector<VectorXi>> selected_span;//[vars,span index]



        span_filter::compute_valid_span(sci_deriv_control_points,b,selected_span,shrink_ratio,2);





        span_filter::compute_boundary_span(b,selected_span,true);

        // for(int i=0;i<selected_span[0].size();i++)
        // {
        //     std::cout<<i<<" "<<selected_span[0][i].transpose()<<std::endl;
        // }

        tbb::enumerable_thread_specific<std::vector<VectorX<double>>> local_root;

            
        auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();

        std::vector<int>multi_root_span; 
        Eigen::VectorXd weights=Eigen::VectorXd::Ones(b->mfa->var(0).tmesh.tensor_prods[0].ctrl_pts.rows());



        tbb::affinity_partitioner ap;

        
        tbb::parallel_for(tbb::blocked_range<size_t>(0,selected_span[0].size()), //
        [&](const tbb::blocked_range<size_t>& range)
        {
            auto& root_thread = local_root.local();

            for(auto i=range.begin();i!=range.end();++i)
            {
                // std::cout<<selected_span[0][i][2]<<std::endl;
                std::vector<VectorX<real_t>> root_block;
                std::vector<real_t> function_value_block;
                root_block.reserve(16);
                function_value_block.reserve(16);
                // std::vector<VectorXd> root_span;
                if(find_boundary_roots::root_finding(b,selected_span,root_block,i,root_finding_grad_epsilon, same_root_epsilon,initial_point_finding_hessian_threshold,max_itr, point_itr_threshold))
                {
                    VectorXi selected_span_index=selected_span[0][i]- b->mfa->var(0).p;
                    size_t index=utility::obtain_index_from_domain_index(selected_span_index,number_in_every_domain);
                    // root_info span_root;
                    // span_root.roots = std::move(root_block);
                    // span_root.function_value = std::move(function_value_block);
                    // span_root.valid_span_index = index;

                    root_thread.insert(root_thread.end(),root_block.begin(),root_block.end());

                }

                if(i%100000==0)
                {
                    std::cout<<"finish find span "<<i<<std::endl;
                }
                // break;
                //multiplicity_root[index].insert(multiplicity_root[index].end(),multi_root_span.begin(),multi_root_span.end());
            }
        },ap               
        );




        size_t total_size = 0;
        for (const auto& thread_indices : local_root)
            total_size += thread_indices.size();

        root.reserve(total_size);

        for (auto& thread_root_block : local_root)
        {
            root.insert(root.end(), thread_root_block.begin(), thread_root_block.end());
        }   


        string test_file=cp_tracing_file+"_test.obj";


        tracking_utility::convert_to_obj(test_file,root);
        // find_boundary_roots::test_root_finding(b,root,root_finding_grad_epsilon);

        std::cout<<"finish finding root "<<root.size()<<std::endl;

        std::vector<CP_Trace<double>> traces;

        traces.resize(root.size());

        xy_cp_tracking::find_trace(step_size.back(),step_size[0],max_step,b,root, traces,hessian_threshold_for_cpt_tracking,root_finding_grad_epsilon,correction_max_itr);

        double max_dis_stop_square = 0.0;
        for(int  i=0;i<step_size.size()-1;++i)
        {
            max_dis_stop_square+=step_size[i]*step_size[i];
        }
        max_dis_stop_square*=25.0;


        degenerate_case_tracing::tracing_from_all_degenerate_points(b, degenerate_points, traces, step_size.back(), step_size[0], 0.1, initial_point_finding_hessian_threshold, root_finding_grad_epsilon, max_itr, correction_max_itr, max_dis_stop_square,point_itr_threshold);

        CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces);



    });

}