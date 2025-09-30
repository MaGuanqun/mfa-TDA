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

#include "find_boundary_roots.h"
#include "boundary_critical_point_tracking.h"

#include "tracking_utility.h"
#include "trace_deduplication.h"

// #include "trace.h"

// #include "../morse_smale/find_isocontour.h"
// #include"../morse_smale/span_filter.h"
// #include "../morse_smale/find_root.h"
// #include "../critical_point/find_all_root.h"
// #include "../morse_smale/ridge_valley_graph.h"
// #include "../morse_smale/find_root_h.h"

// #include "../morse_smale/connect_rv_graph.h"


using namespace std;


// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }

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
    Degenerate_case_tracing<double>::read_degenerate_point(singular_point_file,degenerate_points);


    std::vector<CP_Trace<double>> traces;
    std::vector<double> step_size;

    VectorXd core_mins;

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;
        core_mins=b->core_mins;

        double min_ = local_domain_range.minCoeff();
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        step_size.resize(span_num.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);

       
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

        std::vector<std::vector<VectorXi>> selected_span;//[vars,span index]




        span_filter::compute_valid_span(sci_deriv_control_points,b,selected_span,shrink_ratio,2,true);


        span_filter::compute_boundary_span(b,selected_span,true);

        std::cout<<"valid span num "<<selected_span[0].size()<<std::endl;
            
        auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();


        VectorXi point_num_in_block = b->mfa->var(0).p +2 * VectorXi::Ones(b->mfa->var(0).p.size()); //n
 

        Find_boundary_roots find_boundary_roots(root_finding_grad_epsilon,b->core_mins,b->core_maxs,point_num_in_block,span_num,same_root_epsilon,0,max_itr,point_itr_threshold,b);

        find_boundary_roots.root_finding(selected_span[0], root);


        std::cout<<"find root num before deduplicate between spans "<<root.size()<<std::endl;

        std::vector<VectorX<double>> root_unique;
        spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,step_size[0],step_size.back());


        std::cout<<"finish finding root before deduplicate between spans "<<root.size()<<" after "<<root_unique.size()<<std::endl;

        root.clear();
        root.shrink_to_fit();

        string test_file=cp_tracing_file+"_test.obj";

        tracking_utility::convert_to_obj(test_file,root_unique);





        traces.resize(root_unique.size());
        int function_type=0;


        Boundary_critical_point_tracking boundary_critical_point_tracking(b->core_mins, b->core_maxs,root_finding_grad_epsilon, step_size.back(), step_size[0], function_type, correction_max_itr,b);

        boundary_critical_point_tracking.find_trace(root_unique, traces);

        double max_dis_stop_square = 0.0;
        for(int  i=0;i<step_size.size()-1;++i)
        {
            max_dis_stop_square+=step_size[i]*step_size[i];
        }
        max_dis_stop_square*=25.0;


        Degenerate_case_tracing degenerate_case_tracing(b->core_mins, b->core_maxs, point_num_in_block, &find_boundary_roots, step_size.back(), step_size[0], root_finding_grad_epsilon,correction_max_itr, 0, b);

        degenerate_case_tracing.tracing_from_all_degenerate_points(degenerate_points, traces, 0.1, max_dis_stop_square);


    });

    int trace_size=0;
    for(auto& trace:traces)
    {
        if((!trace.duplicated) && trace.traces.size()>1)
        {            
            trace_size++;
        }
    }
    std::cout<<"trace before splitting "<<traces.size()<<" real traces "<<trace_size<<std::endl;

        //the result will sort degenerate_points by time
    deduplication::deduplicate_traces(traces, degenerate_points, step_size[0], step_size.back(), core_mins);

    std::cout<<"trace_after splitting "<<traces.size()<<std::endl;

    trace_size=0;
    for(auto& trace:traces)
    {
        if((!trace.duplicated) && trace.traces.size()>1)
        {            
            trace_size++;
        }
    }
    std::cout<<"traces after deduplication "<<trace_size<<std::endl;

    CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces);

}