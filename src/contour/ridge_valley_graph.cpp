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

#include "find_isocontour.h"
#include"span_filter.h"
#include "find_root.h"
#include "../critical_point/find_all_root.h"
#include "ridge_valley_graph.h"
#include "find_root_rv_graph.h"

#include "connect_rv_graph.h"

#include <tbb/global_control.h>

// #include "transfer_data.h"
// #include <fstream>
// #include <iostream>
// #include <Eigen/Dense>


using namespace std;

struct root_info
{
    std::vector<VectorX<real_t>> roots;
    std::vector<real_t> function_value;
    size_t valid_span_index;
};


// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }


int main(int argc, char** argv)
{

    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    // string inControlPoint = "derivative_control_point.dat";
    string inCriticalPoint = "default_critical_point.dat";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;

    int max_step = 5000; //max point number in one isocontour

    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";


    string ridge_valley_graph_file = "ridge_valley_graph.dat";
    double connect_traces_threshold = 2.0;


    string edge_type = "";



    double tolerance = 1e-4; // same_root_epsilon
    double step_size = 5.0;
    double epsilon = 1e-3; // decide if the result is close enough to the initail point and stop iteration
    double threshold_correction = 1e-4; //activate the correction of point tracing
    double grad_square_threshold = 1e-10;

    double trace_split_grad_square_threshold = 1e-4; //check if need to split the trace

    int correction_max_itr = 30;

    // double center_threshold_square = epsilon*epsilon;
    double trace_threshold_square = step_size*step_size; //check duplication

    double threshold_connect_traces_square = 1e-4;

    // double same_root_epsilon = SAME_ROOT_EPSILON;


    real_t root_finding_epsilon = 1e-8;

    real_t threshold_distance_stop_tracing = 0.3;
    string input_shrink_ratio = "0-1-0-1";

    string func_error_file = "func_error.dat";

    int z_zero = 1;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('j', "inCriticalPoint",  inCriticalPoint,  " diy input critical point file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "ridge_valley_graph_file", ridge_valley_graph_file, " file name of ridge valley graph");
    ops >> opts::Option('e', "same_root_epsilon", epsilon, " epsilon for reaching start point");
    ops >> opts::Option('z', "step_size",    step_size,       " step size");
    ops >> opts::Option('a', "connect_traces_threshold",    connect_traces_threshold,       " threshold to connect traces");
    ops >> opts::Option('c', "threshold_correction",    threshold_correction,       " threshold activate the correction of point tracing");
    ops >> opts::Option('g', "grad_square_threshold",    grad_square_threshold,       " gradient is smaller enough to change to permutate the point a little bit");
    ops >> opts::Option('p', "trace_split_grad_square_threshold",    trace_split_grad_square_threshold,       "determine split a trace");
    ops >> opts::Option('q', "edge_type",    edge_type,       "edge type file, (pseudo) ridge/valley");
    ops >> opts::Option('x', "root_finding_epsilon",    root_finding_epsilon,       "first root finding epsilon");
    ops >> opts::Option('y', "threshold_distance_stop_tracing",    threshold_distance_stop_tracing,       "when the distance from previous point is smaller than this*step_size, stop tracing");
    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
    ops >> opts::Option('w', "z_zero",    z_zero,       " output z as zero");
    ops >> opts::Option('s', "func_error_file",    func_error_file,       " file name of function error");


    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }






    std::cout << "Max concurrency0: " 
              << tbb::this_task_arena::max_concurrency() 
              << std::endl;



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



    threshold_correction = root_finding_epsilon;



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


 


   




    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        // Eigen::VectorXd weights=Eigen::VectorXd::Ones(b->vars[0].mfa_data->tmesh.tensor_prods[0].ctrl_pts.rows());
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;


        double min_ = local_domain_range.minCoeff();
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        // std::cout<<"span_num "<<span_num.transpose()<<std::endl;


        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());
        double min_span_size = Span_size.minCoeff()/step_size;
        step_size = min_span_size;


        tolerance = 0.5*step_size; // same_root_epsilon
        epsilon = step_size; // decide if the result is close enough to the initail point and stop iteration
        trace_threshold_square = step_size*step_size; //check duplication

        // std::cout<<"step "<<step_size<<std::endl;
        // std::cout<<"trace_threshold_square__"<<trace_threshold_square<<std::endl;


        threshold_connect_traces_square = connect_traces_threshold * connect_traces_threshold* step_size*step_size;

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);
            std::vector<std::vector<VectorXi>> valid_span; //[nvars][]


            std::vector<std::vector<VectorX<real_t>>> root; //the inner vector store the root in a span
            std::vector<std::vector<real_t>> function_value;

            std::vector<VectorX<real_t>> root_record;
            std::vector<real_t> record_function_value;

            std::vector<size_t> valid_span_index;
           
          


                root.clear(); //the inner vector store the root in a span
                function_value.clear();
                root_record.clear();
                record_function_value.clear();
                valid_span_index.clear();
                std::vector<TraceInSpan<double>> ridge_valley_in_span;

                // std::cout<<"start finding root"<<std::endl;
                


                tbb::enumerable_thread_specific<std::vector<root_info>> local_root;

                // std::vector<std::vector<VectorX<real_t>>> root_block; //the inner vector store the root in a span
                // std::vector<std::vector<real_t>> function_value_block;
                // std::vector<size_t> valid_span_index_block;

                auto start_time = std::chrono::high_resolution_clock::now();

                tbb::affinity_partitioner ap;


                // std::vector<VectorX<real_t>> root_block;
                // std::vector<real_t> function_value_block;

                // find_root_rv_graph::span_range(b);



                std::vector<double> shrink_ratio_normalized(shrink_ratio.size());
                if(shrink_ratio[0]==0 && shrink_ratio[2]==0 && shrink_ratio[1]==1 && shrink_ratio[3]==1)
                {
                    shrink_ratio_normalized = shrink_ratio;
                }
                else
                {
                    for(int i=0;i<2;i++)
                    {
                        shrink_ratio_normalized[i]=(shrink_ratio[i]-b->core_mins[0])/local_domain_range[0];
                        shrink_ratio_normalized[i+2]=(shrink_ratio[i+2]-b->core_mins[1])/local_domain_range[1];
                    }
                }               

                // std::cout<<"shrink_ratio_normalized "<<shrink_ratio_normalized[0]<<" "<<shrink_ratio_normalized[1]<<" "<<shrink_ratio_normalized[2]<<" "<<shrink_ratio_normalized[3]<<std::endl;
                

                tbb::parallel_for(tbb::blocked_range<size_t>(0,spanned_block_num),
                [&](const tbb::blocked_range<size_t>& range)
                {

                    auto& root_thread = local_root.local();

                    for(auto i=range.begin();i!=range.end();++i)
                    {
                        std::vector<VectorX<real_t>> root_block;
                        std::vector<real_t> function_value_block;
                        root_block.reserve(16);
                        function_value_block.reserve(16);
                        if(find_root_rv_graph::root_finding(b,root_block,number_in_every_domain,i, root_finding_epsilon,function_value_block,tolerance,shrink_ratio_normalized))
                        {
                            root_info span_root;
                            span_root.roots = std::move(root_block);
                            span_root.function_value = std::move(function_value_block);
                            span_root.valid_span_index = i;

                            root_thread.emplace_back(std::move(span_root));

                        }
                    }
                },ap

                );

                size_t total_size = 0;
                for (const auto& thread_indices : local_root)
                    total_size += thread_indices.size();

                root.reserve(total_size);
                function_value.reserve(total_size);
                valid_span_index.reserve(total_size);

                for (auto& thread_root_block : local_root)
                {
                    for(auto& root_info : thread_root_block)
                    {
                        root.emplace_back(std::move(root_info.roots));
                        valid_span_index.emplace_back(root_info.valid_span_index);
                        function_value.emplace_back(std::move(root_info.function_value));
                    }
                }
                               


                auto root_end_time = std::chrono::high_resolution_clock::now();

                // std::cout<<"finish finding root"<<root.size()<<std::endl;

                int root_num=0;
                for(auto i=0;i<root.size();++i)
                {
                    root_num+=root[i].size();
                }

                std::cout<<"total root number "<< root_num<<std::endl;

                ridge_valley_in_span.resize(valid_span_index.size());

                ridge_valley_graph::find_ridge_valley_graph(step_size,max_step,b,root, valid_span_index,trace_threshold_square,ridge_valley_in_span,threshold_correction,trace_split_grad_square_threshold,correction_max_itr,threshold_distance_stop_tracing*threshold_distance_stop_tracing);

                std::cout<<"finish tracing"<<std::endl;

                // for(int i=0;i<ridge_valley_in_span.size();i++)
                // {
                //    if(ridge_valley_graph::test_correct(ridge_valley_in_span[i]))
                //    {
                //        std::cout<<"occur error +1"<<i<<std::endl;
                //    }
                // }


                // std::cout<<"total points "<<num<<std::endl;

                ridge_valley_graph::getFunctionValue(b,ridge_valley_in_span);
                // std::cout<<"getting function "<<std::endl;

                // for(int i=0;i<ridge_valley_in_span.size();i++)
                // {
                //    if(ridge_valley_graph::test_correct(ridge_valley_in_span[i]))
                //    {
                //        std::cout<<"occur error +2"<<i<<std::endl;
                //    }
                // }


                ridge_valley_graph::determine_ridge_valley(ridge_valley_in_span,b);


                // for(int i=0;i<ridge_valley_in_span.size();i++)
                // {
                //    if(ridge_valley_graph::test_correct(ridge_valley_in_span[i]))
                //    {
                //        std::cout<<"occur error +3"<<i<<std::endl;
                //    }
                // }


                std::vector<VectorX<double>> critical_points;
                connect_rv_graph::read_critical_point(inCriticalPoint,critical_points,false,false,0.0,0.01);


                // for(int i=0;i<ridge_valley_in_span.size();i++)
                // {
                //    if(ridge_valley_graph::test_correct(ridge_valley_in_span[i]))
                //    {
                //        std::cout<<"occur error +4"<<i<<std::endl;
                //    }
                // }

                
                std::cout<<"print z zero "<<z_zero<<std::endl;

                connect_rv_graph::connect_ridge_valley_graph(ridge_valley_in_span,critical_points,b,spanned_block_num, span_num,
                valid_span_index,
                threshold_connect_traces_square, ridge_valley_graph_file,step_size*step_size,z_zero,edge_type);
                

                auto rv_graph_end_time = std::chrono::high_resolution_clock::now();

                ridge_valley_graph::get_rv_error(b,ridge_valley_in_span,func_error_file);



                std::cout<<"root finding running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(root_end_time - start_time).count()/1000<<std::endl;
                std::cout<<"ridge valley graph running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(rv_graph_end_time - root_end_time).count()/1000<<std::endl;
                std::cout<<"total running time, millisecond : "<< std::chrono::duration_cast<std::chrono::microseconds>(rv_graph_end_time - start_time).count()/1000<<std::endl;

            

            

        // }

    });


}