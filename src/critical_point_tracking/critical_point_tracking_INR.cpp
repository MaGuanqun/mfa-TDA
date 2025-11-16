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
#include "connect_trajectory_degenerate_point.h"
#include "critical_point_utility.h"
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

    string input_function_name = "quartic_potential";               // diy input file


    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;

    float shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string cp_tracing_file = "cp_tracing.dat";



    std::vector<float> same_root_epsilon; // same_root_epsilon
    float time_step = 1e-3;
    
    float grad_threshold = 1e-5;

    int correction_max_itr = 30;




    float initial_point_finding_hessian_threshold = 1e-20;
    float root_finding_grad_epsilon = 1e-8;

    string input_shrink_ratio = "0-1-0-1-0-1";
    float dxy_dt_gradient_epsilon = 1e-10;

    string singular_point_file = "singular_point.dat";

    int max_itr=50;

    float point_itr_threshold = 0.5;

    float  spatial_step_size = 1.0;
        string input_model="";
    string edge_type_file ="";

    ops >> opts::Option('f', "input_function_name",  input_function_name,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "cp_tracing_file", cp_tracing_file, " file name of cp_tracing");
    ops >> opts::Option('z', "time_step",    time_step,       " time step size");
    ops >> opts::Option('g', "spatial_step_size",    spatial_step_size,       " spatial step size");

    ops >> opts::Option('x', "root_finding_grad_epsilon",    root_finding_grad_epsilon,       "first root finding epsilon");

    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
    ops >> opts::Option('m', "max_itr", max_itr, " max iteration");
    ops >> opts::Option('s', "singular_point_file", singular_point_file, " singular point file name");

    ops >> opts::Option('p', "point_itr_threshold", point_itr_threshold, " stop iteration when point update is less than point_itr_threshold * step size");

    ops >> opts::Option('i', "input_model", input_model, " input INR model");

    ops >> opts::Option('e', "edge_type_file", edge_type_file, " edge type file name");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }


    std::istringstream iss(input_shrink_ratio);
    std::vector<float> shrink_ratio;
    float number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    }  


    INRModel<float> inr_model(input_function_name,input_model);
    int function_type=-1; 



    Eigen::VectorXf local_domain_range=inr_model.domain_max-inr_model.domain_min;
    VectorXf core_maxs = inr_model.domain_max;
    VectorXf core_mins = inr_model.domain_min;

    VectorXi span_num = inr_model.block_num;

    VectorXf Span_size = local_domain_range.cwiseQuotient(span_num.cast<float>());


    
    std::vector<float> step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);
    step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time

    float d_max_square_= spatial_step_size*spatial_step_size/16* step_size[0]* step_size[0]; 

    std::vector<VectorX<float>> degenerate_points;
    Degenerate_case_tracing<float>::read_degenerate_point(singular_point_file,degenerate_points);


    std::vector<CP_Trace<float>> traces;
   
    same_root_epsilon = step_size; // same_root_epsilon

    std::cout<<Span_size.transpose()<<std::endl;
    std::cout<<"spatial step size "<<step_size[0]<<" "<<"time step " <<step_size.back()<<std::endl;

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);


        std::vector<VectorX<float>> root; //the inner vector store the root in a span

        std::vector<VectorXi> selected_span;
        span_filter::compute_boundary_span(span_num,selected_span,true);
        std::cout<<"valid span num "<<selected_span.size()<<std::endl;
        // std::vector<VectorXi> selected_span;

        auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();


        VectorXi point_num_in_block = inr_model.point_num_in_block; //number of initial points in a block

        Find_boundary_roots find_boundary_roots(root_finding_grad_epsilon,core_mins,core_maxs,point_num_in_block,span_num,same_root_epsilon,function_type,max_itr,point_itr_threshold,static_cast<Block<float>*>(nullptr),&inr_model);

        find_boundary_roots.root_finding(selected_span, root);


        std::cout<<"find root num before deduplicate between spans "<<root.size()<<std::endl;

        std::vector<VectorX<float>> root_unique;
        spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,same_root_epsilon[0],same_root_epsilon.back());


        std::cout<<"finish finding root before deduplicate between spans "<<root.size()<<" after "<<root_unique.size()<<std::endl;

        root.clear();
        root.shrink_to_fit();

        auto finding_end_time = std::chrono::high_resolution_clock::now();

        string test_file=cp_tracing_file+"_test.obj";

        tracking_utility::convert_to_obj(test_file,root_unique);

        auto tracking_start_time = std::chrono::high_resolution_clock::now();

        traces.resize(root_unique.size());

        Boundary_critical_point_tracking boundary_critical_point_tracking(core_mins, core_maxs,root_finding_grad_epsilon, step_size.back(), step_size[0], d_max_square_,function_type, correction_max_itr, static_cast<Block<float>*>(nullptr), &inr_model);

        boundary_critical_point_tracking.find_trace(root_unique, traces);

        std::cout<<"finish boundary critical point tracing "<<std::endl;


        Degenerate_case_tracing degenerate_case_tracing(core_mins, core_maxs, point_num_in_block, &find_boundary_roots, step_size.back(), step_size[0], root_finding_grad_epsilon,correction_max_itr, function_type, static_cast<Block<float>*>(nullptr), &inr_model);
        degenerate_case_tracing.tracing_from_all_degenerate_points(degenerate_points, traces, 0.1, d_max_square_);




    int trace_size=0;
    for(auto& trace:traces)
    {
        if((!trace.duplicated) && trace.traces.size()>=1)
        {            
            trace_size++;
        }
    }
    std::cout<<"trace before splitting "<<traces.size()<<" real traces "<<trace_size<<std::endl;

        //the result will sort degenerate_points by time
    deduplication::deduplicate_traces(traces, degenerate_points, step_size[0], step_size.back(), core_mins);

    std::cout<<"trace_after splitting "<<traces.size()<<std::endl;


    connect_trajectory_degenerate_point::connect_trajectory(traces, degenerate_points, step_size[0], step_size.back(), core_mins);

    trace_size=0;
    for(auto& trace:traces)
    {
        if((!trace.duplicated) && trace.traces.size()>=1)
        {            
            trace_size++;
        }
    }
    std::cout<<"traces after deduplication "<<trace_size<<std::endl;


  auto tracking_end_time = std::chrono::high_resolution_clock::now();

    std::cout<<"overall extraction time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(tracking_end_time - cpt_extract_start_time).count()/1000<<std::endl;

    std::cout<<"overall finding time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(finding_end_time - cpt_extract_start_time).count()/1000<<std::endl;


    std::cout<<"overall tracking time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(tracking_end_time - tracking_start_time).count()/1000<<std::endl;


    std::vector<int> critical_point_types;
    if(edge_type_file!="")
        critical_point_utility::compute_critical_point_type(traces, degenerate_points, critical_point_types,function_type, static_cast<Block<float>*>(nullptr), &inr_model);

    CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces, degenerate_points,&critical_point_types, edge_type_file);

    critical_point_utility::accuracy(traces, degenerate_points, function_type,static_cast<Block<float>*>(nullptr), &inr_model);


    // CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces,degenerate_points);


  
}