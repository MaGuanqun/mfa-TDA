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
    tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 8);
}


void save_root(std::vector<VectorX<double>>& root_unique, string degenerate_point_file, int spatial_step_size)
{
    std::vector<MatrixXd> root_matrix(1);

    root_matrix[0].resize(root_unique.size(),root_unique[0].size());
    for(int j=0;j<root_unique.size();j++)
    {
        root_matrix[0].row(j) = root_unique[j].transpose();
    }

    string degenerate_file_name = degenerate_point_file + std::to_string(spatial_step_size) + ".dat";

    utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);

    std::cout<<"save root with step size "<<degenerate_file_name<<std::endl;

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

    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string cp_tracing_file = "cp_tracing.dat";



    std::vector<double> same_root_epsilon; // same_root_epsilon
    double time_step = 1e-3;
    
    double grad_threshold = 1e-5;

    int correction_max_itr = 30;




    double initial_point_finding_hessian_threshold = 1e-20;
    double root_finding_grad_epsilon = 1e-8;

    string input_shrink_ratio = "0-1-0-1-0-1";
    double dxy_dt_gradient_epsilon = 1e-10;

    string singular_point_file = "singular_point.dat";

    int max_itr=50;

    double point_itr_threshold = 0.5;

    double  spatial_step_size = 1.0;
        string input_model="";
    string edge_type_file ="";

    string boundary_start="";
    int compute_boundary_start=1;

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

    ops >> opts::Option('j', "boundary_start", boundary_start, " boundary_start_file");
    ops >> opts::Option('c', "compute_boundary_start", compute_boundary_start, " compute_boundary_start");

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


    INRModel<double> inr_model(input_function_name,input_model);
    int function_type=-1;

    // Vector of modules, one per TBB worker; each thread uses thread_modules_[its index] with no per-call clone/lock.
    {
        int nw = std::max(1, tbb::this_task_arena::max_concurrency());
        std::cout<<"number of TBB workers "<<nw<<std::endl;
        inr_model.prepare_thread_modules(static_cast<size_t>(nw));
    }

    Eigen::VectorXd local_domain_range=inr_model.domain_max-inr_model.domain_min;
    VectorXd core_maxs = inr_model.domain_max;
    VectorXd core_mins = inr_model.domain_min;

    VectorXi span_num = inr_model.block_num;

    VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());


    
    std::vector<double> step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);
    step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time

    double d_max_square_= spatial_step_size*spatial_step_size/16* step_size[0]* step_size[0]; 

    std::vector<VectorX<double>> degenerate_points;
    Degenerate_case_tracing<double>::read_degenerate_point(singular_point_file,degenerate_points);


    std::vector<CP_Trace<double>> traces;
   
    same_root_epsilon = step_size; // same_root_epsilon

    std::cout<<Span_size.transpose()<<std::endl;
    std::cout<<"spatial step size "<<step_size[0]<<" "<<"time step " <<step_size.back()<<std::endl;

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);


        std::vector<VectorX<double>> root; //the inner vector store the root in a span

        std::vector<VectorXi> selected_span;
        span_filter::compute_boundary_span(span_num,selected_span,true);
        std::cout<<"valid span num "<<selected_span.size()<<std::endl;
        // std::vector<VectorXi> selected_span;

        auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();


        VectorXi point_num_in_block = inr_model.point_num_in_block; //number of initial points in a block
        std::vector<VectorX<double>> root_unique;
        Find_boundary_roots find_boundary_roots(root_finding_grad_epsilon,core_mins,core_maxs,point_num_in_block,span_num,same_root_epsilon,function_type,max_itr,point_itr_threshold,static_cast<Block<double>*>(nullptr),&inr_model);
        if(compute_boundary_start==1){
            find_boundary_roots.root_finding(selected_span, root);

            spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,same_root_epsilon[0],same_root_epsilon.back());

            save_root(root_unique, boundary_start, spatial_step_size);

            auto temp_step_size= step_size;
            auto temp_spatial_ratio =spatial_step_size;
            std::cout<<"find root num before deduplicate between spans "<<root.size()<<std::endl;
            std::cout<<"finish finding root before deduplicate between spans "<<root.size()<<" after "<<root_unique.size()<<std::endl;
            for (int i = spatial_step_size/2; i > 1; i /= 2)
            {
                root_unique.clear();
                temp_step_size[0] *=2;
                temp_step_size.back() *=2;
                spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,temp_step_size[0],temp_step_size.back());
                save_root(root_unique, boundary_start, i);
                std::cout<<"finish finding root before deduplicate between spans "<<root.size()<<" after "<<root_unique.size()<<std::endl;
                
            }
            root.clear();
            root.shrink_to_fit();
        }
        else
        {
            string name  =  boundary_start + std::to_string(int(spatial_step_size)) + ".dat";
            Degenerate_case_tracing<double>::read_degenerate_point(name,root_unique);

            std::cout<<"read root num from file "<<root_unique.size()<<std::endl;
        }

        auto finding_end_time = std::chrono::high_resolution_clock::now();

        // string test_file=cp_tracing_file+"_test.obj";

        // tracking_utility::convert_to_obj(test_file,root_unique);

        auto tracking_start_time = std::chrono::high_resolution_clock::now();

        traces.resize(root_unique.size());

        Boundary_critical_point_tracking boundary_critical_point_tracking(core_mins, core_maxs,root_finding_grad_epsilon, step_size.back(), step_size[0], d_max_square_,function_type, correction_max_itr, static_cast<Block<double>*>(nullptr), &inr_model);

        boundary_critical_point_tracking.find_trace(root_unique, traces);

        std::cout<<"finish boundary critical point tracing "<<std::endl;


        Degenerate_case_tracing degenerate_case_tracing(core_mins, core_maxs, point_num_in_block, &find_boundary_roots, step_size.back(), step_size[0], root_finding_grad_epsilon,correction_max_itr, function_type, static_cast<Block<double>*>(nullptr), &inr_model);
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
        critical_point_utility::compute_critical_point_type(traces, degenerate_points, critical_point_types,function_type, static_cast<Block<double>*>(nullptr), &inr_model);

    CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces, degenerate_points,&critical_point_types, edge_type_file);

    critical_point_utility::accuracy(traces, degenerate_points, function_type,static_cast<Block<double>*>(nullptr), &inr_model);


    // CP_Trace_fuc::convert_to_obj(cp_tracing_file,traces,degenerate_points);


  
}