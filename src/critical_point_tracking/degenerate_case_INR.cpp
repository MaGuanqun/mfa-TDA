// #include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

// #include <diy/master.hpp>
// #include <diy/reduce-operations.hpp>
// #include <diy/decomposition.hpp>
// #include <diy/assigner.hpp>
// #include <diy/io/block.hpp>

#include <chrono>


#include <tbb/tbb.h>
#include <atomic>

// #include "opts.h"

// #include "block.hpp"


#include <fstream>
#include <iostream>
#include <Eigen/Dense>

#include "degenerate_case.h"

// #include "tracking_utility.h"

#include "spatial_hashing_spatial_temporal.h"
#include "closed_form_function.h"
#include "../INRModel.h"

using namespace std;


namespace {
    tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
}

void choose_span(std::vector<VectorXi>& record_span)
{
    std::vector<std::vector<int>> span_index(3);
    span_index[0]={4,5};
    span_index[1]={4,5};
    span_index[2]={2,3};

    std::cout<<"total number of spans to process: "<<span_index[0].size()<<" "<<span_index[1].size()<<" "<<span_index[2].size()<<std::endl;

    for(int i=0;i<span_index[0].size();++i)
    {
        for(int j=0;j<span_index[1].size();++j)
        {
            for(int k=0;k<span_index[2].size();++k)
            {
                VectorXi temp(3);
                temp<<span_index[0][i],span_index[1][j],span_index[2][k];
                record_span.emplace_back(temp);
            }
        }
    }
}

int main(int argc, char** argv)
{

    // setenv("CUDA_VISIBLE_DEVICES", "", /*overwrite=*/1);


    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;


    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string degenerate_point_file = "degenerate_point.dat";
    
    double J_threshold = std::numeric_limits<double>::epsilon();

    int correction_max_itr = 30;
    double time_step = 1e-3;
    double spatial_step_size = 1.0;

    string input_shrink_ratio = "0-1-0-1-0-1";


    double grad_epsilon = std::numeric_limits<double>::epsilon();
    string input_function_name="quartic_potential";
    
    int max_itr=50;
    string input_model="";

    ops >> opts::Option('f', "input_function_name",  input_function_name,  " diy input function name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "degenerate_point_file", degenerate_point_file, " file name of degenerate points");
    ops >> opts::Option('z', "time_step",    time_step,       " time step size");
    ops >> opts::Option('s', "spatial_step_size",    spatial_step_size,       " spatial step size");
    ops >> opts::Option('j', "J_threshold",    J_threshold,       " Determine whether J is a zero vector");
    ops >> opts::Option('m', "input_model", input_model, " input INR model");

    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");


    ops >> opts::Option('g', "grad_epsilon", grad_epsilon, " gradient epsilon for root finding");
  
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

    printf(input_model.c_str());
    std::ifstream model_file(input_model);
    if (!model_file) {
        std::cerr << "Error: Input model file '" << input_model << "' does not exist." << std::endl;
        return 1;
    }

    INRModel<double> inr_model(input_function_name,input_model);
    int function_type=-1; 

    auto start_time = std::chrono::high_resolution_clock::now();
        

        Eigen::VectorXd local_domain_range=inr_model.domain_max-inr_model.domain_min;
        VectorXd core_maxs = inr_model.domain_max;
        VectorXd core_mins = inr_model.domain_min;

        VectorXi span_num = inr_model.block_num;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        std::vector<double> step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);
        step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time



        std::vector<VectorX<double>> root; //the inner vector store the root in a span


        VectorXi point_num_in_block = inr_model.point_num_in_block + VectorXi::Ones(inr_model.point_num_in_block.size()); //number of initial points in a block

        // VectorXd p_test(3);
        // p_test<<0.5,0.5,0.5;
        // VectorXd result;
        // inr_model.query(p_test,result);
        // std::cout<<"test query "<<result.transpose()<<std::endl;
        // VectorXi deriv(3);
        // deriv<<1,1,0;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test second derivative "<<result.transpose()<<std::endl;
        // deriv<<1,0,0;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test first derivative " <<result.transpose()<<std::endl;
        // deriv<<1,1,1;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test third derivative " <<result.transpose()<<std::endl;

        // inr_model.derivative(p_test);

        Tracking_degenerate_case<double> tracking_degenerate_case(core_mins, core_maxs, J_threshold, grad_epsilon,step_size, max_itr, function_type, nullptr, &inr_model);

        // std::vector<VectorXi> record_span;
        // choose_span(record_span);

        tracking_degenerate_case.degenerate_finding(root,point_num_in_block, span_num);

        std::cout<<root.size()<<" roots before deduplication"<<std::endl;
        if(root.empty())
        {
            auto end_time = std::chrono::high_resolution_clock::now();
            std::cout<<"degenerate case extraction time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000<<std::endl;
            return 1;
        }


        std::vector<VectorX<double>> root_unique;
        spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,step_size[0],step_size.back());
        std::cout<<"degenerate case size "<<root_unique.size()<<std::endl;

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout<<"degenerate case extraction time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000<<std::endl;

        //save roots to a file
        std::vector<MatrixXd> root_matrix(1);

        root_matrix[0].resize(root_unique.size(),root_unique[0].size());
        for(int j=0;j<root_unique.size();j++)
        {
            root_matrix[0].row(j) = root_unique[j].transpose();
        }

        string degenerate_file_name = degenerate_point_file + std::to_string(int(spatial_step_size)) + ".dat";

        utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);

        for (int i = spatial_step_size/2; i > 1; i /= 2)
        {
                root_unique.clear();
                step_size[0] *=2;
                step_size.back() *=2;
                spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique, step_size[0], step_size.back());
                root_matrix[0].resize(root_unique.size(),root_unique[0].size());
                for(int j=0;j<root_unique.size();j++)
                {
                    root_matrix[0].row(j) = root_unique[j].transpose();
                }

                degenerate_file_name = degenerate_point_file + std::to_string(i) + ".dat";
                utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);
        }


}