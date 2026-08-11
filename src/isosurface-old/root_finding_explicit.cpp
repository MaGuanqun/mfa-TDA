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

#include "save_control_data.hpp"
#include "find_initial_root.h"
#include "tracking_utility.h"
#include "closed_form_function.h"
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



using namespace std;


// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }

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



    std::vector<double> same_root_epsilon; // same_root_epsilon
    double time_step = 1e-3;
    
    double grad_threshold = 1e-5;

    int correction_max_itr = 30;




    real_t initial_point_finding_hessian_threshold = 1e-20;
    real_t root_finding_grad_epsilon = 1e-8;

    string input_shrink_ratio = "0-1-0-1-0-1";
    real_t dxy_dt_gradient_epsilon = 1e-10;

    string root_file = "root.dat";

    int max_itr=50;


    double  spatial_step_size = 1.0;

    double root_finding_epsilon = 1e-8;

    double function_value = 0.0;

    ops >> opts::Option('f', "input_function_name",  input_function_name,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('g', "spatial_step_size",    spatial_step_size,       " spatial step size");

    ops >> opts::Option('x', "root_finding_epsilon",    root_finding_epsilon,       "root finding epsilon");

    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
    ops >> opts::Option('m', "max_itr", max_itr, " max iteration");
    ops >> opts::Option('s', "root_file", root_file, " save obtained root");
    ops >> opts::Option('v', "function_value", function_value, " function value");



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

    int function_type = closed_form_function::initial_func_type(input_function_name);

    Eigen::VectorXd local_domain_range=closed_form_function::domain_max(function_type)-closed_form_function::domain_min(function_type);
    VectorXd core_maxs = closed_form_function::domain_max(function_type);
    VectorXd core_mins = closed_form_function::domain_min(function_type);

    VectorXi span_num = closed_form_function::block_num(function_type);

    VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());


    double step_size = Span_size.minCoeff()/spatial_step_size;
    
    std::cout<<"span_size " <<Span_size.transpose()<< "step size "<<step_size<<std::endl;

    VectorXi num_of_initial_point = closed_form_function::point_num_in_block(function_type);
    for (int i = 0; i < num_of_initial_point.size(); ++i)
        num_of_initial_point[i] *= span_num[i];
   
    Find_initial_root<double> find_initial_root(core_mins, core_maxs, step_size, root_finding_epsilon, max_itr, function_type, nullptr, nullptr, num_of_initial_point);


    std::vector<VectorX<double>> root;
    
    find_initial_root.root_finding(root, function_value);

    std::cout<<"root num "<<root.size()<<std::endl;

    find_initial_root.test_function_value(root);
    
     tracking_utility::save_root(root, root_file, static_cast<int>(spatial_step_size));

    return 0;
}