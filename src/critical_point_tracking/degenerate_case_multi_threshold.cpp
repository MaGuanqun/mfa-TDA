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


#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include "save_control_data.hpp"
#include "degenerate_case_multi_threshold.h"
#include "critical_point/span_filter.hpp"
#include "tracking_utility.h"
#include "spatial_hashing_spatial_temporal.h"



using namespace std;


// namespace {
//     tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
// }

int main(int argc, char** argv)
{

    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file


    string inControlPoint = "derivative_control_point.dat";

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;


    double shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string degenerate_point_file = "degenerate_point";
    
    double J_threshold = 1e-5;

    int max_itr = 50;
    double time_step = 1e-3;
    double spatial_step_size = 1.0;

    real_t hessian_threshold = 1e-20;

    string input_shrink_ratio = "0-1-0-1-0-1";

    real_t point_itr_threshold = 0.5;

    real_t grad_epsilon = 1e-8;

    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "degenerate_point_file", degenerate_point_file, " file name of degenerate points");
    ops >> opts::Option('z', "time_step",    time_step,       " time step size");
    ops >> opts::Option('s', "spatial_step_size",    spatial_step_size,       " spatial step size");
    ops >> opts::Option('j', "J_threshold",    J_threshold,       " Determine whether J is a zero vector");

    ops >> opts::Option('a', "inControlPoint",  inControlPoint,  " diy input derivative control point file name");
    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");

    ops >> opts::Option('p', "point_itr_threshold", point_itr_threshold, " stop iteration when point is away from block center than point_itr_threshold * block size");

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

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<double> step_size;

    std::vector<double> J_threshold_list;
    J_threshold_list.push_back(J_threshold);
    for(int i=1;i<7;i++)
    {
        J_threshold_list.push_back(J_threshold*std::pow(10,-i));
    }


    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;

        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        std::vector<double> step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);
        step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

        std::vector<std::vector<VectorX<double>>> root; //the inner vector store the root in a span

        std::vector<std::vector<VectorXi>> selected_span;//[vars,span index]

        std::vector<size_t> valid_span_index;


        span_filter::compute_valid_span(sci_deriv_control_points,b,selected_span,shrink_ratio,2,true);

        VectorXi point_num_in_block = b->mfa->var(0).p + VectorXi::Ones(b->mfa->var(0).p.size()); //number of initial points in a block

        selected_span[0].resize(100);


        Degenerate_case_multi_threshold<double> tracking_degenerate_case(b->core_mins, b->core_maxs, J_threshold_list, grad_epsilon,step_size, max_itr,0,b);
        tracking_degenerate_case.degenerate_finding(root,
            point_num_in_block, span_num, selected_span[0]);

        std::cout<<"end finding degenerate case"<<std::endl;

        std::vector<std::vector<VectorX<double>>> root_unique(root.size());
        for(int i=0;i<root.size();i++)
        {
            spatial_hashing_spatial_temporal::find_all_unique_root(root[i], root_unique[i], step_size[0], step_size.back());

            std::cout<<"degenerate case size "<< J_threshold_list[i] << " "<<root_unique[i].size()<<std::endl;
        }
        // spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,step_size[0],step_size.back());
       

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout<<"degenerate case extraction time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000<<std::endl;

        //save roots to a file
        std::vector<MatrixXd> root_matrix(1);

        for (size_t k = 0; k < root_unique.size(); k++)
        {

        root_matrix[0].resize(root_unique[k].size(),root_unique[k][0].size());
        for(int j=0;j<root_unique[k].size();j++)
        {
            root_matrix[0].row(j) = root_unique[k][j].transpose();
        }

        string degenerate_file_name = degenerate_point_file + std::to_string(int(spatial_step_size)) + "_"+ std::to_string(-int(log10(J_threshold_list[k]))) + ".dat";

        // string degenerate_file_name = degenerate_point_file + std::to_string(int(spatial_step_size)) + ".dat";

        utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);



        for (int i = spatial_step_size/2; i > 1; i /= 2)
        {
                root_unique.clear();
                step_size[0] *=2;
                step_size.back() *=2;
                spatial_hashing_spatial_temporal::find_all_unique_root(root[k], root_unique[k], step_size[0], step_size.back());
                root_matrix[0].resize(root_unique[k].size(),root_unique[k][0].size());
                for(int j=0;j<root_unique[k].size();j++)
                {
                    root_matrix[0].row(j) = root_unique[k][j].transpose();
                }

                degenerate_file_name = degenerate_point_file + std::to_string(i) + "_"+ std::to_string(-int(log10(J_threshold_list[k]))) + ".dat";
                utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);
        }

        step_size[0]=Span_size.head(Span_size.size()-1).minCoeff()/64.0;
        step_size.back() = Span_size[Span_size.size()-1]/64.0; 
        for (int i = 64; i > 1; i /= 2)
        {
       
            root_unique.clear();
            step_size[0] *=2;
            step_size.back() *=2;
            spatial_hashing_spatial_temporal::find_all_unique_root(root[k], root_unique[k], step_size[0], step_size.back());
            root_matrix[0].resize(root_unique[k].size(),root_unique[k][0].size());
            for(int j=0;j<root_unique[k].size();j++)
            {
                root_matrix[0].row(j) = root_unique[k][j].transpose();
            }

            degenerate_file_name = degenerate_point_file + std::to_string(i) + "_"+ std::to_string(-int(log10(J_threshold_list[k]))) + ".dat";
            utility::writeMatrixVector(degenerate_file_name.c_str(),root_matrix);

        }
    }

    });

}