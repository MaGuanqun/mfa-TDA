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
#include "degenerate_case.h"
#include "critical_point/span_filter.hpp"
#include "tracking_utility.h"



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

    string degenerate_point_file = "degenerate_point.dat";



    double tolerance = 1e-4; // same_root_epsilon
    
    double J_threshold = 1e-5;

    int correction_max_itr = 30;
    double step_size = 1e-3;

    real_t hessian_threshold = 1e-20;

    string input_shrink_ratio = "0-1-0-1-0-1";



    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "degenerate_point_file", degenerate_point_file, " file name of degenerate points");
    ops >> opts::Option('z', "step_size",    step_size,       " step size");
    ops >> opts::Option('j', "J_threshold",    J_threshold,       " Determine whether J is a zero vector");

    ops >> opts::Option('a', "inControlPoint",  inControlPoint,  " diy input derivative control point file name");
    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");
  
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

    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
    {
        
        Eigen::VectorXd local_domain_range=b->core_maxs-b->core_mins;

        std::vector<VectorX<double>> root_ori;

        double min_ = local_domain_range.minCoeff();
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());
        double min_span_size = Span_size[2]/step_size;
        step_size = min_span_size;


        tolerance = 0.5*step_size; // same_root_epsilon

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

        std::vector<VectorX<real_t>> root; //the inner vector store the root in a span

        std::vector<std::vector<VectorXi>> selected_span;//[vars,span index]

        std::vector<size_t> valid_span_index;


        span_filter::compute_valid_span(sci_deriv_control_points,b,selected_span,shrink_ratio,2);

        tbb::enumerable_thread_specific<std::vector<VectorXd>> local_root;

        std::vector<int>multi_root_span; 
        Eigen::VectorXd weights=Eigen::VectorXd::Ones(b->mfa->var(0).tmesh.tensor_prods[0].ctrl_pts.rows());



        tbb::affinity_partitioner ap;

        std::cout<<"start find degenerate point "<<std::endl;

        tbb::parallel_for(tbb::blocked_range<size_t>(0,selected_span[0].size()), //
        [&](const tbb::blocked_range<size_t>& range)
        {
            auto& root_thread = local_root.local();

            for(auto i=range.begin();i!=range.end();++i)
            {

            // for(auto i=0;i!=selected_span[0].size();++i)
            // {
                if(i%1000==0)
                {
                    std::cout<<"find span "<<i<<std::endl;
                }
                // std::cout<<"find span "<<i<<std::endl;
                // std::cout<<selected_span[0][i][2]<<std::endl;
                std::vector<VectorX<real_t>> root_block;
                root_block.reserve(16);
                // std::vector<VectorXd> root_span;
                if(cp_tracking_degenerate_case::degenerate_finding(b,selected_span,root_block,i,J_threshold, tolerance,hessian_threshold))
                {

                    root_thread.insert(root_thread.end(), root_block.begin(), root_block.end());

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

        
        for (const auto& thread_vec : local_root) {
            root.insert(root.end(), thread_vec.begin(), thread_vec.end());
        }

        std::cout<<"degenerate case size "<<root.size()<<std::endl;

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout<<"degenerate case extraction time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000<<std::endl;

        //save roots to a file
        std::vector<MatrixXd> root_matrix(1);

        root_matrix[0].resize(root.size(),root[0].size());
        for(int j=0;j<root.size();j++)
        {
            root_matrix[0].row(j) = root[j].transpose();
        }

        utility::writeMatrixVector(degenerate_point_file.c_str(),root_matrix);

    });

}