// #include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

// #include <diy/master.hpp>
// #include <diy/reduce-operations.hpp>
// #include <diy/decomposition.hpp>
// #include <diy/assigner.hpp>
// #include <diy/io/block.hpp>

#include <chrono>
#include <mpi.h>

#include <tbb/tbb.h>
#include <atomic>

#include "opts.h"
// #include "block.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Dense>

#include "degenerate_case.h"
#include "../utility/utility_function.h"

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

/// Gather root vectors from all MPI ranks to rank 0 (same as in critical_point_tracking_INR.cpp).
static void mpi_gather_roots(std::vector<VectorX<double>>& local_root, std::vector<VectorX<double>>& gathered,
                             int dim, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    int local_count = static_cast<int>(local_root.size());
    std::vector<int> recv_counts(static_cast<size_t>(size), 0);
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, comm);

    if (rank == 0) {
        int total = 0;
        for (int c : recv_counts) total += c;
        std::vector<double> recvbuf(static_cast<size_t>(total) * dim);
        std::vector<int> recvcounts(static_cast<size_t>(size)), displs(static_cast<size_t>(size));
        int off = 0;
        for (int i = 0; i < size; ++i) {
            recvcounts[i] = recv_counts[i] * dim;
            displs[i] = off;
            off += recvcounts[i];
        }
        std::vector<double> sendbuf(static_cast<size_t>(local_count) * dim);
        for (int i = 0; i < local_count; ++i)
            for (int d = 0; d < dim; ++d)
                sendbuf[static_cast<size_t>(i) * dim + d] = local_root[static_cast<size_t>(i)](d);
        MPI_Gatherv(sendbuf.data(), local_count * dim, MPI_DOUBLE,
                    recvbuf.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, 0, comm);
        gathered.resize(static_cast<size_t>(total));
        for (int i = 0; i < total; ++i) {
            gathered[static_cast<size_t>(i)].resize(dim);
            for (int d = 0; d < dim; ++d)
                gathered[static_cast<size_t>(i)](d) = recvbuf[static_cast<size_t>(i) * dim + d];
        }
    } else {
        std::vector<double> sendbuf(static_cast<size_t>(local_count) * dim);
        for (int i = 0; i < local_count; ++i)
            for (int d = 0; d < dim; ++d)
                sendbuf[static_cast<size_t>(i) * dim + d] = local_root[static_cast<size_t>(i)](d);
        MPI_Gatherv(sendbuf.data(), local_count * dim, MPI_DOUBLE,
                    nullptr, nullptr, nullptr, MPI_DOUBLE, 0, comm);
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
    int num_procs = 8;
    int initial_point_num_in_a_block = -1;

    ops >> opts::Option('f', "input_function_name",  input_function_name,  " diy input function name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('b', "degenerate_point_file", degenerate_point_file, " file name of degenerate points");
    ops >> opts::Option('z', "time_step",    time_step,       " time step size");
    ops >> opts::Option('s', "spatial_step_size",    spatial_step_size,       " spatial step size");
    ops >> opts::Option('j', "J_threshold",    J_threshold,       " Determine whether J is a zero vector");
    ops >> opts::Option('m', "input_model", input_model, " input INR model");
    ops >> opts::Option('n', "num_procs", num_procs, " number of processes when run without mpirun (default 8)");

    ops >> opts::Option('k', "shrink range",    input_shrink_ratio,       " shrink the range of the pointset, by \"x1-x2-y1-y2-...\"");


    ops >> opts::Option('g', "grad_epsilon", grad_epsilon, " gradient epsilon for root finding");
    ops >> opts::Option('o', "initial_point_num_in_a_block", initial_point_num_in_a_block, " initial point number in a block");
  
    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // When run as single process, optionally spawn workers (same pattern as critical_point_tracking_INR).
    MPI_Comm world_comm = MPI_COMM_WORLD;
    int world_rank = world.rank();
    int world_size = world.size();
    MPI_Comm parent_comm;
    MPI_Comm_get_parent(&parent_comm);
    if (parent_comm != MPI_COMM_NULL) {
        MPI_Intercomm_merge(parent_comm, 1, &world_comm);
        MPI_Comm_rank(world_comm, &world_rank);
        MPI_Comm_size(world_comm, &world_size);
    } else if (world_size == 1 && num_procs > 1) {
        int n_spawn = num_procs - 1;
        MPI_Comm child_comm;
        int err = MPI_Comm_spawn(argv[0], argv, n_spawn, MPI_INFO_NULL, 0, MPI_COMM_SELF, &child_comm, MPI_ERRCODES_IGNORE);
        if (err == MPI_SUCCESS) {
            MPI_Intercomm_merge(child_comm, 0, &world_comm);
            world_rank = 0;
            world_size = num_procs;
            if (world_rank == 0)
                std::cout << "Spawned " << n_spawn << " workers (total " << world_size << " processes)" << std::endl;
        }
    } else if (world_size > 1) {
        world_comm = MPI_COMM_WORLD;
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

    std::ifstream model_file(input_model);
    if (!model_file) {
        std::cerr << "Error: Input model file '" << input_model << "' does not exist." << std::endl;
        return 1;
    }

    INRModel<double> inr_model(input_function_name,input_model,initial_point_num_in_a_block);
    int function_type=-1;

    // {
    //     int nw = std::max(1, tbb::this_task_arena::max_concurrency());
    //     nw=1;
    //     inr_model.prepare_thread_modules(static_cast<size_t>(nw));
    // }

    auto start_time = std::chrono::high_resolution_clock::now();
        

        Eigen::VectorXd local_domain_range=inr_model.domain_max-inr_model.domain_min;
        VectorXd core_maxs = inr_model.domain_max;
        VectorXd core_mins = inr_model.domain_min;

        VectorXi span_num = inr_model.block_num;

        VectorXd Span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());

        std::vector<double> step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/spatial_step_size);
        step_size.back() = Span_size[Span_size.size()-1]/time_step; // the last dimension is time



        std::vector<VectorX<double>> root; //the inner vector store the root in a span


        VectorXi point_num_in_block = inr_model.point_num_in_block ;//+ VectorXi::Ones(inr_model.point_num_in_block.size()); //number of initial points in a block

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

        Tracking_degenerate_case<double> tracking_degenerate_case(core_mins, core_maxs, J_threshold, grad_epsilon, step_size, max_itr, function_type, nullptr, &inr_model);

        // Partition spans across MPI ranks (same pattern as critical_point_tracking_INR).
        int total_blocks = span_num.prod();
        VectorXi number_in_every_dim;
        utility::obtain_number_in_every_domain(span_num, number_in_every_dim);
        size_t chunk = (static_cast<size_t>(total_blocks) + static_cast<size_t>(world_size) - 1) / static_cast<size_t>(world_size);
        size_t lo = static_cast<size_t>(world_rank) * chunk;
        size_t hi = std::min(lo + chunk, static_cast<size_t>(total_blocks));
        std::vector<VectorXi> my_spans;
        my_spans.reserve(hi - lo);
        for (size_t i = lo; i < hi; ++i) {
            VectorXi block_index;
            utility::obtainDomainIndex(i, block_index, number_in_every_dim);
            my_spans.push_back(block_index);
        }
        // Empty file name so workers do not write temp files.
        tracking_degenerate_case.degenerate_finding(root, point_num_in_block, span_num, my_spans);

        if (world_size > 1) {
            std::vector<VectorX<double>> gathered;
            int dim = static_cast<int>(core_mins.size());
            mpi_gather_roots(root, gathered, dim, world_comm);
            if (world_rank == 0)
                root = std::move(gathered);
        }

        if (world_rank == 0) {
            std::cout << root.size() << " roots before deduplication" << std::endl;
        }
        int exit_early = 0;
        if (world_rank == 0 && root.empty()) exit_early = 1;
        if (world_size > 1)
            MPI_Bcast(&exit_early, 1, MPI_INT, 0, world_comm);
        if (exit_early) {
            if (world_rank == 0) {
                auto end_time = std::chrono::high_resolution_clock::now();
                std::cout << "degenerate case extraction time, millisecond : " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000 << std::endl;
            }
            MPI_Barrier(world_comm);
            return 1;
        }
        if (world_rank != 0) {
            MPI_Barrier(world_comm);
            return 0;
        }

        std::vector<VectorX<double>> root_unique;
        spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique, step_size[0], step_size.back());
        std::cout << "degenerate case size " << root_unique.size() << std::endl;

        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout << "degenerate case extraction time, millisecond : " << std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count()/1000 << std::endl;

        // Save roots to a file
        std::vector<MatrixXd> root_matrix(1);
        root_matrix[0].resize(root_unique.size(), root_unique[0].size());
        for (size_t j = 0; j < root_unique.size(); j++)
        {
            root_matrix[0].row(static_cast<Eigen::Index>(j)) = root_unique[j].transpose();
        }

        string degenerate_file_name = degenerate_point_file + std::to_string(int(spatial_step_size)) + "_" + std::to_string(initial_point_num_in_a_block) + ".dat";
        utility::writeMatrixVector(degenerate_file_name.c_str(), root_matrix);

        for (int i = spatial_step_size/2; i > 1; i /= 2)
        {
            root_unique.clear();
            step_size[0] *= 2;
            step_size.back() *= 2;
            spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique, step_size[0], step_size.back());
            root_matrix[0].resize(root_unique.size(), root_unique[0].size());
            for (size_t j = 0; j < root_unique.size(); j++)
            {
                root_matrix[0].row(static_cast<Eigen::Index>(j)) = root_unique[j].transpose();
            }
            degenerate_file_name = degenerate_point_file + std::to_string(i) + ".dat";
            utility::writeMatrixVector(degenerate_file_name.c_str(), root_matrix);
        }
    
    MPI_Barrier(world_comm);
}