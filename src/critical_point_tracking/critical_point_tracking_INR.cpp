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
#include <mpi.h>

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


/// Gather root vectors from all MPI ranks to rank 0.
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

/// Scatter root_unique from rank 0 to all ranks (each gets its chunk for parallel find_trace).
static void mpi_scatter_root_unique(std::vector<VectorX<double>>& root_unique, std::vector<VectorX<double>>& my_root_unique,
                                    int dim, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    int total = (rank == 0) ? static_cast<int>(root_unique.size()) : 0;
    MPI_Bcast(&total, 1, MPI_INT, 0, comm);
    int chunk = (total + size - 1) / size;
    std::vector<int> sendcounts(static_cast<size_t>(size)), displs(static_cast<size_t>(size));
    int off = 0;
    for (int r = 0; r < size; ++r) {
        displs[r] = off;
        int n = std::min(chunk, total - r * chunk);
        if (n < 0) n = 0;
        sendcounts[r] = n * dim;
        off += sendcounts[r];
    }
    int my_count = (rank < size) ? std::min(chunk, total - rank * chunk) : 0;
    if (my_count < 0) my_count = 0;
    my_root_unique.resize(static_cast<size_t>(my_count));
    std::vector<double> sendbuf;
    if (rank == 0)
        sendbuf.resize(static_cast<size_t>(total) * dim);
    if (rank == 0)
        for (int i = 0; i < total; ++i)
            for (int d = 0; d < dim; ++d)
                sendbuf[i * dim + d] = root_unique[static_cast<size_t>(i)](d);
    std::vector<double> recvbuf(static_cast<size_t>(my_count) * dim);
    MPI_Scatterv(rank == 0 ? sendbuf.data() : nullptr, sendcounts.data(), displs.data(), MPI_DOUBLE,
                 recvbuf.data(), my_count * dim, MPI_DOUBLE, 0, comm);
    for (int i = 0; i < my_count; ++i) {
        my_root_unique[static_cast<size_t>(i)].resize(dim);
        for (int d = 0; d < dim; ++d)
            my_root_unique[static_cast<size_t>(i)](d) = recvbuf[static_cast<size_t>(i) * dim + d];
    }
}

/// Gather traces from all ranks to rank 0 (serialize CP_Trace to flat buffer: n_pts, pts, connect_info, dup).
static void mpi_gather_traces(std::vector<CP_Trace<double>>& local_traces, std::vector<CP_Trace<double>>& gathered,
                              MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    const int dim = 3;
    std::vector<double> sendbuf;
    size_t pos = 0;
    auto write_int = [&](int x) { sendbuf.resize(pos + 1); sendbuf[pos++] = static_cast<double>(x); };
    auto write_dbl = [&](double x) { sendbuf.resize(pos + 1); sendbuf[pos++] = x; };
    for (const auto& tr : local_traces) {
        int n = static_cast<int>(tr.traces.size());
        write_int(n);
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < dim; ++d)
                write_dbl(tr.traces[static_cast<size_t>(i)](d));
        write_int(tr.connect_info[0]);
        write_int(tr.connect_info[1]);
        write_int(tr.duplicated ? 1 : 0);
    }
    int sendcount = static_cast<int>(sendbuf.size());
    std::vector<int> recvcounts(static_cast<size_t>(size));
    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, comm);
    std::vector<int> displs(static_cast<size_t>(size));
    int total = 0;
    if (rank == 0) {
        for (int r = 0; r < size; ++r) { displs[r] = total; total += recvcounts[r]; }
    }
    std::vector<double> recvbuf(rank == 0 ? static_cast<size_t>(total) : 0);
    MPI_Gatherv(sendbuf.data(), sendcount, MPI_DOUBLE, recvbuf.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, 0, comm);
    if (rank == 0) {
        gathered.clear();
        size_t idx = 0;
        for (int r = 0; r < size; ++r) {
            size_t end = idx + static_cast<size_t>(recvcounts[r]);
            while (idx < end) {
                CP_Trace<double> tr;
                int n = static_cast<int>(recvbuf[idx++]);
                tr.traces.resize(static_cast<size_t>(n));
                for (int i = 0; i < n; ++i) {
                    tr.traces[static_cast<size_t>(i)].resize(dim);
                    for (int d = 0; d < dim; ++d)
                        tr.traces[static_cast<size_t>(i)](d) = recvbuf[idx++];
                }
                tr.connect_info[0] = static_cast<int>(recvbuf[idx++]);
                tr.connect_info[1] = static_cast<int>(recvbuf[idx++]);
                tr.duplicated = (static_cast<int>(recvbuf[idx++]) != 0);
                gathered.push_back(std::move(tr));
            }
        }
    }
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

    int max_itr=100;

    double point_itr_threshold = 0.5;

    double  spatial_step_size = 1.0;
        string input_model="";
    string edge_type_file ="";

    string boundary_start="";
    int compute_boundary_start=1;
    int num_procs = 8;  // number of processes for parallel root_finding (when run as single process)
    int initial_point_num_in_a_block = -1;

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
    ops >> opts::Option('n', "num_procs", num_procs, " number of processes (when run without mpirun; default 8)");
    ops >> opts::Option('o', "initial_point_num_in_a_block", initial_point_num_in_a_block, " initial point number in a block");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // When run as single process (e.g. ./critical_point_tracking_INR ...), spawn workers and merge into one comm.
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


    INRModel<double> inr_model(input_function_name,input_model,initial_point_num_in_a_block);
    int function_type=-1; 

    // // Vector of modules, one per TBB worker; each thread uses thread_modules_[its index] with no per-call clone/lock.
    // {
    //     int nw = std::max(1, tbb::this_task_arena::max_concurrency());
    //     nw=1;
    //     // std::cout<<"number of TBB workers "<<nw<<std::endl;
    //     inr_model.prepare_thread_modules(static_cast<size_t>(nw));
    // }

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

    if (world_rank == 0){
    std::cout<<Span_size.transpose()<<std::endl;
    std::cout<<"spatial step size "<<step_size[0]<<" "<<"time step " <<step_size.back()<<std::endl;
    }

        int spanned_block_num =span_num.prod();

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);


        std::vector<VectorX<double>> root; //the inner vector store the root in a span

        std::vector<VectorXi> selected_span;
        span_filter::compute_boundary_span(span_num,selected_span,true);
        if (world_rank == 0)
            std::cout<<"valid span num "<<selected_span.size()<<std::endl;
        // std::vector<VectorXi> selected_span;

        auto cpt_extract_start_time = std::chrono::high_resolution_clock::now();


        VectorXi point_num_in_block = inr_model.point_num_in_block; //number of initial points in a block
        std::cout<<"point_num_in_block "<<point_num_in_block.transpose()<<std::endl;
        std::vector<VectorX<double>> root_unique;
        Find_boundary_roots find_boundary_roots(root_finding_grad_epsilon,core_mins,core_maxs,point_num_in_block,span_num,same_root_epsilon,function_type,max_itr,point_itr_threshold,static_cast<Block<double>*>(nullptr),&inr_model);

        boundary_start = boundary_start +"_" + std::to_string(initial_point_num_in_a_block) + "_";

        if(compute_boundary_start==1){
            if (world_size > 1) {
                // MPI: partition spans across ranks, each process does its chunk
                size_t chunk = (selected_span.size() + static_cast<size_t>(world_size) - 1) / world_size;
                size_t lo = world_rank * chunk;
                size_t hi = std::min(lo + chunk, selected_span.size());
                std::vector<VectorXi> my_spans;
                if (lo < hi)
                    my_spans.assign(selected_span.begin() + lo, selected_span.begin() + hi);
                find_boundary_roots.root_finding(my_spans, root);
                std::vector<VectorX<double>> gathered;
                int dim = static_cast<int>(core_mins.size());
                mpi_gather_roots(root, gathered, dim, world_comm);
                if (world_rank == 0)
                    root = std::move(gathered);
            } else {
                find_boundary_roots.root_finding(selected_span, root);
            }

            if (world_rank == 0) {
                spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,same_root_epsilon[0],same_root_epsilon.back());

                save_root(root_unique, boundary_start, spatial_step_size);

                auto temp_step_size= step_size;
                std::cout<<"find root num before deduplicate between spans "<<root.size()<<std::endl;
                std::cout<<"finish finding root before deduplicate between spans "<<root.size()<<" after "<<root_unique.size()<<std::endl;
                for (int i = spatial_step_size/2; i > 1; i /= 2)
                {
                    root_unique.clear();
                    temp_step_size[0] *=2;
                    temp_step_size.back() *=2;
                    spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,temp_step_size[0],temp_step_size.back());
                    save_root(root_unique, boundary_start, i);
                    std::cout<<"finish finding root before deduplicate between spans "<<i<<" "<<root.size()<<" after "<<root_unique.size()<<std::endl;
                }

                std::vector<double> new_step_size(Span_size.size(),Span_size.head(Span_size.size()-1).minCoeff()/64.0);
                new_step_size.back() = Span_size[Span_size.size()-1]/64.0; // the last dimension is time
                temp_step_size= new_step_size;
                for (int i = 64; i > 1; i /= 2)
                {
                    root_unique.clear();
                    spatial_hashing_spatial_temporal::find_all_unique_root(root, root_unique,temp_step_size[0],temp_step_size.back());
                    save_root(root_unique, boundary_start, i);
                    std::cout<<"finish finding root before deduplicate between spans "<<i<<" "<<root.size()<<" after "<<root_unique.size()<<std::endl;
                    temp_step_size[0] *=2;
                    temp_step_size.back() *=2;
                }
            }
            root.clear();
            root.shrink_to_fit();
        }
        else
        {
            if (world_rank == 0) {
                string name  =  boundary_start  + std::to_string(int(spatial_step_size)) + ".dat";
                Degenerate_case_tracing<double>::read_degenerate_point(name,root_unique);

            }
        }

        auto finding_end_time = std::chrono::high_resolution_clock::now();

        if (world_rank == 0)
            std::cout<<"finding time, millisecond : "<<std::chrono::duration_cast<std::chrono::microseconds>(finding_end_time - cpt_extract_start_time).count()/1000<<std::endl;


            // MPI_Barrier(world_comm);
            // return 0;
        // string test_file=cp_tracing_file+"_test.obj";

        // tracking_utility::convert_to_obj(test_file,root_unique);


        auto tracking_start_time = std::chrono::high_resolution_clock::now();

        if (world_size > 1) {
            // MPI: scatter root_unique, each rank runs find_trace on its chunk, gather traces
            std::vector<VectorX<double>> my_root_unique;
            int dim = static_cast<int>(core_mins.size());
            mpi_scatter_root_unique(root_unique, my_root_unique, dim, world_comm);
            std::vector<CP_Trace<double>> my_traces;
            Boundary_critical_point_tracking boundary_critical_point_tracking(core_mins, core_maxs, root_finding_grad_epsilon, step_size.back(), step_size[0], d_max_square_, function_type, correction_max_itr, static_cast<Block<double>*>(nullptr), &inr_model);
            my_traces.resize(my_root_unique.size());
            boundary_critical_point_tracking.find_trace(my_root_unique, my_traces);

            // mpi_gather_traces(my_traces, traces, world_comm);
            // if (world_rank == 0)
            // {
            //     std::cout << "finish boundary critical point tracing (MPI)" << std::endl;
            // }

            std::vector<VectorX<double>> my_degenerate_points;
            mpi_scatter_root_unique(degenerate_points, my_degenerate_points, dim, world_comm);
            Degenerate_case_tracing degenerate_case_tracing(core_mins, core_maxs, point_num_in_block, &find_boundary_roots, step_size.back(), step_size[0], root_finding_grad_epsilon, correction_max_itr, function_type, static_cast<Block<double>*>(nullptr), &inr_model);
            degenerate_case_tracing.tracing_from_all_degenerate_points(my_degenerate_points, my_traces, 0.1, d_max_square_);
            mpi_gather_traces(my_traces, traces, world_comm);
            if (world_rank == 0) {
                std::cout << "finish tracing (MPI)" << std::endl;
            }
                
        } else {
            traces.resize(root_unique.size());
            Boundary_critical_point_tracking boundary_critical_point_tracking(core_mins, core_maxs, root_finding_grad_epsilon, step_size.back(), step_size[0], d_max_square_, function_type, correction_max_itr, static_cast<Block<double>*>(nullptr), &inr_model);
            boundary_critical_point_tracking.find_trace(root_unique, traces);

            Degenerate_case_tracing degenerate_case_tracing(core_mins, core_maxs, point_num_in_block, &find_boundary_roots, step_size.back(), step_size[0], root_finding_grad_epsilon, correction_max_itr, function_type, static_cast<Block<double>*>(nullptr), &inr_model);
            degenerate_case_tracing.tracing_from_all_degenerate_points(degenerate_points, traces, 0.1, d_max_square_);
        }


    if (world_rank != 0) {
        MPI_Barrier(world_comm);
        return 0;
    }

    if (world_rank == 0) {
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
        deduplication::deduplicate_traces(traces, degenerate_points, step_size[0], step_size.back(), core_mins, core_maxs);

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
    }
    MPI_Barrier(world_comm);
    return 0;
}