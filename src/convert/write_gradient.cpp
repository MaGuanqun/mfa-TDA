//--------------------------------------------------------------
// write_gradient
//
// Writes the full gradient vector field v = grad(s) on a regular space-time
// grid in the .vff format consumed by
// src/critical_point_tracking/feature_flow_fields (the Feature Flow Fields
// tracker).
//
// Works for both model types:
//   - MFA model: a .mfa diy block file               (-f, default)
//   - INR model: a TorchScript .pt file + the matching analytical "function
//                name" (for the domain bounds), evaluated with autograd        (-i, -n)
//
// All gradient components are stored (C = D), i.e. both the spatial partials
// and the time partial. The grid is the model's domain subdivision multiplied
// by the upsample factor, x fastest, last axis = time:
//   - MFA: ndom_pts[i] = upsample_factor[i] * (#knot spans in dim i)
//   - INR: ndom_pts[i] = upsample_factor[i] * (#domain blocks in dim i)
//
// .vff layout (little-endian): "VFF1", uint32 dtype(0=f64), uint32 D, uint32 C,
//   uint32 n[D], float64 dmin[D], float64 dmax[D], then C*prod(n) float64 values,
//   AoS (components innermost): offset = (i0 + n0*(i1 + n1*(...)))*C + c.
//
// INR evaluation follows the same conventions as
// src/feature_flow_fields/gradient_norm.cpp and src/INRModel.h.
//--------------------------------------------------------------

#include    "mfa/mfa.hpp"

#include    <iostream>
#include    <fstream>
#include    <vector>
#include    <string>
#include    <cstdint>
#include    <chrono>
#include    <algorithm>

#include    <diy/master.hpp>
#include    <diy/io/block.hpp>
#include    <diy/mpi/collectives.hpp>

#include    <tbb/tbb.h>

#include    "opts.h"

#include    "block.hpp"
#include    "mfa_extend.h"
#include    "utility_function.h"
#include    "INRModel.h"

using namespace std;

// ---------------------------------------------------------------------------
// Write a gradient vector field to a .vff file.
//   data: node-major AoS buffer of size prod(ndom_pts) * C, components innermost.
// ---------------------------------------------------------------------------
template<typename T>
void write_vff_file(const string& file_name, const VectorXi& ndom_pts,
                    const VectorX<T>& mins, const VectorX<T>& maxs,
                    int C, const std::vector<double>& data)
{
    const int D = (int)ndom_pts.size();

    std::ofstream out(file_name.c_str(), std::ios::binary);
    if (!out)
    {
        std::cerr << "write_vff_file: cannot open " << file_name << " for writing" << std::endl;
        return;
    }

    out.write("VFF1", 4);
    std::uint32_t dtype = 0;                          // float64
    std::uint32_t Du    = (std::uint32_t)D;
    std::uint32_t Cu    = (std::uint32_t)C;
    out.write(reinterpret_cast<const char*>(&dtype), sizeof(dtype));
    out.write(reinterpret_cast<const char*>(&Du), sizeof(Du));
    out.write(reinterpret_cast<const char*>(&Cu), sizeof(Cu));

    std::vector<std::uint32_t> nn(D);
    for (int a = 0; a < D; a++) nn[a] = (std::uint32_t)ndom_pts(a);
    out.write(reinterpret_cast<const char*>(nn.data()), sizeof(std::uint32_t) * D);

    std::vector<double> mn(D), mx(D);
    for (int a = 0; a < D; a++) { mn[a] = (double)mins(a); mx[a] = (double)maxs(a); }
    out.write(reinterpret_cast<const char*>(mn.data()), sizeof(double) * D);
    out.write(reinterpret_cast<const char*>(mx.data()), sizeof(double) * D);

    out.write(reinterpret_cast<const char*>(data.data()), sizeof(double) * data.size());

    std::cout << "write_gradient: wrote " << file_name << " ("
              << (data.size() / std::max(1, C)) << " nodes, " << C
              << " components, D=" << D << ")" << std::endl;
}

// build the per-dimension vertex coordinates of the sampling grid
template<typename T>
void build_grid(int D, const VectorXi& ndom_pts,
                const VectorX<T>& mins, const VectorX<T>& maxs,
                std::vector<std::vector<T>>& vertex_domain)
{
    vertex_domain.assign(D, std::vector<T>());
    for (int i = 0; i < D; i++)
    {
        vertex_domain[i].resize(ndom_pts(i));
        T span = maxs(i) - mins(i);
        T d    = (ndom_pts(i) > 1) ? span / (ndom_pts(i) - 1) : T(0);
        for (int j = 0; j < ndom_pts(i); j++)
            vertex_domain[i][j] = mins(i) + T(j) * d;
    }
}

// ---------------------------------------------------------------------------
// MFA: gradient vector field on a regular grid.
// ---------------------------------------------------------------------------
template<typename T>
void save_vff_mfa(const string& file_name, Block<real_t>* block, size_t dom_dim,
                  std::vector<int>& upsample_factor)
{
    auto& tc = block->mfa->var(0).tmesh.tensor_prods[0];
    VectorXi span_num = tc.nctrl_pts - block->mfa->var(0).p;     // #knot spans per dim

    VectorXi ndom_pts(dom_dim);
    for (int i = 0; i < (int)dom_dim; i++)
        ndom_pts(i) = upsample_factor[i] * span_num(i);
    long long npts = ndom_pts.prod();

    VectorX<T> mins = block->core_mins.head(dom_dim);
    VectorX<T> maxs = block->core_maxs.head(dom_dim);

    std::vector<std::vector<T>> vertex_domain;
    build_grid<T>((int)dom_dim, ndom_pts, mins, maxs, vertex_domain);

    VectorXi number_in_every_domain(dom_dim);
    utility::obtain_number_in_every_domain(ndom_pts, number_in_every_domain);

    std::cout << "save_vff (MFA) grid (D=" << dom_dim << "): "
              << ndom_pts.transpose() << std::endl;

    // node-major AoS buffer: npts nodes, D components each (components innermost)
    std::vector<double> data(static_cast<size_t>(npts) * dom_dim);

    tbb::affinity_partitioner ap;
    tbb::parallel_for((tbb::blocked_range<size_t>(0, (size_t)npts)),
    [&](const tbb::blocked_range<size_t>& interval)
    {
        VectorXi domain_index_;
        VectorX<T> out(1);
        for (size_t j = interval.begin(); j < interval.end(); ++j)
        {
            utility::obtainDomainIndex(j, domain_index_, number_in_every_domain);

            VectorX<T> coordinate(dom_dim);
            for (int m = 0; m < (int)dom_dim; m++)
            {
                T c = vertex_domain[m][domain_index_(m)];
                if (c < mins(m)) c = mins(m);
                if (c > maxs(m)) c = maxs(m);
                coordinate[m] = c;
            }

            VectorXi deriv = VectorXi::Zero(dom_dim);
            for (int c = 0; c < (int)dom_dim; c++)
            {
                deriv[c] = 1;
                mfa_extend::recover_mfa(block, coordinate, out, deriv);
                data[j * dom_dim + c] = static_cast<double>(out[0]);
                deriv[c] = 0;
            }
        }
    }, ap);

    write_vff_file<T>(file_name, ndom_pts, mins, maxs, (int)dom_dim, data);
}

// ---------------------------------------------------------------------------
// INR: gradient vector field on a regular grid, via batched autograd.
//
// The grid nodes are distributed across MPI ranks (each rank owns a contiguous
// range of the x-fastest flattened index and loads its own copy of the model,
// so the work runs in parallel on CPU). Each rank's slice is gathered to rank 0
// which writes the .vff file. Within a rank, torch autograd is single-threaded
// (INRModel sets at::set_num_threads(1)); launch one MPI rank per core to use
// all cores, e.g.  mpirun -n <ncores> write_gradient -i model.pt -n name ...
// ---------------------------------------------------------------------------
template<typename T>
void save_vff_inr(const string& file_name, INRModel<T>& model,
                  std::vector<int>& upsample_factor,
                  const diy::mpi::communicator& world)
{
    const int  D     = (int)model.domain_min.size();
    const int  rank  = world.rank();
    const int  nproc = world.size();
    const int  root  = 0;

    VectorXi ndom_pts(D);
    for (int i = 0; i < D; i++)
        ndom_pts(i) = upsample_factor[i] * model.block_num(i);   // #domain blocks per dim
    long long npts = ndom_pts.prod();

    const VectorX<T>& mins = model.domain_min;
    const VectorX<T>& maxs = model.domain_max;

    std::vector<std::vector<T>> vertex_domain;
    build_grid<T>(D, ndom_pts, mins, maxs, vertex_domain);

    VectorXi number_in_every_domain(D);
    utility::obtain_number_in_every_domain(ndom_pts, number_in_every_domain);

    if (rank == root)
        std::cout << "save_vff (INR) grid (D=" << D << "): " << ndom_pts.transpose()
                  << " on " << nproc << " MPI rank(s)" << std::endl;

    // distribute the npts nodes contiguously across ranks (block partition)
    const long long base = npts / nproc;
    const long long rem  = npts % nproc;
    const long long my_start = rank * base + std::min<long long>(rank, rem);
    const long long my_count = base + (rank < rem ? 1 : 0);
    const long long my_end   = my_start + my_count;

    // this rank's portion of the AoS buffer: my_count nodes, D components each
    std::vector<double> local(static_cast<size_t>(my_count) * D);

    // chunked batched autograd to bound memory for large grids
    const size_t chunk = static_cast<size_t>(1) << 16;
    for (long long start = my_start; start < my_end; start += (long long)chunk)
    {
        long long end = std::min<long long>(start + (long long)chunk, my_end);

        std::vector<VectorX<T>> pts((size_t)(end - start));
        for (long long j = start; j < end; ++j)
        {
            VectorXi domain_index_;
            utility::obtainDomainIndex((size_t)j, domain_index_, number_in_every_domain);

            VectorX<T> coordinate(D);
            for (int m = 0; m < D; m++)
            {
                T c = vertex_domain[m][domain_index_(m)];
                if (c < mins(m)) c = mins(m);
                if (c > maxs(m)) c = maxs(m);
                coordinate[m] = c;
            }
            pts[(size_t)(j - start)] = coordinate;
        }

        std::vector<VectorX<T>> grads;
        model.eval_grad_batch_in_domain(pts, grads);

        for (long long j = start; j < end; ++j)
        {
            const VectorX<T>& g = grads[(size_t)(j - start)];
            const size_t off = (size_t)(j - my_start) * D;
            for (int c = 0; c < D; c++)
                local[off + c] = static_cast<double>(g[c]);
        }
    }

    // gather every rank's contiguous slice to root, in rank order
    if (nproc == 1)
    {
        write_vff_file<T>(file_name, ndom_pts, mins, maxs, D, local);
        return;
    }

    if (rank == root)
    {
        std::vector<std::vector<double>> gathered;
        diy::mpi::gather(world, local, gathered, root);

        std::vector<double> data(static_cast<size_t>(npts) * D);
        size_t pos = 0;
        for (int r = 0; r < nproc; ++r)
        {
            std::copy(gathered[r].begin(), gathered[r].end(), data.begin() + pos);
            pos += gathered[r].size();
        }
        write_vff_file<T>(file_name, ndom_pts, mins, maxs, D, data);
    }
    else
    {
        diy::mpi::gather(world, local, root);
    }
}

// expand a user upsample-factor list to exactly 'dim' entries
static void expand_upsample(std::vector<int>& uf, int dim)
{
    if (uf.empty())
        uf.push_back(1);
    if ((int)uf.size() == 1)
    {
        for (int i = 1; i < dim; ++i)
            uf.push_back(uf[0]);
    }
    else if ((int)uf.size() < dim)
    {
        for (int i = (int)uf.size(); i < dim; ++i)
            uf.push_back(1);
    }
}

int main(int argc, char** argv)
{
    using T = real_t;

    // initialize MPI
    diy::mpi::environment  env(argc, argv);       // equivalent of MPI_Init/MPI_Finalize
    diy::mpi::communicator world;                 // equivalent of MPI_COMM_WORLD

    string infile           = "approx.mfa";       // MFA diy input file (MFA mode)
    string inr_model        = "";                 // TorchScript .pt (selects INR mode)
    string inr_name         = "";                 // INR analytical function name (domain bounds)
    string output_vff_name  = "output.vff";       // output .vff gradient vector field
    string input_upsample_factor = "1";           // upsample factor "x-y-z-..."
    bool   help             = false;

    opts::Options ops;
    ops >> opts::Option('f', "infile",   infile,                 " input MFA .mfa diy file (MFA mode, default)");
    ops >> opts::Option('i', "inr",      inr_model,              " input INR TorchScript .pt model (selects INR mode)");
    ops >> opts::Option('n', "name",     inr_name,               " INR function name for domain bounds, e.g. vortex_street_3d");
    ops >> opts::Option('w', "output vff", output_vff_name,      " output .vff gradient vector field");
    ops >> opts::Option('u', "upsample", input_upsample_factor,  " upsample factor per dim \"x-y-z-...\"");
    ops >> opts::Option('h', "help",     help,                   " show help");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // parse '-' separated upsample factors
    std::vector<int> upsample_factor;
    {
        std::istringstream iuf(input_upsample_factor);
        std::string token;
        int number;
        while (std::getline(iuf, token, '-'))
        {
            std::istringstream tokenStream(token);
            if (tokenStream >> number)
                upsample_factor.push_back(number);
        }
    }

    const bool use_inr = !inr_model.empty();

    if (use_inr)
    {
        // -------- INR mode --------
        if (inr_name.empty())
        {
            std::cerr << "INR mode requires -n <function name> for the domain bounds" << std::endl;
            return 1;
        }

        INRModel<T> model(inr_name, inr_model);
        if (!model.isLoaded())
        {
            std::cerr << "failed to load INR model " << inr_model << std::endl;
            return 1;
        }

        const int D = (int)model.domain_min.size();
        if (world.rank() == 0)
            std::cout << "INR model '" << inr_name << "' domain dim " << D
                      << " min " << model.domain_min.transpose()
                      << " max " << model.domain_max.transpose() << std::endl;

        expand_upsample(upsample_factor, D);

        auto start_time = std::chrono::high_resolution_clock::now();
        save_vff_inr<T>(output_vff_name, model, upsample_factor, world);
        auto end_time = std::chrono::high_resolution_clock::now();
        if (world.rank() == 0)
            std::cout << "gradient field running time, millisecond : "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
                      << std::endl;
    }
    else
    {
        // -------- MFA mode --------
        diy::FileStorage storage("./DIY.XXXXXX");
        diy::Master      master(world,
                1,
                -1,
                &Block<real_t>::create,
                &Block<real_t>::destroy);
        diy::ContiguousAssigner assigner(world.size(), -1); // number of blocks set by read_blocks()

        diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
        std::cout << master.size() << " blocks read from file " << infile << "\n\n";

        master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
        {
            int dom_dim = b->mfa->dom_dim;
            std::vector<int> uf = upsample_factor;
            expand_upsample(uf, dom_dim);

            auto start_time = std::chrono::high_resolution_clock::now();
            save_vff_mfa<T>(output_vff_name, b, (size_t)dom_dim, uf);
            auto end_time = std::chrono::high_resolution_clock::now();
            std::cout << "gradient field running time, millisecond : "
                      << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
                      << std::endl;
        });
    }

    return 0;
}
