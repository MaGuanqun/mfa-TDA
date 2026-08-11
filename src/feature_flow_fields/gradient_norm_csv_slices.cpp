//
// Spatial gradient-norm statistics over a series of per-time-slice CSVs.
//
// Each CSV stores 3D critical points (x,y,z) at one time index; the time
// coordinate is taken from the file index and mapped onto the MFA time domain
// [core_mins[D-1], core_maxs[D-1]]:
//   t(i) = t_min + i * (t_max - t_min) / (last - first)
// so index `first` -> t_min and index `last` -> t_max.
//
// For every point across all slices, evaluates ||(fx,fy,fz)|| using the MFA
// (same recover_mfa scheme as gradient_norm.cpp) and reports min/max/avg over
// the pooled set of all points.
//
// Inputs:
//   -f  .mfa   diy MFA input file
//   -p  prefix path without "_<index>.csv", e.g. build/.../vortex.mfa_16
//   -a  first slice index (default 0)
//   -b  last  slice index inclusive (required)
//   -v  if 1, also print per-slice averages
//
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <limits>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include <tbb/tbb.h>

#include "opts.h"
#include "block.hpp"
#include "mfa_extend.h"

using namespace std;

// .csv : one point per row (comma separated). A non-numeric first line is
// treated as a header and skipped. Only the first ncol columns are used.
template<typename T>
bool read_points_csv(const std::string& filename, int ncol, std::vector<VectorX<T>>& points)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
        std::cerr << "read_points_csv: cannot open " << filename << std::endl;
        return false;
    }

    points.clear();
    std::string line;
    bool first = true;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        std::vector<T> vals;
        std::stringstream ss(line);
        std::string tok;
        bool numeric = true;
        while (std::getline(ss, tok, ','))
        {
            try { vals.push_back(static_cast<T>(std::stod(tok))); }
            catch (...) { numeric = false; break; }
        }

        if (!numeric)
        {
            if (first) { first = false; continue; }
            continue;
        }
        first = false;

        if (static_cast<int>(vals.size()) < ncol) continue;
        VectorX<T> p(ncol);
        for (int a = 0; a < ncol; ++a) p[a] = vals[a];
        points.emplace_back(std::move(p));
    }
    return true;
}

// MFA: spatial gradient norm ||(fx,...,f_{D-2})|| (drops last/time axis).
template<typename T>
T mfa_spatial_grad_norm(const Block<T>* b, const VectorX<T>& vertex)
{
    const int D = static_cast<int>(b->core_mins.size());

    VectorX<T> p(D);
    for (int d = 0; d < D; ++d)
    {
        T c = (d < vertex.size()) ? vertex[d] : T(0);
        if (c < b->core_mins[d]) c = b->core_mins[d];
        if (c > b->core_maxs[d]) c = b->core_maxs[d];
        p[d] = c;
    }

    VectorX<T> grad(D);
    VectorX<T> out(1);
    VectorXi   deriv = VectorXi::Zero(D);
    for (int d = 0; d < D; ++d)
    {
        deriv[d] = 1;
        mfa_extend::recover_mfa(b, p, out, deriv);
        grad[d] = out[0];
        deriv[d] = 0;
    }

    const int n = std::max(1, D - 1);
    return grad.head(n).norm();
}

int main(int argc, char** argv)
{
    using T = double;

    diy::mpi::environment  env(argc, argv);
    diy::mpi::communicator world;

    string mfa_file   = "approx.mfa";
    string csv_prefix = "";
    int    first_idx  = 0;
    int    last_idx   = -1;
    int    verbose    = 0;
    bool   help       = false;

    opts::Options ops;
    ops >> opts::Option('f', "mfa",     mfa_file,   " input MFA .mfa diy file");
    ops >> opts::Option('p', "prefix",  csv_prefix, " CSV path prefix (files are <prefix>_<i>.csv)");
    ops >> opts::Option('a', "first",   first_idx,  " first slice index (mapped to t_min)");
    ops >> opts::Option('b', "last",    last_idx,   " last slice index inclusive (mapped to t_max)");
    ops >> opts::Option('v', "verbose", verbose,    " 1: also print per-slice averages");
    ops >> opts::Option('h', "help",    help,       " show help");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    if (csv_prefix.empty() || last_idx < first_idx)
    {
        std::cerr << "need -p <prefix> and -b <last> with last >= first (="
                  << first_idx << ")" << std::endl;
        return 1;
    }

    const int num_intervals = last_idx - first_idx; // 0 => all slices at t_min

    diy::FileStorage storage("./DIY.XXXXXX");
    diy::Master      master(world, -1, -1,
                            &Block<T>::create, &Block<T>::destroy,
                            &storage, &Block<T>::save, &Block<T>::load);
    diy::ContiguousAssigner assigner(world.size(), -1);

    diy::io::read_blocks(mfa_file.c_str(), world, assigner, master, &Block<T>::load);
    std::cout << master.size() << " blocks read from file " << mfa_file << std::endl;
    if (master.size() == 0)
    {
        std::cerr << "no blocks in " << mfa_file << std::endl;
        return 1;
    }

    T global_sum = T(0);
    T global_min = std::numeric_limits<T>::infinity();
    T global_max = -std::numeric_limits<T>::infinity();
    size_t global_count = 0;
    size_t files_ok = 0;

    auto t0 = std::chrono::high_resolution_clock::now();

    master.foreach([&](Block<T>* b, const diy::Master::ProxyWithLink&)
    {
        const int D = static_cast<int>(b->core_mins.size());
        if (D < 2)
        {
            std::cerr << "MFA domain dim " << D << " < 2; need a time axis" << std::endl;
            return;
        }

        const T t_min = b->core_mins[D - 1];
        const T t_max = b->core_maxs[D - 1];
        std::cout << "MFA domain min " << b->core_mins.transpose()
                  << " max " << b->core_maxs.transpose() << std::endl;
        std::cout << "time mapping: index " << first_idx << " -> " << t_min
                  << ", index " << last_idx << " -> " << t_max << std::endl;

        const int n_slices = last_idx - first_idx + 1;
        std::vector<std::vector<VectorX<T>>> per_slice(static_cast<size_t>(n_slices));
        std::vector<T>             slice_t(static_cast<size_t>(n_slices), t_min);
        std::vector<std::string>   slice_path(static_cast<size_t>(n_slices));
        std::vector<char>          slice_ok(static_cast<size_t>(n_slices), 0);

        // Parallel CSV I/O + attach time coordinate.
        tbb::parallel_for(0, n_slices, [&](int k)
        {
            const int idx = first_idx + k;
            slice_path[static_cast<size_t>(k)] =
                csv_prefix + "_" + std::to_string(idx) + ".csv";
            slice_t[static_cast<size_t>(k)] = (num_intervals == 0)
                ? t_min
                : t_min + static_cast<T>(idx - first_idx)
                          * (t_max - t_min) / static_cast<T>(num_intervals);

            std::vector<VectorX<T>> pts3;
            if (!read_points_csv<T>(slice_path[static_cast<size_t>(k)], 3, pts3)
                || pts3.empty())
                return;

            auto& out = per_slice[static_cast<size_t>(k)];
            out.resize(pts3.size());
            const T t = slice_t[static_cast<size_t>(k)];
            for (size_t i = 0; i < pts3.size(); ++i)
            {
                VectorX<T> p4(D);
                for (int d = 0; d < D - 1; ++d)
                    p4[d] = (d < pts3[i].size()) ? pts3[i][d] : T(0);
                p4[D - 1] = t;
                out[i] = std::move(p4);
            }
            slice_ok[static_cast<size_t>(k)] = 1;
        });

        // Flatten into one point list; remember exclusive end index per slice.
        std::vector<VectorX<T>> points;
        std::vector<size_t> slice_end(static_cast<size_t>(n_slices), 0);
        size_t total = 0;
        for (int k = 0; k < n_slices; ++k)
        {
            if (slice_ok[static_cast<size_t>(k)])
            {
                ++files_ok;
                total += per_slice[static_cast<size_t>(k)].size();
            }
            slice_end[static_cast<size_t>(k)] = total;
        }
        points.reserve(total);
        for (int k = 0; k < n_slices; ++k)
        {
            auto& s = per_slice[static_cast<size_t>(k)];
            points.insert(points.end(),
                          std::make_move_iterator(s.begin()),
                          std::make_move_iterator(s.end()));
            s.clear();
            s.shrink_to_fit();
        }

        if (points.empty())
            return;

        std::vector<T> norms(points.size(), T(0));

        // Parallel MFA gradient-norm evaluation + reduction.
        struct Stats
        {
            T sum = T(0);
            T minv = std::numeric_limits<T>::infinity();
            T maxv = -std::numeric_limits<T>::infinity();

            Stats() = default;
            Stats(Stats&, tbb::split) {}

            void join(const Stats& o)
            {
                sum += o.sum;
                minv = std::min(minv, o.minv);
                maxv = std::max(maxv, o.maxv);
            }
        };

        Stats stats = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(0, points.size()),
            Stats{},
            [&](const tbb::blocked_range<size_t>& range, Stats local) -> Stats
            {
                for (size_t i = range.begin(); i != range.end(); ++i)
                {
                    const T nrm = mfa_spatial_grad_norm<T>(b, points[i]);
                    norms[i] = nrm;
                    local.sum += nrm;
                    local.minv = std::min(local.minv, nrm);
                    local.maxv = std::max(local.maxv, nrm);
                }
                return local;
            },
            [](Stats a, const Stats& b) -> Stats
            {
                a.join(b);
                return a;
            });

        global_sum = stats.sum;
        global_min = stats.minv;
        global_max = stats.maxv;
        global_count = points.size();

        if (verbose)
        {
            size_t start = 0;
            for (int k = 0; k < n_slices; ++k)
            {
                const size_t end = slice_end[static_cast<size_t>(k)];
                if (!slice_ok[static_cast<size_t>(k)] || end == start)
                {
                    if (!slice_ok[static_cast<size_t>(k)])
                        std::cerr << "skipping missing/unreadable "
                                  << slice_path[static_cast<size_t>(k)] << std::endl;
                    else if (verbose)
                        std::cout << slice_path[static_cast<size_t>(k)]
                                  << ": 0 points" << std::endl;
                    start = end;
                    continue;
                }

                T slice_sum = T(0);
                T slice_min = norms[start];
                T slice_max = norms[start];
                for (size_t i = start; i < end; ++i)
                {
                    slice_sum += norms[i];
                    slice_min = std::min(slice_min, norms[i]);
                    slice_max = std::max(slice_max, norms[i]);
                }
                const size_t n = end - start;
                std::cout << slice_path[static_cast<size_t>(k)]
                          << ": n=" << n
                          << " t=" << slice_t[static_cast<size_t>(k)]
                          << " avg=" << (slice_sum / static_cast<T>(n))
                          << " min=" << slice_min
                          << " max=" << slice_max << std::endl;
                start = end;
            }
        }
    });

    auto t1 = std::chrono::high_resolution_clock::now();

    if (global_count == 0)
    {
        std::cerr << "no points read from any CSV" << std::endl;
        return 1;
    }

    const T avg = global_sum / static_cast<T>(global_count);
    std::cout << "spatial gradient norm over " << global_count
              << " points from " << files_ok << " CSV files:" << std::endl;
    std::cout << "  min = " << global_min << std::endl;
    std::cout << "  max = " << global_max << std::endl;
    std::cout << "  avg = " << avg << std::endl;
    std::cout << "evaluation time (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << std::endl;

    return 0;
}
