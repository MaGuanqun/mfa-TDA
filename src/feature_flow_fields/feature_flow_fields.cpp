//
// Feature Flow Fields critical-point tracking (Theisel & Seidel, VisSym 2003).
//
// Discrete reimplementation: integrate stream lines of the feature flow field f
// (built from a sampled vector field v = grad_space(s) via multilinear
// interpolation and the analytic in-cell derivatives of that interpolant),
// seeded from critical points extracted at the discrete time steps. No Newton
// correction step (as in the original paper).
//
// Inputs:
//   -f  .vff  discrete vector field (gradient field) on a regular space-time grid
//   -s  .csv  seed critical points (header line + rows of x,y,[z,]t)
//   -b        output base name (writes <base>.obj)
//
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <functional>
#include <list>
#include <map>

#include "opts.h"
// degenerate_case_tracing.h pulls in <mfa/mfa.hpp>, which defines the global
// VectorX/MatrixX alias templates that feature_flow_fields.h relies on. It must
// be included BEFORE feature_flow_fields.h (which no longer does
// `using namespace Eigen;`) so those names resolve unambiguously.
#include "degenerate_case_tracing.h"
#include "feature_flow_fields.h"

using namespace std;

int main(int argc, char** argv)
{
    using T = double;

    string vff_file    = "field.vff";
    string seed_file   = "seeds.csv";
    string out_base    = "fff_tracking";

    double step        = -1.0;   // arc-length step; <0 => derive from grid spacing
    int    max_steps   = 5000;
    double spatial_eps = -1.0;   // <0 => derive from grid spacing
    double temporal_eps = -1.0;  // <0 => derive from grid spacing
    bool   help        = false;

    opts::Options ops;
    ops >> opts::Option('f', "vff",        vff_file,    " input .vff vector field (gradient field)");
    ops >> opts::Option('s', "seeds",      seed_file,   " input seed critical-point .dat");
    ops >> opts::Option('b', "out",        out_base,    " output base name (writes <base>.obj)");
    ops >> opts::Option('l', "step",       step,        " RK4 arc-length step (default: min spatial spacing)");
    ops >> opts::Option('m', "max_steps",  max_steps,   " max RK4 steps per direction");
    ops >> opts::Option('e', "spatial_eps", spatial_eps, " spatial epsilon for covered-seed test (default: spatial spacing)");
    ops >> opts::Option('t', "temporal_eps", temporal_eps, " temporal epsilon for covered-seed test (default: time spacing)");
    ops >> opts::Option('h', "help",       help,        " show help");

    if (!ops.parse(argc, argv) || help)
    {
        std::cout << ops;
        return 1;
    }

    fff::VectorFieldGrid<T> grid;
    if (!fff::vff_io::load(vff_file, grid))
        return 1;

    std::cout << "loaded vff: D=" << grid.D << " C=" << grid.C << " n=[";
    for (int a = 0; a < grid.D; ++a) std::cout << grid.n[a] << (a + 1 < grid.D ? "," : "");
    std::cout << "]\n  dmin=" << grid.dmin.transpose() << "\n  dmax=" << grid.dmax.transpose()
              << "\n  spacing=" << grid.spacing.transpose() << std::endl;

    // default step / eps from grid spacing (spatial axes = all but the last)
    double min_spatial = grid.spacing.head(grid.D - 1).minCoeff();
    double time_space  = grid.spacing[grid.D - 1];
    if (step        < 0) step        = min_spatial;
    if (spatial_eps < 0) spatial_eps = min_spatial;
    if (temporal_eps < 0) temporal_eps = time_space;

    std::cout << "step=" << step << " max_steps=" << max_steps
              << " spatial_eps=" << spatial_eps << " temporal_eps=" << temporal_eps << std::endl;

    std::vector<VectorX<T>> seeds;
    // if (!fff::read_seeds_csv(seed_file, grid.D, seeds))
    //     return 1;
    Degenerate_case_tracing<double>::read_degenerate_point(seed_file,seeds);
    std::cout << "read " << seeds.size() << " seed critical points from " << seed_file << std::endl;
    if (seeds.empty())
    {
        std::cerr << "no seeds -- nothing to track" << std::endl;
        return 1;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<CP_Trace<T>> traces;
    fff::track_all(grid, seeds, static_cast<T>(step), max_steps,
                   static_cast<T>(spatial_eps), static_cast<T>(temporal_eps), traces);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "tracking time (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;

    // write output as .obj (reuse the project's trace writer; empty degenerate points)
    std::vector<VectorX<T>> degenerate_points;
    VectorX<T> domain_min   = grid.dmin;
    VectorX<T> domain_range = grid.dmax - grid.dmin;

    string out_file = out_base + ".obj";
    CP_Trace_fuc::convert_to_obj(out_file, traces, degenerate_points, domain_min, domain_range);
    std::cout << "wrote " << out_file << std::endl;

    return 0;
}
