//
// Stable Feature Flow Fields critical-point tracking
// (Weinkauf, Theisel, Van Gelder & Pang, IEEE TVCG 2011).
//
// Discrete reimplementation built on top of the plain Feature Flow Fields
// driver. Integrates stream lines of the *stable* feature flow field
// h = (+/-)f + tau * g, which adds an attracting correction g around the
// feature line so numerical drift is automatically corrected during the
// ongoing integration (see stable_feature_flow_fields.h for the math).
//
// The plain FFF tool (feature_flow_fields) is kept separately and unchanged.
//
// Inputs:
//   -f  .vff  discrete vector field (gradient field) on a regular space-time grid
//   -s  .csv  seed critical points (header line + rows of x,y,[z,]t)
//   -b        output base name (writes <base>.obj)
//   -k        convergence strength (paper's k > 0; 0 => plain FFF behavior)
//
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <functional>
#include <list>
#include <map>

#include "opts.h"
#include "stable_feature_flow_fields.h"

using namespace std;

int main(int argc, char** argv)
{
    using T = double;

    string vff_file    = "field.vff";
    string seed_file   = "seeds.csv";
    string out_base    = "sfff_tracking";

    double step         = -1.0;  // arc-length step; <0 => derive from grid spacing
    int    max_steps    = 5000;
    double spatial_eps  = -1.0;  // <0 => derive from grid spacing
    double temporal_eps = -1.0;  // <0 => derive from grid spacing
    double k_strength   = 1.0;   // paper's k > 0 (convergence strength); 0 => plain FFF
    double tau_max      = 100.0; // upper clamp for the adaptive tau (paper used 100)
    double fd_frac      = 0.5;   // finite-difference step as a fraction of grid spacing
    double loop_eps     = -1.0;  // loop-closure threshold; <0 => derive from step (0 disables)
    int    self_skip    = -1;    // trailing window (steps) for spiral/self-revisit detection; <0 => derive (0 disables)
    bool   help         = false;

    opts::Options ops;
    ops >> opts::Option('f', "vff",          vff_file,     " input .vff vector field (gradient field)");
    ops >> opts::Option('s', "seeds",        seed_file,    " input seed critical-point .csv");
    ops >> opts::Option('b', "out",          out_base,     " output base name (writes <base>.obj)");
    ops >> opts::Option('l', "step",         step,         " RK4 arc-length step (default: min spatial spacing)");
    ops >> opts::Option('m', "max_steps",    max_steps,    " max RK4 steps per direction");
    ops >> opts::Option('e', "spatial_eps",  spatial_eps,  " spatial epsilon for covered-seed test (default: spatial spacing)");
    ops >> opts::Option('t', "temporal_eps", temporal_eps, " temporal epsilon for covered-seed test (default: time spacing)");
    ops >> opts::Option('k', "strength",     k_strength,   " stable-FFF convergence strength k>0 (0 => plain FFF)");
    ops >> opts::Option('x', "tau_max",      tau_max,      " upper clamp for the adaptive correction strength tau");
    ops >> opts::Option('d', "fd",           fd_frac,      " finite-difference step as a fraction of grid spacing");
    ops >> opts::Option('p', "loop_eps",     loop_eps,     " loop-closure threshold; stop+close when a trace returns within this of its seed (default: step; 0 disables)");
    ops >> opts::Option('w', "self_skip",    self_skip,    " trailing window (steps) for spiral self-revisit detection; stop when a trace re-enters an earlier-visited cell (default: derived; 0 disables)");
    ops >> opts::Option('h', "help",         help,         " show help");

    if (!ops.parse(argc, argv) || help)
    {
        std::cout << ops;
        return 1;
    }

    sfff::VectorFieldGrid<T> grid;
    if (!sfff::vff_io::load(vff_file, grid))
        return 1;

    std::cout << "loaded vff: D=" << grid.D << " C=" << grid.C << " n=[";
    for (int a = 0; a < grid.D; ++a) std::cout << grid.n[a] << (a + 1 < grid.D ? "," : "");
    std::cout << "]\n  dmin=" << grid.dmin.transpose() << "\n  dmax=" << grid.dmax.transpose()
              << "\n  spacing=" << grid.spacing.transpose() << std::endl;

    // default step / eps from grid spacing (spatial axes = all but the last)
    double min_spatial = grid.spacing.head(grid.D - 1).minCoeff();
    double time_space  = grid.spacing[grid.D - 1];
    if (step         < 0) step         = min_spatial;
    if (spatial_eps  < 0) spatial_eps  = min_spatial;
    if (temporal_eps < 0) temporal_eps = time_space;
    if (loop_eps     < 0) loop_eps     = step;        // closed-loop detection on by default
    if (self_skip    < 0)                              // spiral self-revisit detection on by default
    {
        // Exclude enough of the trailing path that the locally-adjacent points
        // (which sit within ~1 occupancy cell of the head) cannot self-match;
        // a few cell widths of arc length is sufficient.
        double cell = std::max(spatial_eps, temporal_eps);
        self_skip = std::max(4, static_cast<int>(std::ceil(3.0 * cell / step)));
    }

    std::cout << "step=" << step << " max_steps=" << max_steps
              << " spatial_eps=" << spatial_eps << " temporal_eps=" << temporal_eps
              << " loop_eps=" << loop_eps << (loop_eps > 0 ? "" : " (seed-return detection disabled)")
              << " self_skip=" << self_skip << (self_skip > 0 ? "" : " (spiral detection disabled)") << std::endl;
    std::cout << "stable FFF: k=" << k_strength << " tau_max=" << tau_max
              << " fd_frac=" << fd_frac;
    if (k_strength <= 0.0)
        std::cout << "  (k<=0 -> integrating plain FFF)";
    else if (grid.D != 3)
        std::cout << "  (D!=3 -> stabilization unavailable, integrating plain FFF)";
    std::cout << std::endl;

    std::vector<VectorX<T>> seeds;
    if (!sfff::read_seeds_csv(seed_file, grid.D, seeds))
        return 1;
    std::cout << "read " << seeds.size() << " seed critical points from " << seed_file << std::endl;
    if (seeds.empty())
    {
        std::cerr << "no seeds -- nothing to track" << std::endl;
        return 1;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<CP_Trace<T>> traces;
    sfff::track_all(grid, seeds, static_cast<T>(step), max_steps,
                    static_cast<T>(spatial_eps), static_cast<T>(temporal_eps),
                    static_cast<T>(k_strength), static_cast<T>(tau_max),
                    static_cast<T>(fd_frac), static_cast<T>(loop_eps), self_skip, traces);

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
