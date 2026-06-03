#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include <chrono>
#include <tbb/tbb.h>

#include "opts.h"
#include "block.hpp"
#include "closed_form_function.h"
#include "marching_triangles.h"
#include "../critical_point_tracking/degenerate_case_tracing.h"

#include <fstream>
#include <Eigen/Dense>

using namespace std;

int main(int argc, char** argv)
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    string input_function_name = "ellipsoid";
    bool help = false;
    opts::Options ops;

    string root_filename = "root.dat";
    string singular_point_file = "";
    string output_mesh_prefix = "isosurface_sheet";
    string input_shrink_ratio = "0-1-0-1-0-1";

    double spatial_step_size = 1.0;
    double function_value = 0.0;
    double root_finding_epsilon = 1e-8;
    double same_vertex_epsilon = -1.0;
    double hessian_rank_threshold = 1e-10;
    double stop_curve_distance = 0.05;
    double d_min_ratio = 0.1;
    int max_projection_itr = 50;
    int max_sheets = 100;
    bool enable_crack_closing = true;

    ops >> opts::Option('f', "input_function_name", input_function_name, "closed-form function name");
    ops >> opts::Option('h', "help", help, "show help");
    ops >> opts::Option('g', "spatial_step_size", spatial_step_size, "spatial step size (used in root filename)");
    ops >> opts::Option('s', "root_file", root_filename, "root file from root_finding");
    ops >> opts::Option('d', "singular_point_file", singular_point_file, "degenerate curve points (optional)");
    ops >> opts::Option('o', "output_mesh_prefix", output_mesh_prefix, "output OBJ prefix per sheet");
    ops >> opts::Option('v', "function_value", function_value, "isovalue");
    ops >> opts::Option('x', "root_finding_epsilon", root_finding_epsilon, "surface projection tolerance");
    ops >> opts::Option('H', "hessian_rank_threshold", hessian_rank_threshold, "min |eigenvalue| for full-rank seed");
    ops >> opts::Option('c', "stop_curve_distance", stop_curve_distance, "stop growth near degenerate points");
    ops >> opts::Option('n', "d_min_ratio", d_min_ratio, "minimum marching step ratio of normal step size");
    ops >> opts::Option('m', "max_projection_itr", max_projection_itr, "max surface projection iterations");
    ops >> opts::Option('k', "shrink range", input_shrink_ratio, "unused, kept for CLI compatibility");
    ops >> opts::Option('r', "enable_crack_closing", enable_crack_closing, "run Akkouche crack-fixing after growth");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    int function_type = closed_form_function::initial_func_type(input_function_name);
    VectorXd core_maxs = closed_form_function::domain_max(function_type);
    VectorXd core_mins = closed_form_function::domain_min(function_type);
    VectorXi span_num = closed_form_function::block_num(function_type);
    VectorXd local_domain_range = core_maxs - core_mins;
    VectorXd span_size = local_domain_range.cwiseQuotient(span_num.cast<double>());
    double step_size = span_size.minCoeff() / spatial_step_size;




    std::vector<VectorX<double>> seed_roots;

    Degenerate_case_tracing<double>::read_degenerate_point(root_filename,seed_roots);

    std::cout << "loaded " << seed_roots.size() << " roots from " << root_filename << std::endl;

    if (seed_roots.empty())
    {
        std::cerr << "No roots found. Run root_finding_explicit first." << std::endl;
        return 1;
    }

    double d_min = step_size * d_min_ratio;

    if (same_vertex_epsilon < 0)
        same_vertex_epsilon = d_min * 0.25;

    marching_triangles::MarchingTriangles<double> mesher(
        core_mins, core_maxs, function_type, function_value,
        root_finding_epsilon, max_projection_itr, same_vertex_epsilon,
        hessian_rank_threshold, stop_curve_distance, d_min,step_size,
        1.5, nullptr, nullptr, enable_crack_closing);

    // mesher.set_degenerate_points(degenerate_points);

    std::vector<std::vector<VectorX<double>>> sheet_vertices;
    std::vector<std::vector<marching_triangles::Triangle<double>>> sheet_triangles;

    if (!mesher.extract_all_sheets(seed_roots, sheet_vertices, sheet_triangles))
    {
        std::cerr << "Marching triangles produced no sheets." << std::endl;
        return 1;
    }

    int written = 0;
    for (size_t s = 0; s < sheet_vertices.size() && written < max_sheets; ++s)
    {
        const std::string out_name = output_mesh_prefix + std::to_string(s) + ".obj";
        marching_triangles::MarchingTriangles<double>::save_mesh_obj(
            out_name, sheet_vertices[s], sheet_triangles[s]);
        std::cout << "wrote " << out_name << std::endl;
        ++written;
    }

    return 0;
}
