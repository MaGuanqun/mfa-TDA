// Sample the existing analytical quartic_potential_2, without fitting a model.
// VFF1 stores (df/dx, df/dy, df/dt), x fastest and time slowest. SFFF uses
// the first two components and differentiates the multilinear interpolant.
#include <iostream>
#include <string>
#include <stdexcept>
#include <functional>
#include <list>
#include <map>

#include "opts.h"
#include "feature_flow_fields.h"

using std::string;
using Eigen::VectorXd;
using Eigen::VectorXi;
#include "closed_form_function.h"

int main(int argc, char** argv)
{
    int nx = 100, ny = 100, nt = 100;
    string output = "quartic_potential_2_gradient.vff";
    bool help = false;
    opts::Options options;
    options >> opts::Option('x', "nx", nx, " number of x samples, including endpoints (default 100)");
    options >> opts::Option('y', "ny", ny, " number of y samples, including endpoints (default 100)");
    options >> opts::Option('t', "nt", nt, " number of time samples, including endpoints (default 100)");
    options >> opts::Option('o', "output", output, " output float64 .vff file");
    options >> opts::Option('h', "help", help, " show help");
    if (!options.parse(argc, argv)) { std::cerr << options; return 1; }
    if (help) { std::cout << options; return 0; }
    if (nx < 2 || ny < 2 || nt < 2)
    {
        std::cerr << "Each grid dimension must be at least 2.\n";
        return 1;
    }

    try
    {
        fff::VectorFieldGrid<double> grid;
        grid.D = 3;
        grid.C = 3;
        grid.n = {nx, ny, nt};
        const int function_type = closed_form_function::initial_func_type("quartic_potential_2");
        grid.dmin = closed_form_function::domain_min(function_type);
        grid.dmax = closed_form_function::domain_max(function_type);
        grid.compute_spacing();
        std::size_t count = grid.C;
        for (int n : grid.n)
        {
            if (count > grid.data.max_size() / static_cast<std::size_t>(n))
                throw std::overflow_error("Requested grid is too large");
            count *= static_cast<std::size_t>(n);
        }
        grid.data.resize(count);
        VectorX<double> point(3), value(1);
        VectorXi derivative = VectorXi::Zero(3);
        std::size_t offset = 0;
        for (int it = 0; it < nt; ++it)
        for (int iy = 0; iy < ny; ++iy)
        for (int ix = 0; ix < nx; ++ix)
        {
            point << grid.dmin[0] + ix * grid.spacing[0],
                     grid.dmin[1] + iy * grid.spacing[1],
                     grid.dmin[2] + it * grid.spacing[2];
            for (int component = 0; component < grid.C; ++component)
            {
                derivative.setZero();
                derivative[component] = 1;
                closed_form_function::quartic_potential_2(point, value, derivative);
                grid.data[offset++] = value[0];
            }
        }
        if (!fff::vff_io::save(output, grid)) return 1;
        std::cout << "Sampled quartic_potential_2 on " << nx << " x " << ny << " x " << nt
                  << " nodes\n  domain: " << grid.dmin.transpose() << " -> " << grid.dmax.transpose()
                  << "\n  spacing: " << grid.spacing.transpose()
                  << "\n  components: df/dx, df/dy, df/dt\nWrote " << output << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << "Sampling failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
