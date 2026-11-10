#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace local_reference {
using Triple = std::array<double, 3>;

struct Dataset {
    Triple min, max, base_intervals;
};

// Count nonzero spans and reject nonuniform knots: a VTI needs one spacing
// per axis, so it cannot assign the same point count to unequal spans.
inline std::size_t uniform_span_count(const std::vector<double>& knots) {
    if(knots.size() < 2 || knots.front() != 0. || knots.back() != 1.)
        throw std::invalid_argument("MFA knots must cover [0,1]");
    std::vector<double> widths;
    for(std::size_t i = 1; i < knots.size(); ++i) {
        const double width = knots[i] - knots[i-1];
        if(!std::isfinite(width) || width < 0)
            throw std::invalid_argument("Invalid MFA knot sequence");
        if(width > 0) widths.push_back(width);
    }
    if(widths.empty()) throw std::invalid_argument("MFA has no nonzero knot spans");
    const double expected = 1. / widths.size();
    for(double width : widths)
        if(std::abs(width - expected) > 1e-10 * expected)
            throw std::invalid_argument("--points-per-span requires uniform MFA knot spans; use --ratio for nonuniform models");
    return widths.size();
}

inline Dataset dataset(const std::string& name) {
    if(name == "boussinesq_3d")
        return {{-.5, -.5, 0}, {.5, 2.5, 1.5}, {10, 30, 15}};
    if(name == "vortex_street_3d")
        return {{-.5, -.5, 13.5}, {7.5, .5, 15}, {80, 10, 15}};
    throw std::invalid_argument("--dataset must be boussinesq_3d or vortex_street_3d");
}

struct Grid {
    Triple min, max, spacing;
    std::array<std::size_t, 3> dimensions;
    std::size_t size = 1;

    Grid(const Dataset& domain, const Triple& region_min, const Triple& region_max,
         const Triple& ratios, const Triple& padding, bool align_to_domain = false) {
        for(int a = 0; a < 3; ++a) {
            if(!std::isfinite(region_min[a]) || !std::isfinite(region_max[a]) ||
               region_min[a] >= region_max[a] || region_min[a] < domain.min[a] ||
               region_max[a] > domain.max[a])
                throw std::invalid_argument("Region must have increasing bounds inside the original dataset domain (x, y, t)");
            if(!std::isfinite(ratios[a]) || ratios[a] <= 0 ||
               !std::isfinite(padding[a]) || padding[a] < 0)
                throw std::invalid_argument("Sampling ratios must be positive; padding must be nonnegative");
            min[a] = std::max(domain.min[a], region_min[a] - padding[a]);
            max[a] = std::min(domain.max[a], region_max[a] + padding[a]);
            const double nominal = (domain.max[a] - domain.min[a]) /
                                   domain.base_intervals[a] / ratios[a];
            double intervals = (max[a] - min[a]) / nominal;
            if(align_to_domain) {
                const double total = domain.base_intervals[a] * ratios[a];
                if(!std::isfinite(total) || total > std::numeric_limits<int>::max() - 1. ||
                   total < (a == 2 ? 2 : 1) || total != std::floor(total))
                    throw std::invalid_argument("Invalid number of full-domain sampling intervals");
                const double eps = 64 * std::numeric_limits<double>::epsilon() * std::max(1., total);
                double first = std::max(0., std::floor((min[a]-domain.min[a])/nominal + eps));
                double last = std::min(total, std::ceil((max[a]-domain.min[a])/nominal - eps));
                const double required = a == 2 ? 2 : 1;
                if(last-first < required) {
                    last = std::min(total, first+required);
                    first = std::max(0., last-required);
                }
                min[a] = domain.min[a] + first*nominal;
                max[a] = last == total ? domain.max[a] : domain.min[a] + last*nominal;
                intervals = last-first;
            }
            // Remove roundoff near integer counts, and include both exact endpoints.
            // Legacy persistence tracking needs at least three frames.
            const double count = std::max(a == 2 ? 2. : 1., std::ceil(intervals -
                32 * std::numeric_limits<double>::epsilon() * std::max(1., intervals)));
            if(!std::isfinite(count) || count > std::numeric_limits<int>::max() - 1.)
                throw std::invalid_argument("Sampling grid axis exceeds VTK integer dimensions");
            dimensions[a] = static_cast<std::size_t>(count) + 1;
            spacing[a] = align_to_domain ? nominal : (max[a] - min[a]) / count;
            if(size > std::numeric_limits<std::size_t>::max() / dimensions[a] / sizeof(double))
                throw std::invalid_argument("Sampling grid size overflows addressable storage");
            size *= dimensions[a];
        }
    }

    // File layout is (nt, ny, nx), with x changing fastest.
    Triple point(std::size_t index) const {
        Triple p;
        for(int a = 0; a < 3; ++a) {
            const auto j = index % dimensions[a];
            index /= dimensions[a];
            p[a] = j == dimensions[a] - 1 ? max[a] : min[a] + j * spacing[a];
        }
        return p;
    }
};
} // namespace local_reference
