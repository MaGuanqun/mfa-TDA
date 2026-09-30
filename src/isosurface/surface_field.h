#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace marching_triangles {
template<class T> using Point = Eigen::Matrix<T, 3, 1>;

template<class T> struct SurfaceField {
    std::function<T(const Point<T>&)> value;
    std::function<Point<T>(const Point<T>&)> gradient;
};

template<class T> struct MeshingOptions {
    Point<T> domain_min = Point<T>::Constant(-4);
    Point<T> domain_max = Point<T>::Constant(4);
    T iso_value = 0;
    T step = T(.2);
    T min_step = T(.02);
    T projection_tolerance = T(1e-8);
    T vertex_tolerance = T(1e-6);
    T projection_alpha = T(1.5);
    T stop_distance = 0;
    int projection_iterations = 80;
    int quality_iterations = 8;
    size_t max_triangles = 200000;
    bool close_cracks = true;

    void validate() const {
        if (!domain_min.allFinite() || !domain_max.allFinite() ||
            !(domain_max.array() > domain_min.array()).all() ||
            !std::isfinite(iso_value) || !std::isfinite(step) || step <= 0 ||
            !std::isfinite(min_step) || min_step <= 0 || min_step > step ||
            !std::isfinite(projection_tolerance) || projection_tolerance <= 0 ||
            !std::isfinite(vertex_tolerance) || vertex_tolerance <= 0 || vertex_tolerance >= step ||
            !std::isfinite(projection_alpha) || projection_alpha <= 1 || projection_alpha >= 2 ||
            !std::isfinite(stop_distance) || stop_distance < 0 ||
            projection_iterations <= 0 || quality_iterations < 0 || quality_iterations > 100 || max_triangles < 4)
            throw std::invalid_argument("Invalid meshing bounds, step, tolerance, or iteration limit");
    }
};

// Section 3.2.1: gradient march with alpha > 1 to bracket the root, then bisection.
// Backtracking safeguards the march against singular gradients and distant roots.
template<class T> class SurfaceProjector {
public:
    SurfaceProjector(const SurfaceField<T>& field, const MeshingOptions<T>& options)
        : field_(field), options_(options) {}

    bool in_domain(const Point<T>& p) const {
        return p.allFinite() && (p.array() >= options_.domain_min.array()).all() &&
               (p.array() <= options_.domain_max.array()).all();
    }
    bool regular(const Point<T>& p) const {
        const Point<T> g = field_.gradient(p);
        return g.allFinite() && g.stableNorm() > std::numeric_limits<T>::min();
    }
    bool project(Point<T>& p, T max_displacement = std::numeric_limits<T>::infinity()) const {
        if (!p.allFinite()) return false;
        const Point<T> origin = p;
        Point<T> x = p;
        T fx = field_.value(x) - options_.iso_value;
        for (int i = 0; i < options_.projection_iterations; ++i) {
            if (!std::isfinite(fx)) return false;
            if (std::abs(fx) <= options_.projection_tolerance) {
                if (!in_domain(x) || !regular(x)) return false;
                p = x;
                return true;
            }
            Point<T> g = field_.gradient(x);
            const T norm = g.stableNorm();
            if (!g.allFinite() || norm <= std::numeric_limits<T>::min()) return false;
            g /= norm;
            T step = options_.projection_alpha * fx / norm;
            if (!std::isfinite(step)) return false;
            Point<T> y;
            T fy = std::numeric_limits<T>::quiet_NaN();
            bool accepted = false;
            for (int j = 0; j < 40; ++j) {
                y = x - step * g;
                if (y.allFinite() && (y - origin).norm() <= max_displacement) {
                    fy = field_.value(y) - options_.iso_value;
                    if (std::isfinite(fy) && (std::signbit(fx) != std::signbit(fy) || std::abs(fy) < std::abs(fx))) {
                        accepted = true;
                        break;
                    }
                }
                step *= T(.5);
            }
            if (!accepted) return false;
            if (std::signbit(fx) != std::signbit(fy)) {
                for (++i; i < options_.projection_iterations; ++i) {
                    Point<T> mid = (x + y) * T(.5);
                    const T fm = field_.value(mid) - options_.iso_value;
                    if (!std::isfinite(fm)) return false;
                    if (std::abs(fm) <= options_.projection_tolerance) {
                        if (!in_domain(mid) || !regular(mid)) return false;
                        p = mid;
                        return true;
                    }
                    if (std::signbit(fm) == std::signbit(fx)) { x = mid; fx = fm; }
                    else { y = mid; fy = fm; }
                }
                return false;
            }
            x = y;
            fx = fy;
        }
        return false;
    }
private:
    const SurfaceField<T>& field_;
    const MeshingOptions<T>& options_;
};
} // namespace marching_triangles
