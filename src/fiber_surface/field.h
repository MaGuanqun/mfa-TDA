#pragma once

#include "../isosurface/mesh_geometry.h"
#include <string>

namespace fiber_surface {
using Vec = Eigen::Vector3d;
using Range = Eigen::Vector2d;
using Jacobian = Eigen::Matrix<double, 2, 3>;
using Hessian = Eigen::Matrix3d;
using Mesh = marching_triangles::TriangleMesh<double>;
using marching_triangles::Triangle;
using marching_triangles::Edge;
using marching_triangles::edge_key;

struct Box {
    Vec min = Vec::Constant(-2), max = Vec::Constant(2);
    bool contains(const Vec& x) const {
        return x.allFinite() && (x.array() >= min.array()).all() && (x.array() <= max.array()).all();
    }
    unsigned faces(const Vec& x, double tolerance) const {
        unsigned mask = 0;
        for (int d = 0; d < 3; ++d) {
            if (std::abs(x[d] - min[d]) <= tolerance) mask |= 1u << (2*d);
            if (std::abs(x[d] - max[d]) <= tolerance) mask |= 1u << (2*d+1);
        }
        return mask;
    }
    void fix(Vec& x, unsigned faces) const {
        for (int d = 0; d < 3; ++d) {
            if (faces & (1u << (2*d))) x[d] = min[d];
            if (faces & (1u << (2*d+1))) x[d] = max[d];
        }
    }
};

// Derivatives are with respect to shared physical coordinates, never MFA parameters.
struct Field {
    std::function<Range(const Vec&)> value;
    std::function<Jacobian(const Vec&)> jacobian;
    std::function<std::array<Hessian, 2>(const Vec&)> hessians;
    Box box;
};

struct Segment {
    Range a, b, d, n, w;
    double length;
    Segment(Range first = Range(1, -2), Range last = Range(1, .5)) : a(first), b(last) {
        d = b-a; length = d.norm();
        if (!a.allFinite() || !b.allFinite() || !std::isfinite(length) || length <= 0)
            throw std::invalid_argument("Range segment endpoints must be finite and distinct");
        n = Range(-d[1], d[0])/length; w = d/d.squaredNorm();
    }
    double phi(const Range& f) const { return n.dot(f-a); }
    double parameter(const Range& f) const { return w.dot(f-a); }
    Vec gradient(const Jacobian& j) const { return j.transpose()*n; }
    Vec parameter_gradient(const Jacobian& j) const { return j.transpose()*w; }
    Hessian hessian(const std::array<Hessian, 2>& h) const { return n[0]*h[0]+n[1]*h[1]; }
    Hessian parameter_hessian(const std::array<Hessian, 2>& h) const { return w[0]*h[0]+w[1]*h[1]; }
};

struct Options {
    double step = .2, min_step = .025, fiber_step = .15, spatial_tolerance = .015;
    double residual_tolerance = 1e-9, derivative_tolerance = 1e-10;
    double min_angle_degrees = .1, max_normal_degrees = 55;
    int seed_grid = 7, iterations = 40, max_trace_steps = 10000;
    size_t max_triangles = 30000;
    void validate(const Field& f, const Segment& segment) const {
        for (double x : {step, min_step, fiber_step, spatial_tolerance, residual_tolerance, derivative_tolerance})
            if (!std::isfinite(x) || x <= 0) throw std::invalid_argument("Steps and tolerances must be positive and finite");
        if (min_step > step || seed_grid < 2 || seed_grid > 100 || iterations < 1 || max_trace_steps < 4 ||
            max_triangles < 4 || !std::isfinite(min_angle_degrees) || min_angle_degrees < 0 || min_angle_degrees >= 60 ||
            !std::isfinite(max_normal_degrees) || max_normal_degrees <= 0 || max_normal_degrees >= 90 ||
            residual_tolerance >= .01*segment.length)
            throw std::invalid_argument("Invalid meshing limits (range residual tolerance must be much smaller than segment length)");
        if (!f.value || !f.jacobian || !f.hessians || !f.box.min.allFinite() || !f.box.max.allFinite() ||
            (f.box.max.array() <= f.box.min.array()).any()) throw std::invalid_argument("Invalid field or domain");
    }
};

struct Event { Vec position; std::string reason; };
struct VertexInfo { int endpoint = -1, fiber = -1; unsigned faces = 0; };
struct Fiber { int endpoint = -1; std::vector<int> vertices; bool closed = false, complete = false; };
struct Result {
    Mesh mesh;
    std::vector<VertexInfo> info;
    std::vector<Fiber> fibers;
    std::map<Edge, int> fiber_edges;
    std::set<Edge> box_edges;
    std::vector<Event> unresolved;
    std::array<int, 2> fiber_seed_hits{0,0};
    int surface_seed_hits = 0;
    double seconds = 0;
};

struct Seeds {
    std::vector<Vec> surface;
    std::array<std::vector<Vec>,2> endpoints;
    size_t starting_points = 0;
    double seconds = 0;
};

// The restart state includes the virtual fiber uses as well as the mesh.
// Ordinary exposed halfedges are reconstructed from oriented triangles.
struct Checkpoint {
    Result result;
    std::map<Edge,int> pending, fiber_uses;
};

Seeds find_seeds(const Field& field, const Segment& segment, const Options& options);
Result trace_endpoint_fibers(const Field& field, const Segment& segment,
                             const Options& options, const Seeds& seeds);

inline std::vector<Vec> discovery_grid(const Box& box, int n) {
    std::vector<Vec> points;
    for (int k=0; k<n; ++k) for (int j=0; j<n; ++j) for (int i=0; i<n; ++i)
        points.push_back(box.min+(box.max-box.min).cwiseProduct(Vec(i,j,k)/double(n-1)));
    return points;
}
} // namespace fiber_surface
