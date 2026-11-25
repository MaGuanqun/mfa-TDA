#pragma once

#include "../isosurface/mesh_geometry.h"
#include <optional>
#include <string>

namespace fiber_surface {
using Vec = Eigen::Vector3d;
using Range = Eigen::Vector2d;
using Jacobian = Eigen::Matrix<double, 2, 3>;
using Hessian = Eigen::Matrix3d;
using SeedShape = std::array<int, 3>;
using Mesh = marching_triangles::TriangleMesh<double>;
using marching_triangles::Triangle;
using marching_triangles::Edge;
using marching_triangles::edge_key;

inline Vec edge_midpoint(const Mesh& mesh, int u, int v) {
    const auto [i,j]=edge_key(u,v);
    const Vec& a=mesh.vertices[i];
    return a+(mesh.vertices[j]-a)/2;
}

inline Vec face_centroid(const Mesh& mesh, const Triangle& t) {
    const Vec& a=mesh.vertices[t[0]];
    return a+((mesh.vertices[t[1]]-a)+(mesh.vertices[t[2]]-a))/3;
}

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
// Evaluation callbacks must support concurrent calls during seed discovery.
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
    double parameter(const Range& f) const { return d.dot(f-a)/d.squaredNorm(); }
    static bool contains_parameter(double tau) { return 0 <= tau && tau <= 1; }
    Vec gradient(const Jacobian& j) const { return j.transpose()*n; }
    Vec parameter_gradient(const Jacobian& j) const { return j.transpose()*w; }
    Hessian hessian(const std::array<Hessian, 2>& h) const { return n[0]*h[0]+n[1]*h[1]; }
    Hessian parameter_hessian(const std::array<Hessian, 2>& h) const { return w[0]*h[0]+w[1]*h[1]; }
};

// A simple directed range curve. A repeated first point closes the loop.
struct RangeCurve {
    std::vector<Segment> segments;
    RangeCurve(const Segment& segment = {}) : segments{segment} {}
    explicit RangeCurve(std::vector<Range> points);
    int size() const { return static_cast<int>(segments.size()); }
    bool closed() const { return segments.back().b == segments.front().a; }
    int vertex_count() const { return size() + (closed() ? 0 : 1); }
    int corner(int side, int endpoint) const { return (side + endpoint) % vertex_count(); }
    const Range& point(int i) const { return i == size() ? segments.back().b : segments[i].a; }
};

struct Options {
    static constexpr double fiber_field_tolerance = 1e-9;
    static constexpr int default_max_seed_starts = 1000000;
    double step = .15, min_step = step/8, spatial_tolerance = .075*step;
    std::optional<double> fiber_seed_radius;
    double residual_tolerance = 1e-9, derivative_tolerance = 1e-10;
    // Optional tighter residual tolerance for surface roots on adjacent patches.
    double surface_solve_tolerance = 0;
    double min_angle_degrees = .1, max_normal_degrees = 55;
    // Count on the shortest physical axis for every field; zero is an alias for 7.
    int seed_grid = 7, iterations = 40, max_trace_steps = 10000;
    int max_seed_starts = default_max_seed_starts;
    size_t max_triangles = 30000;
    double seed_coverage_radius() const { return fiber_seed_radius.value_or(spatial_tolerance); }
    double fiber_plane_tolerance() const { return 1e-6 * step; }
    void validate(const Field& f, const Segment& segment) const {
        for (double x : {step, min_step, spatial_tolerance, seed_coverage_radius(), residual_tolerance, derivative_tolerance})
            if (!std::isfinite(x) || x <= 0) throw std::invalid_argument("Steps and tolerances must be positive and finite");
        if (!std::isfinite(surface_solve_tolerance) || surface_solve_tolerance<0 ||
            min_step > step || seed_grid < 0 || seed_grid == 1 || iterations < 1 || max_trace_steps < 4 ||
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
    std::vector<int> fiber_seed_hits{0,0};
    // Assembled curves label each triangle with its zero-based range segment.
    std::vector<int> face_segments;
    int face_segment(int face) const { return face_segments.empty() ? 0 : face_segments[face]; }
    size_t box_contact_triangles = 0;
    int surface_seed_hits = 0;
    double seconds = 0;
};

struct Seeds {
    std::vector<Vec> surface;
    std::array<std::vector<Vec>,2> endpoints;
    size_t starting_points = 0;
    SeedShape grid_shape{0,0,0};
    double seconds = 0;
};

// The restart state includes the virtual fiber uses as well as the mesh.
// Ordinary exposed halfedges are reconstructed from oriented triangles.
struct Checkpoint {
    Result result;
    std::map<Edge,int> pending, fiber_uses;
};

Seeds find_seeds(const Field& field, const Segment& segment, const Options& options,
                 bool discover_last_endpoint = true);
Result trace_endpoint_fibers(const Field& field, const Segment& segment,
                             const Options& options, const Seeds& seeds);

inline size_t discovery_size(const SeedShape& shape, int limit=Options::default_max_seed_starts) {
    size_t size=1;
    for (int n:shape) {
        if (limit<1 || n<1 || size>size_t(limit)/size_t(n))
            throw std::invalid_argument("Seed grid exceeds its allocation limit or has an empty axis");
        size*=size_t(n);
    }
    return size;
}

inline SeedShape discovery_shape(const Field& field, const Options& options) {
    const Vec lengths=field.box.max-field.box.min;
    if (!lengths.allFinite() || (lengths.array()<=0).any()) throw std::invalid_argument("Invalid seed domain");
    SeedShape shape;
    const int shortest=options.seed_grid==0 ? 7 : options.seed_grid;
    if (shortest<2) throw std::invalid_argument("--seed-grid must be 0 (default count 7) or at least 2");
    for (int d=0; d<3; ++d) {
        const double n=std::ceil(shortest*(lengths[d]/lengths.minCoeff()));
        if (!std::isfinite(n) || n>options.max_seed_starts)
            throw std::invalid_argument("Seed axis exceeds the seed grid allocation limit");
        shape[d]=int(n);
    }
    discovery_size(shape,options.max_seed_starts);
    return shape;
}

inline std::vector<Vec> discovery_grid(const Box& box, const SeedShape& shape, int limit=Options::default_max_seed_starts) {
    std::vector<Vec> points;
    points.reserve(discovery_size(shape,limit));
    for (int k=0; k<shape[2]; ++k) for (int j=0; j<shape[1]; ++j) for (int i=0; i<shape[0]; ++i) {
        const SeedShape index{i,j,k}; Vec point;
        for (int d=0; d<3; ++d) {
            const double u=shape[d]==1 ? .5 : double(index[d])/(shape[d]-1);
            point[d]=index[d]==shape[d]-1 && shape[d]>1 ? box.max[d] : box.min[d]+(box.max[d]-box.min[d])*u;
        }
        points.push_back(point);
    }
    return points;
}
inline std::vector<Vec> discovery_grid(const Box& box, int n) {
    return discovery_grid(box,SeedShape{n,n,n});
}
inline std::vector<Vec> discovery_grid(const Field& field, const Options& options) {
    return discovery_grid(field.box,discovery_shape(field,options),options.max_seed_starts);
}
} // namespace fiber_surface
