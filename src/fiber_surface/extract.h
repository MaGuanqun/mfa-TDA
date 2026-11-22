#pragma once
#include "trace.h"

namespace fiber_surface {
class Extractor {
    Field field_;
    Segment segment_;
    Options options_;
    Corrector corrector_;
    Result result_;
    std::map<Edge,int> pending_; // oriented virtual front uses of unfilled fiber edges
    std::map<Edge,int> fiber_uses_;

    bool intended(const Edge& e) const;
    bool active(const Edge& e) const;
    std::vector<Edge> fronts() const;
    Vec gradient(const Vec& x) const;
    double epsilon() const;
    bool geometry(const Triangle& t, bool delaunay, bool overlap=true) const;
    bool insert(const Triangle& t, bool delaunay, bool overlap=true);
    int vertex(const Vec& x, VertexInfo info={});
    void rollback_vertices(size_t n);
    bool new_triangle(const Edge& e, const Vec& x);
    bool clip(const Edge& e, const Vec& outside);
    bool existing(const Edge& e, const Vec& target, bool boundary_only=false);
    bool advance(const Edge& e);
    void grow();
    void stitch();
    void repair();
    void snap_box_slivers();
    bool covered(const Vec& seed) const;
    bool seed_triangle(const Vec& p);
public:
    Extractor(Field field, Segment segment, Options options={});
    Extractor(const Extractor&)=delete;
    Extractor& operator=(const Extractor&)=delete;
    Checkpoint march(const Seeds& seeds, Result fibers);
    Checkpoint finish(Checkpoint state, bool curve_patch=false);
    Result run();
};
// Use a tighter surface root solve where neighboring patches meet.
inline Options curve_options(const RangeCurve& curve, Options options) {
    if (curve.size()>1) {
        const double tolerance=options.surface_solve_tolerance>0 ? options.surface_solve_tolerance : options.residual_tolerance;
        options.surface_solve_tolerance=std::min(1e-10,tolerance);
    }
    return options;
}
std::vector<Seeds> find_curve_seeds(const Field&, const RangeCurve&, const Options&);
std::vector<Result> trace_curve_fibers(const Field&, const RangeCurve&, const Options&, const std::vector<Seeds>&);
Result join_curve_patches(const RangeCurve&, const std::vector<Result>&);
Result extract_curve(const Field&, const RangeCurve&, const Options&);
} // namespace fiber_surface
