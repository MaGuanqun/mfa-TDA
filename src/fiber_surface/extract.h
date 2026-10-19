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
    bool covered(const Vec& seed) const;
    bool seed_triangle(const Vec& p);
public:
    Extractor(Field field, Segment segment, Options options={});
    Extractor(const Extractor&)=delete;
    Extractor& operator=(const Extractor&)=delete;
    Checkpoint march(const Seeds& seeds, Result fibers);
    Checkpoint finish(Checkpoint state);
    Result run();
};
} // namespace fiber_surface
