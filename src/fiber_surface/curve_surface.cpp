#include "extract.h"

namespace fiber_surface {
namespace {
int turn(const Range& a, const Range& b, const Range& c) {
    const auto u=b.cast<long double>()-a.cast<long double>();
    const auto v=c.cast<long double>()-a.cast<long double>();
    const long double cross=u.x()*v.y()-u.y()*v.x();
    return (cross>0)-(cross<0);
}
bool between(const Range& a, const Range& p, const Range& b) {
    return (p.array()>=a.cwiseMin(b).array()).all() && (p.array()<=a.cwiseMax(b).array()).all();
}
bool intersect(const Segment& s, const Segment& t) {
    const int a=turn(s.a,s.b,t.a), b=turn(s.a,s.b,t.b);
    const int c=turn(t.a,t.b,s.a), d=turn(t.a,t.b,s.b);
    return (a*b<0 && c*d<0) || (!a && between(s.a,t.a,s.b)) || (!b && between(s.a,t.b,s.b)) ||
           (!c && between(t.a,s.a,t.b)) || (!d && between(t.a,s.b,t.b));
}

// Copy a frozen component in trace order, relabeling its corner and vertex IDs.
int append_fiber(Result& target, const Result& source, Fiber fiber, int corner) {
    const int component=static_cast<int>(target.fibers.size());
    fiber.endpoint=corner;
    for (int& v:fiber.vertices) {
        target.mesh.vertices.push_back(source.mesh.vertices[v]);
        target.info.push_back({corner,component,source.info[v].faces});
        v=static_cast<int>(target.info.size())-1;
    }
    const auto& vertices=fiber.vertices;
    for (size_t j=1; j<vertices.size(); ++j)
        target.fiber_edges[edge_key(vertices[j-1],vertices[j])]=component;
    if (fiber.closed) target.fiber_edges[edge_key(vertices.back(),vertices.front())]=component;
    target.fibers.push_back(std::move(fiber));
    return component;
}
} // namespace

RangeCurve::RangeCurve(std::vector<Range> points) {
    for (const auto& p:points)
        if (!p.allFinite()) throw std::invalid_argument("Range curve vertices must be finite");
    points.erase(std::unique(points.begin(),points.end()),points.end());
    const bool loop=points.size()>1 && points.back()==points.front();
    if (loop) points.pop_back();
    // Removing a collinear point also exposes its neighbors to the same check.
    for (size_t i=loop?0:1; points.size()>=3 && i<points.size()-(loop?0:1);) {
        const Range& a=points[(i+points.size()-1)%points.size()];
        const Range& b=points[i];
        const Range& c=points[(i+1)%points.size()];
        if (turn(a,b,c)) { ++i; continue; }
        if (!between(a,b,c)) throw std::invalid_argument("Range curve retraces a segment");
        points.erase(points.begin()+i);
        i=loop?0:std::max<size_t>(1,i-1);
    }
    if (points.size()<(loop?3u:2u)) throw std::invalid_argument("Range curve has too few distinct vertices");
    for (size_t i=1; i<points.size(); ++i) segments.emplace_back(points[i-1],points[i]);
    if (loop) segments.emplace_back(points.back(),points.front());
    for (int i=0; i<size(); ++i) for (int j=i+2; j<size(); ++j) {
        if (loop && i==0 && j==size()-1) continue;
        if (intersect(segments[i],segments[j])) throw std::invalid_argument("Range curve must be simple");
    }
}

std::vector<Seeds> find_curve_seeds(const Field& f, const RangeCurve& curve, const Options& o) {
    std::vector<Seeds> seeds;
    for (int i=0; i<curve.size(); ++i)
        seeds.push_back(find_seeds(f,curve.segments[i],o,!curve.closed() && i+1==curve.size()));
    return seeds;
}

std::vector<Result> trace_curve_fibers(const Field& f, const RangeCurve& curve, const Options& o,
                                       const std::vector<Seeds>& seeds) {
    std::vector<Result> patches(curve.size());
    if (curve.size()==1) {
        patches[0]=trace_endpoint_fibers(f,curve.segments[0],o,seeds[0]);
        return patches;
    }
    std::vector<Result> corners;
    for (int k=0; k<curve.vertex_count(); ++k) {
        const int side=std::min(k,curve.size()-1), endpoint=k-side;
        Seeds candidates;
        candidates.endpoints[endpoint]=seeds[side].endpoints[endpoint];
        auto corner=trace_endpoint_fibers(f,curve.segments[side],o,candidates);
        auto& owner=patches[side];
        owner.seconds+=corner.seconds;
        owner.unresolved.insert(owner.unresolved.end(),corner.unresolved.begin(),corner.unresolved.end());
        corners.push_back(std::move(corner));
    }
    for (int i=0; i<curve.size(); ++i) {
        auto& patch=patches[i];
        patch.seconds+=seeds[i].seconds;
        for (int endpoint=0; endpoint<2; ++endpoint) {
            const auto& corner=corners[curve.corner(i,endpoint)];
            for (const Fiber& fiber:corner.fibers) append_fiber(patch,corner,fiber,endpoint);
            patch.fiber_seed_hits[endpoint]=corner.fiber_seed_hits[0]+corner.fiber_seed_hits[1];
        }
    }
    return patches;
}

Result join_curve_patches(const RangeCurve& curve, const std::vector<Result>& patches) {
    Result surface;
    surface.fiber_seed_hits.resize(curve.vertex_count());
    // Stable identity is (range vertex, traced component, sample index).
    std::vector<std::vector<int>> corner_fibers(curve.vertex_count());
    std::vector<bool> seen(curve.vertex_count(),false);
    for (int i=0; i<curve.size(); ++i) {
        const auto& patch=patches.at(i);
        std::vector<int> ids(patch.info.size(),-1);
        std::array<size_t,2> counts{};
        for (const auto& fiber:patch.fibers) {
            const int k=curve.corner(i,fiber.endpoint);
            auto& components=corner_fibers[k];
            const size_t component=counts[fiber.endpoint]++;
            if (!seen[k]) components.push_back(append_fiber(surface,patch,fiber,k));
            const Fiber& shared=surface.fibers.at(components.at(component));
            if (shared.vertices.size()!=fiber.vertices.size() || shared.closed!=fiber.closed)
                throw std::runtime_error("Inconsistent shared corner component");
            for (size_t j=0; j<fiber.vertices.size(); ++j) {
                const int v=fiber.vertices[j];
                ids[v]=shared.vertices[j];
                if (patch.mesh.vertices[v]!=surface.mesh.vertices[ids[v]])
                    throw std::runtime_error("Inconsistent shared corner sample");
            }
        }
        for (int e=0; e<2; ++e) {
            const int k=curve.corner(i,e);
            if (counts[e]!=corner_fibers[k].size()) throw std::runtime_error("Missing shared corner component");
            surface.fiber_seed_hits[k]=patch.fiber_seed_hits[e];
            seen[k]=true;
        }
        for (size_t j=0; j<ids.size(); ++j) if (ids[j]<0) {
            ids[j]=static_cast<int>(surface.info.size());
            surface.info.push_back(patch.info[j]);
            surface.mesh.vertices.push_back(patch.mesh.vertices[j]);
        }
        for (auto t:patch.mesh.triangles) {
            for (int& v:t) v=ids[v];
            if (!surface.mesh.topology_ok(t)) {
                surface.unresolved.push_back({surface.mesh.vertices[t[0]],"Inconsistent triangle incidence at a range corner"});
                continue;
            }
            surface.mesh.add(t);
            surface.face_segments.push_back(i);
        }
        for (auto e:patch.box_edges) surface.box_edges.insert(edge_key(ids[e.first],ids[e.second]));
        surface.unresolved.insert(surface.unresolved.end(),patch.unresolved.begin(),patch.unresolved.end());
        surface.surface_seed_hits+=patch.surface_seed_hits;
        surface.seconds+=patch.seconds;
        surface.box_contact_triangles+=patch.box_contact_triangles;
    }
    return surface;
}

Result extract_curve(const Field& field, const RangeCurve& curve, const Options& requested) {
    const auto options=curve_options(curve,requested);
    if (curve.size()==1) return Extractor(field,curve.segments.front(),options).run();
    auto seeds=find_curve_seeds(field,curve,options);
    auto patches=trace_curve_fibers(field,curve,options,seeds);
    for (int i=0; i<curve.size(); ++i) {
        Extractor extractor(field,curve.segments[i],options);
        auto state=extractor.march(seeds[i],std::move(patches[i]));
        patches[i]=extractor.finish(std::move(state),true).result;
    }
    return join_curve_patches(curve,patches);
}
} // namespace fiber_surface
