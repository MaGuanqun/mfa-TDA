#include "trace.h"
#include <chrono>

namespace fiber_surface {
namespace {
constexpr char stalled[]="Fiber trace made insufficient progress";

template<class PointAt>
bool near_trace(const Vec& seed, size_t count, bool closed, double radius, PointAt point) {
    if (count<2) return false;
    const size_t edges=closed ? count : count-1;
    for (size_t i=0; i<edges; ++i) {
        const Vec closest=marching_triangles::closest_on_segment(seed,point(i),point((i+1)%count));
        if ((seed-closest).norm()<=radius) return true;
    }
    return false;
}
}

bool FiberTracer::contact(const Vec& p, const Vec& trial, const Range& q, double h, Vec& x) const {
    double fraction=1;
    const Vec delta=trial-p;
    for (int d=0; d<3; ++d) {
        if (delta[d]>0) fraction=std::min(fraction,(c_.field.box.max[d]-p[d])/delta[d]);
        if (delta[d]<0) fraction=std::min(fraction,(c_.field.box.min[d]-p[d])/delta[d]);
    }
    x=p+std::max(0.,fraction)*delta;
    const unsigned faces=c_.field.box.faces(x,1e-10*h);
    c_.field.box.fix(x,faces);
    if (!faces || !c_.fiber(x,q,Vec::Zero(),x,2*h,faces)) return false;
    return (x-p).dot(delta)>=-1e-10*h*delta.norm();
}

// A step either reaches a verified contact or returns a corrected, regular sample.
const char* FiberTracer::try_step(const Sample& start, const Range& q, double sign, double h, Sample& next) const {
    const Vec& p=start.point;
    next.on_box=false;
    auto reach_box=[&](const Vec& trial) -> const char* {
        next.on_box=contact(p,trial,q,h,next.point);
        return next.on_box ? nullptr : "Fiber-box contact correction failed";
    };
    std::array<Vec,4> k;
    k[0]=start.tangent;
    for (int i=1; i<4; ++i) {
        const Vec trial=p+(i==3 ? h : .5*h)*k[i-1];
        if (!c_.field.box.contains(trial)) return reach_box(trial);
        if (!c_.tangent(trial,k[i])) return "Fiber RK4 stage lost Jacobian rank";
        k[i]*=sign;
    }
    const Vec predictor=p+h*(k[0]+2*k[1]+2*k[2]+k[3])/6;
    if (!c_.field.box.contains(predictor)) return reach_box(predictor);
    next.point=predictor;
    if (!c_.fiber(next.point,q,start.tangent,predictor,h)) {
        if (!c_.field.box.contains(next.point)) return reach_box(next.point);
        return "Interior fiber Newton correction failed";
    }
    if (!c_.tangent(next.point,next.tangent)) return "Fiber correction reached rank loss";
    next.tangent*=sign;
    return nullptr;
}

const char* FiberTracer::advance(const Sample& start, const Range& q, double sign, double h, Sample& next) const {
    const char* failure=try_step(start,q,sign,h,next);
    if (!failure || h/2<c_.options.min_step) return failure;
    // Keep the nominal output spacing. Both halves must succeed before emitting
    // a sample, unless the first half reaches a verified domain boundary.
    Sample middle;
    if (const char* first=advance(start,q,sign,h/2,middle)) return first;
    if (middle.on_box) { next=middle; return nullptr; }
    return advance(middle,q,sign,h/2,next);
}

FiberTracer::Trace FiberTracer::direction(const Vec& seed, const Range& q, double sign, std::vector<Event>& events) const {
    Trace trace; trace.points.push_back(seed);
    const double h=c_.options.step;
    Sample current{seed,Vec::Zero()};
    if (!c_.tangent(seed,current.tangent)) { events.push_back({seed,"Endpoint fiber lost Jacobian rank"}); return trace; }
    current.tangent*=sign;
    for (int step=0; step<c_.options.max_trace_steps; ++step) {
        Sample next{Vec::Zero(),Vec::Zero()};
        if (const char* reason=advance(current,q,sign,h,next)) { events.push_back({current.point,reason}); return trace; }
        if (next.on_box) {
            if ((next.point-current.point).norm()>h*1e-7) trace.points.push_back(next.point);
            trace.complete=true; return trace;
        }
        if ((next.point-current.point).norm()<h*.1) { events.push_back({current.point,stalled}); return trace; }
        if (trace.points.size()>=4 && (next.point-seed).norm()<=h) {
            trace.closed=true; trace.complete=true; return trace;
        }
        trace.points.push_back(next.point);
        current=next;
    }
    events.push_back({trace.points.back(),"Fiber trace step limit reached"}); return trace;
}

bool FiberTracer::covered(const Vec& seed, int endpoint) const {
    for (const auto& fiber:result_.fibers) if (fiber.endpoint==endpoint) {
        if (near_trace(seed,fiber.vertices.size(),fiber.closed,c_.options.seed_coverage_radius(),
                       [&](size_t i)->const Vec& { return result_.mesh.vertices[fiber.vertices[i]]; })) return true;
    }
    return false;
}

FiberTracer::Fragment FiberTracer::discover(const Vec& seed, const Range& q, size_t index) const {
    Fragment f{index,{}, {},std::nullopt};
    f.trace=direction(seed,q,1,f.events);
    if (f.trace.closed) return f;
    const bool stalled_forward=!f.trace.complete;
    Trace backward=direction(seed,q,-1,f.events);
    std::reverse(backward.points.begin(),backward.points.end());
    backward.points.pop_back();
    backward.points.insert(backward.points.end(),f.trace.points.begin(),f.trace.points.end());
    f.trace.points=std::move(backward.points);
    f.trace.complete=f.trace.complete && backward.complete;
    if (!f.events.empty() && std::all_of(f.events.begin(),f.events.end(),[](const Event& e) { return e.reason==stalled; }))
        f.evidence=witnesses(f.trace,seed,q,!backward.complete,stalled_forward);
    return f;
}

void FiberTracer::run_seeds(const std::array<std::vector<Vec>,2>& endpoints) {
    for (int endpoint=0; endpoint<2; ++endpoint) {
        const Range q=endpoint ? c_.segment.b : c_.segment.a;
        const auto& seeds=endpoints[endpoint];
        std::vector<Fragment> fragments;
        const double radius=c_.options.seed_coverage_radius();
        for (size_t i=0; i<seeds.size(); ++i) {
            const Vec& seed=seeds[i];
            ++result_.fiber_seed_hits[endpoint];
            if (covered(seed,endpoint)) continue;
            const bool on_fragment=std::any_of(fragments.begin(),fragments.end(),[&](const Fragment& f) {
                const auto& points=f.trace.points;
                return near_trace(seed,points.size(),f.trace.closed,radius,[&](size_t k)->const Vec& { return points[k]; });
            });
            if (!on_fragment) fragments.push_back(discover(seed,q,i));
        }
        // Freeze only after normal discovery. A repaired curve must not hide
        // an independent, as-yet unverified fragment under the seed coverage radius.
        for (auto& f:fragments) {
            if (f.evidence) f.recovery=recover(f.trace,q,seeds,f.seed_index,*f.evidence);
        }
        std::vector<const Fragment*> repairs;
        for (const auto& f:fragments) {
            if (f.recovery==Recovery::complete && std::none_of(repairs.begin(),repairs.end(),[&](const Fragment* other) {
                    return same_component(f.trace,other->trace,q);
                })) repairs.push_back(&f);
        }
        for (auto& f:fragments) {
            auto& forward=f.trace;
            const Vec& seed=seeds[f.seed_index];
            if (f.recovery==Recovery::complete) {
                if (std::find(repairs.begin(),repairs.end(),&f)==repairs.end()) continue;
            } else {
                if (forward.complete && std::any_of(repairs.begin(),repairs.end(),[&](const Fragment* other) {
                        return same_component(other->trace,forward,q);
                    })) continue;
                // A different seed may have supplied the unique verified continuation.
                if (f.recovery!=Recovery::ambiguous && f.evidence &&
                    std::count_if(repairs.begin(),repairs.end(),[&](const Fragment* other) {
                        return contains_fragment(other->trace,q,*f.evidence);
                    })==1) continue;
                result_.unresolved.insert(result_.unresolved.end(),f.events.begin(),f.events.end());
            }
            if (forward.points.size()<2) {
                result_.unresolved.push_back({seed,"Fiber component has fewer than two samples"});
                continue;
            }
            Fiber fiber; fiber.endpoint=endpoint; fiber.closed=forward.closed; fiber.complete=forward.complete;
            const int component=int(result_.fibers.size());
            for (const Vec& p:forward.points) {
                fiber.vertices.push_back(int(result_.mesh.vertices.size()));
                result_.mesh.vertices.push_back(p);
                result_.info.push_back({endpoint,component,c_.field.box.faces(p,c_.options.step*1e-8)});
            }
            const size_t n=fiber.vertices.size();
            for (size_t k=0; k+(fiber.closed?0:1)<n; ++k)
                result_.fiber_edges[edge_key(fiber.vertices[k],fiber.vertices[(k+1)%n])]=component;
            result_.fibers.push_back(std::move(fiber));
        }
    }
}

void FiberTracer::run(const std::vector<Vec>& grid) {
    std::array<std::vector<Vec>,2> endpoints;
    for (int endpoint=0; endpoint<2; ++endpoint)
        for (Vec seed:grid)
            if (c_.fiber_seed(seed,endpoint ? c_.segment.b : c_.segment.a)) endpoints[endpoint].push_back(seed);
    run_seeds(endpoints);
}

Result trace_endpoint_fibers(const Field& field, const Segment& segment,
                             const Options& options, const Seeds& seeds) {
    const auto start=std::chrono::steady_clock::now();
    Result result;
    Corrector corrector(field,segment,options);
    FiberTracer(corrector,result).run_seeds(seeds.endpoints);
    result.seconds=seeds.seconds+std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
} // namespace fiber_surface
