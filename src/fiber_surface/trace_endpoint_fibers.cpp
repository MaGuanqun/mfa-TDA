#include "trace.h"
#include <chrono>

namespace fiber_surface {
void FiberTracer::event(const Vec& x, const std::string& why) { result_.unresolved.push_back({x,why}); }

bool FiberTracer::contact(const Vec& p, const Vec& trial, const Range& q, Vec& x) const {
    double fraction=1;
    const Vec delta=trial-p;
    for (int d=0; d<3; ++d) {
        if (delta[d]>0) fraction=std::min(fraction,(c_.field.box.max[d]-p[d])/delta[d]);
        if (delta[d]<0) fraction=std::min(fraction,(c_.field.box.min[d]-p[d])/delta[d]);
    }
    x=p+std::max(0.,fraction)*delta;
    const unsigned faces=c_.field.box.faces(x,1e-10*c_.options.step);
    c_.field.box.fix(x,faces);
    if (!faces || !c_.fiber(x,q,Vec::Zero(),x,2*c_.options.fiber_step,faces)) return false;
    return (x-p).dot(delta)>=-1e-10*c_.options.step*delta.norm();
}

FiberTracer::Trace FiberTracer::direction(const Vec& seed, const Range& q, double sign) {
    Trace trace; trace.points.push_back(seed);
    const double h=c_.options.fiber_step;
    Vec start_t;
    if (!c_.tangent(seed,start_t)) { event(seed,"Endpoint fiber lost Jacobian rank"); return trace; }
    start_t*=sign;
    for (int step=0; step<c_.options.max_trace_steps; ++step) {
        const Vec p=trace.points.back(); Vec t;
        if (!c_.tangent(p,t)) { event(p,"Endpoint fiber lost Jacobian rank"); return trace; }
        t*=sign;
        Vec v2,v3,v4,predictor;
        bool hit=false, failed=false; Vec contact_point;
        auto stage=[&](const Vec& x, Vec& v) {
            if (hit || failed) return;
            if (!c_.field.box.contains(x)) {
                hit=contact(p,x,q,contact_point);
                failed=!hit; return;
            }
            if (!c_.tangent(x,v)) { failed=true; return; }
            v*=sign;
        };
        stage(Vec(p+.5*h*t),v2);
        if (!hit && !failed) stage(Vec(p+.5*h*v2),v3);
        if (!hit && !failed) stage(Vec(p+h*v3),v4);
        if (!hit && !failed) {
            predictor=p+h*(t+2*v2+2*v3+v4)/6;
            if (!c_.field.box.contains(predictor)) {
                hit=contact(p,predictor,q,contact_point); failed=!hit;
            } else {
                Vec corrected=predictor;
                if (!c_.fiber(corrected,q,t,predictor,h)) {
                    if (!c_.field.box.contains(corrected)) hit=contact(p,corrected,q,contact_point);
                    failed=!hit;
                } else predictor=corrected;
            }
        }
        if (failed) { event(p,"RK4 stage, fiber correction, or fiber-box contact failed"); return trace; }
        if (hit) {
            if ((contact_point-p).norm()>h*1e-7) trace.points.push_back(contact_point);
            trace.complete=true; return trace;
        }
        Vec next_t;
        if (!c_.tangent(predictor,next_t)) { event(p,"Fiber correction reached rank loss"); return trace; }
        if (trace.points.size()>=4 && (predictor-seed).norm()<=h && sign*next_t.dot(start_t)>.8) {
            trace.closed=true; trace.complete=true; return trace;
        }
        if ((predictor-p).norm()<h*.1) { event(p,"Fiber trace made insufficient progress"); return trace; }
        trace.points.push_back(predictor);
    }
    event(trace.points.back(),"Fiber trace step limit reached"); return trace;
}

bool FiberTracer::covered(const Vec& seed, int endpoint) const {
    Vec t; if (!c_.tangent(seed,t)) return false;
    for (const auto& fiber:result_.fibers) if (fiber.endpoint==endpoint) {
        const size_t n=fiber.vertices.size();
        for (size_t k=0; k+(fiber.closed?0:1)<n; ++k) {
            const Vec& a=result_.mesh.vertices[fiber.vertices[k]];
            const Vec& b=result_.mesh.vertices[fiber.vertices[(k+1)%n]];
            Vec near=marching_triangles::closest_on_segment(seed,a,b);
            if ((near-seed).norm()>c_.options.fiber_step*.6) continue;
            const Range q=endpoint ? c_.segment.b : c_.segment.a;
            if (!c_.fiber(near,q,t,seed,c_.options.fiber_step*.6)) continue;
            Vec nt;
            if ((near-seed).norm()<c_.options.fiber_step*1e-5 && c_.tangent(near,nt) && std::abs(t.dot(nt))>.8) return true;
        }
    }
    return false;
}

void FiberTracer::run_seeds(const std::array<std::vector<Vec>,2>& endpoints) {
    for (int endpoint=0; endpoint<2; ++endpoint) {
        const Range q=endpoint ? c_.segment.b : c_.segment.a;
        for (const Vec& seed:endpoints[endpoint]) {
            ++result_.fiber_seed_hits[endpoint];
            if (covered(seed,endpoint)) continue;
            Trace forward=direction(seed,q,1);
            if (!forward.closed) {
                Trace backward=direction(seed,q,-1);
                std::reverse(backward.points.begin(),backward.points.end());
                backward.points.pop_back();
                backward.points.insert(backward.points.end(),forward.points.begin(),forward.points.end());
                forward.points=std::move(backward.points);
                forward.complete=forward.complete && backward.complete;
            }
            if (forward.points.size()<2) { event(seed,"Fiber component has fewer than two samples"); continue; }
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
