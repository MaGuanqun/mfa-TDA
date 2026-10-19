#include "extract.h"
#include "../isosurface/crack_closing.h"
#include <chrono>
#include <deque>

namespace fiber_surface {
Extractor::Extractor(Field field, Segment segment, Options options)
    : field_(std::move(field)), segment_(std::move(segment)), options_(options), corrector_(field_,segment_,options_) {
    options_.validate(field_,segment_);
}

Checkpoint Extractor::march(const Seeds& seeds, Result fibers) {
    const auto start=std::chrono::steady_clock::now();
    result_=std::move(fibers); pending_.clear(); fiber_uses_.clear();
    for (const auto& entry:result_.fiber_edges) {
        if (!result_.fibers[entry.second].complete) continue;
        int a=entry.first.first,b=entry.first.second;
        Vec mid=(result_.mesh.vertices[a]+result_.mesh.vertices[b])/2,t;
        const int endpoint=result_.fibers[entry.second].endpoint;
        if (!corrector_.tangent(result_.mesh.vertices[a],t) ||
            !corrector_.fiber(mid,endpoint ? segment_.b : segment_.a,t,mid,options_.step)) continue;
        const Jacobian j=field_.jacobian(mid);
        const Range inward=(endpoint ? -1. : 1.)*segment_.d/segment_.length;
        const Vec lift=j.transpose()*(j*j.transpose()).ldlt().solve(inward);
        if ((result_.mesh.vertices[b]-result_.mesh.vertices[a]).cross(gradient(mid)).dot(lift)<0) std::swap(a,b);
        pending_[{a,b}]=entry.second;
    }
    fiber_uses_=pending_;
    grow(); stitch();
    for (const Vec& seed:seeds.surface) {
        if (result_.mesh.triangles.size()>=options_.max_triangles) break;
        ++result_.surface_seed_hits;
        if (covered(seed)) continue;
        if (seed_triangle(seed)) { grow(); stitch(); }
    }
    result_.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return {std::move(result_),std::move(pending_),std::move(fiber_uses_)};
}

Result Extractor::run() {
    const auto seeds=find_seeds(field_,segment_,options_);
    auto fibers=trace_endpoint_fibers(field_,segment_,options_,seeds);
    auto state=march(seeds,std::move(fibers));
    return finish(std::move(state)).result;
}

bool Extractor::intended(const Edge& e) const {
    return result_.fiber_edges.count(edge_key(e.first,e.second)) || result_.box_edges.count(edge_key(e.first,e.second));
}

bool Extractor::active(const Edge& e) const {
    return pending_.count(e) || (result_.mesh.boundary.count(e) && !intended(e));
}

std::vector<Edge> Extractor::fronts() const {
    std::vector<Edge> out;
    for (const auto& e:pending_) out.push_back(e.first);
    for (const auto& e:result_.mesh.boundary) if (!intended(e.first)) out.push_back(e.first);
    return out;
}

Vec Extractor::gradient(const Vec& x) const { return segment_.gradient(field_.jacobian(x)); }

double Extractor::epsilon() const { return options_.step*1e-9; }

bool Extractor::geometry(const Triangle& t, bool delaunay, bool overlap) const {
    const auto& mesh=result_.mesh;
    if (!mesh.topology_ok(t)) return false;
    const Vec &a=mesh.vertices[t[0]], &b=mesh.vertices[t[1]], &c=mesh.vertices[t[2]];
    const double longest=std::max({(b-a).norm(),(c-a).norm(),(c-b).norm()});
    const double shortest=std::min({(b-a).norm(),(c-a).norm(),(c-b).norm()});
    if (longest>2*options_.step || shortest<=epsilon()) return false;
    const Vec cross=(b-a).cross(c-a);
    if (cross.norm()<1e-8*longest*longest) return false;
    const Vec normal=cross.normalized();
    for (int k=0; k<3; ++k) {
        const int u=t[k],v=t[(k+1)%3]; const Edge key=edge_key(u,v);
        if (intended(key) && mesh.edges.count(key)) return false;
        auto fiber=result_.fiber_edges.find(key);
        if (fiber!=result_.fiber_edges.end() && !pending_.count({v,u})) return false;
        // Traced boundary edges may only be consumed in their prescribed order.
        if (result_.info[u].fiber>=0 && result_.info[u].fiber==result_.info[v].fiber && fiber==result_.fiber_edges.end()) return false;
        const Vec e1=mesh.vertices[t[(k+1)%3]]-mesh.vertices[t[k]], e2=mesh.vertices[t[(k+2)%3]]-mesh.vertices[t[k]];
        const double angle=std::acos(std::clamp(e1.normalized().dot(e2.normalized()),-1.,1.))*180/std::acos(-1.);
        if (angle<options_.min_angle_degrees) return false;
    }
    const double min_cos=std::cos(options_.max_normal_degrees*std::acos(-1.)/180);
    for (const Vec& x:std::array<Vec,3>{a,b,c}) {
        if (!corrector_.member(x) || std::abs(segment_.phi(field_.value(x)))>options_.residual_tolerance*1.01) return false;
        const Vec g=gradient(x);
        if (!g.allFinite() || g.norm()<=options_.derivative_tolerance || normal.dot(g.normalized())<min_cos) return false;
    }
    // Check chords as well as corrected vertices: nonlinear segment bounds
    // can fail inside a triangle whose three corners satisfy the constraints.
    const std::array<Vec,4> samples{(a+b)/2,(b+c)/2,(c+a)/2,(a+b+c)/3};
    for (int sample_id=0; sample_id<4; ++sample_id) {
        Vec x=samples[sample_id];
        if (!corrector_.member(x)) return false;
        const Vec sample=x;
        bool corrected=false;
        if (sample_id<3 && result_.fiber_edges.count(edge_key(t[sample_id],t[(sample_id+1)%3]))) {
            const auto& info=result_.info[t[sample_id]]; Vec tangent;
            corrected=corrector_.tangent(mesh.vertices[t[sample_id]],tangent) &&
                corrector_.fiber(x,info.endpoint ? segment_.b : segment_.a,tangent,sample,options_.spatial_tolerance);
        } else {
            const unsigned faces=sample_id<3 ? result_.info[t[sample_id]].faces & result_.info[t[(sample_id+1)%3]].faces : 0;
            corrected=corrector_.surface(x,options_.spatial_tolerance,faces);
        }
        if (!corrected || (x-sample).norm()>options_.spatial_tolerance) return false;
        const Vec g=gradient(x);
        if (normal.dot(g.normalized())<min_cos) return false;
    }
    // Intended curved boundary chords have their own spatial error check.
    for (int k=0; k<3; ++k) {
        const int u=t[k],v=t[(k+1)%3]; const auto& iu=result_.info[u]; const auto& iv=result_.info[v];
        Vec x=(mesh.vertices[u]+mesh.vertices[v])/2, original=x;
        if (result_.fiber_edges.count(edge_key(u,v))) {
            Vec tangent; if (!corrector_.tangent(mesh.vertices[u],tangent)) return false;
            if (!corrector_.fiber(x,iu.endpoint ? segment_.b : segment_.a,tangent,original,options_.spatial_tolerance)) return false;
        } else if (iu.faces & iv.faces) {
            if (!corrector_.surface(x,options_.spatial_tolerance,iu.faces & iv.faces)) return false;
        }
    }
    Vec center; double radius2=0;
    if (delaunay && !marching_triangles::circumcircle(a,b,c,center,radius2)) return false;
    for (const auto& old:mesh.triangles) {
        if (marching_triangles::triangles_conflict(mesh.vertices,t,old,epsilon()) ||
            marching_triangles::triangles_conflict(mesh.vertices,old,t,epsilon())) return false;
        if (overlap && marching_triangles::surface_overlap(mesh.vertices,t,old,epsilon(),normal)) return false;
        if (!delaunay) continue;
        const Vec old_normal=(mesh.vertices[old[1]]-mesh.vertices[old[0]]).cross(mesh.vertices[old[2]]-mesh.vertices[old[0]]);
        if (old_normal.dot(normal)<=0) continue;
        bool shared=false;
        for (int i:old) {
            if (std::find(t.begin(),t.end(),i)!=t.end()) shared=true;
            else if ((mesh.vertices[i]-center).squaredNorm()<radius2*(1-1e-9)) return false;
        }
        if (!shared &&
            (marching_triangles::closest_on_triangle(center,mesh.vertices[old[0]],mesh.vertices[old[1]],mesh.vertices[old[2]])-center).squaredNorm()<radius2*(1-1e-9)) return false;
    }
    return true;
}

bool Extractor::insert(const Triangle& t, bool delaunay, bool overlap) {
    if (result_.mesh.triangles.size()>=options_.max_triangles || !geometry(t,delaunay,overlap)) return false;
    result_.mesh.add(t);
    for (int k=0; k<3; ++k) {
        const int a=t[k],b=t[(k+1)%3]; pending_.erase({b,a});
        if ((result_.info[a].faces & result_.info[b].faces) && !result_.fiber_edges.count(edge_key(a,b)))
            result_.box_edges.insert(edge_key(a,b));
    }
    return true;
}

int Extractor::vertex(const Vec& x, VertexInfo info) {
    const int id=int(result_.mesh.vertices.size()); result_.mesh.vertices.push_back(x); result_.info.push_back(info); return id;
}

void Extractor::rollback_vertices(size_t n) { result_.mesh.vertices.resize(n); result_.info.resize(n); }

bool Extractor::new_triangle(const Edge& e, const Vec& x) {
    // Join a nearby compatible front instead of creating almost coincident
    // vertices and leaving a strip too thin to triangulate later.
    bool too_close=false;
    std::set<int> nearby;
    for (const auto& f:fronts()) for (int id:{f.first,f.second})
        if (id!=e.first && id!=e.second && (result_.mesh.vertices[id]-x).norm()<options_.step*.45) nearby.insert(id);
    std::vector<std::pair<double,int>> ordered;
    for (int id:nearby) ordered.emplace_back((result_.mesh.vertices[id]-x).squaredNorm(),id);
    std::sort(ordered.begin(),ordered.end());
    for (const auto& candidate:ordered) {
        if (insert({e.second,e.first,candidate.second},false)) return true;
        if (candidate.first<.04*options_.step*options_.step) too_close=true;
    }
    if (too_close) return false;
    for (const auto& old:result_.mesh.vertices) if ((old-x).norm()<epsilon()) return false;
    const size_t n=result_.mesh.vertices.size();
    const int id=vertex(x,{-1,-1,field_.box.faces(x,epsilon())});
    if (insert({e.second,e.first,id},true)) return true;
    rollback_vertices(n); return false;
}

bool Extractor::clip(const Edge& e, const Vec& outside) {
    // Clip in oriented triangle order, then correct all new boundary points
    // simultaneously with their active box coordinates held fixed.
    struct ClipPoint { Vec x; int id=-1; unsigned faces=0; };
    std::vector<ClipPoint> polygon{{result_.mesh.vertices[e.second],e.second,result_.info[e.second].faces},
        {result_.mesh.vertices[e.first],e.first,result_.info[e.first].faces},{outside,-1,0}};
    for (int face=0; face<6; ++face) {
        const int d=face/2; const double limit=face%2 ? field_.box.max[d] : field_.box.min[d];
        auto inside=[&](const Vec& x) { return face%2 ? x[d]<=limit : x[d]>=limit; };
        std::vector<ClipPoint> next;
        for (size_t i=0; i<polygon.size(); ++i) {
            const auto& a=polygon[i]; const auto& b=polygon[(i+1)%polygon.size()];
            if (inside(a.x)) next.push_back(a);
            if (inside(a.x)!=inside(b.x)) {
                const double u=(limit-a.x[d])/(b.x[d]-a.x[d]);
                Vec x=a.x+u*(b.x-a.x); x[d]=limit;
                next.push_back({x,-1,(a.faces & b.faces) | (1u<<face)});
            }
        }
        polygon=std::move(next);
    }
    if (polygon.size()<3) return false;
    const size_t n=result_.mesh.vertices.size();
    std::vector<int> ids;
    for (auto p:polygon) {
        if (p.id<0) {
            p.faces |= field_.box.faces(p.x,epsilon()); field_.box.fix(p.x,p.faces);
            // Combined range/box events reuse the pretraced contact only.
            int existing=-1;
            for (size_t i=0; i<n; ++i) {
                const auto& info=result_.info[i];
                if ((info.faces & p.faces)!=p.faces || (result_.mesh.vertices[i]-p.x).norm()>options_.step*.4) continue;
                Vec test=p.x;
                if (info.endpoint>=0) {
                    const Vec delta=result_.mesh.vertices[i]-p.x;
                    if (delta.norm()>options_.step*.25) continue;
                    existing=int(i); break;
                }
                if (corrector_.surface(test,options_.step,p.faces) && (test-result_.mesh.vertices[i]).norm()<epsilon()) { existing=int(i); break; }
            }
            if (existing>=0) p.id=existing;
            else {
                if (!corrector_.surface(p.x,options_.step,p.faces)) { rollback_vertices(n); return false; }
                const double s=segment_.parameter(field_.value(p.x));
                if (s<=corrector_.slack() || s>=1-corrector_.slack()) { rollback_vertices(n); return false; }
                p.id=vertex(p.x,{-1,-1,p.faces});
            }
        }
        if (ids.empty() || ids.back()!=p.id) ids.push_back(p.id);
    }
    if (ids.size()>1 && ids.front()==ids.back()) ids.pop_back();
    if (ids.size()<3) { rollback_vertices(n); return false; }
    // Atomic polygon commit: rejected triangulations cannot leave half a clip.
    const Mesh original=result_.mesh; const auto pending=pending_; const auto box=result_.box_edges;
    for (size_t i=1; i+1<ids.size(); ++i) {
        if (!insert({ids[0],ids[i],ids[i+1]},false)) {
            result_.mesh=original; pending_=pending; result_.box_edges=box; rollback_vertices(n); return false;
        }
    }
    return true;
}

bool Extractor::existing(const Edge& e, const Vec& target, bool boundary_only) {
    std::set<int> candidates;
    for (const auto& f:fronts()) { candidates.insert(f.first); candidates.insert(f.second); }
    for (const auto& f:result_.fiber_edges) { candidates.insert(f.first.first); candidates.insert(f.first.second); }
    std::vector<std::pair<double,int>> ordered;
    for (int id:candidates) if (id!=e.first && id!=e.second && (!boundary_only || result_.info[id].endpoint>=0))
        ordered.emplace_back((result_.mesh.vertices[id]-target).squaredNorm(),id);
    std::sort(ordered.begin(),ordered.end());
    for (const auto& pair:ordered) {
        if (pair.first>4*options_.step*options_.step) break;
        if (insert({e.second,e.first,pair.second},false)) return true;
    }
    return false;
}

bool Extractor::advance(const Edge& e) {
    const Vec a=result_.mesh.vertices[e.first],b=result_.mesh.vertices[e.second], edge=b-a;
    Vec mid=(a+b)/2;
    auto pending=pending_.find(e);
    if (pending!=pending_.end()) {
        Vec tangent; if (!corrector_.tangent(a,tangent)) return false;
        const auto& fiber=result_.fibers[pending->second];
        if (!corrector_.fiber(mid,fiber.endpoint ? segment_.b : segment_.a,tangent,mid,options_.step)) return false;
    } else if (!corrector_.surface(mid,options_.step)) return false;
    const Vec g=gradient(mid); Vec tangent=edge.cross(g);
    if (tangent.norm()<options_.derivative_tolerance*edge.norm()) return false;
    tangent.normalize();
    double previous=edge.norm(),next=edge.norm();
    for (const auto& f:fronts()) {
        if (f.second==e.first) previous=(result_.mesh.vertices[f.first]-a).norm();
        if (f.first==e.second) next=(result_.mesh.vertices[f.second]-b).norm();
    }
    double height=std::sqrt(3.)/2*(previous+edge.norm()+next)/3;
    if (height<options_.min_step) height=.75*height+.25*options_.min_step;
    // Curvature limits prediction; chord checks remain the acceptance test.
    const Hessian h=segment_.hessian(field_.hessians(mid));
    if (!h.allFinite()) return false;
    const double curvature=std::abs(tangent.dot(h*tangent))/g.norm();
    if (curvature>0) height=std::min(height,std::sqrt(6*options_.spatial_tolerance/curvature));
    height=std::min(height,options_.step*1.3);
    // Give existing feature samples first use before growing into their strip.
    if (existing(e,Vec(mid+height*tangent),true)) return true;
    for (int attempt=0; attempt<1 && height>=options_.min_step*.5; ++attempt,height*=.65) {
        Vec trial=mid+height*tangent;
        if (!field_.box.contains(trial)) { if (clip(e,trial)) return true; }
        else if (corrector_.surface(trial,height)) {
            if (new_triangle(e,trial)) return true;
        } else if (!field_.box.contains(trial)) { if (clip(e,trial)) return true; }
    }
    for (const auto& f:fronts()) {
        int id=-1;
        if (f.second==e.first) id=f.first;
        if (f.first==e.second) id=f.second;
        if (id>=0 && id!=e.first && id!=e.second && insert({e.second,e.first,id},true)) return true;
    }
    return false;
}

void Extractor::grow() {
    std::deque<Edge> queue; for (const auto& e:fronts()) queue.push_back(e);
    while (!queue.empty() && result_.mesh.triangles.size()<options_.max_triangles) {
        const Edge e=queue.front(); queue.pop_front(); if (!active(e)) continue;
        const size_t before=result_.mesh.triangles.size();
        if (!advance(e)) continue;
        for (size_t i=before; i<result_.mesh.triangles.size(); ++i) {
            const auto t=result_.mesh.triangles[i];
            for (int k=0; k<3; ++k) if (active({t[k],t[(k+1)%3]})) queue.emplace_back(t[k],t[(k+1)%3]);
        }
    }
}

bool Extractor::covered(const Vec& seed) const {
    const Vec g=gradient(seed);
    for (const auto& t:result_.mesh.triangles) {
        const Vec &a=result_.mesh.vertices[t[0]],&b=result_.mesh.vertices[t[1]],&c=result_.mesh.vertices[t[2]];
        if ((b-a).cross(c-a).dot(g)<=0) continue;
        Vec nearest=marching_triangles::closest_on_triangle(seed,a,b,c);
        if ((nearest-seed).norm()>options_.step*.6) continue;
        if (corrector_.surface(nearest,options_.spatial_tolerance*2) && (nearest-seed).norm()<options_.step*.2) return true;
    }
    return false;
}

bool Extractor::seed_triangle(const Vec& p) {
    const Vec normal=gradient(p).normalized(); Eigen::Index axis; normal.cwiseAbs().minCoeff(&axis);
    const Vec u=normal.cross(Vec::Unit(axis)).normalized(),v=normal.cross(u);
    for (double h=options_.step; h>=options_.min_step; h*=.5) {
        Vec a=p+h*u,b=p+h*(.5*u+std::sqrt(3.)/2*v);
        if (!corrector_.surface(a,h) || !corrector_.surface(b,h)) continue;
        const size_t n=result_.mesh.vertices.size();
        const int ia=vertex(p),ib=vertex(a),ic=vertex(b);
        if (insert({ia,ib,ic},false)) return true;
        rollback_vertices(n);
    }
    return false;
}
} // namespace fiber_surface
