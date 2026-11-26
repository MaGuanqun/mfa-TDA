#include "extract.h"
#include "../isosurface/crack_closing.h"
#include <chrono>
#include <deque>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace fiber_surface {
namespace {
bool all_vertices_used(const Mesh& mesh) {
    std::vector<bool> used(mesh.vertices.size(),false);
    for (const auto& triangle:mesh.triangles)
        for (int id:triangle) used[id]=true;
    return std::find(used.begin(),used.end(),false)==used.end();
}

// A last short front can connect contacts on two different box faces without
// reaching their shared box edge. Close only that explicit transverse contact,
// solving the scalar constraint with BOTH box coordinates fixed. This is not an
// arbitrary cap of an unfinished loop, and never changes frozen corner fibers.
size_t finish_box_contacts(Result& r,const Field& f,const Segment& s,const Options& o) {
    size_t count=0;Corrector corrector(f,s,o);
    const auto fronts=r.mesh.boundary;const double eps=o.step*1e-9;
    for(const auto& entry:fronts) {
        if(r.mesh.triangles.size()>=o.max_triangles)break;
        const auto e=entry.first;const int a=e.first,b=e.second;
        if(!r.mesh.boundary.count(e) || r.fiber_edges.count(edge_key(a,b)))continue;
        const unsigned fa=r.info[a].faces,fb=r.info[b].faces;
        if(!fa || !fb || (fa&fb))continue;
        const Vec mid=r.mesh.vertices[a]+(r.mesh.vertices[b]-r.mesh.vertices[a])/2;Vec x=mid;
        if(!corrector.surface(x,o.step,fa|fb) || (x-mid).norm()>o.step)continue;
        bool nearby=false;for(const Vec& p:r.mesh.vertices)if((x-p).norm()<eps)nearby=true;
        if(nearby)continue;
        const int id=int(r.mesh.vertices.size());const Triangle t{b,a,id};r.mesh.vertices.push_back(x);
        const Vec cross=(r.mesh.vertices[a]-r.mesh.vertices[b]).cross(x-r.mesh.vertices[b]);
        bool valid=cross.norm()>eps*eps && cross.dot(s.gradient(f.jacobian(mid)))>0 && r.mesh.topology_ok(t);
        for(int k=0;k<3 && valid;++k) {
            const Vec p=r.mesh.vertices[t[k]],q=r.mesh.vertices[t[(k+1)%3]],z=r.mesh.vertices[t[(k+2)%3]];
            const double angle=std::acos(std::clamp((q-p).normalized().dot((z-p).normalized()),-1.,1.))*180/std::acos(-1.);
            if((q-p).norm()>2*o.step || angle<o.min_angle_degrees)valid=false;
            if(!corrector.member(p,t[k]==id ? -1 : r.info[t[k]].endpoint))valid=false;
        }
        if(valid && !corrector.member(face_centroid(r.mesh,t)))valid=false;
        for(const auto& old:r.mesh.triangles)if(valid && marching_triangles::triangles_conflict(r.mesh.vertices,t,old,eps))valid=false;
        if(!valid){r.mesh.vertices.pop_back();continue;}
        r.info.push_back({-1,-1,fa|fb});r.mesh.add(t);
        r.box_edges.insert(edge_key(a,id));r.box_edges.insert(edge_key(b,id));
        ++count;
    }
    return count;
}
} // namespace

Checkpoint Extractor::finish(Checkpoint state, bool curve_patch) {
    const auto start=std::chrono::steady_clock::now();
    index_growth_=false;face_bins_.clear();indexed_faces_=0;
    result_=std::move(state.result); pending_=std::move(state.pending); fiber_uses_=std::move(state.fiber_uses);
    index_stitch_faces_=((field_.box.max-field_.box.min)/options_.step).maxCoeff()<1e8 &&
        std::all_of(result_.mesh.triangles.begin(),result_.mesh.triangles.end(),[&](const Triangle& t) {
            const auto& v=result_.mesh.vertices;
            return (v[t[1]]-v[t[0]]).norm()<=2*options_.step && (v[t[2]]-v[t[0]]).norm()<=2*options_.step;
        });
    stitch(); repair();
    index_stitch_faces_=false;face_bins_.clear();indexed_faces_=0;
    if (curve_patch) result_.box_contact_triangles=finish_box_contacts(result_,field_,segment_,options_);
    snap_box_slivers();
    for (const auto& e:fronts()) result_.unresolved.push_back({(result_.mesh.vertices[e.first]+result_.mesh.vertices[e.second])/2,
        result_.mesh.triangles.size()>=options_.max_triangles ? "Triangle limit reached; unfinished front" : "Unfinished front: no valid growth or existing-vertex connection"});
    for (const auto& [edge,component]:result_.fiber_edges) {
        if (result_.mesh.edges.count(edge)) continue;
        const Vec mid=edge_midpoint(result_.mesh,edge.first,edge.second);
        const char* reason="Traced fiber edge has no surface attachment";
        if (!corrector_.member(result_.mesh.vertices[edge.first],result_.info[edge.first].endpoint) ||
            !corrector_.member(result_.mesh.vertices[edge.second],result_.info[edge.second].endpoint))
            reason="Traced fiber vertex exceeds endpoint residual tolerance";
        result_.unresolved.push_back({mid,reason});
    }
    if (result_.mesh.triangles.empty() && (!curve_patch || result_.surface_seed_hits)) result_.unresolved.push_back({(field_.box.min+field_.box.max)/2,"No regular surface discovered; finite search does not prove emptiness"});
    result_.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return {std::move(result_),std::move(pending_),std::move(fiber_uses_)};
}

void Extractor::stitch() {
    for (const auto& e:fronts()) if (pending_.count(e))
        existing(e,Vec((result_.mesh.vertices[e.first]+result_.mesh.vertices[e.second])/2));
    crack_closing::CrackCloser<double>::close_cracks(result_.mesh,[&](const Triangle& t) {
        return insert(t,false);
    },options_.max_triangles,[&](int a,int b){return !intended({a,b});},2*options_.step);
}

void Extractor::snap_box_slivers() {
    auto& mesh=result_.mesh;
    const double repair_distance=0.01*options_.step;
    std::map<int,std::vector<int>> neighbors;
    for (const auto& [edge,face]:mesh.boundary) {
        neighbors[edge.first].push_back(edge.second);
        neighbors[edge.second].push_back(edge.first);
    }
    for (const auto& [tip,ends]:neighbors) {
        auto& info=result_.info[tip];
        if (ends.size()!=2 || info.fiber>=0 || info.endpoint>=0 || info.faces) continue;
        const int a=ends[0], b=ends[1];
        if (intended({tip,a}) || intended({tip,b}) || mesh.edges.count(edge_key(a,b))) continue;
        const unsigned face=result_.info[a].faces & result_.info[b].faces;
        if (!face || (face & (face-1))) continue;

        const Vec original=mesh.vertices[tip];
        const Vec& start=mesh.vertices[a];
        const Vec& end=mesh.vertices[b];
        const Vec chord=end-start;
        auto near_chord=[&](const Vec& x) {
            const double along=(x-start).dot(chord);
            return along>0 && along<chord.squaredNorm() &&
                   (x-marching_triangles::closest_on_segment(x,start,end)).norm()<=repair_distance;
        };
        if (!near_chord(original)) continue;
        Vec point=original;
        if (!corrector_.surface(point,repair_distance,face) || !near_chord(point) ||
            field_.box.faces(point,epsilon())!=face) continue;
        const double tau=segment_.parameter(field_.value(point));
        if (!(tau>0 && tau<1)) continue; // Segment ends belong to the traced fibers.

        // Keep the vertex ID and its incident triangles. Commit both seam edges
        // together only after the changed triangles and both curved chords pass.
        mesh.vertices[tip]=point;
        bool valid=true;
        for (int id:ends) {
            Vec mid=edge_midpoint(mesh,tip,id);
            if (!corrector_.surface(mid,options_.spatial_tolerance,face)) { valid=false; break; }
        }
        for (const auto& t:mesh.triangles) {
            if (!valid) break;
            const auto corner=std::find(t.begin(),t.end(),tip);
            if (corner==t.end()) continue;
            const auto k=corner-t.begin();
            const int u=t[(k+1)%3], v=t[(k+2)%3];
            const Vec base=mesh.vertices[v]-mesh.vertices[u];
            const Vec old_normal=base.cross(original-mesh.vertices[u]);
            const Vec new_normal=base.cross(point-mesh.vertices[u]);
            valid=!(result_.info[u].faces & result_.info[v].faces & face) &&
                  old_normal.dot(new_normal)>0 && geometry(t,false,false);
        }
        if (!valid) { mesh.vertices[tip]=original; continue; }
        info.faces=face;
        result_.box_edges.insert(edge_key(tip,a));
        result_.box_edges.insert(edge_key(tip,b));
    }
}

std::optional<Extractor::CavityRepair> Extractor::try_cavity(const std::set<int>& removed,size_t before) const {
    Extractor trial(field_,segment_,options_);
    trial.index_stitch_faces_=index_stitch_faces_;
    auto& r=trial.result_;
    r.info=result_.info;
    r.fiber_edges=result_.fiber_edges;
    r.mesh.vertices=result_.mesh.vertices;
    std::set<Triangle> forbidden;
    for (size_t i=0;i<result_.mesh.triangles.size();++i) {
        const auto& t=result_.mesh.triangles[i];
        if (!removed.count(int(i))) r.mesh.add(t);
        else forbidden.insert(marching_triangles::make_triangle_key(t[0],t[1],t[2]));
    }
    // Crack closure can only use vertices still present on the mesh.
    if (!all_vertices_used(r.mesh)) return {};
    for (const auto& e:r.mesh.boundary)
        if ((r.info[e.first.first].faces & r.info[e.first.second].faces) && !r.fiber_edges.count(edge_key(e.first.first,e.first.second)))
            r.box_edges.insert(edge_key(e.first.first,e.first.second));
    trial.pending_=pending_;
    for (const auto& e:fiber_uses_) if (!r.mesh.edges.count(edge_key(e.first.first,e.first.second))) trial.pending_[e.first]=e.second;
    crack_closing::CrackCloser<double>::close_cracks(r.mesh,[&](const Triangle& t) {
        if (forbidden.count(marching_triangles::make_triangle_key(t[0],t[1],t[2]))) return false;
        return trial.insert(t,false);
    },options_.max_triangles,[&](int a,int b){return !trial.intended({a,b});},2*options_.step);
    if (trial.fronts().size()>=before) return {};
    return CavityRepair{std::move(r.mesh),std::move(r.box_edges),std::move(trial.pending_)};
}

void Extractor::repair() {
    // Removing faces and closing their boundary loops can only reuse vertices of
    // existing faces. An initially unused vertex would invalidate every trial.
    const bool can_repair_cavities=all_vertices_used(result_.mesh);
    while (can_repair_cavities && !fronts().empty() && result_.mesh.triangles.size()<options_.max_triangles) {
        std::set<int> adjacent;
        for (const auto& e:result_.mesh.boundary) if (!intended(e.first)) adjacent.insert(e.second);
        std::vector<std::set<int>> cavities;
        for (int f:adjacent) cavities.push_back({f});
        for (const auto& loop:result_.mesh.boundary_loops()) if (loop.size()<=8) {
            std::set<int> ring;
            for (size_t i=0; i<loop.size(); ++i) {
                auto found=result_.mesh.boundary.find({loop[i],loop[(i+1)%loop.size()]});
                if (found!=result_.mesh.boundary.end() && !intended(found->first)) ring.insert(found->second);
            }
            for (int a:ring) for (int b:ring) if (a<b) cavities.push_back({a,b});
            if (ring.size()>2) cavities.push_back(ring);
        }
        const size_t before=fronts().size();
        // Each concurrent trial owns a mesh. Bound the total face count and
        // commit the first improvement in the original cavity order.
        if (cavities.empty()) break;
        const size_t face_budget=std::max<size_t>(1,2000000 / result_.mesh.triangles.size());
        const size_t workers=field_.value_jacobian ? std::min<size_t>(tbb::this_task_arena::max_concurrency(),face_budget) : 1;
        std::optional<CavityRepair> winner;
        for (size_t first=0;first<cavities.size() && !winner;first+=workers) {
            std::vector<std::optional<CavityRepair>> trials(std::min(workers,cavities.size()-first));
            auto attempt=[&](size_t i) { trials[i]=try_cavity(cavities[first+i],before); };
            if (trials.size()==1) attempt(0);
            else tbb::parallel_for(size_t(0),trials.size(),attempt);
            for (auto& trial:trials) if (trial) { winner=std::move(trial);break; }
        }
        if (!winner) break;
        result_.mesh=std::move(winner->mesh);
        result_.box_edges=std::move(winner->box_edges);
        pending_=std::move(winner->pending);
        face_bins_.clear();indexed_faces_=0;
    }
    // A remaining triangular contour has a uniquely prescribed cap. Keep
    // the true 3D intersection and function checks; avoid a redundant plane
    // overlap rejection on strongly curved, already fixed neighboring faces.
    for (const auto& loop:result_.mesh.boundary_loops()) if (loop.size()==3)
        insert({loop[0],loop[2],loop[1]},false,false);
}
} // namespace fiber_surface
