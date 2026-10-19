#include "extract.h"
#include "../isosurface/crack_closing.h"
#include <chrono>
#include <deque>

namespace fiber_surface {
Checkpoint Extractor::finish(Checkpoint state) {
    const auto start=std::chrono::steady_clock::now();
    result_=std::move(state.result); pending_=std::move(state.pending); fiber_uses_=std::move(state.fiber_uses);
    stitch(); repair();
    for (const auto& e:fronts()) result_.unresolved.push_back({(result_.mesh.vertices[e.first]+result_.mesh.vertices[e.second])/2,
        result_.mesh.triangles.size()>=options_.max_triangles ? "Triangle limit reached; unfinished front" : "Unfinished front: no valid growth or existing-vertex connection"});
    for (const auto& e:result_.fiber_edges) if (!result_.mesh.edges.count(e.first))
        result_.unresolved.push_back({result_.mesh.vertices[e.first.first],"Traced fiber edge has no surface attachment"});
    if (result_.mesh.triangles.empty()) result_.unresolved.push_back({(field_.box.min+field_.box.max)/2,"No regular surface discovered; finite search does not prove emptiness"});
    result_.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return {std::move(result_),std::move(pending_),std::move(fiber_uses_)};
}

void Extractor::stitch() {
    for (const auto& e:fronts()) if (pending_.count(e))
        existing(e,Vec((result_.mesh.vertices[e.first]+result_.mesh.vertices[e.second])/2));
    crack_closing::CrackCloser<double>::close_cracks(result_.mesh,[&](const Triangle& t) {
        return insert(t,false);
    },options_.max_triangles);
}

void Extractor::repair() {
    while (!fronts().empty() && result_.mesh.triangles.size()<options_.max_triangles) {
        const auto saved=result_; const auto saved_pending=pending_;
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
        const size_t before=fronts().size(); bool improved=false;
        for (const auto& removed:cavities) {
            result_=saved; pending_=saved_pending;
            Mesh trial; trial.vertices=saved.mesh.vertices;
            std::set<Triangle> forbidden;
            for (size_t i=0; i<saved.mesh.triangles.size(); ++i) {
                const auto t=saved.mesh.triangles[i];
                if (!removed.count(int(i))) trial.add(t);
                else forbidden.insert(marching_triangles::make_triangle_key(t[0],t[1],t[2]));
            }
            result_.mesh=std::move(trial); result_.box_edges.clear();
            for (const auto& e:result_.mesh.boundary) if ((result_.info[e.first.first].faces & result_.info[e.first.second].faces) && !result_.fiber_edges.count(edge_key(e.first.first,e.first.second)))
                result_.box_edges.insert(edge_key(e.first.first,e.first.second));
            for (const auto& e:fiber_uses_) if (!result_.mesh.edges.count(edge_key(e.first.first,e.first.second))) pending_[e.first]=e.second;
            crack_closing::CrackCloser<double>::close_cracks(result_.mesh,[&](const Triangle& t) {
                if (forbidden.count(marching_triangles::make_triangle_key(t[0],t[1],t[2]))) return false;
                return insert(t,false);
            },options_.max_triangles);
            if (fronts().size()<before) {
                std::vector<bool> used(result_.mesh.vertices.size(),false);
                for (const auto& t:result_.mesh.triangles) for (int id:t) used[id]=true;
                if (std::find(used.begin(),used.end(),false)==used.end()) { improved=true; break; }
            }
        }
        if (!improved) { result_=saved; pending_=saved_pending; break; }
    }
    // A remaining triangular contour has a uniquely prescribed cap. Keep
    // the true 3D intersection and function checks; avoid a redundant plane
    // overlap rejection on strongly curved, already fixed neighboring faces.
    for (const auto& loop:result_.mesh.boundary_loops()) if (loop.size()==3)
        insert({loop[0],loop[2],loop[1]},false,false);
}
} // namespace fiber_surface
