#pragma once

#include "crack_closing.h"
#include "mesh_quality.h"
#include <deque>
#include <fstream>
#include <iomanip>
#include <string>

namespace marching_triangles {
template<class T> class MarchingTriangles {
public:
    using Mesh = TriangleMesh<T>;
    MarchingTriangles(SurfaceField<T> field, MeshingOptions<T> options)
        : field_(std::move(field)), options_(std::move(options)), projector_(field_,options_) {
        options_.validate();
        if(!field_.value||!field_.gradient) throw std::invalid_argument("A scalar field and its gradient are required");
    }
    // The projector holds references to this object's field/options.
    MarchingTriangles(const MarchingTriangles&) = delete;
    MarchingTriangles& operator=(const MarchingTriangles&) = delete;
    void set_degenerate_points(const std::vector<Point<T>>& points) { stop_points_=points; }
    const std::string& diagnostic() const {return diagnostic_;}

    bool extract_all_sheets(const std::vector<Point<T>>& seeds,
                            std::vector<std::vector<Point<T>>>& sheet_vertices,
                            std::vector<std::vector<Triangle>>& sheet_triangles) {
        sheet_vertices.clear(); sheet_triangles.clear(); diagnostic_.clear();
        std::vector<Mesh> sheets;
        for(auto seed:seeds) {
            if(!projector_.project(seed)||near_stop(seed)) continue;
            bool covered=false;
            for(const Mesh& sheet:sheets) if(covers(sheet,seed)) {covered=true;break;}
            if(covered) continue;
            Mesh mesh;
            if(!seed_triangle(mesh,seed)) continue;
            grow(mesh);
            if(options_.close_cracks && mesh.triangles.size()<options_.max_triangles) {
                crack_closing::CrackCloser<T>::close_cracks(mesh,[&](const Triangle& t) {
                    return insert(mesh,t,false);
                },options_.max_triangles);
                repair_stalled_cracks(mesh);
                // The tangent-plane guard is deliberately conservative. A
                // final triangular hole already has all three edges fixed;
                // admit its cap using actual 3D intersection, orientation and
                // projection tests (the small-hole rule in section 3.2.2).
                for(const auto& loop:mesh.boundary_loops())
                    if(loop.size()==3) insert(mesh,{loop[0],loop[2],loop[1]},false,false);
            }
            // Closure alone can leave valid but arbitrarily thin triangles.
            // Improve the closed sheet while preserving its topology and
            // projecting every relocated vertex back onto the level set.
            MeshQualityOptimizer<T> optimizer(field_,options_,stop_points_);
            optimizer.improve(mesh);
            if(!mesh.boundary.empty()) {
                diagnostic_="Surface remains open ("+std::to_string(mesh.boundary.size())+
                    " boundary edges). Check domain clipping, singularities, step size, and triangle limit.";
            }
            if(!mesh.boundary.empty() && mesh.triangles.size()>=options_.max_triangles)
                diagnostic_="Triangle limit reached before extraction completed.";
            sheets.push_back(std::move(mesh));
        }
        for(auto& sheet:sheets) {
            sheet_vertices.push_back(std::move(sheet.vertices));
            sheet_triangles.push_back(std::move(sheet.triangles));
        }
        return !sheets.empty();
    }
    static void save_mesh_obj(const std::string& filename,const std::vector<Point<T>>& vertices,
                              const std::vector<Triangle>& triangles) {
        std::ofstream out(filename);
        if(!out) throw std::runtime_error("Cannot write "+filename);
        out<<std::setprecision(std::numeric_limits<T>::max_digits10);
        for(const auto& v:vertices) out<<"v "<<v[0]<<' '<<v[1]<<' '<<v[2]<<'\n';
        for(const auto& t:triangles) out<<"f "<<t[0]+1<<' '<<t[1]+1<<' '<<t[2]+1<<'\n';
        out.close();
        if(!out) throw std::runtime_error("Failed writing "+filename);
    }
private:
    SurfaceField<T> field_;
    MeshingOptions<T> options_;
    SurfaceProjector<T> projector_;
    std::vector<Point<T>> stop_points_;
    std::string diagnostic_;
    T geometric_epsilon() const {return options_.step*T(1e-9);}
    bool near_stop(const Point<T>& p) const {
        for(const auto& q:stop_points_) if((p-q).norm()<options_.stop_distance) return true;
        return false;
    }
    bool covers(const Mesh& mesh,const Point<T>& p) const {
        const Point<T> g=field_.gradient(p);
        for(const auto& t:mesh.triangles) {
            const auto& a=mesh.vertices[t[0]];const auto& b=mesh.vertices[t[1]];const auto& c=mesh.vertices[t[2]];
            const Point<T> n=(b-a).cross(c-a);
            if(n.dot(g)<=0) continue;
            Point<T> q=closest_on_triangle(p,a,b,c);
            const T h=std::max({(a-b).norm(),(b-c).norm(),(c-a).norm()});
            if((q-p).norm()>h*T(.5)) continue;
            const Point<T> chord=q;
            if(!projector_.project(q,h*T(.5))) continue;
            // A seed belongs to the triangle's curved surface patch. Comparing
            // projection displacement distinguishes nearby disconnected sheets.
            const T tolerance=std::max(options_.vertex_tolerance,T(.15)*h);
            if((q-p).norm()<tolerance && (q-chord).norm()<=h*T(.5)) return true;
        }
        return false;
    }
    bool geometry_ok(const Mesh& mesh,const Triangle& t,bool delaunay,bool check_surface_overlap=true) const {
        if(!mesh.topology_ok(t)) return false;
        const Point<T>& a=mesh.vertices[t[0]];const Point<T>& b=mesh.vertices[t[1]];const Point<T>& c=mesh.vertices[t[2]];
        const T h=std::max({(a-b).norm(),(b-c).norm(),(c-a).norm()});
        const T shortest=std::min({(a-b).norm(),(b-c).norm(),(c-a).norm()});
        if(h>2*options_.step || shortest<=options_.vertex_tolerance) return false;
        const Point<T> n=(b-a).cross(c-a);
        if(n.norm()<=T(1e-8)*h*h) return false;
        // Respect the gradient orientation, including during crack closing.
        for(const Point<T>& p:std::array<Point<T>,4>{a,b,c,(a+b+c)/3}) {
            const Point<T> g=field_.gradient(p);
            if(!g.allFinite() || n.dot(g)<=T(1e-8)*n.norm()*g.norm()) return false;
            if(near_stop(p)) return false;
        }
        // Do not bridge across a thin cavity or cap an unrelated tubular end.
        // Midpoints and centroid must project to the same nearby surface patch.
        for(Point<T> p:std::array<Point<T>,4>{(a+b)/2,(b+c)/2,(c+a)/2,(a+b+c)/3}) {
            const Point<T> original=p;
            if(!projector_.project(p,h*T(.35)) || (p-original).norm()>h*T(.35)) return false;
        }
        Point<T> center; T r2=0;
        if(delaunay && !circumcircle(a,b,c,center,r2)) return false;
        const Point<T> surface_normal=(field_.gradient(a).normalized()+field_.gradient(b).normalized()+field_.gradient(c).normalized()).normalized();
        for(const auto& old:mesh.triangles) {
            if(triangles_conflict(mesh.vertices,t,old,geometric_epsilon())) return false;
            if(check_surface_overlap && surface_overlap(mesh.vertices,t,old,geometric_epsilon(),surface_normal)) return false;
            if(!delaunay) continue;
            const Point<T> on=(mesh.vertices[old[1]]-mesh.vertices[old[0]]).cross(mesh.vertices[old[2]]-mesh.vertices[old[0]]);
            if(on.dot(n)<=0) continue;
            int shared=0;
            for(int v:old) {
                if(std::find(t.begin(),t.end(),v)!=t.end()) {++shared;continue;}
                if((mesh.vertices[v]-center).squaredNorm()<r2*(1-T(1e-9))) return false;
            }
            // Global surface test, not just incident vertices. For adjacent
            // faces, shared vertices/edges necessarily touch the circumsphere;
            // use their non-shared vertices plus the exact overlap test above.
            if(!shared && (closest_on_triangle(center,mesh.vertices[old[0]],mesh.vertices[old[1]],mesh.vertices[old[2]])-center).squaredNorm()<r2*(1-T(1e-9))) return false;
        }
        return true;
    }
    bool insert(Mesh& mesh,const Triangle& t,bool delaunay,bool check_surface_overlap=true) const {
        if(mesh.triangles.size()>=options_.max_triangles||!geometry_ok(mesh,t,delaunay,check_surface_overlap)) return false;
        mesh.add(t);
        return true;
    }
    // A greedy split can leave an unmeshable sliver. Re-triangulate its local
    // cavity with a different diagonal instead of inserting an inverted face.
    // Each accepted trial strictly reduces the number of boundary edges.
    void repair_stalled_cracks(Mesh& mesh) const {
        while(!mesh.boundary.empty()) {
            std::set<int> adjacent;
            for(const auto& e:mesh.boundary) adjacent.insert(e.second);
            std::vector<std::set<int>> cavities;
            for(int face:adjacent) cavities.push_back({face});
            for(const auto& loop:mesh.boundary_loops()) {
                if(loop.size()>6) continue;
                std::set<int> ring;
                for(size_t i=0;i<loop.size();++i)
                    ring.insert(mesh.boundary.at({loop[i],loop[(i+1)%loop.size()]}));
                for(int a:ring)for(int b:ring)if(a<b)cavities.push_back({a,b});
                if(ring.size()>2)cavities.push_back(ring);
            }
            bool improved=false;
            for(const auto& removed:cavities) {
                Mesh trial;trial.vertices=mesh.vertices;
                std::set<Triangle> forbidden;
                for(size_t i=0;i<mesh.triangles.size();++i)
                    if(!removed.count(static_cast<int>(i))) trial.add(mesh.triangles[i]);
                    else {
                        const auto& old=mesh.triangles[i];
                        forbidden.insert(make_triangle_key(old[0],old[1],old[2]));
                    }
                crack_closing::CrackCloser<T>::close_cracks(trial,[&](const Triangle& t) {
                    if(forbidden.count(make_triangle_key(t[0],t[1],t[2]))) return false;
                    return insert(trial,t,false);
                },options_.max_triangles);
                if(trial.boundary.size()<mesh.boundary.size()) {
                    std::vector<bool> used(trial.vertices.size(),false);
                    for(const auto& t:trial.triangles)for(int i:t)used[i]=true;
                    if(std::find(used.begin(),used.end(),false)!=used.end()) continue;
                    mesh=std::move(trial);improved=true;break;
                }
            }
            if(!improved) break;
        }
    }
    bool seed_triangle(Mesh& mesh,const Point<T>& seed) const {
        const Point<T> g=field_.gradient(seed).normalized();
        if(!g.allFinite()) return false;
        Eigen::Index axis; g.cwiseAbs().minCoeff(&axis);
        const Point<T> u=g.cross(Point<T>::Unit(axis)).normalized(),v=g.cross(u);
        T step=options_.step;
        for(int attempt=0;attempt<12;++attempt,step*=T(.5)) {
            Point<T> p=seed+step*u,q=seed+step*(T(.5)*u+T(std::sqrt(3)/2)*v);
            if(!projector_.project(p,step)||!projector_.project(q,step)) continue;
            mesh.vertices={seed,p,q};
            if(insert(mesh,{0,1,2},false)) return true;
        }
        mesh.vertices.clear(); return false;
    }
    void grow(Mesh& mesh) const {
        std::deque<Edge> queue;
        for(const auto& e:mesh.boundary) queue.push_back(e.first);
        auto enqueue=[&](const Triangle& t) {
            for(int i=0;i<3;++i) if(mesh.boundary.count({t[i],t[(i+1)%3]})) queue.emplace_back(t[i],t[(i+1)%3]);
        };
        while(!queue.empty() && mesh.triangles.size()<options_.max_triangles) {
            const Edge e=queue.front();queue.pop_front();
            if(!mesh.boundary.count(e)) continue;
            const int a=e.first,b=e.second;
            int prev=-1,next=-1;
            for(const auto& f:mesh.boundary) {
                if(f.first.second==a) prev=f.first.first;
                if(f.first.first==b) next=f.first.second;
            }
            const Point<T> edge=mesh.vertices[b]-mesh.vertices[a];
            Point<T> mid=(mesh.vertices[a]+mesh.vertices[b])/2;
            if(!projector_.project(mid,edge.norm())) continue;
            Point<T> tangent=edge.cross(field_.gradient(mid));
            const T tn=tangent.stableNorm();
            if(!tangent.allFinite()||tn<=std::numeric_limits<T>::min()) continue;
            tangent/=tn;
            const T h0=prev<0?edge.norm():(mesh.vertices[a]-mesh.vertices[prev]).norm();
            const T h1=next<0?edge.norm():(mesh.vertices[next]-mesh.vertices[b]).norm();
            T d=T(std::sqrt(3)/2)*(h0+edge.norm()+h1)/3;
            if(d<options_.min_step) d=T(.75)*d+T(.25)*options_.min_step;
            d=std::min(d,options_.step*T(1.5));
            Point<T> p=mid+d*tangent;
            bool added=false;
            if(projector_.project(p,d)) {
                bool duplicate=false;
                for(const auto& v:mesh.vertices) if((v-p).norm()<=options_.vertex_tolerance) {duplicate=true;break;}
                if(!duplicate) {
                    const int id=static_cast<int>(mesh.vertices.size());
                    mesh.vertices.push_back(p);
                    const Triangle t{a,id,b};
                    if(insert(mesh,t,true)) {enqueue(t);added=true;}
                    else mesh.vertices.pop_back();
                }
            }
            if(added) continue;
            for(int c:std::array<int,2>{prev,next}) {
                if(c<0||c==a||c==b) continue;
                const Triangle t{a,c,b};
                if(insert(mesh,t,true)) {enqueue(t);break;}
            }
        }
    }
};
} // namespace marching_triangles
