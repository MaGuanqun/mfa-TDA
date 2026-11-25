// Diagnose and repair local defects on an assembled native range-triangle boundary.
#include "../pipeline.h"
#include "../corrector.h"
#include "../volume.h"
#include <fstream>
#include <iostream>

using namespace fiber_surface;
namespace {
using Pair=std::pair<size_t,size_t>;

std::vector<Pair> contacts(const VolumeBoundary& boundary,double eps) {
    std::vector<Pair> pairs;
    marching_triangles::visit_triangle_intersections(boundary.mesh.vertices,boundary.mesh.triangles,eps,
        [&](size_t a,size_t b){pairs.emplace_back(a,b);});
    std::sort(pairs.begin(),pairs.end());
    return pairs;
}
void write_contacts(const VolumeBoundary& boundary,const std::vector<Pair>& pairs,const std::string& path) {
    std::ofstream out(path);out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<"face_a,face_b,patch_a,patch_b,shared_vertices,x,y,z\n"<<std::setprecision(16);
    const auto& mesh=boundary.mesh;
    for(auto [a,b]:pairs) {
        const auto& ta=mesh.triangles[a];const auto& tb=mesh.triangles[b];
        int shared=0; Vec center=Vec::Zero();
        for(int id:ta) {shared+=std::count(tb.begin(),tb.end(),id);center+=mesh.vertices[id];}
        for(int id:tb) center+=mesh.vertices[id];
        center/=6;
        out<<a<<','<<b<<','<<boundary.patches[a]<<','<<boundary.patches[b]<<','<<shared<<','
           <<center.x()<<','<<center.y()<<','<<center.z()<<'\n';
    }
}
// All collision queries use all patch labels, including when only one patch is edited.
class JointRepair {
    VolumeBoundary& b;
    const Field& field;
    const RangeTriangle& range;
    const Options& options;
    double eps, normal_cos;
    using Cell=std::array<int,3>;
    std::map<Cell,std::vector<int>> bins;
    std::vector<std::set<int>> incident;
    std::vector<Vec> starting_positions;
    std::set<int> removed;

    template<class Visit> void cells(const Triangle& face,Visit visit) const {
        Vec lo=b.mesh.vertices[face[0]], hi=lo;
        for(int v:face) {lo=lo.cwiseMin(b.mesh.vertices[v]);hi=hi.cwiseMax(b.mesh.vertices[v]);}
        Cell first,last;
        for(int d=0;d<3;++d) {
            first[d]=int(std::floor((lo[d]-eps)/(4*options.step)));
            last[d]=int(std::floor((hi[d]+eps)/(4*options.step)));
        }
        for(int x=first[0];x<=last[0];++x) for(int y=first[1];y<=last[1];++y)
            for(int z=first[2];z<=last[2];++z) visit(Cell{x,y,z});
    }
    void index(int face) {
        const auto& t=b.mesh.triangles[face];
        for(int v:t) incident[v].insert(face);
        // Older bin entries may remain after a flip. Exact tests use the current triangle.
        cells(t,[&](Cell cell){bins[cell].push_back(face);});
    }
    std::vector<int> owners(int a,int c) const {
        std::vector<int> ids;
        std::set_intersection(incident[a].begin(),incident[a].end(),incident[c].begin(),incident[c].end(),
                              std::back_inserter(ids));
        return ids;
    }
    bool geometry(const Triangle& face,int label) const {
        const Vec &a=b.mesh.vertices[face[0]], &c=b.mesh.vertices[face[1]], &d=b.mesh.vertices[face[2]];
        const Vec normal=(c-a).cross(d-a);
        const double longest=std::max({(a-c).norm(),(c-d).norm(),(d-a).norm()});
        if(longest>2*options.step+eps || normal.norm()<=eps*longest) return false;
        for(int v:face) {
            const auto& x=b.mesh.vertices[v];
            if(!field.box.contains(x) || std::abs(range.sides[label-1].phi(field.value(x)))>1.01*options.residual_tolerance) return false;
        }
        for(const Vec& x:volume_triangle_samples(a,c,d)) {
            if(!field.box.contains(x)) return false;
            const Vec expected=-range.sides[label-1].gradient(field.jacobian(x));
            if(!expected.allFinite() || expected.norm()<=options.derivative_tolerance ||
               normal.normalized().dot(expected.normalized())<normal_cos) return false;
            Vec p=x;
            if(!correct_volume_sample(p,label,field,range,options)) return false;
        }
        return true;
    }
    size_t clashes(const std::vector<Triangle>& proposed,const std::vector<int>& replaced) const {
        size_t count=0;
        for(size_t i=0;i<proposed.size();++i) {
            const auto& t=proposed[i];
            std::set<int> nearby;
            cells(t,[&](const Cell& cell) {
                const auto it=bins.find(cell);
                if(it!=bins.end()) nearby.insert(it->second.begin(),it->second.end());
            });
            for(int id:replaced) nearby.erase(id);
            for(int id:nearby) if(!removed.count(id)) count+=marching_triangles::triangles_conflict(b.mesh.vertices,t,b.mesh.triangles[id],eps);
            for(size_t j=0;j<i;++j) count+=marching_triangles::triangles_conflict(b.mesh.vertices,t,proposed[j],eps);
        }
        return count;
    }
    void install(const std::vector<Triangle>& faces,const std::vector<int>& ids,int label) {
        for(size_t i=0;i<faces.size();++i) {
            int id;
            if(i<ids.size()) {
                id=ids[i];
                for(int v:b.mesh.triangles[id]) incident[v].erase(id);
                b.mesh.triangles[id]=faces[i];
            } else {
                id=int(b.mesh.triangles.size());
                b.mesh.triangles.push_back(faces[i]);b.patches.push_back(label);
            }
            index(id);
        }
    }
    bool try_edge(int a,int c) {
        auto ids=owners(a,c);
        if(ids.size()!=2) return false;
        int label=b.patches[ids[0]];
        const int other=b.patches[ids[1]];
        if(label<1 || label>3 || other<1 || other>3) return false;
        const auto f=b.mesh.triangles[ids[0]],g=b.mesh.triangles[ids[1]];
        auto start=std::find(f.begin(),f.end(),a)-f.begin();
        if(f[(start+1)%3]!=c) std::swap(a,c);
        auto pos=std::find(g.begin(),g.end(),c)-g.begin();
        if(g[(pos+1)%3]!=a) return false;
        int v=opposite(f,a,c), w=opposite(g,a,c);
        const size_t before=clashes({f,g},ids);
        if(!before) return false;
        auto improves=[&](const std::vector<Triangle>& trial) {
            return std::all_of(trial.begin(),trial.end(),[&](const Triangle& t){return geometry(t,label);})
                    && clashes(trial,ids)<before;
        };
        if(label!=other) {
            const auto& q=shared_corner(label,other);
            const Vec x=b.mesh.vertices[a], y=b.mesh.vertices[c];
            if((field.value(x)-q).norm()>Options::fiber_field_tolerance ||
               (field.value(y)-q).norm()>Options::fiber_field_tolerance) return false;
            Corrector solver(field,range.sides[label-1],options);
            for(double alpha:{.5,.25,.75}) {
                Vec mpos=x+alpha*(y-x);
                if(!solver.fiber(mpos,q,(y-x).normalized(),mpos,.5*options.step)) continue;
                const int m=int(b.mesh.vertices.size());
                b.mesh.vertices.push_back(mpos);
                const std::vector<Triangle> left{{a,m,v},{m,c,v}},right{{c,m,w},{m,a,w}};
                if(geometry(left[0],label) && geometry(left[1],label) && geometry(right[0],other) && geometry(right[1],other)
                   && clashes({left[0],left[1],right[0],right[1]},ids)<before) {
                    incident.emplace_back();starting_positions.push_back(mpos);
                    install(left,{ids[0]},label);install(right,{ids[1]},other);++seam_splits;return true;
                }
                b.mesh.vertices.pop_back();
            }
            return false;
        }
        if(v!=w && owners(v,w).empty()) {
            std::vector<Triangle> trial{{v,a,w},{v,w,c}};
            if(improves(trial)) {install(trial,ids,label);++flips;return true;}
        }
        // Bisect the common edge on both faces. The new point must belong to this
        // side's interior; an endpoint or a corner requires a shared seam edit.
        Vec mid=volume_midpoint(b.mesh.vertices[a],b.mesh.vertices[c]);
        if(!correct_volume_sample(mid,label,field,range,options)) return false;
        const Range value=field.value(mid);
        for(int k=0;k<3;++k) if(k!=label-1 && range.sides[k].phi(value)<=options.residual_tolerance) return false;
        const int m=int(b.mesh.vertices.size());
        b.mesh.vertices.push_back(mid);
        std::vector<Triangle> trial{{a,m,v},{c,m,w},{m,c,v},{m,a,w}};
        if(improves(trial)) {
            starting_positions.push_back(mid);incident.emplace_back();install(trial,ids,label);++splits;return true;
        }
        b.mesh.vertices.pop_back();return false;
    }
    const Range& shared_corner(int left,int right) const {
        return range.corners[(left%3+1==right) ? right-1 : left-1];
    }
    bool refine_face(int face) {
        const auto t=b.mesh.triangles[face];
        const int label=b.patches[face];
        if(label<1 || label>3) return false;
        const size_t before=clashes({t},{face});
        if(!before) return false;
        const Vec center=volume_centroid(b.mesh.vertices[t[0]],b.mesh.vertices[t[1]],b.mesh.vertices[t[2]]);
        std::vector<Vec> candidates{center};
        for(int v:t) candidates.push_back(volume_midpoint(center,b.mesh.vertices[v]));
        for(Vec point:candidates) {
            if(!correct_volume_sample(point,label,field,range,options)) continue;
            const auto value=field.value(point);
            bool interior=true;
            for(int k=0;k<3;++k) if(k!=label-1 && range.sides[k].phi(value)<=Options::fiber_field_tolerance) interior=false;
            if(!interior) continue;
            const int m=int(b.mesh.vertices.size());
            b.mesh.vertices.push_back(point);
            std::vector<Triangle> trial{{t[0],t[1],m},{t[1],t[2],m},{t[2],t[0],m}};
            if(std::all_of(trial.begin(),trial.end(),[&](const Triangle& tri){return geometry(tri,label);})
               && clashes(trial,{face})<before) {
                incident.emplace_back();starting_positions.push_back(point);
                install(trial,{face},label);++face_splits;return true;
            }
            b.mesh.vertices.pop_back();
        }
        return false;
    }
    bool seam(int a,int c,int left,int right) const {
        if(left==right) return true;
        if(left<1 || left>3 || right<1 || right>3) return false;
        const Range& q=shared_corner(left,right);
        const auto &x=b.mesh.vertices[a], &y=b.mesh.vertices[c];
        const double tolerance=Options::fiber_field_tolerance;
        if((field.value(x)-q).norm()>tolerance || (field.value(y)-q).norm()>tolerance) return false;
        Corrector solver(field,range.sides[left-1],options);
        const Vec dir=(y-x).normalized();
        for(double t:{.25,.5,.75}) {
            Vec p=x+t*(y-x);
            if(!solver.fiber(p,q,dir,p,options.spatial_tolerance)) return false;
        }
        return true;
    }
    bool mixed_fan(const std::vector<int>& loop,size_t apex,const std::set<int>& labels) {
        const size_t n=loop.size();
        auto at=[&](size_t i){return loop[(apex+i)%n];};
        auto outside=[&](int u,int v){return b.patches[owners(u,v).front()];};
        std::vector<Triangle> faces;
        std::vector<std::vector<int>> choices;
        for(size_t i=1;i+1<n;++i) {
            int a=at(0),c=at(i),d=at(i+1);
            if(i>1 && !owners(a,c).empty()) return false;
            faces.push_back({a,c,d});choices.emplace_back();
            for(int label:labels) {
                if(!geometry(faces.back(),label) || !seam(c,d,label,outside(c,d))) continue;
                if(i==1 && !seam(a,c,label,outside(a,c))) continue;
                if(i+2==n && !seam(d,a,label,outside(d,a))) continue;
                choices.back().push_back(label);
            }
            if(choices.back().empty()) return false;
        }
        if(clashes(faces,{})) return false;
        // Dynamic programming over the triangles of the fan. Neighboring side
        // assignments may differ only across a model-confirmed corner fiber.
        std::vector<std::map<int,int>> parent(faces.size());
        for(int first:choices.front()) parent.front()[first]=0;
        for(size_t i=1;i<faces.size();++i) for(int next:choices[i]) {
            for(auto [prev,unused]:parent[i-1]) if(seam(at(0),at(i+1),prev,next)) {
                parent[i][next]=prev;break;
            }
        }
        if(parent.back().empty()) return false;
        int label=parent.back().begin()->first;
        for(size_t i=faces.size();i-->0;) {
            install({faces[i]},{},label);
            label=parent[i].at(label);
        }
        return true;
    }
    bool joins_fans(int v,int prev,int next) const {
        std::map<int,std::vector<int>> link;
        for(int id:incident[v]) {
            const auto& t=b.mesh.triangles[id];
            int i=int(std::find(t.begin(),t.end(),v)-t.begin());
            int a=t[(i+1)%3],c=t[(i+2)%3];
            link[a].push_back(c);link[c].push_back(a);
        }
        std::set<int> seen{prev};std::vector<int> stack{prev};
        while(!stack.empty()) {
            int a=stack.back();stack.pop_back();
            for(int c:link.at(a)) if(seen.insert(c).second) stack.push_back(c);
        }
        // Closing one fan at a vertex shared by several would permanently pinch it.
        return !seen.count(next) || seen.size()==link.size();
    }
    bool cap_loop(const std::vector<int>& loop) {
        std::vector<int> sides;
        std::set<int> labels;
        Vec center=Vec::Zero();
        const size_t n=loop.size();
        for(size_t i=0;i<n;++i) {
            const int a=loop[i],c=loop[(i+1)%n];
            const auto use=owners(a,c);
            if(use.size()!=1 || !joins_fans(a,loop[(i+n-1)%n],c)) return false;
            sides.push_back(b.patches[use.front()]);
            labels.insert(sides.back());
            center+=b.mesh.vertices[a];
        }
        if(*labels.begin()<1 || *labels.rbegin()>3) return false;
        for(size_t apex=0;apex<n;++apex) if(mixed_fan(loop,apex,labels)) return true;
        center/=n;
        if(labels.size()==1) {
            if(!correct_volume_sample(center,sides[0],field,range,options)) return false;
        } else {
            if(labels.size()!=2) return false;
            std::vector<int> changes;
            for(size_t i=0;i<n;++i) if(sides[i]!=sides[(i+n-1)%n]) changes.push_back(loop[i]);
            if(changes.size()!=2) return false;
            const Vec &a=b.mesh.vertices[changes[0]],&c=b.mesh.vertices[changes[1]];
            center=volume_midpoint(a,c);
            Corrector solver(field,range.sides[sides[0]-1],options);
            if(!solver.fiber(center,shared_corner(*labels.begin(),*labels.rbegin()),
                             (c-a).normalized(),center,options.spatial_tolerance)) return false;
        }
        const int m=int(b.mesh.vertices.size());
        b.mesh.vertices.push_back(center);
        std::vector<Triangle> trial;
        for(size_t i=0;i<n;++i) {
            const Triangle t{loop[i],loop[(i+1)%n],m};
            if(!geometry(t,sides[i]) || !seam(loop[i],m,sides[i],sides[(i+n-1)%n])) break;
            trial.push_back(t);
        }
        if(trial.size()==n && !clashes(trial,{})) {
            incident.emplace_back();starting_positions.push_back(center);
            for(size_t i=0;i<n;++i) install({trial[i]},{},sides[i]);
            return true;
        }
        b.mesh.vertices.pop_back();return false;
    }
    bool close_after_move(const std::vector<int>& loop,int depth,size_t first=0) {
        if(cap_loop(loop)) return true;
        if(!depth) return false;
        for(size_t i=first;i<loop.size();++i) {
            int v=loop[i];
            int label=b.patches[*incident[v].begin()];
            if(label<1 || label>3) continue;
            std::vector<int> ids(incident[v].begin(),incident[v].end());
            if(std::any_of(ids.begin(),ids.end(),[&](int id){return b.patches[id]!=label;})) continue;
            const Vec original=b.mesh.vertices[v];
            const auto val=field.value(original);
            bool corner=false;
            for(const auto& q:range.corners) corner|=(val-q).norm()<=Options::fiber_field_tolerance;
            if(corner || field.box.faces(original,eps)) continue;
            Vec g=range.sides[label-1].gradient(field.jacobian(original));
            if(g.norm()<=options.derivative_tolerance) continue;
            const Vec t1=g.unitOrthogonal(),t2=g.normalized().cross(t1);
            std::vector<Triangle> star;
            for(int id:ids) star.push_back(b.mesh.triangles[id]);
            for(double radius:{.25,.5}) for(int direction=0;direction<8;++direction) {
                double theta=direction*std::acos(-1.)/4;
                Vec p=original+radius*options.step*(t1*std::cos(theta)+t2*std::sin(theta));
                if(!correct_volume_sample(p,label,field,range,options) || (p-starting_positions[v]).norm()>.5*options.step) continue;
                bool interior=true;
                for(int k=0;k<3;++k) if(k!=label-1 && range.sides[k].phi(field.value(p))<=Options::fiber_field_tolerance) interior=false;
                if(!interior) continue;
                b.mesh.vertices[v]=p;
                if(std::all_of(star.begin(),star.end(),[&](const Triangle& t){return geometry(t,label);}) && !clashes(star,ids)) {
                    // Make the moved star visible to the cap and to a second move.
                    for(int id:ids) index(id);
                    if(close_after_move(loop,depth-1,i+1)) {++rim_moves;return true;}
                }
                b.mesh.vertices[v]=original;
            }
        }
        return false;
    }
    bool expand_rim(const std::vector<int>& loop,int depth) {
        if(cap_loop(loop)) return true;
        if(!depth) return false;
        // Replace a collar face and the hole in the same transaction. Both
        // remaining edges of that face must have a second owner.
        for(size_t i=0;i<loop.size();++i) {
            int a=loop[i],c=loop[(i+1)%loop.size()];
            auto use=owners(a,c);
            if(use.size()!=1) continue;
            int id=use.front();
            const auto face=b.mesh.triangles[id];
            int w=opposite(face,a,c);
            if(std::find(loop.begin(),loop.end(),w)!=loop.end() ||
               owners(a,w).size()!=2 || owners(c,w).size()!=2 || field.box.faces(b.mesh.vertices[w],eps)) continue;
            std::vector<int> expanded=loop;
            expanded.insert(expanded.begin()+i+1,w);
            for(int v:face) incident[v].erase(id);
            removed.insert(id);
            if(expand_rim(expanded,depth-1)) return true;
            removed.erase(id);
            for(int v:face) incident[v].insert(id);
        }
        return false;
    }
    bool move_vertex(int v) {
        const auto& star=incident[v];
        if(star.empty()) return false;
        const int label=b.patches[*star.begin()];
        if(label<1 || label>3) return false;
        const Vec original=b.mesh.vertices[v];
        if(field.box.faces(original,eps)) return false;
        for(const auto& q:range.corners) if((field.value(original)-q).norm()<=Options::fiber_field_tolerance) return false;
        std::vector<int> ids(star.begin(),star.end());
        std::set<int> neighbors;
        std::vector<Triangle> faces;
        for(int face:ids) {
            if(b.patches[face]!=label) return false;
            faces.push_back(b.mesh.triangles[face]);
            for(int u:faces.back()) if(u!=v) neighbors.insert(u);
        }
        Vec centroid=Vec::Zero();
        for(int u:neighbors) {
            if(owners(v,u).size()!=2) return false;
            centroid+=b.mesh.vertices[u];
        }
        const auto original_count=clashes(faces,ids);
        if(!original_count) return false;
        std::vector<Vec> targets{centroid/neighbors.size()};
        for(int u:neighbors) targets.push_back(b.mesh.vertices[u]);
        for(const Vec& target:targets) for(double t:{.5,.25,.75,-.25}) {
            Vec p=original+t*(target-original);
            if(!correct_volume_sample(p,label,field,range,options) || (p-starting_positions[v]).norm()>.5*options.step) continue;
            if(!range.sides[label-1].contains_parameter(range.sides[label-1].parameter(field.value(p)))) continue;
            bool on_corner=false;
            for(const auto& q:range.corners) on_corner|=(field.value(p)-q).norm()<=Options::fiber_field_tolerance;
            if(on_corner) continue;
            b.mesh.vertices[v]=p;
            if(std::all_of(faces.begin(),faces.end(),[&](const Triangle& face){return geometry(face,label);})
               && clashes(faces,ids)<original_count) {
                for(int id:ids) index(id);
                ++moves;return true;
            }
            b.mesh.vertices[v]=original;
        }
        return false;
    }
    static int opposite(const Triangle& t,int a,int c) {
        for(int v:t) if(v!=a && v!=c) return v;
        throw std::runtime_error("Degenerate mesh edge");
    }
public:
    size_t flips=0,splits=0,seam_splits=0,face_splits=0,moves=0,rim_moves=0,fills=0;
    JointRepair(VolumeBoundary& boundary,const Field& f,const RangeTriangle& r,const Options& o):
        b(boundary),field(f),range(r),options(o),eps(o.step*1e-8),
        normal_cos(std::cos(o.max_normal_degrees*std::acos(-1.)/180)),incident(b.mesh.vertices.size()),starting_positions(b.mesh.vertices) {
        for(size_t i=0;i<b.mesh.triangles.size();++i) index(int(i));
    }
    void untangle(std::vector<Pair> pairs) {
        while(!pairs.empty()) {
            std::set<int> faces;
            for(auto [a,c]:pairs) {faces.insert(int(a));faces.insert(int(c));}
            bool changed=false;
            for(int face:faces) {
                const auto t=b.mesh.triangles[face];
                bool edited=false;
                for(int i=0;i<3 && !edited;++i) edited=try_edge(t[i],t[(i+1)%3]);
                for(int v:t) if(!edited) edited=move_vertex(v);
                if(!edited) edited=refine_face(face);
                changed|=edited;
            }
            if(!changed) break;
            auto next=contacts(b,eps);
            if(next.size()>=pairs.size()) throw std::runtime_error("Local and global contact checks disagree");
            pairs=std::move(next);
        }
    }
    void fill_rims() {
        // Every proposed disk must follow oriented open edges. Check each
        // cross-label join against its corner fiber and every face against all patches.
        std::map<int,std::vector<int>> rim,out;
        for(size_t i=0;i<b.mesh.triangles.size();++i) {
            const auto& t=b.mesh.triangles[i];
            for(int j=0;j<3;++j) {
                int a=t[j],c=t[(j+1)%3];
                if(owners(a,c).size()!=1) continue;
                rim[a].push_back(c);rim[c].push_back(a);out[c].push_back(a);
            }
        }
        std::set<int> seen;
        for(const auto& [start,neighbors]:rim) {
            if(seen.count(start)) continue;
            std::vector<int> component{start};seen.insert(start);
            for(size_t i=0;i<component.size();++i) for(int next:rim.at(component[i]))
                if(seen.insert(next).second) component.push_back(next);
            if(component.size()<3 || std::any_of(component.begin(),component.end(),[&](int v){
                return field.box.faces(b.mesh.vertices[v],eps);
            })) continue;
            if(std::all_of(component.begin(),component.end(),[&](int v){return rim[v].size()==2 && out[v].size()==1;})) {
                std::vector<int> loop{start};
                while(out[loop.back()][0]!=start) loop.push_back(out[loop.back()][0]);
                if(expand_rim(loop,2)) ++fills;
                else {
                    bool closable=true;
                    for(size_t i=0;i<loop.size();++i)
                        closable &= joins_fans(loop[i],loop[(i+loop.size()-1)%loop.size()],loop[(i+1)%loop.size()]);
                    if(closable) fills+=close_after_move(loop,2);
                }
                continue;
            }
            // At a branch try oriented lobes returning to the same junction.
            // Paths meeting another junction have ambiguous connectivity.
            for(int v:component) if(rim[v].size()>2) for(int next:out[v]) {
                std::vector<int> loop{v};
                while(next!=v && rim[next].size()==2 && out[next].size()==1) {
                    loop.push_back(next);next=out[next].front();
                }
                if(next==v && loop.size()>=3) fills+=expand_rim(loop,2);
            }
        }
        size_t dest=0;
        for(size_t src=0;src<b.mesh.triangles.size();++src) if(!removed.count(int(src))) {
            b.mesh.triangles[dest]=b.mesh.triangles[src];b.patches[dest]=b.patches[src];++dest;
        }
        b.mesh.triangles.resize(dest);b.patches.resize(dest);
    }
};

}
int main(int argc,char** argv) {
    try {
        if(argc!=4 && argc!=5) throw std::invalid_argument("Use: fiber_repair_boundary STAGES INPUT.ply OUTPUT_PREFIX [--repair]");
        const bool repair=argc==5 && std::string(argv[4])=="--repair";
        if(argc==5 && !repair) throw std::invalid_argument("Expected --repair");
        const auto cfg=read_config(argv[1]);
        if(cfg.curve.size()!=3 || !cfg.curve.closed()) throw std::invalid_argument("Expected a closed range triangle");
        auto b=read_boundary_ply(argv[2]);
        const double eps=cfg.options.step*1e-8;
        auto topology=[&](){return marching_triangles::validate_mesh(b.mesh.vertices,b.mesh.triangles,eps,false);};
        const auto input=topology();
        if(input.degenerate_faces || input.duplicate_faces || input.nonmanifold_edges || input.inconsistent_edges)
            throw std::invalid_argument("Expected nondegenerate triangles with consistent oriented edge incidence");
        const std::string prefix=argv[3];
        auto pairs=contacts(b,eps);
        write_contacts(b,pairs,prefix+"_pairs.csv");
        std::cout<<"Input intersecting triangle pairs: "<<pairs.size()<<std::endl;
        if(repair) {
            const auto model=load_field(cfg,true);
            const RangeTriangle range(cfg.curve.point(0),cfg.curve.point(1),cfg.curve.point(2));
            JointRepair edits(b,model.field,range,cfg.options);
            edits.untangle(pairs);
            edits.fill_rims();
            const auto output=topology();
            const auto remaining=contacts(b,eps);
            if(output.invalid_indices || output.degenerate_faces || output.duplicate_faces ||
               output.nonmanifold_edges || output.inconsistent_edges || output.unused_vertices>input.unused_vertices ||
               output.nonmanifold_vertices>input.nonmanifold_vertices || output.boundary_edges>input.boundary_edges ||
               remaining.size()>pairs.size()) throw std::runtime_error("Repair failed the global topology check");
            write_contacts(b,remaining,prefix+"_remaining_pairs.csv");
            write_boundary_ply(prefix+".ply",b,false);
            std::cout<<"Edge flips "<<edits.flips<<", projected edge bisections "<<edits.splits<<", shared seam bisections "<<edits.seam_splits<<", face subdivisions "<<edits.face_splits<<", projected vertex moves "<<edits.moves<<", rim vertex moves "<<edits.rim_moves<<", filled rims "<<edits.fills<<'\n';
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
