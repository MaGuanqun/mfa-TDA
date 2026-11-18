#include "volume.h"
#include "extract.h"
#include <mapbox/earcut.hpp>
#include <deque>

namespace fiber_surface {
RangeTriangle::RangeTriangle(Range a, Range b, Range c)
    : corners{a,b,c}, sides{Segment(a,b),Segment(b,c),Segment(c,a)} {
    const Range u=b-a,v=c-a;
    const double cross=u[0]*v[1]-u[1]*v[0];
    if (!std::isfinite(cross) || cross<=1e-12*u.norm()*v.norm())
        throw std::invalid_argument("Range triangle must have three noncollinear, counterclockwise vertices");
}
bool RangeTriangle::contains(const Range& f,double tolerance) const {
    if (!f.allFinite()) return false;
    for (const auto& s:sides) if (s.phi(f)<-tolerance) return false;
    return true;
}
namespace {
void event(VolumeBoundary& b,const Vec& x,const std::string& why) { b.unresolved.push_back({x,why}); }
int vertex(VolumeBoundary& b,const Vec& x,double eps) {
    for (size_t i=0;i<b.mesh.vertices.size();++i) if ((b.mesh.vertices[i]-x).norm()<eps) return int(i);
    b.mesh.vertices.push_back(x); return int(b.mesh.vertices.size())-1;
}
void facet(VolumeBoundary& b,Triangle t,int label) {
    b.mesh.add(t); b.patches.push_back(label);
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
        const Vec mid=volume_midpoint(r.mesh.vertices[a],r.mesh.vertices[b]);Vec x=mid;
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
            Vec sample=volume_midpoint(p,q);const Vec original=sample;
            const unsigned fp=t[k]==id ? fa|fb : r.info[t[k]].faces;
            const unsigned fq=t[(k+1)%3]==id ? fa|fb : r.info[t[(k+1)%3]].faces;
            if(!corrector.member(sample) || !corrector.surface(sample,o.spatial_tolerance,fp&fq) || (sample-original).norm()>o.spatial_tolerance)valid=false;
        }
        Vec centroid=volume_centroid(r.mesh.vertices[a],r.mesh.vertices[b],x);
        if(valid && (!corrector.member(centroid) || !corrector.surface(centroid,o.spatial_tolerance)))valid=false;
        for(const auto& old:r.mesh.triangles)if(valid && marching_triangles::triangles_conflict(r.mesh.vertices,t,old,eps))valid=false;
        if(!valid){r.mesh.vertices.pop_back();continue;}
        r.info.push_back({-1,-1,fa|fb});r.mesh.add(t);
        r.box_edges.insert(edge_key(a,id));r.box_edges.insert(edge_key(b,id));
        r.unresolved.erase(std::remove_if(r.unresolved.begin(),r.unresolved.end(),[&](const Event& issue){
            return issue.reason=="Unfinished front: no valid growth or existing-vertex connection" && (issue.position-mid).norm()<10*eps;
        }),r.unresolved.end());
        ++count;
    }
    return count;
}
Range project(const Vec& x,int face) {
    const int d=face/2,u=(d+1)%3,v=(d+2)%3;
    return Range(x[u],(face%2?1.:-1.)*x[v]);
}
double area(const std::vector<int>& loop,const Mesh& mesh,int face) {
    double a=0;
    const Range origin=project(mesh.vertices[loop.front()],face);
    for (size_t i=0;i<loop.size();++i) {
        const Range p=project(mesh.vertices[loop[i]],face)-origin,q=project(mesh.vertices[loop[(i+1)%loop.size()]],face)-origin;
        a+=p[0]*q[1]-p[1]*q[0];
    }
    return a/2;
}
bool in_ring(const Range& p,const std::vector<int>& ring,const Mesh& mesh,int face) {
    bool inside=false;
    for (size_t i=0,j=ring.size()-1;i<ring.size();j=i++) {
        const Range a=project(mesh.vertices[ring[i]],face)-p,b=project(mesh.vertices[ring[j]],face)-p;
        if ((a[1]>0)!=(b[1]>0) && 0<(b[0]-a[0])*(-a[1])/(b[1]-a[1])+a[0]) inside=!inside;
    }
    return inside;
}
void cap(VolumeBoundary& b,const Field& f,const RangeTriangle& range,const Options& o,int face) {
    const double eps=o.step*1e-7;
    const unsigned mask=1u<<face;
    std::set<Edge> graph;
    // Opposite halfedge winding closes the already oriented side patches.
    for (const auto& entry:b.mesh.boundary) {
        const auto e=entry.first;
        if ((f.box.faces(b.mesh.vertices[e.first],eps)&f.box.faces(b.mesh.vertices[e.second],eps)&mask)!=0)
            graph.insert({e.second,e.first});
    }
    const int d=face/2,u=(d+1)%3,v=(d+2)%3;
    std::array<Vec,4> corners;
    for (int k=0;k<4;++k) {
        Vec x=f.box.min; x[d]=face%2?f.box.max[d]:f.box.min[d];
        x[u]=(k==1||k==2)?f.box.max[u]:f.box.min[u]; x[v]=k>=2?f.box.max[v]:f.box.min[v]; corners[k]=x;
    }
    if (!(face%2)) std::reverse(corners.begin(),corners.end());
    for (int k=0;k<4;++k) {
        const Vec a=corners[k],z=corners[(k+1)%4],delta=z-a;
        const int ia=vertex(b,a,eps),iz=vertex(b,z,eps);
        std::vector<std::pair<double,int>> line{{0,ia},{1,iz}};
        for (size_t i=0;i<b.mesh.vertices.size();++i) {
            if (int(i)==ia || int(i)==iz) continue;
            const Vec& x=b.mesh.vertices[i]; const double t=(x-a).dot(delta)/delta.squaredNorm();
            if (t>0 && t<1 && (x-(a+t*delta)).norm()<eps) line.emplace_back(t,int(i));
        }
        std::sort(line.begin(),line.end());
        for (size_t i=1;i<line.size();++i) {
            const int first=line[i-1].second,last=line[i].second;
            const Vec mid=volume_midpoint(b.mesh.vertices[first],b.mesh.vertices[last]);
            if (range.contains(f.value(mid))) graph.insert({first,last});
        }
    }
    // A separate constrained search detects face-contact components missed by
    // side growth, rather than interpreting a missing graph as an empty cap.
    const auto shape=discovery_shape(f,o);
    for (const auto& s:range.sides) {
        Corrector c(f,s,o);
        for (int i=0;i<shape[u];++i) for (int j=0;j<shape[v];++j) {
            Vec p=f.box.min; p[d]=face%2?f.box.max[d]:f.box.min[d];
            p[u]+=double(i)/(shape[u]-1)*(f.box.max[u]-f.box.min[u]);
            p[v]+=double(j)/(shape[v]-1)*(f.box.max[v]-f.box.min[v]);
            if (!c.surface(p,(f.box.max-f.box.min).norm(),mask) || !range.contains(f.value(p),o.residual_tolerance)) continue;
            double distance=std::numeric_limits<double>::infinity();
            for (auto e:graph) distance=std::min(distance,(p-marching_triangles::closest_on_segment(p,b.mesh.vertices[e.first],b.mesh.vertices[e.second])).norm());
            if (distance>2*o.spatial_tolerance) { event(b,p,"Independent box-face search found an unassembled contact"); return; }
        }
    }
    std::map<int,int> next,previous;
    for (auto e:graph) if (!next.emplace(e.first,e.second).second || !previous.emplace(e.second,e.first).second) {
        event(b,b.mesh.vertices[e.first],"Ambiguous box-face boundary graph"); return;
    }
    if (next.size()!=previous.size()) { event(b,corners[0],"Open box-face boundary graph"); return; }
    std::vector<std::vector<int>> loops;
    while (!graph.empty()) {
        std::vector<int> loop; int first=graph.begin()->first,current=first;
        do {
            auto it=next.find(current);
            if (it==next.end() || !graph.erase({current,it->second})) { event(b,b.mesh.vertices[current],"Open box-face boundary loop"); return; }
            loop.push_back(current); current=it->second;
        } while (current!=first);
        if (loop.size()<3 || std::abs(area(loop,b.mesh,face))<eps*eps) { event(b,b.mesh.vertices[first],"Degenerate box-face region"); return; }
        loops.push_back(std::move(loop));
    }
    std::vector<std::vector<int>> children(loops.size());
    for (size_t i=0;i<loops.size();++i) if (area(loops[i],b.mesh,face)<0) {
        int parent=-1; double smallest=std::numeric_limits<double>::infinity();
        const Range p=project(b.mesh.vertices[loops[i][0]],face);
        for (size_t j=0;j<loops.size();++j) {
            const double a=area(loops[j],b.mesh,face);
            if (a>0 && a<smallest && in_ring(p,loops[j],b.mesh,face)) {parent=int(j);smallest=a;}
        }
        if (parent<0) { event(b,b.mesh.vertices[loops[i][0]],"Cap hole has no enclosing selected region"); return; }
        children[parent].push_back(int(i));
    }
    for (size_t i=0;i<loops.size();++i) if (area(loops[i],b.mesh,face)>0) {
        std::vector<std::vector<std::array<double,2>>> polygon;
        std::vector<int> ids,constrained;
        const Range origin=project(b.mesh.vertices[loops[i].front()],face);
        auto ring=[&](int j) {
            auto clean=loops[j];
            constrained.insert(constrained.end(),clean.begin(),clean.end());
            bool changed=true;
            while(changed && clean.size()>3) {
                changed=false;
                for(size_t k=0;k<clean.size();++k) {
                    const Vec& p=b.mesh.vertices[clean[k]];
                    const Vec& a=b.mesh.vertices[clean[(k+clean.size()-1)%clean.size()]];
                    const Vec& z=b.mesh.vertices[clean[(k+1)%clean.size()]];
                    if((p-marching_triangles::closest_on_segment(p,a,z)).norm()<eps) {
                        clean.erase(clean.begin()+k);changed=true;break;
                    }
                }
            }
            polygon.emplace_back();
            for (int id:clean) { const Range p=project(b.mesh.vertices[id],face)-origin; polygon.back().push_back({p[0],p[1]}); ids.push_back(id); }
        };
        ring(int(i)); for (int j:children[i]) ring(j);
        const auto indices=mapbox::earcut<uint32_t>(polygon);
        std::vector<Triangle> triangles;
        for (size_t k=0;k<indices.size();k+=3) triangles.push_back({ids[indices[k]],ids[indices[k+1]],ids[indices[k+2]]});
        // Earcut may omit collinear boundary samples. Reinsert each one on its
        // containing edge so cap/side and cap/cap seams retain identical IDs.
        for (int id:constrained) {
            bool used=false; for (const auto& t:triangles) if (std::find(t.begin(),t.end(),id)!=t.end()) used=true;
            if (used) continue;
            bool split=false;
            const size_t n=triangles.size();
            for (size_t k=0;k<n;++k) for (int j=0;j<3;++j) {
                const Triangle t=triangles[k]; const int a=t[j],z=t[(j+1)%3],c=t[(j+2)%3];
                if ((b.mesh.vertices[id]-marching_triangles::closest_on_segment(b.mesh.vertices[id],b.mesh.vertices[a],b.mesh.vertices[z])).norm()<eps) {
                    triangles[k]={a,id,c}; triangles.push_back({id,z,c}); split=true; break;
                }
            }
            if (!split) { event(b,b.mesh.vertices[id],"Cap triangulation omitted a constrained vertex"); return; }
        }
        for (auto t:triangles) {
            const Vec normal=(b.mesh.vertices[t[1]]-b.mesh.vertices[t[0]]).cross(b.mesh.vertices[t[2]]-b.mesh.vertices[t[0]]);
            if (normal[d]*(face%2?1:-1)<0) std::swap(t[1],t[2]);
            facet(b,t,4+face);
        }
    }
}
void compact(VolumeBoundary& b) {
    Mesh mesh; std::vector<int> ids(b.mesh.vertices.size(),-1);
    for (auto t:b.mesh.triangles) {
        for (int& id:t) { const int old=id; if (ids[old]<0) {ids[old]=int(mesh.vertices.size());mesh.vertices.push_back(b.mesh.vertices[old]);} id=ids[old]; }
        mesh.add(t);
    }
    b.mesh=std::move(mesh);
}
} // namespace
VolumeBoundary assemble_volume_boundary(const Field& f,const RangeTriangle& range,const Options& requested) {
    Options o=requested;o.surface_solve_tolerance=std::min(1e-10,o.residual_tolerance);
    VolumeBoundary b;
    std::array<Seeds,3> seeds;
    std::array<Result,3> corners;
    std::array<std::vector<int>,3> corner_ids;
    for (int i=0;i<3;++i) {
        seeds[i]=find_seeds(f,range.sides[i],o);
        Seeds only; only.endpoints[0]=seeds[i].endpoints[0];
        corners[i]=trace_endpoint_fibers(f,range.sides[i],o,only);
        b.corner_components[i]=corners[i].fibers.size();
        b.unresolved.insert(b.unresolved.end(),corners[i].unresolved.begin(),corners[i].unresolved.end());
        for (const Vec& x:corners[i].mesh.vertices) { corner_ids[i].push_back(int(b.mesh.vertices.size())); b.mesh.vertices.push_back(x); }
    }
    for (int i=0;i<3;++i) {
        Result shared; std::vector<int> ids;
        for (int endpoint=0;endpoint<2;++endpoint) {
            const int k=(i+endpoint)%3,offset=int(shared.mesh.vertices.size()),component=int(shared.fibers.size());
            for (size_t j=0;j<corners[k].mesh.vertices.size();++j) {
                shared.mesh.vertices.push_back(corners[k].mesh.vertices[j]);
                auto info=corners[k].info[j]; info.endpoint=endpoint; info.fiber+=component; shared.info.push_back(info);
                ids.push_back(corner_ids[k][j]);
            }
            for (auto fiber:corners[k].fibers) {fiber.endpoint=endpoint;for(int& v:fiber.vertices)v+=offset;shared.fibers.push_back(fiber);}
            for (auto entry:corners[k].fiber_edges) shared.fiber_edges[{entry.first.first+offset,entry.first.second+offset}]=entry.second+component;
        }
        Extractor extractor(f,range.sides[i],o);
        auto state=extractor.finish(extractor.march(seeds[i],std::move(shared)));
        auto& side=state.result;
        b.box_contact_triangles+=finish_box_contacts(side,f,range.sides[i],o);
        for(const auto& issue:side.unresolved) {
            // An inactive range edge is normal for a filled selection (even a
            // wholly selected box). Coverage is checked independently below.
            if(side.mesh.triangles.empty() && seeds[i].surface.empty() &&
               issue.reason=="No regular surface discovered; finite search does not prove emptiness") continue;
            b.unresolved.push_back(issue);
        }
        for (size_t j=ids.size();j<side.mesh.vertices.size();++j) {
            const Vec& x=side.mesh.vertices[j];
            if (side.info[j].faces) ids.push_back(vertex(b,x,o.step*1e-7));
            else { ids.push_back(int(b.mesh.vertices.size())); b.mesh.vertices.push_back(x); }
        }
        for (const auto& t:side.mesh.triangles) facet(b,{ids[t[0]],ids[t[2]],ids[t[1]]},i+1);
    }
    for (int face=0;face<6;++face) cap(b,f,range,o,face);
    compact(b);
    return b;
}
} // namespace fiber_surface
