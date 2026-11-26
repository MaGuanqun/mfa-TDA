#include "volume.h"
#include "extract.h"
#include <mapbox/earcut.hpp>
#include <tbb/parallel_for.h>
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
// Improve the planar triangulation while keeping every boundary edge fixed.
// Near-collinear rim samples can otherwise form a thin or inverted ear when
// restored to an edge of the simplified polygon.
void improve_cap(std::vector<Triangle>& cells,const Mesh& mesh,int face) {
    std::map<Edge,std::set<size_t>> incident;
    std::set<Edge> pending;
    auto edges_of=[](const Triangle& t) {
        return std::array<Edge,3>{edge_key(t[0],t[1]),edge_key(t[1],t[2]),edge_key(t[2],t[0])};
    };
    for (size_t i=0;i<cells.size();++i) for (auto e:edges_of(cells[i])) incident[e].insert(i);
    for (const auto& [e,faces]:incident) if (faces.size()==2) pending.insert(e);
    auto quality=[&](const Triangle& t) {
        const Range a=project(mesh.vertices[t[1]]-mesh.vertices[t[0]],face);
        const Range b=project(mesh.vertices[t[2]]-mesh.vertices[t[0]],face);
        return (a.x()*b.y()-a.y()*b.x())/(a.squaredNorm()+b.squaredNorm()+(a-b).squaredNorm());
    };
    while (!pending.empty()) {
        const auto e=*pending.begin(); pending.erase(pending.begin());
        const auto it=incident.find(e);
        if (it==incident.end() || it->second.size()!=2) continue;
        const size_t i=*it->second.begin(),j=*it->second.rbegin();
        const auto first=cells[i],second=cells[j];
        int k=0; while (edge_key(first[k],first[(k+1)%3])!=e) ++k;
        const int a=first[k],z=first[(k+1)%3],c=first[(k+2)%3];
        int d=second[0]; for (int id:second) if (id!=a && id!=z) d=id;
        if (c==d || incident.count(edge_key(c,d))) continue;
        const Triangle t{c,d,z},s{d,c,a};
        const double before=std::min(quality(first),quality(second));
        const double after=std::min(quality(t),quality(s));
        if (!(after>0 && after>before+1e-13)) continue;
        for (size_t id:{i,j}) for (auto old:edges_of(cells[id])) {
            auto& list=incident[old]; list.erase(id);
            if (list.empty()) incident.erase(old);
        }
        cells[i]=t; cells[j]=s;
        for (size_t id:{i,j}) for (auto added:edges_of(cells[id])) {
            incident[added].insert(id); pending.insert(added);
        }
    }
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
    const double radius=2*o.spatial_tolerance;
    // Most probes converge onto an assembled outline. Test only segments
    // whose expanded bounds contain the projected point, stopping at a hit.
    struct OutlineSegment {Vec a,z; double umin,umax,vmin,vmax;};
    std::vector<OutlineSegment> outline;
    for (const auto& e:graph) {
        const Vec a=b.mesh.vertices[e.first],z=b.mesh.vertices[e.second];
        outline.push_back({a,z,std::min(a[u],z[u])-radius,std::max(a[u],z[u])+radius,
                               std::min(a[v],z[v])-radius,std::max(a[v],z[v])+radius});
    }
    for (const auto& s:range.sides) {
        Corrector c(f,s,o);
        std::vector<std::optional<Vec>> missed(shape[u]);
        tbb::parallel_for(0,shape[u],[&](int i) {
          for (int j=0;j<shape[v];++j) {
            Vec p=f.box.min; p[d]=face%2?f.box.max[d]:f.box.min[d];
            p[u]+=double(i)/(shape[u]-1)*(f.box.max[u]-f.box.min[u]);
            p[v]+=double(j)/(shape[v]-1)*(f.box.max[v]-f.box.min[v]);
            if (!c.surface(p,(f.box.max-f.box.min).norm(),mask) || !range.contains(f.value(p),o.residual_tolerance)) continue;
            const bool covered=std::any_of(outline.begin(),outline.end(),[&](const OutlineSegment& edge) {
                if (p[u]<edge.umin || p[u]>edge.umax || p[v]<edge.vmin || p[v]>edge.vmax) return false;
                return (p-marching_triangles::closest_on_segment(p,edge.a,edge.z)).squaredNorm()<=radius*radius;
            });
            if (!covered) { missed[i]=p; break; }
          }
        });
        for (const auto& point:missed) if (point) {
            event(b,*point,"Independent box-face search found an unassembled contact"); return;
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
    std::vector<Triangle> cap_faces;
    for (size_t i=0;i<loops.size();++i) if (area(loops[i],b.mesh,face)>0) {
        std::vector<std::vector<std::array<double,2>>> polygon;
        std::vector<int> ids;
        std::vector<int> rings{int(i)};
        rings.insert(rings.end(),children[i].begin(),children[i].end());
        const Range origin=project(b.mesh.vertices[loops[i].front()],face);
        auto ring=[&](int j) {
            auto clean=loops[j];
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
        for (int j:rings) ring(j);
        const auto indices=mapbox::earcut<uint32_t>(polygon);
        std::vector<Triangle> triangles;
        for (size_t k=0;k<indices.size();k+=3) triangles.push_back({ids[indices[k]],ids[indices[k+1]],ids[indices[k+2]]});
        // Recover omitted samples on their own ring edge. A distance test
        // alone can match several nearly collinear interior diagonals.
        std::set<int> used;
        for (const auto& t:triangles) used.insert(t.begin(),t.end());
        for (int j:rings) {
            const auto& loop=loops[j];
            auto anchor=std::find_if(loop.begin(),loop.end(),[&](int id){return used.count(id);});
            if (anchor==loop.end()) {event(b,b.mesh.vertices[loop[0]],"Cap triangulation omitted a ring");return;}
            size_t start=size_t(anchor-loop.begin());
            for (size_t walk=1;walk<=loop.size();++walk) {
                const size_t stop=(size_t(anchor-loop.begin())+walk)%loop.size();
                if (!used.count(loop[stop])) continue;
                if ((start+1)%loop.size()!=stop) {
                    const int first=loop[start],last=loop[stop];
                    size_t face_id=0;
                    int edge=-1;
                    for (;face_id<triangles.size() && edge<0;++face_id) {
                        const auto& t=triangles[face_id];
                        for (int e=0;e<3;++e) if (edge_key(t[e],t[(e+1)%3])==edge_key(first,last)) {edge=e;break;}
                    }
                    if (edge<0) {event(b,b.mesh.vertices[first],"Cap triangulation omitted a constrained edge");return;}
                    --face_id;
                    const Triangle t=triangles[face_id];
                    const int third=t[(edge+2)%3];
                    int prev=first;
                    for (size_t k=(start+1)%loop.size();;k=(k+1)%loop.size()) {
                        Triangle split{prev,loop[k],third};
                        if (t[edge]!=first) std::swap(split[0],split[1]);
                        if (prev==first) triangles[face_id]=split;
                        else triangles.push_back(split);
                        if (k==stop) break;
                        prev=loop[k];
                    }
                }
                start=stop;
            }
        }
        improve_cap(triangles,b.mesh,face);
        for (auto t:triangles) {
            const Vec normal=(b.mesh.vertices[t[1]]-b.mesh.vertices[t[0]]).cross(b.mesh.vertices[t[2]]-b.mesh.vertices[t[0]]);
            if (normal[d]*(face%2?1:-1)<0) std::swap(t[1],t[2]);
            cap_faces.push_back(t);
        }
    }
    // Commit a complete face only after the oriented boundary of its
    // triangulation agrees exactly with the assembled cap and side rims.
    std::set<Edge> expected;
    for (const auto& loop:loops) for (size_t i=0;i<loop.size();++i)
        expected.insert({loop[i],loop[(i+1)%loop.size()]});
    std::map<Edge,int> directed;
    for (const auto& t:cap_faces) {
        const Vec a=b.mesh.vertices[t[1]]-b.mesh.vertices[t[0]],z=b.mesh.vertices[t[2]]-b.mesh.vertices[t[0]];
        if (a.cross(z).norm()<=o.step*1e-8*std::max(a.norm(),z.norm())) {
            event(b,b.mesh.vertices[t[0]],"Degenerate cap triangle");return;
        }
        for (int k=0;k<3;++k) {
            const Edge e{t[k],t[(k+1)%3]};
            if (++directed[e]>1) {event(b,b.mesh.vertices[e.first],"Repeated oriented cap edge");return;}
        }
    }
    for (const auto& [e,count]:directed) if (!directed.count({e.second,e.first}) && !expected.erase(e)) {
        event(b,b.mesh.vertices[e.first],"Unexpected exposed cap edge");return;
    }
    if (!expected.empty()) {event(b,b.mesh.vertices[expected.begin()->first],"Cap omitted a rim edge");return;}
    for (const auto& t:cap_faces) facet(b,t,4+face);
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
    const RangeCurve loop({range.corners[0],range.corners[1],range.corners[2],range.corners[0]});
    const Options o=curve_options(loop,requested);
    return assemble_volume_boundary(f,range,o,extract_curve(f,loop,o));
}
VolumeBoundary assemble_volume_boundary(const Field& f,const RangeTriangle& range,const Options& o,Result sides) {
    VolumeBoundary b;
    b.mesh.vertices=std::move(sides.mesh.vertices);
    b.unresolved=std::move(sides.unresolved);
    b.box_contact_triangles=sides.box_contact_triangles;
    for (const auto& fiber:sides.fibers) ++b.corner_components[fiber.endpoint];
    for (size_t i=0; i<sides.mesh.triangles.size(); ++i) {
        auto t=sides.mesh.triangles[i];
        std::swap(t[1],t[2]); // The CCW curve normals point into the selected region.
        facet(b,t,sides.face_segments[i]+1);
    }
    for (int face=0;face<6;++face) cap(b,f,range,o,face);
    compact(b);
    return b;
}
} // namespace fiber_surface
