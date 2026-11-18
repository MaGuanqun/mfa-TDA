#include "volume.h"
#include "corrector.h"
#include <tetgen.h>
#include <numeric>

namespace fiber_surface {
VolumeMesh fill_volume(const VolumeBoundary& boundary,const BoundaryAssessment& report) {
    if (!report.valid || !boundary.unresolved.empty()) throw std::invalid_argument("TetGen requires a validated closed boundary");
    tetgenio in,out; in.firstnumber=0;
    in.numberofpoints=int(boundary.mesh.vertices.size()); in.pointlist=new REAL[3*in.numberofpoints];
    for (int i=0;i<in.numberofpoints;++i) for(int d=0;d<3;++d) in.pointlist[3*i+d]=boundary.mesh.vertices[i][d];
    in.numberoffacets=int(boundary.mesh.triangles.size()); in.facetlist=new tetgenio::facet[in.numberoffacets]();
    in.facetmarkerlist=new int[in.numberoffacets];
    for (int i=0;i<in.numberoffacets;++i) {
        auto& facet=in.facetlist[i]; facet.numberofpolygons=1; facet.polygonlist=new tetgenio::polygon[1]();
        facet.polygonlist[0].numberofvertices=3; facet.polygonlist[0].vertexlist=new int[3];
        for(int k=0;k<3;++k)facet.polygonlist[0].vertexlist[k]=boundary.mesh.triangles[i][k];
        // Unique parent IDs let assessment prove coverage of EVERY input facet
        // after Steiner subdivision. Semantic patch labels are kept separately.
        in.facetmarkerlist[i]=i+1;
    }
    in.numberofholes=int(report.holes.size());
    if (in.numberofholes) {in.holelist=new REAL[3*in.numberofholes];for(int i=0;i<in.numberofholes;++i)for(int d=0;d<3;++d)in.holelist[3*i+d]=report.holes[i][d];}
    char switches[]="pV"; // No background mesh, interior samples, -q/-a/-m/-i/-c.
    try {tetrahedralize(switches,&in,&out);} catch(int code) {throw std::runtime_error("TetGen failed with code "+std::to_string(code));}
    VolumeMesh mesh;
    for(int i=0;i<out.numberofpoints;++i)mesh.vertices.emplace_back(out.pointlist[3*i],out.pointlist[3*i+1],out.pointlist[3*i+2]);
    if(out.numberofcorners!=4)throw std::runtime_error("TetGen returned non-linear cells");
    for(int i=0;i<out.numberoftetrahedra;++i) {
        Tetrahedron t;for(int k=0;k<4;++k)t[k]=out.tetrahedronlist[4*i+k];
        const Vec a=mesh.vertices[t[0]],b=mesh.vertices[t[1]],c=mesh.vertices[t[2]],d=mesh.vertices[t[3]];
        if((b-a).dot((c-a).cross(d-a))<0)std::swap(t[0],t[1]);
        mesh.cells.push_back(t);
    }
    for(int i=0;i<out.numberoftrifaces;++i) {
        if(!out.trifacemarkerlist || out.trifacemarkerlist[i]<=0)continue;
        Triangle t{out.trifacelist[3*i],out.trifacelist[3*i+1],out.trifacelist[3*i+2]};
        const int parent=out.trifacemarkerlist[i]-1;
        if(parent>=0 && size_t(parent)<boundary.mesh.triangles.size()) {
            const auto p=boundary.mesh.triangles[parent];
            const Vec normal=(boundary.mesh.vertices[p[1]]-boundary.mesh.vertices[p[0]]).cross(boundary.mesh.vertices[p[2]]-boundary.mesh.vertices[p[0]]);
            if((mesh.vertices[t[1]]-mesh.vertices[t[0]]).cross(mesh.vertices[t[2]]-mesh.vertices[t[0]]).dot(normal)<0)std::swap(t[1],t[2]);
        }
        mesh.boundary_faces.push_back(t);
        mesh.parent_facets.push_back(parent);
    }
    return mesh;
}
namespace {
// Exact-sign orientation supplied by TetGen's adaptive predicates. Its sign is
// opposite the usual (b-a)x(c-a).(d-a) convention used in exported cells.
double orientation(const Vec& a,const Vec& b,const Vec& c,const Vec& d) {
    double aa[3]{a[0],a[1],a[2]},bb[3]{b[0],b[1],b[2]},cc[3]{c[0],c[1],c[2]},dd[3]{d[0],d[1],d[2]};
    return -orient3d(aa,bb,cc,dd);
}
bool overlapping(const VolumeMesh& m,const Tetrahedron& a,const Tetrahedron& b,double eps) {
    auto separated=[&](const Vec& axis) {
        // Axis lengths carry spatial units. An absolute cutoff silently drops
        // every separating axis when a valid mesh uses small coordinate units.
        const double length=axis.stableNorm();
        if(length==0 || !std::isfinite(length))return false;
        const Vec n=axis/length;
        double alo=1e300,ahi=-1e300,blo=1e300,bhi=-1e300;
        const Vec origin=m.vertices[a[0]];
        for(int i:a){double x=n.dot(m.vertices[i]-origin);alo=std::min(alo,x);ahi=std::max(ahi,x);}
        for(int i:b){double x=n.dot(m.vertices[i]-origin);blo=std::min(blo,x);bhi=std::max(bhi,x);}
        return ahi<=blo+eps || bhi<=alo+eps;
    };
    for(const auto& t:{a,b})for(int i=0;i<4;++i)
        if(separated(Vec((m.vertices[t[(i+1)%4]]-m.vertices[t[i]]).cross(m.vertices[t[(i+2)%4]]-m.vertices[t[i]]))))return false;
    for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)for(int k=0;k<4;++k)for(int l=k+1;l<4;++l)
        if(separated(Vec((m.vertices[a[j]]-m.vertices[a[i]]).cross(m.vertices[b[l]]-m.vertices[b[k]]))))return false;
    return true;
}
}
VolumeAssessment assess_volume(const VolumeMesh& m,const VolumeBoundary& boundary,const BoundaryAssessment& br,
                               const Field& f,const RangeTriangle& range,const Options& o) {
    VolumeAssessment r; const double eps=o.step*1e-8;
    r.boundary_volume=std::accumulate(br.signed_volumes.begin(),br.signed_volumes.end(),0.);
    std::set<Tetrahedron> unique;
    struct Use {int cell,opposite;};
    std::map<Triangle,std::vector<Use>> faces;
    struct Bounds {Vec lo,hi;int cell;}; std::vector<Bounds> bounds;
    for(size_t i=0;i<m.cells.size();++i) {
        const auto t=m.cells[i];bool valid=true;
        for(int id:t)if(id<0 || size_t(id)>=m.vertices.size() || !m.vertices[id].allFinite())valid=false;
        auto key=t;std::sort(key.begin(),key.end());
        if(!valid || !unique.insert(key).second || std::adjacent_find(key.begin(),key.end())!=key.end()) {++r.invalid_cells;continue;}
        const double volume=orientation(m.vertices[t[0]],m.vertices[t[1]],m.vertices[t[2]],m.vertices[t[3]])/6;
        if(!std::isfinite(volume) || volume<=eps*eps*eps || !std::isfinite(r.volume+volume)) {++r.invalid_cells;continue;}
        r.volume+=volume;r.min_volume=r.min_volume==0?volume:std::min(r.min_volume,volume);r.max_volume=std::max(r.max_volume,volume);
        Bounds box{m.vertices[t[0]],m.vertices[t[0]],int(i)};
        for(int id:t){box.lo=box.lo.cwiseMin(m.vertices[id]);box.hi=box.hi.cwiseMax(m.vertices[id]);} bounds.push_back(box);
        for(int k=0;k<4;++k)faces[marching_triangles::make_triangle_key(t[(k+1)%4],t[(k+2)%4],t[(k+3)%4])].push_back({int(i),t[k]});
    }
    std::map<Triangle,int> parents;
    if(m.boundary_faces.size()!=m.parent_facets.size()){++r.conformity_failures;return r;}
    for(size_t i=0;i<m.boundary_faces.size();++i) {
        auto t=m.boundary_faces[i];std::sort(t.begin(),t.end());
        if(!parents.emplace(t,m.parent_facets[i]).second)++r.conformity_failures;
    }
    std::vector<double> covered(boundary.mesh.triangles.size(),0);
    Mesh output_boundary;
    // Keep a compact vertex array for independent surface topology checks.
    std::map<int,int> boundary_ids;
    for(const auto& entry:faces) {
        auto t=entry.first; const auto& use=entry.second;
        const Vec a=m.vertices[t[0]],b=m.vertices[t[1]],c=m.vertices[t[2]];
        if(use.size()==2) {
            const double first=orientation(a,b,c,m.vertices[use[0].opposite]);
            const double second=orientation(a,b,c,m.vertices[use[1].opposite]);
            // Compare signs directly: multiplying two valid small determinants
            // can underflow to zero and falsely mark an interior face invalid.
            if(parents.count(t) || !std::isfinite(first) || !std::isfinite(second) ||
               first==0 || second==0 || std::signbit(first)==std::signbit(second))++r.invalid_faces;
        } else if(use.size()==1) {
            const auto it=parents.find(t);
            if(it==parents.end() || it->second<0 || size_t(it->second)>=covered.size()) {++r.conformity_failures;continue;}
            const int parent=it->second;const auto input=boundary.mesh.triangles[parent];
            const Vec pa=boundary.mesh.vertices[input[0]],pb=boundary.mesh.vertices[input[1]],pc=boundary.mesh.vertices[input[2]];
            const Vec normal=(pb-pa).cross(pc-pa).normalized();
            for(const Vec& p:{a,b,c}) if((p-marching_triangles::closest_on_triangle(p,pa,pb,pc)).norm()>eps)++r.conformity_failures;
            if(orientation(pa,pb,pc,m.vertices[use[0].opposite])>=0)++r.invalid_faces;
            covered[parent]+=.5*(b-a).cross(c-a).norm();
            if((b-a).cross(c-a).dot(normal)<0)std::swap(t[1],t[2]);
            for(int& id:t) {auto inserted=boundary_ids.emplace(id,int(boundary_ids.size()));if(inserted.second)output_boundary.vertices.push_back(m.vertices[id]);id=inserted.first->second;}
            output_boundary.add(t);
            const int patch=boundary.patches[parent];
            for(const Vec& p:volume_triangle_samples(a,b,c)) {
                if(!f.box.contains(p)){++r.conformity_failures;continue;}
                const Range value=f.value(p);
                if(!value.allFinite()){++r.conformity_failures;continue;}
                double violation=0;for(const auto& s:range.sides)violation=std::max(violation,-s.phi(value));
                r.inequality_violation=std::max(r.inequality_violation,violation);
                Vec corrected=p;
                if(!correct_volume_sample(corrected,patch,f,range,o))++r.conformity_failures;
                else r.spatial_error=std::max(r.spatial_error,(corrected-p).norm());
                if(patch<=3) {
                    const auto& s=range.sides[patch-1];r.side_residual=std::max(r.side_residual,std::abs(s.phi(value)));
                }
            }
        } else ++r.invalid_faces;
    }
    for(const auto& entry:parents)if(!faces.count(entry.first))++r.conformity_failures;
    for(size_t i=0;i<covered.size();++i) {
        const auto t=boundary.mesh.triangles[i];
        const double area=.5*(boundary.mesh.vertices[t[1]]-boundary.mesh.vertices[t[0]]).cross(boundary.mesh.vertices[t[2]]-boundary.mesh.vertices[t[0]]).norm();
        if(std::abs(area-covered[i])>std::max(eps*eps,area*1e-8))++r.conformity_failures;
    }
    r.output_topology=marching_triangles::validate_mesh(output_boundary.vertices,output_boundary.triangles,eps);
    if(!r.output_topology.closed())++r.conformity_failures;
    // Broad phase followed by tetrahedron separating-axis tests, including
    // edge/edge axes. Face incidences alone do not rule out interior overlaps.
    std::sort(bounds.begin(),bounds.end(),[](const Bounds& a,const Bounds& b){return a.lo[0]<b.lo[0];});
    for(size_t i=0;i<bounds.size();++i)for(size_t j=i+1;j<bounds.size() && bounds[j].lo[0]<bounds[i].hi[0]-eps;++j) {
        if((bounds[i].lo.array()>=bounds[j].hi.array()-eps).any() || (bounds[j].lo.array()>=bounds[i].hi.array()-eps).any())continue;
        if(overlapping(m,m.cells[bounds[i].cell],m.cells[bounds[j].cell],eps))++r.intersections;
    }
    if(m.cells.empty()){r.min_volume=0;++r.invalid_cells;}
    if(std::abs(r.volume-r.boundary_volume)>std::max(eps*eps*eps,1e-8*std::abs(r.boundary_volume)))++r.conformity_failures;
    r.valid=br.valid && !r.invalid_cells && !r.invalid_faces && !r.conformity_failures && !r.intersections;
    return r;
}
} // namespace fiber_surface
