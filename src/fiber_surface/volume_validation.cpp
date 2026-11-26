#include "volume.h"
#include "corrector.h"
#include <numeric>

namespace fiber_surface {
bool correct_volume_sample(Vec& x,int patch,const Field& f,const RangeTriangle& range,const Options& o) {
    if(patch<1 || patch>9 || !f.box.contains(x))return false;
    const Vec origin=x;
    unsigned faces=f.box.faces(x,o.step*1e-9);
    if(patch>=4)faces|=1u<<(patch-4);
    f.box.fix(x,faces);
    std::array<bool,3> active{};if(patch<=3)active[patch-1]=true;
    const double tolerance=o.patch_solve_tolerance();
    for(int it=0;it<o.iterations;++it) {
        if(!f.box.contains(x) || (x-origin).norm()>o.spatial_tolerance)return false;
        const Range value=f.value(x);const Jacobian jac=f.jacobian(x);
        if(!value.allFinite() || !jac.allFinite())return false;
        bool converged=true;
        for(int k=0;k<3;++k) {
            const double phi=range.sides[k].phi(value);
            if(phi < -tolerance)active[k]=true;
            if(active[k] && std::abs(phi)>tolerance)converged=false;
        }
        if(converged)return true;
        std::vector<int> free;for(int d=0;d<3;++d)if(!(faces&(3u<<(2*d))))free.push_back(d);
        int rows=int(std::count(active.begin(),active.end(),true));
        if(free.empty() || rows>int(free.size()))return false;
        Eigen::MatrixXd a(rows,free.size());Eigen::VectorXd residual(rows);int row=0;
        for(int k=0;k<3;++k)if(active[k]) {
            const Vec g=range.sides[k].gradient(jac);
            for(size_t d=0;d<free.size();++d)a(row,d)=g[free[d]];
            residual[row++]=-range.sides[k].phi(value);
        }
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> solve(a);solve.setThreshold(o.derivative_tolerance);
        if(solve.rank()!=rows)return false;
        const Eigen::VectorXd delta=solve.solve(residual);
        if(!delta.allFinite())return false;
        Vec step=Vec::Zero();
        for(size_t d=0;d<free.size();++d)step[free[d]]=delta[d];
        const Vec proposed=x+step;
        if(f.box.contains(proposed)) {
            x=proposed;
        } else {
            // Follow the Newton direction to its first box contact. Recompute
            // the residual and solve in that plane on the next iteration.
            double fraction=1.;
            unsigned contact=0;
            for(int d:free) {
                const bool lower=proposed[d]<f.box.min[d];
                if(!lower && proposed[d]<=f.box.max[d])continue;
                const double wall=lower?f.box.min[d]:f.box.max[d];
                const double t=(wall-x[d])/step[d];
                const unsigned face=1u<<(2*d+(lower?0:1));
                if(t<fraction) {fraction=t;contact=face;}
                else if(t==fraction)contact|=face;
            }
            if(!contact || fraction<0)return false;
            x=(x+fraction*step).cwiseMax(f.box.min).cwiseMin(f.box.max);
            faces|=contact;
        }
        f.box.fix(x,faces);
    }
    return false;
}
bool inside_shell(const Vec& p,const Mesh& mesh,const std::vector<int>& faces) {
    double angle=0;
    for (int id:faces) {
        const auto t=mesh.triangles[id];
        const Vec a=mesh.vertices[t[0]]-p,b=mesh.vertices[t[1]]-p,c=mesh.vertices[t[2]]-p;
        angle+=2*std::atan2(a.dot(b.cross(c)),a.norm()*b.norm()*c.norm()+a.dot(b)*c.norm()+b.dot(c)*a.norm()+c.dot(a)*b.norm());
    }
    return std::abs(angle)>2*std::acos(-1.);
}
namespace {
double clearance(const Vec& p,const Mesh& m) {
    double d=std::numeric_limits<double>::infinity();
    for (auto t:m.triangles) d=std::min(d,(p-marching_triangles::closest_on_triangle(p,m.vertices[t[0]],m.vertices[t[1]],m.vertices[t[2]])).norm());
    return d;
}
double ray_hit(const Vec& p,const Vec& n,const Mesh& m,int skip,double eps) {
    double best=std::numeric_limits<double>::infinity();
    for (size_t i=0;i<m.triangles.size();++i) if (int(i)!=skip) {
        const auto t=m.triangles[i]; const Vec a=m.vertices[t[0]],u=m.vertices[t[1]]-a,v=m.vertices[t[2]]-a;
        const Vec h=n.cross(v); const double det=u.dot(h);
        if (std::abs(det)<1e-12*u.norm()*v.norm()) continue;
        const Vec s=p-a,q=s.cross(u); const double x=s.dot(h)/det,y=n.dot(q)/det,d=v.dot(q)/det;
        if (x>=-1e-10 && y>=-1e-10 && x+y<=1+1e-10 && d>eps) best=std::min(best,d);
    }
    return best;
}
// Along an x row, shell membership can change only on a segment touching a
// boundary triangle. Mark a conservative superset using triangle bounding boxes.
// A marked segment (including either endpoint) still uses the solid-angle test.
std::vector<char> possible_crossings(const Mesh& m,const std::vector<Vec>& grid,
                                     const SeedShape& shape,double eps) {
    std::array<std::vector<double>,3> axes;
    size_t stride=1;
    for (int d=0;d<3;++d) {
        for (int i=0;i<shape[d];++i) axes[d].push_back(grid[stride*i][d]);
        stride*=shape[d];
    }
    std::vector<char> marked(grid.size(),0);
    for (const auto& t:m.triangles) {
        const Vec lo=m.vertices[t[0]].cwiseMin(m.vertices[t[1]]).cwiseMin(m.vertices[t[2]]);
        const Vec hi=m.vertices[t[0]].cwiseMax(m.vertices[t[1]]).cwiseMax(m.vertices[t[2]]);
        SeedShape first,last;
        for (int d=0;d<3;++d) {
            const auto& a=axes[d];
            first[d]=int(std::lower_bound(a.begin(),a.end(),lo[d]-eps)-a.begin());
            last[d]=int(std::upper_bound(a.begin(),a.end(),hi[d]+eps)-a.begin());
        }
        // Segment i runs from node i-1 to node i. Its left endpoint may
        // touch the upper bound even if its right endpoint lies beyond it.
        last[0]=std::min(last[0]+1,shape[0]);
        for (int k=first[2];k<last[2];++k) for (int j=first[1];j<last[1];++j) {
            const size_t row=(size_t(k)*shape[1]+j)*shape[0];
            std::fill(marked.begin()+row+first[0],marked.begin()+row+last[0],1);
        }
    }
    return marked;
}
}
BoundaryAssessment validate_volume_boundary(VolumeBoundary& b,const Field& f,const RangeTriangle& range,const Options& o) {
    BoundaryAssessment r;
    const double eps=o.step*1e-8;
    r.topology=marching_triangles::validate_mesh(b.mesh.vertices,b.mesh.triangles,eps);
    if (b.mesh.triangles.empty() || !r.topology.closed() || b.patches.size()!=b.mesh.triangles.size() ||
        std::any_of(b.patches.begin(),b.patches.end(),[](int patch){return patch<1 || patch>9;})) {
        b.unresolved.push_back({(f.box.min+f.box.max)/2,"Boundary is empty, open, nonmanifold, intersecting, or inconsistently labeled: "+r.topology.summary()}); return r;
    }
    // Reject before invoking field callbacks, which may only be defined on the
    // model domain. Every convex facet sample then remains in the box as well.
    for(const Vec& x:b.mesh.vertices) if(!f.box.contains(x)) {
        ++r.constraint_failures;b.unresolved.push_back({x,"Boundary vertex lies outside the model domain"});return r;
    }
    // Growth caches are not part of the boundary's geometric definition and
    // may be absent or stale after importing or editing its vertex/facet arrays.
    Mesh rebuilt;rebuilt.vertices=b.mesh.vertices;
    for(const auto& t:b.mesh.triangles)rebuilt.add(t);
    b.mesh=std::move(rebuilt);
    const auto& m=b.mesh;
    auto member=[&](const Vec& p) { return f.box.contains(p) && range.contains(f.value(p)); };
    for (size_t i=0;i<m.triangles.size();++i) {
        const auto t=m.triangles[i]; const Vec a=m.vertices[t[0]],z=m.vertices[t[1]],c=m.vertices[t[2]],center=volume_centroid(a,z,c);
        const Vec normal=(z-a).cross(c-a).normalized(); const int patch=b.patches[i];
        if (patch<1 || patch>9) {++r.constraint_failures;continue;}
        const auto samples=volume_triangle_samples(a,z,c);
        for (const Vec& x:samples) {
            if (!f.box.contains(x)) {++r.constraint_failures;continue;}
            const Range value=f.value(x);
            if (!value.allFinite()) {++r.constraint_failures;continue;}
            double violation=0; for (const auto& s:range.sides) violation=std::max(violation,-s.phi(value));
            r.inequality_violation=std::max(r.inequality_violation,violation);
            Vec corrected=x;
            if(!correct_volume_sample(corrected,patch,f,range,o))++r.constraint_failures;
            else r.spatial_error=std::max(r.spatial_error,(x-corrected).norm());
            if (patch<=3) {
                const auto& s=range.sides[patch-1];
                r.side_residual=std::max(r.side_residual,std::abs(s.phi(value)));
            } else {
                const int face=patch-4,d=face/2; const double limit=face%2?f.box.max[d]:f.box.min[d];
                if (std::abs(x[d]-limit)>eps) ++r.constraint_failures;
            }
        }
        Vec expected;
        if (patch<=3) expected=-range.sides[patch-1].gradient(f.jacobian(center));
        else { expected=Vec::Zero();expected[(patch-4)/2]=(patch-4)%2?1:-1; }
        if (!expected.allFinite() || expected.norm()<=o.derivative_tolerance || normal.dot(expected)<=0) ++r.orientation_failures;
    }
    // Shared labels identify corner fibers and active box constraints. Test the
    // regular-volume assumption explicitly, including range/box combinations.
    std::vector<std::set<int>> labels(m.vertices.size());
    for (size_t i=0;i<m.triangles.size();++i) for (int id:m.triangles[i]) labels[id].insert(b.patches[i]);
    for (size_t i=0;i<m.vertices.size();++i) {
        const Vec& x=m.vertices[i]; const auto j=f.jacobian(x); const auto value=f.value(x);
        if(!j.allFinite() || !value.allFinite()) {++r.constraint_failures;continue;}
        std::vector<Vec> gradients;
        for (int label:labels[i]) {
            if (label<=3) gradients.push_back(range.sides[label-1].gradient(j));
            else gradients.push_back(Vec::Unit((label-4)/2));
        }
        for (int k=0;k<3;++k) if (labels[i].count(k+1) && labels[i].count((k+2)%3+1)) {
            const double residual=(value-range.corners[k]).norm();
            r.corner_residual=std::max(r.corner_residual,residual);
            if (residual>o.corner_residual_tolerance) ++r.constraint_failures;
        }
        Eigen::MatrixXd active(3,gradients.size());
        for (size_t k=0;k<gradients.size();++k) active.col(k)=gradients[k].normalized();
        bool regular=false;
        if(active.allFinite()) {
            Eigen::FullPivLU<Eigen::MatrixXd> lu(active);lu.setThreshold(o.derivative_tolerance);
            regular=lu.rank()==int(gradients.size());
        }
        if (!regular) {
            ++r.constraint_failures; b.unresolved.push_back({x,"Unsupported dependent active range/box constraints"});
        }
    }
    for(const auto& entry:m.edges) {
        const auto& use=entry.second;
        const int first=b.patches[use.face],second=b.patches[use.second];
        if(first>3 || second>3 || first==second)continue;
        int corner=-1;
        for(int k=0;k<3;++k)if((first==k+1 && second==(k+2)%3+1) || (second==k+1 && first==(k+2)%3+1))corner=k;
        if(corner<0)continue;
        const Vec p=volume_midpoint(m.vertices[entry.first.first],m.vertices[entry.first.second]);
        r.corner_residual=std::max(r.corner_residual,(f.value(p)-range.corners[corner]).norm());
    }
    // Connected boundary shells, traversing side/cap seams as ordinary edges.
    r.shells.assign(m.triangles.size(),-1); std::vector<std::vector<int>> components;
    for (size_t seed=0;seed<m.triangles.size();++seed) if (r.shells[seed]<0) {
        const int id=int(components.size()); components.emplace_back(); std::vector<int> stack{int(seed)};
        r.shells[seed]=id;
        while (!stack.empty()) {
            const int face=stack.back();stack.pop_back(); components.back().push_back(face);
            const auto t=m.triangles[face];
            for (int k=0;k<3;++k) {
                const auto& e=m.edges.at(edge_key(t[k],t[(k+1)%3])); const int other=e.face==face?e.second:e.face;
                if (other>=0 && r.shells[other]<0) {r.shells[other]=id;stack.push_back(other);}
            }
        }
        double volume=0; const Vec origin=m.vertices[m.triangles[seed][0]];
        for (int face:components.back()) {const auto t=m.triangles[face];volume+=(m.vertices[t[0]]-origin).dot((m.vertices[t[1]]-origin).cross(m.vertices[t[2]]-origin))/6;}
        r.signed_volumes.push_back(volume);
    }
    auto material=[&](const Vec& p) {int winding=0;for(size_t k=0;k<components.size();++k) if(inside_shell(p,m,components[k]))winding+=r.signed_volumes[k]>0?1:-1;return winding>0;};
    for (size_t k=0;k<components.size();++k) {
        const auto t=m.triangles[components[k][0]];
        const Vec center=volume_centroid(m.vertices[t[0]],m.vertices[t[1]],m.vertices[t[2]]);
        // Any point on this shell is strictly inside or outside every other
        // shell because embeddedness has already excluded shell contacts.
        int depth=0; for(size_t j=0;j<components.size();++j) if(j!=k && inside_shell(center,m,components[j])) ++depth;
        if (std::abs(r.signed_volumes[k])<=eps*eps*eps || (r.signed_volumes[k]>0)!=(depth%2==0)) {
            ++r.orientation_failures; b.unresolved.push_back({center,"Shell containment disagrees with retained-material orientation"});
        }
        if (r.signed_volumes[k]>=0) continue;
        double best=0; Vec hole=Vec::Zero();
        for (int face:components[k]) {
            const auto q=m.triangles[face]; const Vec a=m.vertices[q[0]],z=m.vertices[q[1]],c=m.vertices[q[2]],p=volume_centroid(a,z,c),n=(z-a).cross(c-a).normalized();
            const double hit=ray_hit(p,n,m,face,eps);
            if (!std::isfinite(hit)) continue;
            const Vec candidate=p+.5*hit*n;
            if (!inside_shell(candidate,m,components[k]) || material(candidate) || member(candidate)) continue;
            const double distance=clearance(candidate,m);
            if (distance>std::max(best,10*eps)) {best=distance;hole=candidate;}
        }
        if (best>0) r.holes.push_back(hole);
        else b.unresolved.push_back({center,"No verified ray-midpoint seed for cavity"});
    }
    // Independent empirical coverage, not a certification of absent components.
    const auto shape=discovery_shape(f,o);
    const auto grid=discovery_grid(f.box,shape,o.max_seed_starts);
    const auto crossings=possible_crossings(m,grid,shape,eps);
    bool retained=false;
    for (size_t i=0;i<grid.size();++i) {
        if (i%shape[0]==0) retained=false; // start the ray outside the domain
        if (crossings[i]) retained=material(grid[i]);
        if (member(grid[i])!=retained && clearance(grid[i],m)>2*o.spatial_tolerance) ++r.coverage_misses;
    }
    if (r.constraint_failures || r.orientation_failures || r.coverage_misses)
        b.unresolved.push_back({(f.box.min+f.box.max)/2,"Boundary model/orientation/coverage validation failed"});
    r.valid=b.unresolved.empty() && !r.constraint_failures && !r.orientation_failures && !r.coverage_misses;
    return r;
}
} // namespace fiber_surface
