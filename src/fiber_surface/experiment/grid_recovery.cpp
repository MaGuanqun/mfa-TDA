// Experimental side-surface reconstruction on a regular lattice of fitted-field samples.
// Two adjacent tetrahedra refer to a crossing by its grid edge/face and range side/corner.
#include "../pipeline.h"
#include "../volume.h"
#include "../cli.h"
#include "../corrector.h"
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>

using namespace fiber_surface;
namespace {
using Key=std::array<std::int64_t,4>;
struct Node {std::int64_t id; Vec x; Range f;};
struct Refinement {Vec lo,hi; double step;};
double positive_spacing(const std::string& arg) {
    size_t end=0;
    const double value=std::stod(arg,&end);
    if(end!=arg.size() || !std::isfinite(value) || value<=0) throw std::invalid_argument("Invalid grid spacing");
    return value;
}
double determinant(const Range& a,const Range& b) {return a.x()*b.y()-a.y()*b.x();}

class GridSurface {
    const RangeTriangle& query;
    const Corrector& solver;
    VolumeBoundary result;
    std::map<Key,int> crossings;
    const bool fit;
    const bool carrier_fit;
    const bool guard_fit;
    size_t failed_carriers=0;
    std::vector<Range> targets;
    std::vector<Vec> linear_vertices;
    std::pair<int,bool> vertex(const Key& key,const Vec& point,const Range& value) {
        if(result.mesh.vertices.size()>=size_t(std::numeric_limits<int>::max()))
            throw std::runtime_error("Surface exceeds mesh vertex index range");
        const auto [it,added]=crossings.emplace(key,int(result.mesh.vertices.size()));
        if(added) {
            result.mesh.vertices.push_back(point);
            if(fit) targets.push_back(value);
            if(carrier_fit) linear_vertices.push_back(point);
        }
        return {it->second,added};
    }
    // Bracket the side root on the same grid edge that the linear cut used.
    bool fit_edge(Vec& point,const Node& a,const Node& b,const Segment& side,double t) const {
        const Vec edge=b.x-a.x;
        const bool positive=side.phi(a.f)>0;
        double left=0,right=1;
        for(int k=0;k<64;++k) {
            const Vec x=a.x+t*edge;
            const Range value=solver.field.value(x);
            const double residual=side.phi(value);
            if(std::abs(residual)<=1e-11) {
                if(!Segment::contains_parameter(side.parameter(value))) return false;
                point=x;return true;
            }
            if((residual>0)==positive) left=t; else right=t;
            const double derivative=side.gradient(solver.field.jacobian(x)).dot(edge);
            const double next=derivative!=0 ? t-residual/derivative : left;
            t=next>left && next<right ? next : (left+right)/2;
        }
        return false;
    }
    // Holding the grid face fixed gives a two-dimensional corner-root solve.
    bool fit_corner(Vec& point,const Node& a,const Node& b,const Node& c,const Range& target) const {
        const Vec ab=b.x-a.x,ac=c.x-a.x,n=ab.cross(ac);
        Vec trial=point;
        const double limit=std::max({ab.norm(),ac.norm(),(c.x-b.x).norm()});
        const unsigned fixed=solver.field.box.faces(a.x,0)&solver.field.box.faces(b.x,0)&solver.field.box.faces(c.x,0);
        if(!solver.fiber(trial,target,n.normalized(),point,limit,fixed)) return false;
        const Vec v=trial-a.x;
        if(std::abs(n.normalized().dot(v))>1e-10*limit) return false;
        const double u=v.cross(ac).dot(n)/n.squaredNorm(),w=ab.cross(v).dot(n)/n.squaredNorm();
        if(!(u>0 && w>0 && u+w<1)) return false;
        point=trial;return true;
    }
    void restore_conflicts(double eps) {
        auto& points=result.mesh.vertices;
        const auto& faces=result.mesh.triangles;
        size_t total_restored=0;
        for(;;) {
            std::set<int> restore;
            auto mark=[&](const Triangle& face) {
                for(int id:face) if(points[id]!=linear_vertices[id]) restore.insert(id);
            };
            for(const auto& t:faces) {
                const Vec a=points[t[1]]-points[t[0]],b=points[t[2]]-points[t[0]],normal=a.cross(b);
                const Vec original=(linear_vertices[t[1]]-linear_vertices[t[0]]).cross(linear_vertices[t[2]]-linear_vertices[t[0]]);
                if(normal.norm()<=eps*std::max(a.norm(),b.norm()) || normal.dot(original)<=0) {
                    if(points[t[0]]==linear_vertices[t[0]] && points[t[1]]==linear_vertices[t[1]] && points[t[2]]==linear_vertices[t[2]])
                        throw std::runtime_error("Degenerate linear grid face at the audit tolerance");
                    mark(t);
                }
            }
            if(restore.empty()) {
                marching_triangles::visit_triangle_intersections(points,faces,eps,[&](size_t a,size_t b) {
                    mark(faces[a]);mark(faces[b]);
                });
            }
            if(restore.empty()) break;
            // A vertex returns to its linear position only once. Recheck nearby
            // triangles through the next global pass before accepting the edit.
            for(int id:restore) points[id]=linear_vertices[id];
            total_restored+=restore.size();
        }
        std::cout<<"Fitted vertices restored after face checks: "<<total_restored<<std::endl;
    }
    void polygon(std::vector<int>& ids,const Vec& outward,int label) {
        if(ids.empty()) return;
        if(ids.size()<3 || outward.norm()==0)
            throw std::runtime_error("Nontransverse cut in a grid tetrahedron");
        const auto& points=carrier_fit ? linear_vertices : result.mesh.vertices;
        Vec center=Vec::Zero();
        for(int i:ids) center+=points[i];
        center/=double(ids.size());
        const Vec u=(points[ids[0]]-center).normalized(), v=outward.normalized().cross(u);
        auto angle=[&](int i) {const Vec p=points[i]-center; return std::atan2(p.dot(v),p.dot(u));};
        std::sort(ids.begin(),ids.end(),[&](int a,int b){return angle(a)<angle(b);});
        for(size_t j=1;j+1<ids.size();++j) {
            result.mesh.triangles.push_back({ids[0],ids[j],ids[j+1]});
            result.patches.push_back(label+1);
        }
    }
    void tet(const std::array<Node,4>& n,const std::array<int,3>& axes,const Vec& step) {
        std::array<std::vector<int>,3> cuts;
        for(int side=0;side<3;++side) {
            const auto& s=query.sides[side];
            std::array<double,4> phi;
            for(int i=0;i<4;++i) phi[i]=s.phi(n[i].f);
            if(*std::min_element(phi.begin(),phi.end())>0 || *std::max_element(phi.begin(),phi.end())<0) continue;
            for(int i=0;i<4;++i) for(int j=i+1;j<4;++j) {
                if(phi[i]*phi[j]>0) continue;
                if(phi[i]==0 && phi[j]==0) throw std::runtime_error("Range side coincides with a grid edge");
                const double t=phi[i]/(phi[i]-phi[j]);
                const double parameter=s.parameter(n[i].f+t*(n[j].f-n[i].f));
                if(parameter<0 || parameter>1) continue;
                if(t==0 || t==1 || parameter==0 || parameter==1)
                    throw std::runtime_error("Range query meets a grid vertex or a corner meets a grid edge");
                const auto [id,added]=vertex({side,n[i].id,n[j].id,-1},n[i].x+t*(n[j].x-n[i].x),s.a+parameter*s.d);
                if(carrier_fit && added && !fit_edge(result.mesh.vertices[id],n[i],n[j],s,t)) ++failed_carriers;
                cuts[side].push_back(id);
            }
        }
        for(int i=0;i<2;++i) for(int j=i+1;j<3;++j) for(int k=j+1;k<4;++k) {
            const Range a=n[j].f-n[i].f,b=n[k].f-n[i].f;
            const double det=determinant(a,b);
            for(int corner=0;corner<3;++corner) {
                const Range d=query.corners[corner]-n[i].f;
                if(det==0) {
                    if(determinant(a,d)==0 && determinant(b,d)==0 &&
                       (query.corners[corner].array()>=n[i].f.cwiseMin(n[j].f).cwiseMin(n[k].f).array()).all() &&
                       (query.corners[corner].array()<=n[i].f.cwiseMax(n[j].f).cwiseMax(n[k].f).array()).all())
                        throw std::runtime_error("Corner has a nontransverse preimage on a grid face");
                    continue;
                }
                const double u=determinant(d,b)/det, v=determinant(a,d)/det;
                if(u<0 || v<0 || u+v>1) continue;
                if(u==0 || v==0 || u+v==1) throw std::runtime_error("Corner meets a grid edge");
                const auto [id,added]=vertex({3+corner,n[i].id,n[j].id,n[k].id},n[i].x+u*(n[j].x-n[i].x)+v*(n[k].x-n[i].x),query.corners[corner]);
                if(carrier_fit && added && !fit_corner(result.mesh.vertices[id],n[i],n[j],n[k],query.corners[corner])) {
                    ++failed_carriers;
                    // A corner fiber can cross a different grid face in the
                    // fitted model. Try its interior root before the global
                    // orientation/intersection pass decides which moves to keep.
                    Vec& x=result.mesh.vertices[id];
                    const unsigned fixed=solver.field.box.faces(x,solver.options.step*1e-10);
                    Vec trial=x,tangent;
                    if(guard_fit && solver.tangent(x,tangent) && solver.fiber(trial,query.corners[corner],tangent,x,
                            std::min(step.norm(),solver.options.spatial_tolerance),fixed)) x=trial;
                }
                cuts[corner].push_back(id);
                cuts[(corner+2)%3].push_back(id);
            }
        }
        for(int side=0;side<3;++side) {
            if(cuts[side].empty()) continue;
            Vec outward;
            for(int j=0;j<3;++j) outward[axes[j]]=-query.sides[side].n.dot(n[j+1].f-n[j].f)/step[axes[j]];
            polygon(cuts[side],outward,side);
        }
    }
public:
    GridSurface(const RangeTriangle& range,const Corrector& corrector,const std::string& mode)
        :query(range),solver(corrector),fit(mode=="fit-grid"),
         carrier_fit(mode=="carrier-grid" || mode=="recover-grid"),guard_fit(mode=="recover-grid") {}
    VolumeBoundary run(const Field& field,const Options& options,double spacing,const std::optional<Refinement>& refined) {
        const Vec extent=field.box.max-field.box.min;
        const Vec counts=(extent/spacing).array().ceil();
        if(!counts.allFinite() || counts.maxCoeff()>1000000) throw std::runtime_error("Requested grid is too large");
        if(refined && (!field.box.contains(refined->lo) || !field.box.contains(refined->hi) ||
           (refined->lo.array()>=refined->hi.array()).any() || refined->step>=spacing))
            throw std::invalid_argument("Refinement requires an interior interval and a finer spacing");
        // Subdivide full coordinate planes. Shared faces on the refined interval's
        // boundary retain identical nodes in both cells, including away from the box.
        std::array<std::vector<double>,3> positions;
        Eigen::Vector3i cells;
        for(int d=0;d<3;++d) {
            auto& axis=positions[d];
            axis.push_back(field.box.min[d]);
            const double base_step=extent[d]/counts[d];
            for(int i=1;i<=int(counts[d]);++i) {
                const double left=axis.back(),right=i==int(counts[d])?field.box.max[d]:field.box.min[d]+i*base_step;
                const double pieces=refined && left<refined->hi[d] && right>refined->lo[d] ?
                    std::ceil((right-left)/refined->step) : 1.;
                if(axis.size()+pieces>1000001) throw std::runtime_error("Requested refined axis is too large");
                for(int j=1;j<int(pieces);++j) axis.push_back(left+(right-left)*(j/pieces));
                axis.push_back(right);
            }
            cells[d]=int(axis.size())-1;
        }
        if((cells.cast<double>().array()+1).prod()>double(std::numeric_limits<std::int64_t>::max()/2))
            throw std::runtime_error("Requested grid exceeds addressable dimensions");
        const int nx=cells.x()+1,ny=cells.y()+1;
        const std::int64_t plane=std::int64_t(nx)*ny;
        std::vector<Range> lower(plane),upper(plane);
        auto coord=[&](int d,int i){return positions[d][i];};
        std::cout<<"Grid nodes: "<<nx<<" x "<<ny<<" x "<<(cells.z()+1)<<std::endl;
        auto sample=[&](int z,std::vector<Range>& slab) {
            tbb::parallel_for(0,ny,[&](int y) {
                for(int x=0;x<nx;++x) {
                    Range f=field.value(Vec(coord(0,x),coord(1,y),coord(2,z)));
                    if(!f.allFinite()) throw std::runtime_error("Nonfinite field sample");
                    slab[std::int64_t(y)*nx+x]=f;
                }
            });
        };
        sample(0,lower);
        for(int z=0;z<cells.z();++z) {
            sample(z+1,upper);
            for(int y=0;y<cells.y();++y) for(int x=0;x<cells.x();++x) {
                std::array<Node,8> cube;
                for(int m=0;m<8;++m) {
                    int xx=x+(m&1),yy=y+((m>>1)&1),zz=z+(m>>2);
                    const auto index=std::int64_t(yy)*nx+xx;
                    cube[m]={plane*zz+index,Vec(coord(0,xx),coord(1,yy),coord(2,zz)),(m&4?upper:lower)[index]};
                }
                bool possible=true,crossed=false;
                for(const auto& side:query.sides) {
                    double low=side.phi(cube[0].f),high=low;
                    for(int m=1;m<8;++m) {const double v=side.phi(cube[m].f);low=std::min(low,v);high=std::max(high,v);}
                    if(high<0) {possible=false;break;}
                    if(low<=0) crossed=true;
                }
                if(!possible || !crossed) continue;
                std::array<int,3> axes{0,1,2};
                do {
                    const int a=1<<axes[0],b=a|(1<<axes[1]);
                    tet({cube[0],cube[a],cube[b],cube[7]},axes,cube[7].x-cube[0].x);
                } while(std::next_permutation(axes.begin(),axes.end()));
            }
            lower.swap(upper);
        }
        crossings.clear();
        if(fit) {
            // Fit both range coordinates, including the along-side coordinate.
            // Shared seams have one vertex and one common corner target. Validate
            // the deformed triangles separately, since a Newton step can fold a sheet.
            size_t failed=0;
            for(size_t i=0;i<targets.size();++i) {
                const Vec start=result.mesh.vertices[i];
                Vec point=start,tangent;
                const unsigned box_faces=field.box.faces(start,options.step*1e-10);
                if(solver.tangent(start,tangent) && solver.fiber(point,targets[i],tangent,start,
                        std::min(spacing,options.spatial_tolerance),box_faces))
                    result.mesh.vertices[i]=point;
                else ++failed;
            }
            std::cout<<"Range fits failed: "<<failed<<" of "<<targets.size()<<std::endl;
        }
        if(carrier_fit) std::cout<<"Carrier fits failed: "<<failed_carriers<<" of "<<result.mesh.vertices.size()<<std::endl;
        if(guard_fit) restore_conflicts(options.step*1e-8);
        return std::move(result);
    }
};

void audit(const VolumeBoundary& surface,const Field& field,const RangeTriangle& query,const Options& options,
           const std::string& prefix,bool all_samples) {
    const auto& vertices=surface.mesh.vertices; const auto& faces=surface.mesh.triangles;
    const auto topology=marching_triangles::validate_mesh(vertices,faces,options.step*1e-8);
    std::cout<<topology.summary()<<std::endl;
    std::vector<double> cumulative;
    std::vector<unsigned> side_mask(vertices.size(),0);
    double area=0;
    for(size_t i=0;i<faces.size();++i) {
        const auto& f=faces[i];
        if(surface.patches[i]<=3) {
            area+=(vertices[f[1]]-vertices[f[0]]).cross(vertices[f[2]]-vertices[f[0]]).norm()/2;
            for(int id:f) side_mask[id]|=1u<<(surface.patches[i]-1);
        }
        cumulative.push_back(area);
    }
    if(!(area>0)) throw std::runtime_error("Empty side surface");
    size_t side_vertices=0,seam_tests=0,seam_failures=0;
    double max_vertex_residual=0,max_seam_residual=0,max_vertex_overshoot=0;
    for(size_t i=0;i<vertices.size();++i) {
        const unsigned mask=side_mask[i];
        if(!mask) continue;  // Cap-only and unused vertices have no range-side target.
        ++side_vertices;
        const Range value=field.value(vertices[i]);
        if(!value.allFinite()) throw std::runtime_error("Nonfinite field at a side vertex");
        for(int k=0;k<3;++k) {
            const unsigned side_bit=1u<<k;
            if(!(mask&side_bit)) continue;
            const auto& side=query.sides[k];
            const double tau=side.parameter(value);
            const double bounded=std::clamp(tau,0.,1.);
            max_vertex_overshoot=std::max(max_vertex_overshoot,std::abs(tau-bounded));
            if(mask==side_bit) {
                const Range q=side.a+bounded*side.d;
                max_vertex_residual=std::max(max_vertex_residual,(value-q).norm());
                break;
            }
            const unsigned corner_mask=(1u<<k)|(1u<<((k+2)%3));
            if((mask&corner_mask)!=corner_mask) continue;
            const double residual=(value-query.corners[k]).norm();
            ++seam_tests;
            if(residual>options.corner_residual_tolerance) ++seam_failures;
            max_seam_residual=std::max(max_seam_residual,residual);
        }
    }
    max_vertex_residual=std::max(max_vertex_residual,max_seam_residual);
    std::vector<double> distance,angle;
    std::mt19937 rng(331);
    std::uniform_real_distribution<double> uniform(0,1);
    int failed=0,overturned=0;
    constexpr int samples=20000;
    for(int k=0;k<samples;++k) {
        const size_t i=std::upper_bound(cumulative.begin(),cumulative.end(),uniform(rng)*area)-cumulative.begin();
        const auto& f=faces[i]; const int label=surface.patches[i];
        const Vec a=vertices[f[0]],ab=vertices[f[1]]-a,ac=vertices[f[2]]-a;
        const double u=std::sqrt(uniform(rng)),v=uniform(rng);
        const Vec x=a+u*((1-v)*ab+v*ac);
        const Vec g=-query.sides[label-1].gradient(field.jacobian(x));
        const double cosine=ab.cross(ac).normalized().dot(g.normalized());
        if(!std::isfinite(cosine)) throw std::runtime_error("Undefined sampled normal");
        angle.push_back(std::acos(std::clamp(cosine,-1.,1.))*180/std::acos(-1.));
        if(cosine<=0) ++overturned;
        Vec projected=x;
        if(correct_volume_sample(projected,label,field,query,options)) distance.push_back((projected-x).norm());
        else ++failed;
    }
    std::sort(distance.begin(),distance.end());std::sort(angle.begin(),angle.end());
    auto percentile=[](const std::vector<double>& values,double p) {
        if(values.empty())return std::string("null");
        std::ostringstream s;s<<std::setprecision(12)<<values[size_t(std::ceil(values.size()*p))-1];return s.str();
    };
    std::ofstream report(prefix+".json");
    report.exceptions(std::ios::failbit|std::ios::badbit);
    report<<std::setprecision(15)<<"{\n";
    if(all_samples) {
        // Only edges shared by different labels represent a corner fiber.
        // Near a fold, accurate endpoints alone can still span the wrong arc.
        std::map<Edge,unsigned> edge_labels;
        for(size_t i=0;i<faces.size();++i) {
            const int label=surface.patches[i];
            if(label>3) continue;
            for(int j=0;j<3;++j) {
                const int a=faces[i][j],b=faces[i][(j+1)%3];
                const unsigned common=side_mask[a]&side_mask[b];
                if(common & (common-1)) edge_labels[edge_key(a,b)]|=1u<<(label-1);
            }
        }
        size_t chord_sites=0,chord_failures=0,chord_free_projections=0;
        double max_chord_projection=0;
        for(const auto& [edge,labels]:edge_labels) for(int k=0;k<3;++k) {
            const unsigned mask=(1u<<k)|(1u<<((k+2)%3));
            if((labels&mask)!=mask) continue;
            const Vec &a=vertices[edge.first],&b=vertices[edge.second];
            const Vec direction=(b-a).normalized();
            Corrector solver(field,query.sides[k],options);
            for(double fraction:{.25,.5,.75}) {
                const Vec start=a+fraction*(b-a);
                Vec p=start;
                ++chord_sites;
                if(!solver.fiber(p,query.corners[k],direction,start,options.spatial_tolerance)) {
                    // The phase plane can be tangent to the fiber at a turn.
                    // A bounded minimum-norm solve still measures spatial distance.
                    p=start;
                    if(!solver.fiber_seed(p,query.corners[k]) || (p-start).norm()>options.spatial_tolerance) {
                        ++chord_failures;
                        continue;
                    }
                    ++chord_free_projections;
                }
                max_chord_projection=std::max(max_chord_projection,(p-start).norm());
            }
        }
        report<<"  \"seam_chord_sites\": "<<chord_sites<<", \"seam_chord_failures\": "<<chord_failures
              <<", \"seam_chord_free_projections\": "<<chord_free_projections
              <<", \"max_seam_chord_projection\": "<<max_chord_projection<<",\n";
        struct FaceChecks {
            size_t checked=0,failed=0,over_55=0;
            double max_angle=0,max_distance=0;
            size_t cap_sites=0,cap_failed=0,cap_wrong_plane=0;
            double cap_distance=0;
        };
        tbb::enumerable_thread_specific<FaceChecks> workers;
        tbb::parallel_for(size_t(0),faces.size(),[&](size_t i) {
            const int label=surface.patches[i];
            auto& checks=workers.local();
            const auto& t=faces[i];
            const Vec a=vertices[t[0]],b=vertices[t[1]],c=vertices[t[2]];
            if(label>3) {
                const int face=label-4,axis=face/2;
                const double plane=face%2?field.box.max[axis]:field.box.min[axis];
                if(a[axis]!=plane || b[axis]!=plane || c[axis]!=plane) ++checks.cap_wrong_plane;
                for(const Vec& x:volume_triangle_samples(a,b,c)) {
                    ++checks.cap_sites;
                    Vec p=x;
                    if(correct_volume_sample(p,label,field,query,options))
                        checks.cap_distance=std::max(checks.cap_distance,(p-x).norm());
                    else ++checks.cap_failed;
                }
                return;
            }
            const Vec normal=(b-a).cross(c-a).normalized();
            double worst=0;
            for(const Vec& x:volume_triangle_samples(a,b,c)) {
                const Vec g=-query.sides[label-1].gradient(field.jacobian(x));
                const double cosine=normal.dot(g.normalized());
                if(!std::isfinite(cosine)) throw std::runtime_error("Undefined field normal");
                ++checks.checked;
                worst=std::max(worst,std::acos(std::clamp(cosine,-1.,1.))*180/std::acos(-1.));
                Vec point=x;
                if(correct_volume_sample(point,label,field,query,options))
                    checks.max_distance=std::max(checks.max_distance,(point-x).norm());
                else ++checks.failed;
            }
            checks.max_angle=std::max(worst,checks.max_angle);
            if(worst>55) ++checks.over_55;
        });
        FaceChecks total;
        for(const auto& checks:workers) {
            total.checked+=checks.checked;
            total.failed+=checks.failed;
            total.over_55+=checks.over_55;
            total.max_angle=std::max(total.max_angle,checks.max_angle);
            total.max_distance=std::max(total.max_distance,checks.max_distance);
            total.cap_sites+=checks.cap_sites;
            total.cap_failed+=checks.cap_failed;
            total.cap_wrong_plane+=checks.cap_wrong_plane;
            total.cap_distance=std::max(total.cap_distance,checks.cap_distance);
        }
        report<<"  \"face_sites_checked\": "<<total.checked<<", \"face_sites_projection_failed\": "<<total.failed
              <<", \"faces_over_55_degrees\": "<<total.over_55<<", \"face_sites_max_angle\": "<<total.max_angle
              <<", \"face_sites_max_projection\": "<<total.max_distance<<",\n";
        report<<"  \"cap_sites_checked\": "<<total.cap_sites<<", \"cap_sites_projection_failed\": "<<total.cap_failed
              <<", \"cap_faces_off_plane\": "<<total.cap_wrong_plane<<", \"cap_max_projection\": "<<total.cap_distance<<",\n";
    }
    report<<"  \"triangles\": "<<faces.size()<<", \"vertices\": "<<vertices.size()
          <<", \"topology\": \""<<topology.summary()<<"\",\n  \"intersections_checked\": "
          <<(!topology.invalid_indices && !topology.degenerate_faces?"true":"false")
          <<", \"side_vertices_checked\": "<<side_vertices<<", \"max_vertex_range_residual\": "<<max_vertex_residual
          <<", \"max_vertex_overshoot\": "<<max_vertex_overshoot
          <<", \"seam_tests\": "<<seam_tests<<", \"seam_failures\": "<<seam_failures
          <<", \"max_seam_range_residual\": "<<max_seam_residual
          <<", \"side_area\": "<<area<<", \"samples\": "<<samples<<", \"failed_projections\": "<<failed
          <<", \"reversed_sample_normals\": "<<overturned<<",\n  \"projection_p95\": "<<percentile(distance,.95)
          <<", \"projection_max\": "<<percentile(distance,1)<<", \"normal_degrees_p95\": "<<percentile(angle,.95)
          <<", \"normal_degrees_max\": "<<percentile(angle,1)<<"\n}\n";
}
}

int main(int argc,char** argv) {
    try {
        if(argc!=5 && argc!=6 && argc!=9)
            throw std::invalid_argument("Use: grid|fit-grid|carrier-grid|recover-grid STAGES SPACING PREFIX [--refine LO HI SPACING], cap-sides STAGES SIDES.ply PREFIX, project-samples STAGES SAMPLES.txt PREFIX, or audit STAGES INPUT.ply PREFIX [--all-face-samples]");
        const auto config=read_config(argv[2]);
        if(config.curve.size()!=3 || !config.curve.closed()) throw std::invalid_argument("Expected a closed range triangle");
        const auto field=load_field(config,true);
        RangeTriangle range(config.curve.point(0),config.curve.point(1),config.curve.point(2));
        const std::string mode=argv[1],prefix=argv[4];
        VolumeBoundary boundary;
        if(mode=="project-samples" && argc==5) {
            std::ifstream input(argv[3]);
            if(!input) throw std::runtime_error("Cannot open sample list");
            std::ofstream out(prefix+".csv");
            out.exceptions(std::ios::badbit | std::ios::failbit);
            out<<"patch,x,y,z,success,corrected_x,corrected_y,corrected_z,distance\n"<<std::setprecision(17);
            size_t count=0,failed=0;
            while(input>>std::ws && !input.eof()) {
                int patch;
                Vec point;
                if(!(input>>patch>>point[0]>>point[1]>>point[2]) || patch<1 || patch>9 || !point.allFinite())
                    throw std::runtime_error("Expected patch and three finite sample coordinates");
                Vec projected=point;
                const bool ok=correct_volume_sample(projected,patch,field.field,range,config.options);
                out<<patch<<','<<point[0]<<','<<point[1]<<','<<point[2]<<','<<ok<<','
                   <<projected[0]<<','<<projected[1]<<','<<projected[2]<<','<<(projected-point).norm()<<'\n';
                ++count;
                if(!ok) ++failed;
            }
            std::cout<<"Checked "<<count<<" supplied samples, failed "<<failed<<std::endl;
            return failed?2:0;
        }
        const bool all_samples=argc==6 && std::string(argv[5])=="--all-face-samples";
        if((mode=="grid" || mode=="fit-grid" || mode=="carrier-grid" || mode=="recover-grid") && (argc==5 || argc==9)) {
            std::optional<Refinement> refinement;
            if(argc==9) {
                if(std::string(argv[5])!="--refine") throw std::invalid_argument("Expected --refine LO HI SPACING");
                refinement=Refinement{coordinates<3>(argv[6]),coordinates<3>(argv[7]),positive_spacing(argv[8])};
            }
            const Corrector solver{field.field,range.sides[0],config.options};
            boundary=GridSurface(range,solver,mode).run(field.field,config.options,positive_spacing(argv[3]),refinement);
            // Side patches only. A closed surrogate alone does not validate a fitted-field volume.
            write_boundary_ply(prefix+".ply",boundary,false);
        } else if(mode=="cap-sides" && argc==5) {
            auto input=read_boundary_ply(argv[3],true);
            Result sides;
            sides.box_contact_triangles=input.box_contact_triangles;
            sides.mesh.vertices=std::move(input.mesh.vertices);
            sides.mesh.triangles=std::move(input.mesh.triangles);
            for(size_t i=0;i<sides.mesh.triangles.size();++i) {
                const int label=input.patches.at(i);
                if(label<1 || label>3) throw std::invalid_argument("cap-sides expects only side patches");
                sides.face_segments.push_back(label-1);
                std::swap(sides.mesh.triangles[i][1],sides.mesh.triangles[i][2]);
            }
            boundary=assemble_volume_boundary(field.field,range,config.options,std::move(sides));
            boundary.corner_components=input.corner_components;
            for(const auto& event:boundary.unresolved)
                std::cerr<<event.position.transpose()<<"  "<<event.reason<<'\n';
            write_boundary_ply(prefix+".ply",boundary,false);
            std::cout<<"Assembly events: "<<boundary.unresolved.size()<<std::endl;
            return boundary.unresolved.empty()?0:2;
        } else if(mode=="audit" && (argc==5 || all_samples)) boundary=read_boundary_ply(argv[3]);
        else throw std::invalid_argument("Unknown mode: "+mode);
        audit(boundary,field.field,range,config.options,prefix,all_samples);
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
