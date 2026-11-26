#include "volume_examples.h"
#include <filesystem>
#include <limits>

using namespace fiber_surface;
static void require(bool ok,const std::string& why) {if(!ok)throw std::runtime_error(why);}

int main(int argc,char** argv) {
    try {
        const std::string name=argc>1?argv[1]:"invalid-patches";
        auto example=volume_example("prism");Options options;
        if(name=="translated-box") {
            // A wholly selected box has no side surfaces. Translating it must
            // preserve its six caps and volume, even far from the origin.
            Field field;
            field.box.min=Vec::Constant(1e9);field.box.max=field.box.min+Vec(3,3,1);
            field.value=[](const Vec&) {return Range(.25,.25);};
            field.jacobian=[](const Vec&) {return Jacobian::Zero().eval();};
            field.hessians=[](const Vec&) {return std::array<Hessian,2>{Hessian::Zero(),Hessian::Zero()};};
            const auto result=extract_volume(field,example.triangle,options);
            write_volume("volume_review_translated/volume",result,field,example.triangle,options);
            require(result.complete(),"Translation lost or invalidated selected box caps");
            require(std::abs(result.report.volume-9)<1e-9,"Translated box volume changed");
            return 0;
        }
        auto result=extract_volume(example.field,example.triangle,options);
        require(result.complete(),"Prism fixture failed");
        if(name=="invalid-patches") {
            for(int patch:{0,-1,10,std::numeric_limits<int>::min(),std::numeric_limits<int>::max()}) {
                Vec point(.25,.25,.5);
                require(!correct_volume_sample(point,patch,example.field,example.triangle,options),"Invalid patch accepted by sample corrector");
                auto bad=result.boundary;bad.patches[0]=patch;
                const auto report=validate_volume_boundary(bad,example.field,example.triangle,options);
                require(!report.valid && !bad.unresolved.empty(),"Invalid boundary patch was not rejected");
            }
        } else if(name=="domain-guard") {
            auto outside=result.boundary;
            for(Vec& x:outside.mesh.vertices)x[0]+=10;
            Field guarded=example.field;
            guarded.value=[original=example.field](const Vec& x) {
                if(!original.box.contains(x))throw std::runtime_error("Out-of-domain value evaluation");
                return original.value(x);
            };
            guarded.jacobian=[original=example.field](const Vec& x) {
                if(!original.box.contains(x))throw std::runtime_error("Out-of-domain derivative evaluation");
                return original.jacobian(x);
            };
            require(!validate_volume_boundary(outside,guarded,example.triangle,options).valid,"Out-of-domain boundary accepted");
        } else if(name=="uncached-boundary") {
            // Imported meshes can provide vertices/facets without growth caches.
            auto imported=result.boundary;
            imported.mesh.edges.clear();imported.mesh.boundary.clear();imported.mesh.faces.clear();
            const auto report=validate_volume_boundary(imported,example.field,example.triangle,options);
            require(report.valid && report.signed_volumes.size()==1,"Boundary validation depends on stale adjacency caches");
            require(std::abs(report.signed_volumes[0]-.5)<1e-9,"Reconstructed shell volume changed");
        } else if(name=="invalid-cells") {
            for(const Tetrahedron& cell:std::vector<Tetrahedron>{{-1,0,1,2},{0,0,0,0}}) {
                auto bad=result.mesh;bad.cells={cell};bad.boundary_faces.clear();bad.parent_facets.clear();
                const auto report=assess_volume(bad,result.boundary,result.boundary_report,example.field,example.triangle,options);
                require(!report.valid && report.invalid_cells>0,"Malformed cell accepted");
                require(std::isfinite(report.min_volume) && std::isfinite(report.max_volume) && std::isfinite(report.volume),
                        "Invalid cells produce non-finite report statistics (invalid JSON)");
            }
        } else if(name=="cap-constraints") {
            Field field;
            field.box.min=Vec::Zero();field.box.max=Vec::Ones();
            field.value=[](const Vec& x) {return Range(x[0]+x[1],x[2]);};
            field.jacobian=[](const Vec&) {Jacobian j;j<<1,1,0,0,0,1;return j;};
            field.hessians=example.field.hessians;
            const RangeTriangle selected(Range(.01,0),Range(1.01,0),Range(.01,1));
            Vec near_wall(.001,.013,.5);
            require(correct_volume_sample(near_wall,3,field,selected,options) &&
                    (near_wall-Vec(0,.01,.5)).norm()<1e-10,
                    "A bounded side projection did not recover the root on the first box wall");
            Vec beyond(.001,.013,.5);
            const RangeTriangle outside(Range(-.01,0),Range(.99,0),Range(-.01,1));
            require(!correct_volume_sample(beyond,3,field,outside,options),
                    "An unattainable side equality was accepted at a box corner");
            Vec edge(0,0,.5);
            require(!correct_volume_sample(edge,4,field,selected,options),
                    "Cap correction left a second active box face to satisfy a range constraint");
            field.value=[](const Vec&) {return Range(.25,.25);};
            field.jacobian=[](const Vec&) {return Jacobian::Zero().eval();};
            Vec nearby(.01,.25,.5);
            require(correct_volume_sample(nearby,4,field,example.triangle,options) && nearby[0]==0,
                    "Cap correction failed to enforce its assigned box plane");
            Vec distant(.1,.25,.5);
            require(!correct_volume_sample(distant,4,field,example.triangle,options),
                    "Cap projection exceeded the spatial tolerance");
        } else if(name=="scale-small" || name=="scale-tiny") {
            const double scale=name=="scale-small"?1e-8:1e-55;
            auto boundary=result.boundary;auto mesh=result.mesh;auto boundary_report=result.boundary_report;
            for(Vec& x:boundary.mesh.vertices)x*=scale;
            for(Vec& x:mesh.vertices)x*=scale;
            for(double& volume:boundary_report.signed_volumes)volume*=scale*scale*scale;
            Field field=example.field;
            field.box.min*=scale;field.box.max*=scale;
            field.value=[scale](const Vec& x) {return Range(x[0]/scale,x[1]/scale);};
            field.jacobian=[scale](const Vec&) {Jacobian j=Jacobian::Zero();j(0,0)=j(1,1)=1/scale;return j;};
            options.step*=scale;options.min_step*=scale;options.spatial_tolerance*=scale;
            const auto report=assess_volume(mesh,boundary,boundary_report,field,example.triangle,options);
            std::cout<<"scale="<<scale<<" invalid_cells="<<report.invalid_cells<<" invalid_faces="<<report.invalid_faces
                <<" intersections="<<report.intersections<<" conformity="<<report.conformity_failures<<'\n';
            require(report.valid,"A change of coordinate units invalidated an otherwise valid mesh");
            require(std::abs(report.volume/(scale*scale*scale)-.5)<1e-9,"Scaled volume is inconsistent");
        } else throw std::invalid_argument("Unknown review test: "+name);
        std::cout<<name<<" passed\n";
        return 0;
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
