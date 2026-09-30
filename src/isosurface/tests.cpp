#include "explicit_field.h"
#include "find_initial_root.h"
#include "marching_triangles.h"
#include "mesh_validation.h"
#include "root_io.h"
#include <filesystem>
#include <iostream>
#include <random>

using namespace marching_triangles;
using P = Point<double>;
void require(bool pass,const std::string& message) {if(!pass)throw std::runtime_error(message);}

SurfaceField<double> sphere(double radius=1,P center=P::Zero(),double sign=1) {
    return {[=](const P& p){return sign*((p-center).squaredNorm()/(radius*radius)-1);},
            [=](const P& p)->P{return sign*2*(p-center)/(radius*radius);}};
}

void check_surface(const std::string& name,SurfaceField<double> field,MeshingOptions<double> options,
                   std::vector<P> seeds,size_t components,long euler) {
    MarchingTriangles<double> mesher(field,options);
    std::vector<std::vector<P>> vertices;
    std::vector<std::vector<Triangle>> triangles;
    require(mesher.extract_all_sheets(seeds,vertices,triangles),name+": no sheets");
    require(vertices.size()==components,name+": wrong component count "+std::to_string(vertices.size()));
    for(size_t s=0;s<vertices.size();++s) {
        auto report=validate_mesh(vertices[s],triangles[s],options.step*1e-9);
        require(report.closed(),name+": "+report.summary()+" "+mesher.diagnostic());
        require(report.euler_characteristic==euler,name+": wrong Euler characteristic");
        for(const auto& p:vertices[s])require(std::abs(field.value(p)-options.iso_value)<=options.projection_tolerance*1.01,name+": off-surface vertex");
        double area=0;
        for(const auto& t:triangles[s]) {
            const P a=vertices[s][t[0]],b=vertices[s][t[1]],c=vertices[s][t[2]],n=(b-a).cross(c-a);
            area+=n.norm()/2;
            require(n.dot(field.gradient((a+b+c)/3))>0,name+": inverted triangle");
        }
        require(area>options.step*options.step*5,name+": incomplete small patch");
        const auto quality=measure_mesh_quality(vertices[s],triangles[s],field);
        require(quality.reversed_faces==0,name+": incorrect normal orientation at a vertex or centroid");
        require(quality.min_angle_degrees>5&&quality.max_normal_error_degrees<45,
                name+": poor triangle quality: "+quality.summary());
        std::cout<<name<<" sheet "<<s<<": V="<<vertices[s].size()<<" F="<<triangles[s].size()<<" "<<report.summary()<<" "<<quality.summary()<<std::endl;
    }
}

void geometry_tests() {
    require(make_triangle_key(1,2,3)==make_triangle_key(2,1,3),"reversed duplicate faces must match");
    std::vector<P> v{{0,0,0},{1,0,0},{0,1,0},{1,1,0},{.2,.2,0},{.7,.2,0},{.2,.7,0},
                     {.25,.25,-1},{.25,.25,1},{.75,.25,0},{-1,0,0},{0,-1,0}};
    require(!triangles_conflict(v,{0,1,2},{1,3,2},1e-10),"legal shared edge rejected");
    require(!triangles_conflict(v,{0,1,2},{0,10,11},1e-10),"legal shared vertex rejected");
    require(triangles_conflict(v,{0,1,2},{4,5,6},1e-10),"coplanar overlap missed");
    require(triangles_conflict(v,{0,1,2},{7,8,9},1e-10),"3D intersection missed");
    require(triangles_conflict(v,{0,1,2},{0,1,4},1e-10),"shared-edge overlap missed");
    require(triangles_conflict(v,{0,1,2},{0,5,6},1e-10),"shared-vertex overlap missed");
    TriangleMesh<double> mesh;mesh.vertices=v;mesh.add({0,1,2});
    require(!mesh.topology_ok({2,1,0}),"reversed duplicate accepted");
    require(!mesh.topology_ok({1,2,3}),"same directed edge reused");
    require(mesh.topology_ok({1,3,2}),"valid wedge rejected");mesh.add({1,3,2});
    require(mesh.boundary.size()==4,"wedge did not remove shared edge");
    require(!mesh.topology_ok({1,2,4}),"third face accepted on edge");
    require(mesh.boundary_loops().size()==1&&mesh.boundary_loops()[0].size()==4,"wrong wedge contour");
    auto field=sphere();MeshingOptions<double> options;
    SurfaceProjector<double> projector(field,options);
    P p(2,0,0);require(projector.project(p)&&std::abs(p.norm()-1)<1e-8,"bracket projection failed");
    p.setZero();require(!projector.project(p),"zero-gradient projection accepted");
    p=P(std::numeric_limits<double>::quiet_NaN(),0,0);require(!projector.project(p),"NaN projection accepted");
    bool threw=false;options.step=0;try{options.validate();}catch(const std::invalid_argument&){threw=true;}require(threw,"invalid step accepted");
    const auto file=std::filesystem::temp_directory_path()/"isosurface-regression-roots.dat";
    write_roots(file.string(),{P(1,2,3),P(-4,5,6)});
    auto roots=read_roots(file.string());require(roots.size()==2&&roots[1]==P(-4,5,6),"binary root round trip failed");
    {std::ofstream out(file,std::ios::binary);int count=1;out.write(reinterpret_cast<char*>(&count),sizeof(count));}
    threw=false;try{read_roots(file.string());}catch(const std::runtime_error&){threw=true;}
    std::filesystem::remove(file);require(threw,"truncated root file accepted");
    for(const std::string name:{"ellipsoid","quartic_potential","quartic_potential_2","rotating_quartic_multiwell"}) {
        const auto setup=explicit_field(name,name=="ellipsoid"?1:0);
        const P sample(.4,.6,1.7),gradient=setup.field.gradient(sample);
        for(int axis=0;axis<3;++axis) {
            P plus=sample,minus=sample;plus[axis]+=1e-5;minus[axis]-=1e-5;
            const double finite_difference=(setup.field.value(plus)-setup.field.value(minus))/2e-5;
            require(std::abs(finite_difference-gradient[axis])<1e-7,name+": incorrect gradient callback");
        }
    }
    std::cout<<"geometry, projection, topology, root I/O: passed"<<std::endl;
}

void quality_tests() {
    auto setup=explicit_field("ellipsoid",1);
    setup.options.step=.4;setup.options.min_step=.04;
    setup.options.quality_iterations=0;
    MarchingTriangles<double> mesher(setup.field,setup.options);
    std::vector<std::vector<P>> vertices;
    std::vector<std::vector<Triangle>> triangles;
    require(mesher.extract_all_sheets({P(1,0,0)},vertices,triangles)&&vertices.size()==1,
            "quality: could not reproduce the unoptimized ellipsoid");
    const auto before=measure_mesh_quality(vertices[0],triangles[0],setup.field);
    require(before.triangles_below_5_degrees>0,"quality: fixture no longer exercises slivers");
    TriangleMesh<double> mesh;mesh.vertices=vertices[0];
    for(const auto& t:triangles[0])mesh.add(t);
    const std::vector<P> stops;
    MeshQualityOptimizer<double> optimizer(setup.field,setup.options,stops);
    optimizer.improve(mesh);
    require(mesh.vertices==vertices[0]&&mesh.triangles==triangles[0],"quality: disabled optimizer changed the mesh");
    setup.options.quality_iterations=8;
    optimizer.improve(mesh);
    const auto after=measure_mesh_quality(mesh.vertices,mesh.triangles,setup.field);
    const auto report=validate_mesh(mesh.vertices,mesh.triangles,setup.options.step*1e-9);
    require(report.closed()&&report.euler_characteristic==2,"quality: "+report.summary());
    require(mesh.vertices.size()==vertices[0].size()&&mesh.triangles.size()==triangles[0].size(),
            "quality: optimization changed vertex or face counts");
    require(after.min_angle_degrees>10&&after.triangles_below_5_degrees==0,
            "quality: thin triangles remain: "+after.summary());
    require(after.max_normal_error_degrees<35&&after.reversed_faces==0,
            "quality: poor normal alignment: "+after.summary());
    for(const auto& p:mesh.vertices)
        require(std::abs(setup.field.value(p)-setup.options.iso_value)<=setup.options.projection_tolerance*1.01,
                "quality: relocation left the implicit surface");
    // Reconstruct connectivity independently of replace_pair's cached maps.
    TriangleMesh<double> rebuilt;rebuilt.vertices=mesh.vertices;
    for(const auto& t:mesh.triangles) {
        require(rebuilt.topology_ok(t),"quality: edge flip corrupted topology");
        rebuilt.add(t);
    }
    require(mesh.boundary==rebuilt.boundary&&mesh.faces==rebuilt.faces&&mesh.edges.size()==rebuilt.edges.size(),
            "quality: stale cached connectivity");
    for(const auto& entry:rebuilt.edges) {
        const auto& actual=mesh.edges.at(entry.first);
        const auto& expected=entry.second;
        require((actual.face==expected.face&&actual.second==expected.second&&actual.from==expected.from&&actual.to==expected.to)||
                (actual.face==expected.second&&actual.second==expected.face&&actual.from==expected.to&&actual.to==expected.from),
                "quality: stale edge incidence");
    }
    std::cout<<"quality before: "<<before.summary()<<"\nquality after: "<<after.summary()<<std::endl;
}

int main(int argc,char** argv) {
    try {
        const std::string filter=argc>1?argv[1]:"all";
        if(filter=="all"||filter=="geometry")geometry_tests();
        if(filter=="all"||filter=="quality")quality_tests();
        MeshingOptions<double> options;options.step=.3;options.min_step=.03;
        if(filter=="all"||filter=="sphere") {
            check_surface("sphere",sphere(),options,{P::Zero(),P(1,0,0),P(0,1,0),P(0,0,-1)},1,2);
            check_surface("negative field",sphere(1,P::Zero(),-1),options,{P(1,0,0)},1,2);
        }
        if(filter=="all"||filter=="ellipsoid") {
            for(double step:{.6,.4,.25}) {
                auto setup=explicit_field("ellipsoid",1);setup.options.step=step;setup.options.min_step=step*.1;
                check_surface("ellipsoid step="+std::to_string(step),setup.field,setup.options,{P(1,0,0),P(0,2,0),P(0,0,-3)},1,2);
            }
        }
        if(filter=="all"||filter=="torus") {
            SurfaceField<double> torus{
                [](const P& p){double q=p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+2;return q*q-9*(p[0]*p[0]+p[1]*p[1]);},
                [](const P& p)->P{double q=p.squaredNorm()+2;return {4*p[0]*(q-4.5),4*p[1]*(q-4.5),4*p[2]*q};}};
            // A tubular surface repeatedly reduces the adaptive edge length.
            // Keep d_min comparable to the seed edge, as intended in 3.2.1.
            options.step=.25;options.min_step=.15;
            check_surface("torus",torus,options,{P(2,0,0),P(1,0,0),P(0,1.5,.5)},1,0);
        }
        if(filter=="all"||filter=="components") {
            SurfaceField<double> pair{
                [](const P& p){return ((p-P(-1.3,0,0)).squaredNorm()-1)*((p-P(1.3,0,0)).squaredNorm()-1);},
                [](const P& p)->P{P a=p-P(-1.3,0,0),b=p-P(1.3,0,0);return 2*a*(b.squaredNorm()-1)+2*b*(a.squaredNorm()-1);}};
            check_surface("two spheres",pair,options,{P(-2.3,0,0),P(.3,0,0),P(-1.3,1,0),P(2.3,0,0)},2,2);
            SurfaceField<double> shells{
                [](const P& p){return (p.squaredNorm()-1)*(p.squaredNorm()-1.21);},
                [](const P& p)->P{return 2*p*(2*p.squaredNorm()-2.21);}};
            check_surface("thin concentric shells",shells,options,{P(1,0,0),P(1.1,0,0),P(0,1,0)},2,2);
        }
        if(filter=="all"||filter=="scale") {
            for(double scale:{.001,1000.}) {
                options.domain_min=P::Constant(-2*scale);options.domain_max=P::Constant(2*scale);
                options.step=.3*scale;options.min_step=.03*scale;options.vertex_tolerance=scale*1e-6;
                check_surface("scaled sphere "+std::to_string(scale),sphere(scale),options,{P(scale,0,0)},1,2);
            }
        }
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<std::endl;return 1;}
}
