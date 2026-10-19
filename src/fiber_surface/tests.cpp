#include "analytical_field.h"
#include "output.h"
#include "extract.h"

using namespace fiber_surface;
static void require(bool ok,const std::string& message) { if (!ok) throw std::runtime_error(message); }

static void derivatives() {
    const auto e=analytical_field("sphere"); const Vec p(.2,-.4,.7); const Segment s(Range(.1,.2),Range(2,.6));
    const auto j=e.field.jacobian(p); const auto h=e.field.hessians(p); const double delta=1e-5;
    for (int d=0; d<3; ++d) {
        const Vec a=p-delta*Vec::Unit(d),b=p+delta*Vec::Unit(d);
        require(((e.field.value(b)-e.field.value(a))/(2*delta)-j.col(d)).norm()<1e-8,"Analytical Jacobian");
        require(std::abs((s.phi(e.field.value(b))-s.phi(e.field.value(a)))/(2*delta)-s.gradient(j)[d])<1e-8,"Segment gradient");
        require((s.gradient(e.field.jacobian(b))-s.gradient(e.field.jacobian(a))-2*delta*s.hessian(h).col(d)).norm()<1e-8,"Segment Hessian");
    }
    require((e.field.value(p)-s.a-s.parameter(e.field.value(p))*s.d-s.phi(e.field.value(p))*s.n).norm()<1e-12,"Range decomposition");
    Options options; Corrector corrector(e.field,e.segment,options);
    Vec pole(0,0,-1); require(corrector.surface(pole,.1),"Rank-one J must not reject regular sphere poles");
    Vec tangent; require(!corrector.tangent(pole,tangent),"Rank-one J cannot define a regular endpoint fiber");
    Field guarded=e.field;
    guarded.value=[e](const Vec& x) { if (!e.field.box.contains(x)) throw std::runtime_error("Out-of-domain evaluation"); return e.field.value(x); };
    Corrector bounded(guarded,e.segment,options); Vec outside(3,0,0);
    require(!bounded.surface(outside,1),"Outside correction should fail before model evaluation");
    auto strip=analytical_field("strip"); Corrector planar(strip.field,strip.segment,options);
    Vec contact(0.01,.01,2);
    require(planar.fiber(contact,Range(0,0),Vec::Zero(),contact,.1,1u<<5) && (contact-Vec(0,0,2)).norm()<1e-10,"Simultaneous fiber/face contact");
    strip.field.box.max[1]=0;
    Corrector coincident(strip.field,strip.segment,options); Vec on_face(.5,0,0);
    require(!coincident.surface(on_face,.1,1u<<3),"A coincident box face is not a transverse boundary");
    Options loose; loose.residual_tolerance=.005;
    Corrector loose_corrector(e.field,e.segment,loose); Result traced;
    FiberTracer(loose_corrector,traced).run(discovery_grid(e.field.box,7));
    require(traced.fibers.size()==1 && traced.fibers.front().closed && traced.unresolved.empty(),
            "Loose sampling slack must not create duplicate endpoint components");
}
int main(int argc,char** argv) {
    try {
        const std::string name=argc>1?argv[1]:"sphere";
        if (name=="derivatives") { derivatives(); return 0; }
        auto example=analytical_field(name); Options options;
        const Result r=Extractor(example.field,example.segment,options).run();
        const auto report=evaluate(r,example.field,example.segment,options);
        write_outputs("test_"+name,r,report,example.field,example.segment,options);
        std::cout<<name<<": triangles="<<r.mesh.triangles.size()<<" components="<<report.components<<" unresolved="<<r.unresolved.size()<<" "<<report.topology.summary()<<std::endl;
        require(report.complete(r),"Extraction is incomplete");
        require(report.vertices.maximum<options.residual_tolerance*1.01 && report.endpoint_residual<options.residual_tolerance*1.01,"Constraint residual");
        require(report.components==(name=="disconnected"?2u:1u),"Incorrect connected components");
        require(report.topology.euler_characteristic==(name=="closed-sphere" || name=="disconnected" ? 2 : 1),"Incorrect Euler characteristic");
        if (name=="sphere" || name=="closed-sphere") {
            const double expected_area=(name=="sphere" ? 3 : 4)*std::acos(-1.);
            require(std::abs(report.area-expected_area)/expected_area<.025,"Sphere coverage/area");
            require(r.fibers.size()==(name=="sphere"?1u:0u),"Sphere endpoint fibers");
            if (name=="sphere") {
                require(r.fibers[0].closed,"Circular endpoint fiber must close");
                require(report.topology.boundary_edges==r.fiber_edges.size(),"Circular opening was changed");
            } else require(report.topology.closed(),"Closed sphere must have no boundary");
        } else {
            require(std::abs(report.area-4)<1e-7,"Strip area/box coverage");
            require(r.fibers.size()==(name=="strip"?2u:4u),"Disconnected endpoint fibers");
            for (const auto& fiber:r.fibers) {
                require(!fiber.closed && fiber.complete,"Strip fiber must terminate on the box");
                require(r.info[fiber.vertices.front()].faces && r.info[fiber.vertices.back()].faces,"Fiber-box contact IDs");
            }
        }
        return 0;
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
