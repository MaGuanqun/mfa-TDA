#include "analytical_field.h"
#include "output.h"
#include "extract.h"
#include <tbb/global_control.h>
#include <tbb/task_arena.h>

using namespace fiber_surface;
static void require(bool ok,const std::string& message) { if (!ok) throw std::runtime_error(message); }

static void box_sliver() {
    auto example=analytical_field("strip");
    example.field.box.min.z()=0;
    example.field.value=[](const Vec& x) { return Range(x.x(),x.y()-x.z()); };
    example.field.jacobian=[](const Vec&) { Jacobian j; j << 1,0,0, 0,1,-1; return j; };
    Options options;
    options.min_angle_degrees=5; // The narrow gap stays open until the box-face repair.
    const unsigned lower_z=1u<<4;
    auto finish=[&](double height,const Vec& apex) {
        Checkpoint state;
        auto& r=state.result;
        r.mesh.vertices={Vec(.2,0,0),Vec(.25,height,height),Vec(.3,0,0),apex};
        r.info.resize(4);
        r.info[0].faces=r.info[2].faces=lower_z;
        r.mesh.add({0,3,1});
        r.mesh.add({1,3,2});
        return Extractor(example.field,example.segment,options).finish(std::move(state)).result;
    };
    const Vec apex(.25,.07,.07);
    // In the plane y=z the Euclidean distance to the chord is sqrt(2)*height.
    const double cutoff_height=0.01*options.step/std::sqrt(2.);
    const auto snapped=finish(.99*cutoff_height,apex);
    require(snapped.mesh.vertices==std::vector<Vec>{Vec(.2,0,0),Vec(.25,0,0),Vec(.3,0,0),apex} &&
            snapped.mesh.triangles==std::vector<Triangle>{{0,3,1},{1,3,2}},
            "Correct the sliver tip on both the box face and the side surface, preserving connectivity");
    require(snapped.info[1].faces==lower_z && snapped.box_edges.count({0,1}) && snapped.box_edges.count({1,2}),
            "Both sides of the tip must be available to the cap");
    const double above=1.01*cutoff_height;
    const auto interior=finish(above,apex);
    require(interior.mesh.vertices[1]==Vec(.25,above,above) && !interior.info[1].faces &&
            !interior.box_edges.count({0,1}) && !interior.box_edges.count({1,2}),
            "A tip beyond one percent of the nominal step must retain its position and edge labels");
    // Moving the tip onto this apex would collapse both existing facets.
    const auto blocked=finish(1e-5,Vec(.25,0,0));
    require(blocked.mesh.vertices[1]==Vec(.25,1e-5,1e-5) && !blocked.info[1].faces &&
            !blocked.box_edges.count({0,1}) && !blocked.box_edges.count({1,2}),
            "A rejected move must restore the vertex and leave both edges unresolved");
}

static void contacts() {
    // Two parallel faces lie in disjoint x intervals within
    // the predicate's contact tolerance. A third face is outside it.
    std::vector<Vec> nearby;
    for(double x:{0.,0.75e-9,4e-9})
        for(const auto& yz:std::array<Range,3>{Range(0,0),Range(1,0),Range(0,1)})
            nearby.emplace_back(x,yz.x(),yz.y());
    const std::vector<Triangle> planes{{6,7,8},{3,4,5},{0,1,2}};
    require(marching_triangles::validate_mesh(nearby,planes,1e-9).intersections==1,
            "Intersection bounds must include near-coplanar contacts within epsilon");
    auto check=[](const std::vector<Vec>& vertices,Triangle a,Triangle b,bool expected) {
        std::sort(a.begin(),a.end());
        do {
            std::sort(b.begin(),b.end());
            do {
                require(marching_triangles::triangles_conflict(vertices,a,b,2e-9)==expected,
                        "Triangle contact classification depends on vertex order or misses an intersection");
                require(marching_triangles::triangles_conflict(vertices,b,a,2e-9)==expected,
                        "Triangle contact classification depends on argument order");
            } while(std::next_permutation(b.begin(),b.end()));
        } while(std::next_permutation(a.begin(),a.end()));
    };
    // Actual disconnected-volume facets. Each pair has a separating plane
    // through its single shared vertex; there is no additional intersection.
    const std::vector<Vec> vertices{
        Vec(.66000537067963705,.092803544662036164,-.34999999996676262),
        Vec(.61306546313914601,.062924631044974588,-.42500000000000043),
        Vec(.6130654630691017,.062924631002032758,-.27500000000000041),
        Vec(.66000535895896895,.092803536926332308,-.19999999999172802),
        Vec(.76554159169867964,.16802696430942007,-.27500000000000041),
        Vec(.6130654630517568,.062924630991399194,-.12500000000000039)};
    check(vertices,{0,1,2},{3,4,0},false);
    check(vertices,{3,2,5},{3,4,0},false);
    // Sharing a vertex must not hide a real intersection beyond that vertex.
    check({Vec(0,0,0),Vec(1,0,0),Vec(0,1,0),Vec(.25,.25,-1),Vec(.25,.25,1)},
          {0,1,2},{0,3,4},true);
    check({Vec(0,0,0),Vec(1,0,0),Vec(0,1,0),Vec(.6,.1,0),Vec(.1,.6,0)},
          {0,1,2},{0,3,4},true);
}

static void seeds() {
    auto e=analytical_field("strip"); Options options;
    e.field.box.min=Vec(-.5,-1,-2); e.field.box.max=Vec(.5,1,2);
    require(discovery_shape(e.field,options)==SeedShape{7,14,28},"Proportional default seed counts");
    options.seed_grid=0;
    require(discovery_shape(e.field,options)==SeedShape{7,14,28},"Zero selects the default shortest-axis count");
    options.seed_grid=3;
    const auto shape=discovery_shape(e.field,options);
    require(shape==SeedShape{3,6,12},"Shortest-axis seed count override");
    const auto grid=discovery_grid(e.field.box,shape);
    require(grid.size()==216 && grid.front()==e.field.box.min && grid.back()==e.field.box.max,"Anisotropic grid bounds/count");
    require(discovery_grid(e.field.box,SeedShape{1,1,1}).front()==Vec::Zero(),"Singleton axes use midpoint");
    e.field.box.max[2]=1.1;
    require(discovery_shape(e.field,options)==SeedShape{3,6,10},"Noninteger aspect ratio rounds up");
    e.field.box.max[2]=2;
    Corrector corrector(e.field,e.segment,options);
    // A linear strip reaches its fiber in one update, even for a step across the box.
    Options one_update=options;
    one_update.iterations=2; // Initial point and the point after one Newton step.
    Corrector full_step(e.field,e.segment,one_update);
    Vec corner(-.5,1,0);
    require(full_step.fiber_seed(corner,Range(.5,-1)) && corner==Vec(.5,-1,0),
            "Seed discovery must take the full in-domain Newton step");
    Seeds reference;
    for (int endpoint=0; endpoint<2; ++endpoint)
        for (Vec p:grid) if (corrector.fiber_seed(p,endpoint ? e.segment.b : e.segment.a)) reference.endpoints[endpoint].push_back(p);
    for (Vec p:grid) if (corrector.surface(p,(e.field.box.max-e.field.box.min).norm())) reference.surface.push_back(p);
    for (int threads:{1,4}) {
        tbb::global_control control(tbb::global_control::max_allowed_parallelism,threads);
        const auto actual=find_seeds(e.field,e.segment,options);
        require(actual.grid_shape==shape && actual.starting_points==grid.size(),"Saved seed grid metadata");
        require(actual.surface==reference.surface && actual.endpoints==reference.endpoints,"Parallel corrections must preserve serial values and ordering");
    }
    options.seed_grid=1000000;
    bool rejected=false;
    try { discovery_shape(e.field,options); } catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"Excessive grid allocation must be rejected");
    e.field.box={Vec::Zero(),Vec(499,499,99)};
    options.seed_grid=40;
    rejected=false;
    try { discovery_shape(e.field,options); } catch (const std::invalid_argument&) { rejected=true; }
    require(rejected,"Isabel seed grid needs an explicit larger allocation limit");
    options.max_seed_starts=1632160;
    require(discovery_shape(e.field,options)==SeedShape{202,202,40},"Full Isabel discovery grid");
}

static void polylines() {
    const Field field=analytical_field("strip").field;
    const Range a(0,0), b(1,0), c(0,1), d(1,1);
    const Options options;
    for (const auto& points:std::vector<std::vector<Range>>{{b,a,c},{a,b,c,a},{a,c,d,b,a}}) {
        const RangeCurve curve(points);
        const Result surface=extract_curve(field,curve,options);
        const auto report=evaluate(surface,field,curve,options);
        require(report.complete(surface) && report.components==1,"Disconnected or incomplete curve surface");
        require(report.topology.euler_characteristic==(curve.closed()?0:1),"Incorrect open/closed surface topology");
        require(surface.fibers.size()==size_t(curve.vertex_count()),"Duplicate or missing corner fibers");
        double length=0;
        for (const auto& segment:curve.segments) length+=segment.length;
        require(std::abs(report.area-4*length)<1e-9,"Expected one full-height strip per segment");
        for (const auto& [edge,component]:surface.fiber_edges) {
            const int k=surface.fibers[component].endpoint;
            const auto& use=surface.mesh.edges.at(edge);
            const bool crease=curve.closed() || (k>0 && k<curve.size());
            require((use.second>=0)==crease,"Incorrect number of triangles at a shared fiber");
            if (crease) {
                const int left=(k+curve.size()-1)%curve.size();
                require(std::set<int>{surface.face_segment(use.face),surface.face_segment(use.second)}==
                        std::set<int>{left,k},"Incorrect patches joined at crease");
            }
        }
    }
    const RangeCurve simplified({a,a,Range(.5,0),b,c,a});
    require(simplified.closed() && simplified.size()==3,"Failed to remove redundant vertices");
    for (const auto& invalid:std::vector<std::vector<Range>>{{a,a},{a,b,a},{a,b,Range(.5,0),c},{a,d,b,c},{a,b,c,a,d}}) {
        bool rejected=false;
        try { RangeCurve curve(invalid); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"Degenerate or self-intersecting range curve accepted");
    }
}

static void segment_sampling() {
    const auto strip=analytical_field("strip");
    Options options;
    options.residual_tolerance=.005;
    Corrector corrector(strip.field,strip.segment,options);
    Samples samples;
    for (double tau : {std::nextafter(0.,-1.), 0., .5, 1., std::nextafter(1.,2.)}) {
        const Vec x(tau,0,0);
        const bool inside=0<=tau && tau<=1;
        require(corrector.member(x)==inside,"Exact closed segment interval in marching");
        const auto before=samples.violations;
        samples.add(x,strip.field,strip.segment);
        require(samples.violations-before==size_t(!inside),"Exact closed segment interval in reporting");
    }
    const Segment oblique(Range(0,0),Range(1,58));
    for (const int endpoint : {0,1}) {
        const double sign=endpoint==0 ? -1. : 1.;
        const Vec near(endpoint+sign*5e-10,0,0);
        require(corrector.member(near,endpoint) && !corrector.member(near),
                "Only a traced vertex may use the endpoint field residual");
        require(!corrector.member(Vec(endpoint+sign*2e-9,0,0),endpoint),
                "Traced vertices must satisfy the fiber field tolerance");
    }
    require(oblique.parameter(oblique.b)==1,"The exact far endpoint must have parameter one");

    // On z=y² the seed triangle's chord samples lie off the surface, but x stays in range.
    Field field=strip.field;
    field.value=[](const Vec& p) { return Range(p.x(),p.z()-p.y()*p.y()); };
    field.jacobian=[](const Vec& p) { Jacobian j; j<<1,0,0,0,-2*p.y(),1; return j; };
    field.hessians=[](const Vec&) {
        Hessian h=Hessian::Zero(); h(1,1)=-2;
        return std::array<Hessian,2>{Hessian::Zero(),h};
    };
    options.spatial_tolerance=1e-12;
    options.residual_tolerance=1e-9;
    options.max_triangles=4;
    Seeds seeds;
    seeds.surface.emplace_back(.5,0,0);
    const Result patch=Extractor(field,strip.segment,options).march(seeds,{}).result;
    const Report report=evaluate(patch,field,strip.segment,options);
    require(!patch.mesh.triangles.empty() && report.midpoints.maximum>1e-4 && report.centroids.maximum>1e-4,
            "Curved patch must accept in-range chord samples without projecting them");
    require(!report.vertices.violations && !report.midpoints.violations && !report.centroids.violations,
            "All accepted triangle samples must be in range");

    // In the plane z=0, exclude a small disk around the origin. Position the
    // equilateral seed triangle so only one midpoint or only its centroid is in the disk.
    field=analytical_field("sphere").field;
    options.min_step=options.step; // Try exactly one triangle size.
    options.residual_tolerance=.005;
    const Segment annulus(Range(.001,0),Range(1,0));
    for (const double x : {0., std::sqrt(3.)*options.step/6}) {
        seeds.surface={Vec(x,-options.step/2,0)};
        const auto state=Extractor(field,annulus,options).march(seeds,{});
        if (x==0) {
            require(!state.result.mesh.triangles.empty(),"An out-of-range midpoint alone must permit the seed triangle");
            const auto sampled=evaluate(state.result,field,annulus,options);
            require(sampled.midpoints.violations && !sampled.vertices.violations && !sampled.centroids.violations,
                    "Midpoint violations should be measured on accepted triangles");
        } else {
            require(state.result.mesh.triangles.empty(),"An out-of-range centroid must reject the seed triangle");
        }
    }
}

static void fiber_tolerances() {
    const auto strip=analytical_field("strip");
    const Range q=Range::Zero();
    const Vec normal=Vec::UnitZ(), anchor=Vec::Zero();
    Options options;
    options.iterations=1; // Observe the acceptance test before a second Newton iterate.
    Corrector corrector(strip.field,strip.segment,options);
    for (double surface_tolerance : {1e-12, .005}) {
        options.residual_tolerance=surface_tolerance;
        for (double error : {1e-9, 1.1e-9}) {
            const bool accept=error<=1e-9;
            Vec seed(error,0,0), sample=seed, contact(error,0,strip.field.box.max.z());
            require(corrector.fiber_seed(seed,q)==accept,"Fiber seed field tolerance");
            require(corrector.fiber(sample,q,normal,anchor,1)==accept,"Interior fiber field tolerance");
            require(corrector.fiber(contact,q,normal,anchor,1,1u<<5)==accept,"Fiber contact field tolerance");
        }
    }
    options.corner_residual_tolerance=1e-11;
    Vec off_corner(5e-10,0,0);
    require(!corrector.fiber_seed(off_corner,q),"The staged corner tolerance must govern seed acceptance");
    for (double h : {.1, .5}) {
        options.step=h;
        for (double distance : {1e-6*h, 1.1e-6*h}) {
            Vec sample(0,0,distance);
            require(corrector.fiber(sample,q,normal,anchor,1)==(distance<=1e-6*h),"Fiber plane tolerance must scale with the trace step");
        }
    }

    // On the unit equator, the sixth step returns within h with tangent dot product about 0.72.
    const auto sphere=analytical_field("sphere");
    const Segment segment(Range(1,0),Range(1,.5));
    Options coarse;
    coarse.step=.925;
    coarse.max_trace_steps=6;
    Seeds seeds;
    seeds.endpoints[0].push_back(Vec::UnitX());
    const Result loop=trace_endpoint_fibers(sphere.field,segment,coarse,seeds);
    require(loop.unresolved.empty() && loop.fibers.size()==1 && loop.fibers[0].closed &&
            loop.fibers[0].vertices.size()==6,"A fiber must close on its first eligible return to the seed");

    // The unit equator is regular, but a full RK4 midpoint from (1,0,0)
    // lands exactly on the rank-deficient cylinder r^2=1+h^2/4.
    Field cylinder=sphere.field;
    constexpr double h=.5, critical_radius2=1+h*h/4;
    cylinder.value=[](const Vec& p) {
        const double r2=p.head<2>().squaredNorm();
        return Range(r2,p.z()*(r2-critical_radius2));
    };
    cylinder.jacobian=[](const Vec& p) {
        Jacobian j;
        j << 2*p.x(),2*p.y(),0,2*p.x()*p.z(),2*p.y()*p.z(),p.head<2>().squaredNorm()-critical_radius2;
        return j;
    };
    Options retry;
    retry.step=retry.min_step=h;
    const Result rejected=trace_endpoint_fibers(cylinder,segment,retry,seeds);
    require(!rejected.unresolved.empty() && rejected.fibers.empty(),"Full RK4 step should hit the singular cylinder");
    retry.min_step=h/8;
    const Result recovered=trace_endpoint_fibers(cylinder,segment,retry,seeds);
    require(recovered.unresolved.empty() && recovered.fibers.size()==1 && recovered.fibers[0].closed,
            "Subdividing rejected steps must recover the regular closed fiber");
    const auto& ids=recovered.fibers[0].vertices;
    for (size_t i=1; i<ids.size(); ++i) {
        const Vec& p=recovered.mesh.vertices[ids[i]];
        const double length=(p-recovered.mesh.vertices[ids[i-1]]).norm();
        require((cylinder.value(p)-segment.a).norm()<=Options::fiber_field_tolerance &&
                length>.75*h && length<1.05*h,"Only nominal-spaced fiber samples should be emitted");
    }
}

static void fiber_seed_coverage() {
    const Field field=analytical_field("disconnected").field;
    const double gap=1./64;
    const Vec left(-gap/2,0,0), right(gap/2,0,.125), other(-gap/2,gap/2,0);
    const Segment segment(field.value(left),field.value(other));
    Seeds seeds;
    seeds.endpoints[0]={left,right};
    seeds.endpoints[1]={other};
    Options options;
    options.step=.25; // The closest point to right is halfway along an edge of the left trace.
    for (double tolerance : {gap, gap/2}) {
        options.spatial_tolerance=tolerance;
        const Result result=trace_endpoint_fibers(field,segment,options,seeds);
        const size_t expected=tolerance==gap ? 2 : 3;
        require(result.unresolved.empty() && result.fibers.size()==expected,
                "Distance to the segment must cover a seed up to and including the spatial tolerance");
        require(result.fibers.back().endpoint==1,"A nearby trace for a different range point must not cover a seed");
    }
}

static void stalled_fiber_recovery() {
    const Vec origin=Vec::Zero(), unit=Vec::UnitX();
    const std::vector<Vec> reference={origin,.4*unit,unit};
    require(agrees_over_arclength({-.2*unit,.7*unit,unit},0,origin,reference,1,.01),
            "Exactly one step of shared arclength must qualify with different sample spacing");
    require(!agrees_over_arclength({origin,.999*unit},0,origin,reference,1,.01) &&
            !agrees_over_arclength(reference,0,origin,{origin,.999*unit},1,.01),
            "A short overlap must not qualify in either direction");
    require(!agrees_over_arclength({origin,Vec(.5,.02,0),unit},0,origin,reference,1,.01),
            "The entire overlap must agree, including a deviation between aligned endpoints");
    // At z=0 the fiber is the unit circle. Away from it the tangent field
    // turns sharply enough for a nominal RK4 advance to return to its start.
    Field field=analytical_field("sphere").field;
    const double frequency=2*std::acos(-1.)/.005;
    field.value=[frequency](const Vec& p) {
        const double r2=p.head<2>().squaredNorm();
        return Range(r2,std::sin(frequency*p.z())/frequency*(1+200*p.y()*(r2-1)));
    };
    field.jacobian=[frequency](const Vec& p) {
        const double r2=p.head<2>().squaredNorm();
        const double z=std::sin(frequency*p.z())/frequency;
        Jacobian j;
        j << 2*p.x(),2*p.y(),0, 400*p.x()*p.y()*z,200*(r2-1+2*p.y()*p.y())*z,
             std::cos(frequency*p.z())*(1+200*p.y()*(r2-1));
        return j;
    };
    Options options; options.step=.4; options.min_step=options.step/8;
    options.fiber_seed_radius=options.step; options.max_trace_steps=400;
    const Segment segment(Range(1,0),Range(2,0));
    const Vec seed(1,0,0), later(std::cos(.3),-std::sin(.3),0);
    Seeds seeds; seeds.endpoints[0]={seed};
    auto stalled=trace_endpoint_fibers(field,segment,options,seeds);
    require(!stalled.unresolved.empty() && stalled.unresolved.front().reason=="Fiber trace made insufficient progress",
            "Fixture must stall at the nominal step");
    seeds.endpoints[0].push_back(later+Vec(0,0,.005));
    auto other=trace_endpoint_fibers(field,segment,options,seeds);
    for (const auto& event:stalled.unresolved)
        require(std::any_of(other.unresolved.begin(),other.unresolved.end(),[&](const Event& e) {
                    return e.position==event.position && e.reason==event.reason;
                }),"A nearby separate circle cannot resolve the original stall");
    seeds.endpoints[0]={seed,later};
    auto fixed=trace_endpoint_fibers(field,segment,options,seeds);
    require(fixed.unresolved.empty() && fixed.fibers.size()==1 && fixed.fibers.front().closed,
            "A later seed on the stalled circle must produce one complete fiber");
    double length=0;
    for (const auto& e:fixed.fiber_edges) {
        const auto [a,b]=e.first;
        require(std::abs(fixed.mesh.vertices[a].z())<1e-9,"Recovery switched between close fibers");
        length+=(fixed.mesh.vertices[a]-fixed.mesh.vertices[b]).norm();
    }
    require(std::abs(length-2*std::acos(-1.))<.01,"Recovery must keep one complete lap");
}

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
    Vec near_axis(.1,0,0); // The full Newton step toward the unit equator lands at x=5.05.
    require(!bounded.fiber_seed(near_axis,Range(1,0)),"An out-of-box Newton step must reject the seed candidate");
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
            "A loose surface root tolerance must not create duplicate endpoint components");
}
static void parallel_march() {
    auto example=analytical_field("strip");
    Options options;
    const auto seeds=find_seeds(example.field,example.segment,options);
    const auto fibers=trace_endpoint_fibers(example.field,example.segment,options,seeds);
    const auto serial=Extractor(example.field,example.segment,options).march(seeds,fibers);
    // Exercise the same lookahead as the CoordNet adapter with nonempty
    // endpoint fibers, including edges clipped to the domain boundary.
    example.field.value_jacobian=[field=example.field](const Vec& x) {
        return FirstOrderSample{field.value(x),field.jacobian(x)};
    };
    tbb::task_arena arena(4);
    const auto parallel=arena.execute([&] {
        return Extractor(example.field,example.segment,options).march(seeds,fibers);
    });
    require(!fibers.fibers.empty() && !serial.result.box_edges.empty(),"Expected fibers and box contacts");
    require(serial.result.mesh.vertices==parallel.result.mesh.vertices &&
            serial.result.mesh.triangles==parallel.result.mesh.triangles &&
            serial.result.box_edges==parallel.result.box_edges && serial.pending==parallel.pending &&
            serial.fiber_uses==parallel.fiber_uses,"Parallel edge proposals changed the march");
}

int main(int argc,char** argv) {
    try {
        const std::string name=argc>1?argv[1]:"sphere";
        if (name=="derivatives") { derivatives(); return 0; }
        if (name=="segment-sampling") { segment_sampling(); return 0; }
        if (name=="parallel-march") { parallel_march(); return 0; }
        if (name=="fiber-tolerances") { fiber_tolerances(); return 0; }
        if (name=="seed-coverage") { fiber_seed_coverage(); stalled_fiber_recovery(); return 0; }
        if (name=="seeds") { seeds(); return 0; }
        if (name=="polylines") { polylines(); return 0; }
        if (name=="contacts") { contacts(); return 0; }
        if (name=="box-sliver") { box_sliver(); return 0; }
        auto example=analytical_field(name); Options options;
        const Result r=Extractor(example.field,example.segment,options).run();
        const auto report=evaluate(r,example.field,example.segment,options);
        write_outputs("test_"+name,r,report,example.field,example.segment,options);
        std::cout<<name<<": triangles="<<r.mesh.triangles.size()<<" components="<<report.components<<" unresolved="<<r.unresolved.size()<<" "<<report.topology.summary()<<std::endl;
        require(report.complete(r),"Extraction is incomplete");
        require(report.vertices.maximum<options.residual_tolerance*1.01,"Surface vertex residual");
        require(report.endpoint_residual<=Options::fiber_field_tolerance,"Endpoint fiber residual");
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
