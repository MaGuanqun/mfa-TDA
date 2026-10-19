#include "pipeline.h"
#include "corrector.h"
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace fiber_surface {
namespace {
constexpr int format_version=1;
fs::path stage_artifact(const fs::path& run,Stage stage,const std::string& name) {
    static const std::array<const char*,5> prefixes{"seeds","endpoint_fibers","marching","stitching","final"};
    return run/(std::string(prefixes.at(static_cast<size_t>(stage)-1))+"_"+name);
}
std::string digest(const std::string& text) {
    uint64_t hash=14695981039346656037ull;
    for (unsigned char c:text) { hash^=c; hash*=1099511628211ull; }
    std::ostringstream out; out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash; return out.str();
}
void atomic_write(const fs::path& path, const std::function<void(std::ostream&)>& write) {
    fs::create_directories(path.parent_path());
    const fs::path temporary=path.string()+".tmp";
    try {
        std::ofstream out(temporary); if (!out) throw std::runtime_error("Cannot create "+temporary.string());
        out<<std::setprecision(std::numeric_limits<double>::max_digits10);
        write(out); out.close(); if (!out) throw std::runtime_error("Failed writing "+temporary.string());
        fs::rename(temporary,path);
    } catch (...) { std::error_code error; fs::remove(temporary,error); throw; }
}
void config_body(std::ostream& out, const RunConfig& c) {
    out<<"FIBER_SURFACE_CONFIG "<<format_version<<"\nmodel "<<std::quoted(c.model)<<"\nexample "<<std::quoted(c.example)
        <<"\nfirst "<<std::quoted(c.first)<<' '<<std::quoted(c.first_digest)
        <<"\nsecond "<<std::quoted(c.second)<<' '<<std::quoted(c.second_digest)
        <<"\nvariables "<<c.variable1<<' '<<c.variable2<<"\nblocks "<<c.block1<<' '<<c.block2
        <<"\nbox_min "<<c.box.min.transpose()<<"\nbox_max "<<c.box.max.transpose()
        <<"\nsegment_a "<<c.segment.a.transpose()<<"\nsegment_b "<<c.segment.b.transpose();
    const auto& o=c.options;
    out<<"\noptions "<<o.step<<' '<<o.min_step<<' '<<o.fiber_step<<' '<<o.spatial_tolerance<<' '
        <<o.residual_tolerance<<' '<<o.derivative_tolerance<<' '<<o.min_angle_degrees<<' '<<o.max_normal_degrees<<' '
        <<o.seed_grid<<' '<<o.iterations<<' '<<o.max_trace_steps<<' '<<o.max_triangles<<"\nEND\n";
}
class Reader {
    fs::path path_;
    std::ifstream in_;
public:
    explicit Reader(fs::path path) : path_(std::move(path)),in_(path_) {
        if (!in_) throw std::runtime_error("Missing/unreadable stage input: "+path_.string());
    }
    [[noreturn]] void fail(const std::string& why) const { throw std::runtime_error("Invalid "+path_.string()+": "+why); }
    template<class T> T get() { T x{}; if (!(in_>>x)) fail("truncated or malformed value"); return x; }
    double real() { double x=get<double>(); if (!std::isfinite(x)) fail("nonfinite number"); return x; }
    std::string quoted() { std::string x; if (!(in_>>std::quoted(x))) fail("malformed string"); return x; }
    void word(const std::string& expected) { if (get<std::string>()!=expected) fail("expected "+expected); }
    void header(const std::string& magic) { word(magic); if (get<int>()!=format_version) fail("unsupported format version"); }
    size_t count() {
        const auto n=get<long long>();
        if (n<0 || n>10000000) fail("invalid count (maximum 10000000)");
        return static_cast<size_t>(n);
    }
    Vec point() { Vec p; for (int d=0; d<3; ++d) p[d]=real(); return p; }
    Range range() { Range p; for (int d=0; d<2; ++d) p[d]=real(); return p; }
    int index(size_t n) { const int id=get<int>(); if (id<0 || size_t(id)>=n) fail("out-of-range vertex/component ID"); return id; }
    bool boolean() { const int x=get<int>(); if (x!=0 && x!=1) fail("invalid boolean"); return x!=0; }
    void end() { word("END"); std::string extra; if (in_>>extra) fail("unexpected trailing data"); }
};
std::string dependencies(const fs::path& run, Stage stage) {
    if (stage==Stage::fibers)
        return digest(file_digest(run/"surface_seeds.dat")+file_digest(run/"endpoint_seeds.dat"));
    if (stage==Stage::marching)
        return digest(file_digest(state_path(run,Stage::fibers))+file_digest(run/"surface_seeds.dat"));
    if (stage==Stage::stitching) return file_digest(state_path(run,Stage::marching));
    throw std::logic_error("This stage has no checkpoint dependencies");
}
void points_header(std::ostream& out,const std::vector<Vec>& points) {
    out<<"# vtk DataFile Version 3.0\nFiber surface intermediate stage\nASCII\nDATASET POLYDATA\nPOINTS "<<points.size()<<" double\n";
    for (const auto& p:points) out<<p.transpose()<<'\n';
}
std::vector<Edge> unfinished(const Checkpoint& state) {
    std::vector<Edge> edges;
    for (const auto& e:state.pending) edges.push_back(e.first);
    for (const auto& e:state.result.mesh.boundary) {
        const Edge key=edge_key(e.first.first,e.first.second);
        if (!state.result.fiber_edges.count(key) && !state.result.box_edges.count(key)) edges.push_back(e.first);
    }
    return edges;
}
void events(std::ostream& out,const std::vector<Event>& unresolved) {
    out<<'[';
    for (size_t i=0; i<unresolved.size(); ++i) {
        const auto& e=unresolved[i]; if (i) out<<',';
        out<<"\n    {\"position\": ["<<e.position[0]<<','<<e.position[1]<<','<<e.position[2]<<"], \"reason\": "<<std::quoted(e.reason)<<'}';
    }
    out<<"\n  ]";
}
}

std::string file_digest(const fs::path& path) {
    std::ifstream in(path,std::ios::binary); if (!in) throw std::runtime_error("Cannot read "+path.string());
    uint64_t hash=14695981039346656037ull; char buffer[65536];
    while (in) {
        in.read(buffer,sizeof(buffer));
        for (std::streamsize i=0; i<in.gcount(); ++i) { hash^=static_cast<unsigned char>(buffer[i]); hash*=1099511628211ull; }
    }
    if (!in.eof()) throw std::runtime_error("Failed reading "+path.string());
    std::ostringstream out; out<<std::hex<<std::setw(16)<<std::setfill('0')<<hash; return out.str();
}
std::string config_signature(const RunConfig& config) {
    std::ostringstream out; out<<std::setprecision(17); config_body(out,config); return digest(out.str());
}
void write_config(const fs::path& run,const RunConfig& config) {
    atomic_write(run/"run.conf",[&](std::ostream& out) { config_body(out,config); });
}
RunConfig read_config(const fs::path& run) {
    Reader in(run/"run.conf"); RunConfig c;
    in.header("FIBER_SURFACE_CONFIG");
    in.word("model"); c.model=in.quoted(); in.word("example"); c.example=in.quoted();
    in.word("first"); c.first=in.quoted(); c.first_digest=in.quoted();
    in.word("second"); c.second=in.quoted(); c.second_digest=in.quoted();
    in.word("variables"); c.variable1=in.get<int>(); c.variable2=in.get<int>();
    in.word("blocks"); c.block1=in.get<int>(); c.block2=in.get<int>();
    in.word("box_min"); c.box.min=in.point(); in.word("box_max"); c.box.max=in.point();
    in.word("segment_a"); const Range a=in.range(); in.word("segment_b"); const Range b=in.range(); c.segment=Segment(a,b);
    in.word("options"); auto& o=c.options;
    o.step=in.real(); o.min_step=in.real(); o.fiber_step=in.real(); o.spatial_tolerance=in.real();
    o.residual_tolerance=in.real(); o.derivative_tolerance=in.real(); o.min_angle_degrees=in.real(); o.max_normal_degrees=in.real();
    o.seed_grid=in.get<int>(); o.iterations=in.get<int>(); o.max_trace_steps=in.get<int>(); o.max_triangles=in.count();
    in.end(); return c;
}
fs::path state_path(const fs::path& run,Stage stage) {
    const char* name=stage==Stage::fibers ? "fibers.state" : stage==Stage::marching ? "marching.state" : "stitched.state";
    return run/name;
}
void write_seeds(const fs::path& run,const RunConfig& config,const Seeds& seeds) {
    const auto signature=config_signature(config);
    atomic_write(run/"surface_seeds.dat",[&](std::ostream& out) {
        out<<"FIBER_SURFACE_SEEDS 1 surface "<<signature<<' '<<seeds.seconds<<' '<<seeds.starting_points<<' '<<seeds.surface.size()<<'\n';
        for (size_t i=0; i<seeds.surface.size(); ++i) out<<i<<' '<<seeds.surface[i].transpose()<<'\n';
        out<<"END\n";
    });
    atomic_write(run/"endpoint_seeds.dat",[&](std::ostream& out) {
        out<<"FIBER_SURFACE_SEEDS 1 endpoints "<<signature<<' '<<seeds.seconds<<' '<<seeds.starting_points<<' '
            <<seeds.endpoints[0].size()+seeds.endpoints[1].size()<<'\n';
        size_t id=0;
        for (int endpoint=0; endpoint<2; ++endpoint) for (const auto& p:seeds.endpoints[endpoint]) out<<id++<<' '<<endpoint<<' '<<p.transpose()<<'\n';
        out<<"END\n";
    });
    for (bool endpoints:{false,true}) atomic_write(run/(endpoints?"endpoint_seeds.vtk":"surface_seeds.vtk"),[&](std::ostream& out) {
        std::vector<Vec> points=seeds.surface; std::vector<int> labels(points.size(),-1);
        if (endpoints) { points.clear(); labels.clear(); for (int e=0; e<2; ++e) for (const auto& p:seeds.endpoints[e]) { points.push_back(p); labels.push_back(e); } }
        points_header(out,points); out<<"VERTICES "<<points.size()<<' '<<points.size()*2<<'\n';
        for (size_t i=0; i<points.size(); ++i) out<<"1 "<<i<<'\n';
        out<<"POINT_DATA "<<points.size()<<"\nSCALARS endpoint int 1\nLOOKUP_TABLE default\n";
        for (int label:labels) out<<label<<'\n';
    });
    atomic_write(stage_artifact(run,Stage::seeds,"report.json"),[&](std::ostream& out) {
        out<<"{\n  \"stage\": \"seeds\",\n  \"starting_points\": "<<seeds.starting_points
            <<",\n  \"surface_candidates\": "<<seeds.surface.size()<<",\n  \"endpoint_candidates\": ["
            <<seeds.endpoints[0].size()<<','<<seeds.endpoints[1].size()<<"],\n  \"stage_seconds\": "<<seeds.seconds<<"\n}\n";
    });
}
Seeds read_seeds(const fs::path& run,const RunConfig& config) {
    Seeds seeds; const auto signature=config_signature(config);
    for (bool endpoints:{false,true}) {
        Reader in(run/(endpoints?"endpoint_seeds.dat":"surface_seeds.dat"));
        in.header("FIBER_SURFACE_SEEDS"); in.word(endpoints?"endpoints":"surface"); in.word(signature);
        const double seconds=in.real(); const size_t starts=in.count(),count=in.count();
        if (seconds<0) in.fail("negative duration");
        if (!endpoints) { seeds.seconds=seconds; seeds.starting_points=starts; }
        else if (seeds.seconds!=seconds || seeds.starting_points!=starts) in.fail("seed files are from different discoveries");
        for (size_t i=0; i<count; ++i) {
            if (in.count()!=i) in.fail("nonsequential seed ID");
            const int endpoint=endpoints?in.index(2):-1; const Vec p=in.point();
            if (!config.box.contains(p)) in.fail("seed outside saved box");
            if (endpoints) seeds.endpoints[endpoint].push_back(p); else seeds.surface.push_back(p);
        }
        in.end();
    }
    return seeds;
}
void validate_seeds(const Seeds& seeds,const Field& field,const Segment& segment,const Options& options) {
    Corrector corrector(field,segment,options);
    for (const auto& p:seeds.surface)
        if (!corrector.member(p) || std::abs(segment.phi(field.value(p)))>options.residual_tolerance*1.01 ||
            segment.gradient(field.jacobian(p)).norm()<=options.derivative_tolerance) throw std::runtime_error("Invalid saved surface seed");
    for (int endpoint=0; endpoint<2; ++endpoint) for (const auto& p:seeds.endpoints[endpoint])
        if (!field.box.contains(p) || !corrector.regular(field.jacobian(p)) ||
            (field.value(p)-(endpoint?segment.b:segment.a)).norm()>options.residual_tolerance*1.01)
            throw std::runtime_error("Invalid saved endpoint seed");
}

void write_checkpoint(const fs::path& run,const RunConfig& config,Stage stage,const Checkpoint& state) {
    const auto& r=state.result;
    atomic_write(state_path(run,stage),[&](std::ostream& out) {
        out<<"FIBER_SURFACE_STATE 1 "<<int(stage)<<' '<<config_signature(config)<<' '<<dependencies(run,stage)<<'\n';
        out<<"vertices "<<r.mesh.vertices.size()<<'\n';
        for (size_t i=0; i<r.mesh.vertices.size(); ++i) out<<r.mesh.vertices[i].transpose()<<' '<<r.info[i].endpoint<<' '<<r.info[i].fiber<<' '<<r.info[i].faces<<'\n';
        out<<"triangles "<<r.mesh.triangles.size()<<'\n';
        for (const auto& t:r.mesh.triangles) out<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
        out<<"fibers "<<r.fibers.size()<<'\n';
        for (const auto& f:r.fibers) { out<<f.endpoint<<' '<<int(f.closed)<<' '<<int(f.complete)<<' '<<f.vertices.size(); for (int id:f.vertices) out<<' '<<id; out<<'\n'; }
        for (const auto& entry:std::array<std::pair<const char*,const std::map<Edge,int>*>,3>{{{"fiber_edges",&r.fiber_edges},{"pending",&state.pending},{"fiber_uses",&state.fiber_uses}}}) {
            out<<entry.first<<' '<<entry.second->size()<<'\n';
            for (const auto& e:*entry.second) out<<e.first.first<<' '<<e.first.second<<' '<<e.second<<'\n';
        }
        out<<"box_edges "<<r.box_edges.size()<<'\n'; for (const auto& e:r.box_edges) out<<e.first<<' '<<e.second<<'\n';
        out<<"events "<<r.unresolved.size()<<'\n'; for (const auto& e:r.unresolved) out<<e.position.transpose()<<' '<<std::quoted(e.reason)<<'\n';
        out<<"discovery "<<r.fiber_seed_hits[0]<<' '<<r.fiber_seed_hits[1]<<' '<<r.surface_seed_hits<<"\nseconds "<<r.seconds<<"\nEND\n";
    });
}
Checkpoint read_checkpoint(const fs::path& run,const RunConfig& config,Stage stage) {
    // Check the whole dependency chain, not just the immediately preceding file.
    if (stage==Stage::marching) read_checkpoint(run,config,Stage::fibers);
    if (stage==Stage::stitching) read_checkpoint(run,config,Stage::marching);
    Reader in(state_path(run,stage)); in.header("FIBER_SURFACE_STATE");
    if (in.get<int>()!=int(stage)) in.fail("wrong pipeline stage");
    if (in.get<std::string>()!=config_signature(config)) in.fail("configuration mismatch; rerun this stage");
    if (in.get<std::string>()!=dependencies(run,stage)) in.fail("upstream output changed; rerun this stage");
    Checkpoint state; auto& r=state.result;
    in.word("vertices"); const size_t nv=in.count();
    for (size_t i=0; i<nv; ++i) {
        const Vec p=in.point(); VertexInfo info; info.endpoint=in.get<int>(); info.fiber=in.get<int>(); info.faces=in.get<unsigned>();
        if (!config.box.contains(p) || info.endpoint < -1 || info.endpoint>1 || info.fiber < -1 || info.faces>63) in.fail("invalid vertex/constraints");
        Vec fixed=p; config.box.fix(fixed,info.faces);
        if ((fixed-p).norm()>config.options.step*1e-8) in.fail("vertex does not satisfy its box-face label");
        r.mesh.vertices.push_back(p); r.info.push_back(info);
    }
    in.word("triangles"); const size_t nt=in.count();
    if (nt>config.options.max_triangles) in.fail("triangle count exceeds saved limit");
    for (size_t i=0; i<nt; ++i) {
        Triangle t; for (int& id:t) id=in.index(nv);
        if (!r.mesh.topology_ok(t)) in.fail("invalid triangle incidence/orientation"); r.mesh.add(t);
    }
    in.word("fibers"); const size_t nf=in.count();
    std::map<Edge,int> expected;
    for (size_t i=0; i<nf; ++i) {
        Fiber f; f.endpoint=in.index(2); f.closed=in.boolean(); f.complete=in.boolean(); const size_t count=in.count();
        if (count<2 || (f.closed && count<3)) in.fail("fiber has too few vertices");
        std::set<int> unique;
        for (size_t k=0; k<count; ++k) {
            const int id=in.index(nv); f.vertices.push_back(id);
            if (!unique.insert(id).second || r.info[id].fiber!=int(i) || r.info[id].endpoint!=f.endpoint) in.fail("inconsistent fiber vertex ID");
        }
        for (size_t k=0; k+(f.closed?0:1)<count; ++k) expected[edge_key(f.vertices[k],f.vertices[(k+1)%count])]=int(i);
        r.fibers.push_back(std::move(f));
    }
    for (const auto& info:r.info) if ((info.fiber>=0 && size_t(info.fiber)>=nf) || ((info.fiber<0)!=(info.endpoint<0))) in.fail("invalid fiber label");
    for (const auto& entry:std::array<std::pair<const char*,std::map<Edge,int>*>,3>{{{"fiber_edges",&r.fiber_edges},{"pending",&state.pending},{"fiber_uses",&state.fiber_uses}}}) {
        in.word(entry.first); const size_t count=in.count();
        for (size_t i=0; i<count; ++i) {
            const int a=in.index(nv),b=in.index(nv),fiber=in.index(nf); const Edge e{a,b};
            if (a==b || !entry.second->emplace(e,fiber).second || !expected.count(edge_key(a,b)) || expected.at(edge_key(a,b))!=fiber)
                in.fail("invalid fiber-edge use");
        }
    }
    if (r.fiber_edges!=expected) in.fail("fiber graph does not match the ordered curves");
    for (const auto& e:state.pending)
        if (!state.fiber_uses.count(e.first) || state.fiber_uses.at(e.first)!=e.second || r.mesh.edges.count(edge_key(e.first.first,e.first.second))) in.fail("invalid pending fiber use");
    in.word("box_edges"); const size_t nb=in.count();
    for (size_t i=0; i<nb; ++i) {
        const int a=in.index(nv),b=in.index(nv); const Edge key=edge_key(a,b);
        auto found=r.mesh.edges.find(key);
        if (a>=b || !r.box_edges.insert(key).second || found==r.mesh.edges.end() || found->second.second>=0 || !(r.info[a].faces & r.info[b].faces)) in.fail("invalid box boundary edge");
    }
    in.word("events"); const size_t ne=in.count();
    for (size_t i=0; i<ne; ++i) { const Vec p=in.point(); r.unresolved.push_back({p,in.quoted()}); }
    in.word("discovery"); r.fiber_seed_hits[0]=int(in.count()); r.fiber_seed_hits[1]=int(in.count()); r.surface_seed_hits=int(in.count());
    in.word("seconds"); r.seconds=in.real(); if (r.seconds<0) in.fail("negative duration"); in.end();
    if (stage==Stage::fibers && (!r.mesh.triangles.empty() || !state.pending.empty() || !state.fiber_uses.empty())) in.fail("fiber checkpoint contains surface-growth state");
    return state;
}

void write_stage_visuals(const fs::path& run,Stage stage,const Checkpoint& state) {
    const auto& r=state.result;
    atomic_write(stage_artifact(run,stage,"mesh.vtk"),[&](std::ostream& out) {
        points_header(out,r.mesh.vertices); out<<"POLYGONS "<<r.mesh.triangles.size()<<' '<<4*r.mesh.triangles.size()<<'\n';
        for (const auto& t:r.mesh.triangles) out<<"3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
        out<<"POINT_DATA "<<r.info.size()<<"\nSCALARS endpoint int 1\nLOOKUP_TABLE default\n";
        for (const auto& info:r.info) out<<info.endpoint<<'\n';
    });
    const auto front=unfinished(state);
    for (bool fibers:{false,true}) atomic_write(stage_artifact(run,stage,fibers?"fibers.vtk":"fronts.vtk"),[&](std::ostream& out) {
        points_header(out,r.mesh.vertices); std::vector<Edge> lines=front;
        if (fibers) { lines.clear(); for (const auto& e:r.fiber_edges) lines.push_back(e.first); }
        out<<"LINES "<<lines.size()<<' '<<3*lines.size()<<'\n'; for (const auto& e:lines) out<<"2 "<<e.first<<' '<<e.second<<'\n';
        out<<"CELL_DATA "<<lines.size()<<"\nSCALARS fiber_component int 1\nLOOKUP_TABLE default\n";
        for (const auto& e:lines) { auto found=r.fiber_edges.find(edge_key(e.first,e.second)); out<<(found==r.fiber_edges.end()?-1:found->second)<<'\n'; }
    });
}
void write_stage_report(const fs::path& run,Stage stage,const Checkpoint& state,double stage_seconds) {
    const auto& r=state.result;
    atomic_write(stage_artifact(run,stage,"report.json"),[&](std::ostream& out) {
        out<<"{\n  \"stage\": "<<int(stage)<<",\n  \"vertices\": "<<r.mesh.vertices.size()<<",\n  \"triangles\": "<<r.mesh.triangles.size()
            <<",\n  \"fiber_components\": "<<r.fibers.size()<<",\n  \"unfinished_fronts\": "<<unfinished(state).size()
            <<",\n  \"stage_seconds\": "<<stage_seconds<<",\n  \"cumulative_extraction_seconds\": "<<r.seconds<<",\n  \"unresolved\": ";
        events(out,r.unresolved); out<<"\n}\n";
    });
}
} // namespace fiber_surface
