#include "output.h"
#include "extract.h"

namespace fiber_surface {
Report evaluate(const Result& r, const Field& field, const Segment& segment, const Options& options) {
    Report report;
    const auto& mesh=r.mesh;
    report.topology=marching_triangles::validate_mesh(mesh.vertices,mesh.triangles,options.step*1e-9);
    for (size_t i=0; i<mesh.vertices.size(); ++i) {
        report.vertices.add(mesh.vertices[i],field,segment,options);
        if (r.info[i].endpoint>=0) report.endpoint_residual=std::max(report.endpoint_residual,
            (field.value(mesh.vertices[i])-(r.info[i].endpoint ? segment.b : segment.a)).norm());
    }
    for (const auto& e:mesh.edges) {
        report.midpoints.add(Vec((mesh.vertices[e.first.first]+mesh.vertices[e.first.second])/2),field,segment,options);
        if (e.second.second<0 && !r.fiber_edges.count(e.first) && !r.box_edges.count(e.first)) ++report.unexpected_boundary;
    }
    for (const auto& e:r.fiber_edges) {
        auto found=mesh.edges.find(e.first);
        if (found==mesh.edges.end() || found->second.second>=0) ++report.incorrect_feature_incidence;
    }
    for (const auto& e:r.box_edges) {
        auto found=mesh.edges.find(e);
        if (found==mesh.edges.end() || found->second.second>=0) ++report.incorrect_feature_incidence;
    }
    std::vector<double> angles;
    for (const auto& t:mesh.triangles) {
        const Vec &a=mesh.vertices[t[0]],&b=mesh.vertices[t[1]],&c=mesh.vertices[t[2]];
        report.centroids.add(Vec((a+b+c)/3),field,segment,options);
        report.area+=(b-a).cross(c-a).norm()/2;
        double minimum=180;
        for (int i=0; i<3; ++i) {
            const Vec u=mesh.vertices[t[(i+1)%3]]-mesh.vertices[t[i]],v=mesh.vertices[t[(i+2)%3]]-mesh.vertices[t[i]];
            minimum=std::min(minimum,std::atan2(u.cross(v).norm(),u.dot(v))*180/std::acos(-1.));
        }
        angles.push_back(minimum);
    }
    std::sort(angles.begin(),angles.end());
    if (!angles.empty()) { report.min_angle=angles.front(); report.median_angle=angles[angles.size()/2]; report.p05_angle=angles[angles.size()/20]; }
    std::vector<std::vector<int>> adjacency(mesh.triangles.size());
    for (const auto& e:mesh.edges) if (e.second.second>=0) {
        adjacency[e.second.face].push_back(e.second.second); adjacency[e.second.second].push_back(e.second.face);
    }
    report.face_components.assign(mesh.triangles.size(),-1);
    for (size_t i=0; i<mesh.triangles.size(); ++i) if (report.face_components[i]<0) {
        std::vector<int> stack{int(i)}; report.face_components[i]=int(report.components);
        while (!stack.empty()) {
            const int j=stack.back(); stack.pop_back();
            for (int k:adjacency[j]) if (report.face_components[k]<0) { report.face_components[k]=int(report.components); stack.push_back(k); }
        }
        ++report.components;
    }
    return report;
}

void write_outputs(const std::string& prefix, const Result& r, const Report& report,
                          const Field& field, const Segment& segment, const Options& options) {
    std::ofstream obj(prefix+".obj"),vtk(prefix+".vtk"),fibers(prefix+"_fibers.vtk"),json(prefix+"_report.json");
    if (!obj || !vtk || !fibers || !json) throw std::runtime_error("Cannot create output files at "+prefix);
    obj<<std::setprecision(17); vtk<<std::setprecision(17); fibers<<std::setprecision(17); json<<std::setprecision(17);
    for (const Vec& p:r.mesh.vertices) obj<<"v "<<p.transpose()<<'\n';
    for (const auto& t:r.mesh.triangles) obj<<"f "<<t[0]+1<<' '<<t[1]+1<<' '<<t[2]+1<<'\n';
    auto points=[&](std::ostream& out) {
        out<<"# vtk DataFile Version 3.0\nContinuous-model fiber surface\nASCII\nDATASET POLYDATA\nPOINTS "<<r.mesh.vertices.size()<<" double\n";
        for (const auto& p:r.mesh.vertices) out<<p.transpose()<<'\n';
    };
    points(vtk);
    vtk<<"POLYGONS "<<r.mesh.triangles.size()<<' '<<4*r.mesh.triangles.size()<<'\n';
    for (const auto& t:r.mesh.triangles) vtk<<"3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
    vtk<<"POINT_DATA "<<r.mesh.vertices.size()<<"\nSCALARS segment_parameter double 1\nLOOKUP_TABLE default\n";
    for (const auto& p:r.mesh.vertices) vtk<<segment.parameter(field.value(p))<<'\n';
    for (const std::string name:{"endpoint","fiber_component","box_faces"}) {
        vtk<<"SCALARS "<<name<<" int 1\nLOOKUP_TABLE default\n";
        for (const auto& info:r.info) vtk<<(name=="endpoint" ? info.endpoint : name=="fiber_component" ? info.fiber : int(info.faces))<<'\n';
    }
    vtk<<"CELL_DATA "<<r.mesh.triangles.size()<<"\nSCALARS component int 1\nLOOKUP_TABLE default\n";
    for (int c:report.face_components) vtk<<c<<'\n';
    points(fibers);
    fibers<<"LINES "<<r.fiber_edges.size()<<' '<<3*r.fiber_edges.size()<<'\n';
    for (const auto& e:r.fiber_edges) fibers<<"2 "<<e.first.first<<' '<<e.first.second<<'\n';
    fibers<<"CELL_DATA "<<r.fiber_edges.size()<<"\nSCALARS endpoint int 1\nLOOKUP_TABLE default\n";
    for (const auto& e:r.fiber_edges) fibers<<r.fibers[e.second].endpoint<<'\n';
    json<<"{\n  \"complete\": "<<(report.complete(r)?"true":"false")
        <<",\n  \"segment_a\": ["<<segment.a[0]<<','<<segment.a[1]<<"],\n  \"segment_b\": ["<<segment.b[0]<<','<<segment.b[1]<<"],"
        <<"\n  \"domain_min\": ["<<field.box.min[0]<<','<<field.box.min[1]<<','<<field.box.min[2]<<"],"
        <<"\n  \"domain_max\": ["<<field.box.max[0]<<','<<field.box.max[1]<<','<<field.box.max[2]<<"],"
        <<"\n  \"step\": "<<options.step<<",\n  \"fiber_step\": "<<options.fiber_step
        <<",\n  \"spatial_tolerance\": "<<options.spatial_tolerance<<",\n  \"range_tolerance\": "<<options.residual_tolerance
        <<",\n  \"discovery_grid\": "<<options.seed_grid<<",\n  \"fiber_seed_hits\": ["<<r.fiber_seed_hits[0]<<','<<r.fiber_seed_hits[1]<<']'
        <<",\n  \"surface_seed_hits\": "<<r.surface_seed_hits
        <<",\n  \"vertices\": "<<r.mesh.vertices.size()<<",\n  \"triangles\": "<<r.mesh.triangles.size()
        <<",\n  \"components\": "<<report.components<<",\n  \"fiber_components\": "<<r.fibers.size()
        <<",\n  \"boundary_edges\": "<<report.topology.boundary_edges<<",\n  \"unexpected_boundary_edges\": "<<report.unexpected_boundary
        <<",\n  \"incorrect_feature_incidence\": "<<report.incorrect_feature_incidence
        <<",\n  \"euler_characteristic\": "<<report.topology.euler_characteristic
        <<",\n  \"nonmanifold_vertices\": "<<report.topology.nonmanifold_vertices
        <<",\n  \"intersections\": "<<report.topology.intersections
        <<",\n  \"topology_valid\": "<<(report.topology.valid()?"true":"false")
        <<",\n  \"endpoint_max_residual\": "<<report.endpoint_residual<<",\n  \"area\": "<<report.area
        <<",\n  \"minimum_angle_degrees\": "<<report.min_angle<<",\n  \"p05_minimum_angle_degrees\": "<<report.p05_angle
        <<",\n  \"median_minimum_angle_degrees\": "<<report.median_angle<<",\n  \"extraction_seconds\": "<<r.seconds;
    auto samples=[&](const char* name,const Samples& s) {
        json<<",\n  \""<<name<<"\": {\"count\": "<<s.count<<", \"mean_abs_phi\": "<<s.mean()<<", \"max_abs_phi\": "<<s.maximum
            <<", \"segment_violations\": "<<s.violations<<", \"min_s\": "<<(s.count?s.min_s:0)<<", \"max_s\": "<<(s.count?s.max_s:0)<<'}';
    };
    samples("vertex_samples",report.vertices); samples("edge_midpoint_samples",report.midpoints); samples("centroid_samples",report.centroids);
    json<<",\n  \"unresolved\": [";
    for (size_t i=0; i<r.unresolved.size(); ++i) {
        const auto& e=r.unresolved[i]; if (i) json<<',';
        json<<"\n    {\"position\": ["<<e.position[0]<<','<<e.position[1]<<','<<e.position[2]<<"], \"reason\": "<<std::quoted(e.reason)<<'}';
    }
    json<<"\n  ]\n}\n";
    obj.close(); vtk.close(); fibers.close(); json.close();
    if (!obj || !vtk || !fibers || !json) throw std::runtime_error("Failed writing outputs at "+prefix);
}
int evaluate_and_write(Result result, const Field& field, const Segment& segment, const Options& options, const std::string& prefix) {
    Report report=evaluate(result,field,segment,options);
    if (!report.topology.valid()) result.unresolved.push_back({(field.box.min+field.box.max)/2,
        "Final mesh connectivity/intersection validation failed: "+report.topology.summary()});
    if (report.vertices.violations || report.midpoints.violations || report.centroids.violations)
        result.unresolved.push_back({(field.box.min+field.box.max)/2,"Final sampling found finite-segment violations"});
    write_outputs(prefix,result,report,field,segment,options);
    std::cout<<"vertices="<<result.mesh.vertices.size()<<" triangles="<<result.mesh.triangles.size()<<" components="<<report.components
        <<" fibers="<<result.fibers.size()<<" unresolved="<<result.unresolved.size()<<'\n'<<report.topology.summary()
        <<"\nvertex max |phi|="<<report.vertices.maximum<<" endpoint max residual="<<report.endpoint_residual
        <<" min angle="<<report.min_angle<<" extraction seconds="<<result.seconds
        <<"\nWrote "<<prefix<<".obj, .vtk, _fibers.vtk, _report.json\n";
    if (!report.complete(result)) std::cerr<<"Incomplete extraction; inspect the report and partial mesh.\n";
    return report.complete(result) ? 0 : 2;
}
int run_and_write(const Field& field, const Segment& segment, const Options& options, const std::string& prefix) {
    return evaluate_and_write(Extractor(field,segment,options).run(),field,segment,options,prefix);
}
} // namespace fiber_surface
