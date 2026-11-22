#include "output.h"
#include "extract.h"

namespace fiber_surface {
namespace {
// Corner samples use the outgoing segment, except at the end of an open curve.
std::vector<int> vertex_segments(const Result& r, const RangeCurve& curve) {
    std::vector<int> labels(r.info.size(),0);
    for (size_t i=0; i<r.mesh.triangles.size(); ++i)
        for (int v:r.mesh.triangles[i]) labels[v]=r.face_segment(i);
    for (size_t v=0; v<r.info.size(); ++v)
        if (r.info[v].endpoint>=0) labels[v]=std::min(r.info[v].endpoint,curve.size()-1);
    return labels;
}
}

Report evaluate(const Result& r, const Field& field, const RangeCurve& curve, const Options& options) {
    Report report;
    const auto& mesh=r.mesh;
    report.topology=marching_triangles::validate_mesh(mesh.vertices,mesh.triangles,options.step*1e-9);
    const auto labels=vertex_segments(r,curve);
    for (size_t i=0; i<mesh.vertices.size(); ++i) {
        const int k=r.info[i].endpoint;
        const bool inside=report.vertices.add(mesh.vertices[i],field,curve.segments[labels[i]]);
        if (k<0 && !inside) ++report.ordinary_vertex_violations;
        if (k>=0) {
            const double residual=(field.value(mesh.vertices[i])-curve.point(k)).norm();
            report.endpoint_residual=std::max(report.endpoint_residual,
                std::isfinite(residual) ? residual : std::numeric_limits<double>::infinity());
            if (curve.closed() || (k>0 && k<curve.size()))
                report.vertices.add(mesh.vertices[i],field,curve.segments[(k+curve.size()-1)%curve.size()]);
        }
    }
    for (const auto& [edge,use]:mesh.edges) {
        const Vec mid=edge_midpoint(mesh,edge.first,edge.second);
        const int first=r.face_segment(use.face);
        report.midpoints.add(mid,field,curve.segments[first]);
        if (use.second>=0 && r.face_segment(use.second)!=first) {
            report.midpoints.add(mid,field,curve.segments[r.face_segment(use.second)]);
            if (!r.fiber_edges.count(edge)) ++report.incorrect_feature_incidence;
        }
        if (use.second<0 && !r.fiber_edges.count(edge) && !r.box_edges.count(edge)) ++report.unexpected_boundary;
    }
    for (const auto& [edge,component]:r.fiber_edges) {
        const int k=r.fibers[component].endpoint;
        std::set<int> expected, actual;
        if (curve.closed() || k>0) expected.insert((k+curve.size()-1)%curve.size());
        if (k<curve.size()) expected.insert(k);
        const auto found=mesh.edges.find(edge);
        if (found!=mesh.edges.end()) {
            actual.insert(r.face_segment(found->second.face));
            if (found->second.second>=0) actual.insert(r.face_segment(found->second.second));
        }
        if (actual!=expected || (found!=mesh.edges.end() && (found->second.second>=0)!=(expected.size()==2)))
            ++report.incorrect_feature_incidence;
    }
    for (const auto& e:r.box_edges) {
        auto found=mesh.edges.find(e);
        if (found==mesh.edges.end() || found->second.second>=0) ++report.incorrect_feature_incidence;
    }
    std::vector<double> angles;
    for (size_t face=0; face<mesh.triangles.size(); ++face) {
        const auto& t=mesh.triangles[face];
        const Vec &a=mesh.vertices[t[0]],&b=mesh.vertices[t[1]],&c=mesh.vertices[t[2]];
        report.centroids.add(face_centroid(mesh,t),field,curve.segments[r.face_segment(face)]);
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
                          const Field& field, const RangeCurve& curve, const Options& options) {
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
    const auto labels=vertex_segments(r,curve);
    vtk<<"POINT_DATA "<<r.mesh.vertices.size()<<"\nSCALARS "<<(curve.size()==1?"segment_parameter":"curve_parameter")
       <<" double 1\nLOOKUP_TABLE default\n";
    for (size_t i=0; i<labels.size(); ++i)
        vtk<<labels[i]+curve.segments[labels[i]].parameter(field.value(r.mesh.vertices[i]))<<'\n';
    for (const std::string name:{"endpoint","fiber_component","box_faces"}) {
        vtk<<"SCALARS "<<name<<" int 1\nLOOKUP_TABLE default\n";
        for (const auto& info:r.info) vtk<<(name=="endpoint" ? info.endpoint : name=="fiber_component" ? info.fiber : int(info.faces))<<'\n';
    }
    vtk<<"CELL_DATA "<<r.mesh.triangles.size()<<"\nSCALARS component int 1\nLOOKUP_TABLE default\n";
    for (int c:report.face_components) vtk<<c<<'\n';
    if (curve.size()>1) {
        vtk<<"SCALARS range_segment int 1\nLOOKUP_TABLE default\n";
        for (int side:r.face_segments) vtk<<side<<'\n';
    }
    points(fibers);
    fibers<<"LINES "<<r.fiber_edges.size()<<' '<<3*r.fiber_edges.size()<<'\n';
    for (const auto& e:r.fiber_edges) fibers<<"2 "<<e.first.first<<' '<<e.first.second<<'\n';
    fibers<<"CELL_DATA "<<r.fiber_edges.size()<<"\nSCALARS endpoint int 1\nLOOKUP_TABLE default\n";
    for (const auto& e:r.fiber_edges) fibers<<r.fibers[e.second].endpoint<<'\n';
    const auto shape=discovery_shape(field,options);
    json<<"{\n  \"complete\": "<<(report.complete(r)?"true":"false");
    if (curve.size()==1) {
        const auto& segment=curve.segments.front();
        json<<",\n  \"segment_a\": ["<<segment.a[0]<<','<<segment.a[1]<<"],\n  \"segment_b\": ["<<segment.b[0]<<','<<segment.b[1]<<']';
    } else {
        json<<",\n  \"closed\": "<<(curve.closed()?"true":"false")<<",\n  \"range_vertices\": [";
        for (int i=0; i<curve.vertex_count(); ++i) {
            const Range& p=curve.point(i);
            if (i) json<<',';
            json<<'['<<p.x()<<','<<p.y()<<']';
        }
        json<<']';
    }
    json<<','<<"\n  \"domain_min\": ["<<field.box.min[0]<<','<<field.box.min[1]<<','<<field.box.min[2]<<"],"
        <<"\n  \"domain_max\": ["<<field.box.max[0]<<','<<field.box.max[1]<<','<<field.box.max[2]<<"],"
        <<"\n  \"step\": "<<options.step<<",\n  \"fiber_step\": "<<options.fiber_step
        <<",\n  \"fiber_field_tolerance\": "<<Options::fiber_field_tolerance
        <<",\n  \"fiber_plane_tolerance\": "<<options.fiber_plane_tolerance()
        <<",\n  \"spatial_tolerance\": "<<options.spatial_tolerance<<",\n  \"range_tolerance\": "<<options.residual_tolerance
        <<",\n  \"discovery_grid\": "<<options.seed_grid
        <<",\n  \"discovery_grid_shape\": ["<<shape[0]<<','<<shape[1]<<','<<shape[2]<<']'
        <<",\n  \"fiber_seed_hits\": [";
    for (size_t i=0; i<r.fiber_seed_hits.size(); ++i) {
        if (i) json<<',';
        json<<r.fiber_seed_hits[i];
    }
    json<<']'
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
    json<<",\n  \"ordinary_vertex_violations\": "<<report.ordinary_vertex_violations;
    json<<",\n  \"unresolved\": [";
    for (size_t i=0; i<r.unresolved.size(); ++i) {
        const auto& e=r.unresolved[i]; if (i) json<<',';
        json<<"\n    {\"position\": ["<<e.position[0]<<','<<e.position[1]<<','<<e.position[2]<<"], \"reason\": "<<std::quoted(e.reason)<<'}';
    }
    json<<"\n  ]\n}\n";
    obj.close(); vtk.close(); fibers.close(); json.close();
    if (!obj || !vtk || !fibers || !json) throw std::runtime_error("Failed writing outputs at "+prefix);
}
int evaluate_and_write(Result result, const Field& field, const RangeCurve& curve, const Options& options, const std::string& prefix) {
    Report report=evaluate(result,field,curve,options);
    if (result.mesh.triangles.empty() && result.unresolved.empty())
        result.unresolved.push_back({field.box.min+(field.box.max-field.box.min)/2,"No regular surface discovered along the range curve"});
    if (!report.topology.valid()) result.unresolved.push_back({(field.box.min+field.box.max)/2,
        "Final mesh connectivity/intersection validation failed: "+report.topology.summary()});
    if (report.ordinary_vertex_violations || report.centroids.violations)
        result.unresolved.push_back({(field.box.min+field.box.max)/2,"Final sampling found finite-segment violations"});
    write_outputs(prefix,result,report,field,curve,options);
    std::cout<<"vertices="<<result.mesh.vertices.size()<<" triangles="<<result.mesh.triangles.size()<<" components="<<report.components
        <<" fibers="<<result.fibers.size()<<" unresolved="<<result.unresolved.size()<<'\n'<<report.topology.summary()
        <<"\nvertex max |phi|="<<report.vertices.maximum<<" endpoint max residual="<<report.endpoint_residual
        <<" min angle="<<report.min_angle<<" extraction seconds="<<result.seconds
        <<"\nWrote "<<prefix<<".obj, .vtk, _fibers.vtk, _report.json\n";
    if (!report.complete(result)) std::cerr<<"Incomplete extraction; inspect the report and partial mesh.\n";
    return report.complete(result) ? 0 : 2;
}
int run_and_write(const Field& field, const RangeCurve& curve, const Options& options, const std::string& prefix) {
    return evaluate_and_write(extract_curve(field,curve,options),field,curve,options,prefix);
}
} // namespace fiber_surface
