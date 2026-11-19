#include "volume.h"
#include <chrono>
#include <filesystem>

namespace fiber_surface {
VolumeResult extract_volume(const Field& f,const RangeTriangle& range,const Options& o,
                            VolumeStage stage,const std::string& boundary_input) {
    for(const auto& s:range.sides)o.validate(f,s);
    VolumeResult r;
    auto timed=[&](int stage,auto work) {const auto start=std::chrono::steady_clock::now();work();r.seconds[stage]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    if (stage==VolumeStage::fill)
        r.boundary=read_boundary_ply(boundary_input);
    else
        timed(0,[&]{r.boundary=assemble_volume_boundary(f,range,o);});
    timed(1,[&]{r.boundary_report=validate_volume_boundary(r.boundary,f,range,o);});
    if(!r.boundary_report.valid){r.failed_stage="boundary_validation";return r;}
    if(stage==VolumeStage::boundary)return r;
    try {timed(2,[&]{r.mesh=fill_volume(r.boundary,r.boundary_report);});}
    catch(const std::exception& e){r.failed_stage="volume_filling";r.boundary.unresolved.push_back({(f.box.min+f.box.max)/2,e.what()});return r;}
    timed(3,[&]{r.report=assess_volume(r.mesh,r.boundary,r.boundary_report,f,range,o);});
    if(!r.report.valid) {
        r.failed_stage="volume_assessment";
        r.boundary.unresolved.push_back({(f.box.min+f.box.max)/2,
            "Tetrahedral validity or boundary conformity failed; output boundary: "+r.report.output_topology.summary()});
    }
    return r;
}
namespace {
std::string json(const std::string& s) {
    std::ostringstream out;
    for(unsigned char c:s) {if(c=='"'||c=='\\')out<<'\\'<<c;else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else out<<c;}
    return out.str();
}
std::ofstream output(const std::string& path) {std::ofstream out(path);if(!out)throw std::runtime_error("Cannot write "+path);out<<std::setprecision(17);out.exceptions(std::ios::badbit|std::ios::failbit);return out;}
void points(std::ostream& out,const std::vector<Vec>& vertices) {out<<"POINTS "<<vertices.size()<<" double\n";for(const Vec& p:vertices)out<<p.transpose()<<'\n';}
}
void write_volume(const std::string& prefix,const VolumeResult& r,const Field& f,const RangeTriangle& range,const Options& o,
                  VolumeStage stage) {
    const auto parent=std::filesystem::path(prefix).parent_path();if(!parent.empty())std::filesystem::create_directories(parent);
    const auto& b=r.boundary;const auto& br=r.boundary_report;const auto& vr=r.report;
    if (stage!=VolumeStage::fill) write_boundary_ply(prefix+"_boundary.ply", b, br.valid);
    {
        auto out=output(prefix+"_boundary.obj");for(const Vec& p:b.mesh.vertices)out<<"v "<<p.transpose()<<'\n';
        for(size_t i=0;i<b.mesh.triangles.size();++i){auto t=b.mesh.triangles[i];out<<"g patch_"<<b.patches[i]<<"\nf "<<t[0]+1<<' '<<t[1]+1<<' '<<t[2]+1<<'\n';}
    }
    {
        auto out=output(prefix+"_boundary.vtk");out<<"# vtk DataFile Version 3.0\nFilled range triangle boundary\nASCII\nDATASET POLYDATA\n";points(out,b.mesh.vertices);
        out<<"POLYGONS "<<b.mesh.triangles.size()<<' '<<4*b.mesh.triangles.size()<<'\n';for(auto t:b.mesh.triangles)out<<"3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
        out<<"CELL_DATA "<<b.mesh.triangles.size()<<"\nSCALARS patch int 1\nLOOKUP_TABLE default\n";for(int p:b.patches)out<<p<<'\n';
        out<<"SCALARS shell int 1\nLOOKUP_TABLE default\n";for(size_t i=0;i<b.mesh.triangles.size();++i)out<<(i<br.shells.size()?br.shells[i]:-1)<<'\n';
        std::vector<unsigned> masks(b.mesh.vertices.size(),0);
        for(size_t i=0;i<b.mesh.triangles.size();++i)for(int id:b.mesh.triangles[i])masks[id]|=1u<<(b.patches[i]-1);
        out<<"POINT_DATA "<<masks.size()<<"\nSCALARS incident_patch_mask unsigned_int 1\nLOOKUP_TABLE default\n";for(auto mask:masks)out<<mask<<'\n';
    }
    {
        auto out=output(prefix+"_boundary.poly");out<<"# Only validated boundaries may be filled; see report.json.\n"<<b.mesh.vertices.size()<<" 3 0 0\n";
        for(size_t i=0;i<b.mesh.vertices.size();++i)out<<i<<' '<<b.mesh.vertices[i].transpose()<<'\n';
        out<<b.mesh.triangles.size()<<" 1\n";
        for(size_t i=0;i<b.mesh.triangles.size();++i){auto t=b.mesh.triangles[i];out<<"1 0 "<<i+1<<" # patch "<<b.patches[i]<<"\n3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';}
        out<<br.holes.size()<<'\n';for(size_t i=0;i<br.holes.size();++i)out<<i<<' '<<br.holes[i].transpose()<<'\n';out<<"0\n";
    }
    {
        auto out=output(prefix+"_volume.vtk");out<<"# vtk DataFile Version 3.0\nFilled range triangle tetrahedra\nASCII\nDATASET UNSTRUCTURED_GRID\n";points(out,r.mesh.vertices);
        out<<"CELLS "<<r.mesh.cells.size()<<' '<<5*r.mesh.cells.size()<<'\n';for(auto t:r.mesh.cells)out<<"4 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<' '<<t[3]<<'\n';
        out<<"CELL_TYPES "<<r.mesh.cells.size()<<'\n';for(size_t i=0;i<r.mesh.cells.size();++i)out<<"10\n";
    }
    {
        auto out=output(prefix+"_tet_boundary.vtk");out<<"# vtk DataFile Version 3.0\nTetGen boundary subdivisions\nASCII\nDATASET POLYDATA\n";points(out,r.mesh.vertices);
        out<<"POLYGONS "<<r.mesh.boundary_faces.size()<<' '<<4*r.mesh.boundary_faces.size()<<'\n';for(auto t:r.mesh.boundary_faces)out<<"3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
        out<<"CELL_DATA "<<r.mesh.boundary_faces.size()<<"\nSCALARS parent_facet int 1\nLOOKUP_TABLE default\n";for(int id:r.mesh.parent_facets)out<<id<<'\n';
        out<<"SCALARS patch int 1\nLOOKUP_TABLE default\n";for(int id:r.mesh.parent_facets)out<<(id>=0&&size_t(id)<b.patches.size()?b.patches[id]:-1)<<'\n';
        std::vector<unsigned> masks(r.mesh.vertices.size(),0);
        for(size_t i=0;i<r.mesh.boundary_faces.size();++i) {const int id=r.mesh.parent_facets[i];if(id>=0&&size_t(id)<b.patches.size())for(int v:r.mesh.boundary_faces[i])masks[v]|=1u<<(b.patches[id]-1);}
        out<<"POINT_DATA "<<masks.size()<<"\nSCALARS incident_patch_mask unsigned_int 1\nLOOKUP_TABLE default\n";for(auto mask:masks)out<<mask<<'\n';
    }
    {
        auto out=output(prefix+"_report.json");
        const char* stage_name=stage==VolumeStage::boundary?"boundary":stage==VolumeStage::fill?"fill":"all";
        out<<"{\n  \"stage\": \""<<stage_name<<"\",\n  \"complete\": "<<(r.complete()?"true":"false")<<",\n  \"failed_stage\": \""<<r.failed_stage<<"\",\n  \"tetgen_switches\": \"pV\",\n  \"coverage\": \"empirical, not certified\",\n  \"triangle\": [";
        for(int i=0;i<3;++i){if(i)out<<',';out<<'['<<range.corners[i][0]<<','<<range.corners[i][1]<<']';}
        out<<"],\n  \"box_min\": ["<<f.box.min[0]<<','<<f.box.min[1]<<','<<f.box.min[2]<<"],\n  \"box_max\": ["<<f.box.max[0]<<','<<f.box.max[1]<<','<<f.box.max[2]<<"],\n";
        out<<"  \"step\": "<<o.step<<", \"fiber_step\": "<<o.fiber_step<<", \"spatial_tolerance\": "<<o.spatial_tolerance<<", \"residual_tolerance\": "<<o.residual_tolerance<<",\n";
        out<<"  \"boundary_vertices\": "<<b.mesh.vertices.size()<<", \"boundary_triangles\": "<<b.mesh.triangles.size()<<", \"box_contact_triangles\": "<<b.box_contact_triangles<<",\n  \"tet_vertices\": "<<r.mesh.vertices.size()<<", \"tetrahedra\": "<<r.mesh.cells.size()<<",\n";
        out<<"  \"corner_fiber_components\": ["<<b.corner_components[0]<<','<<b.corner_components[1]<<','<<b.corner_components[2]<<"],\n  \"shell_signed_volumes\": [";
        for(size_t i=0;i<br.signed_volumes.size();++i){if(i)out<<',';out<<br.signed_volumes[i];}out<<"],\n  \"hole_seeds\": [";
        for(size_t i=0;i<br.holes.size();++i){if(i)out<<',';const Vec p=br.holes[i];out<<'['<<p[0]<<','<<p[1]<<','<<p[2]<<']';}out<<"],\n";
        out<<"  \"boundary_valid\": "<<(br.valid?"true":"false")<<", \"topology\": \""<<br.topology.summary()<<"\",\n";
        out<<"  \"constraint_failures\": "<<br.constraint_failures<<", \"orientation_failures\": "<<br.orientation_failures<<", \"coverage_misses\": "<<br.coverage_misses<<",\n";
        out<<"  \"max_side_residual\": "<<br.side_residual<<", \"max_corner_residual\": "<<br.corner_residual<<", \"max_inequality_violation\": "<<br.inequality_violation<<", \"max_spatial_error\": "<<br.spatial_error<<",\n";
        out<<"  \"tet_max_side_residual\": "<<vr.side_residual<<", \"tet_max_inequality_violation\": "<<vr.inequality_violation<<", \"tet_max_spatial_error\": "<<vr.spatial_error<<",\n";
        out<<"  \"volume\": "<<vr.volume<<", \"boundary_volume\": "<<vr.boundary_volume<<", \"min_cell_volume\": "<<vr.min_volume<<", \"max_cell_volume\": "<<vr.max_volume<<",\n";
        out<<"  \"invalid_cells\": "<<vr.invalid_cells<<", \"invalid_faces\": "<<vr.invalid_faces<<", \"conformity_failures\": "<<vr.conformity_failures<<", \"tet_intersections\": "<<vr.intersections<<",\n";
        out<<"  \"tet_boundary_topology\": \""<<vr.output_topology.summary()<<"\",\n";
        out<<"  \"stage_seconds\": ["<<r.seconds[0]<<','<<r.seconds[1]<<','<<r.seconds[2]<<','<<r.seconds[3]<<"],\n  \"unresolved\": [";
        for(size_t i=0;i<b.unresolved.size();++i){if(i)out<<',';const auto& e=b.unresolved[i];out<<"{\"position\":["<<e.position[0]<<','<<e.position[1]<<','<<e.position[2]<<"],\"reason\":\""<<json(e.reason)<<"\"}";}out<<"]\n}\n";
    }
}
int run_volume(const Field& f,const RangeTriangle& range,const Options& o,const std::string& prefix,
               VolumeStage stage,const std::string& boundary_input) {
    auto result=extract_volume(f,range,o,stage,boundary_input);
    write_volume(prefix,result,f,range,o,stage);
    if (stage==VolumeStage::boundary) {
        std::cout<<"Boundary valid: "<<result.boundary_report.valid<<", triangles: "<<result.boundary.mesh.triangles.size()
                 <<"\nBoundary: "<<prefix<<"_boundary.ply\nReport: "<<prefix<<"_report.json\n";
        return result.boundary_report.valid?0:2;
    }
    std::cout<<"Volume complete: "<<result.complete()<<", boundary triangles: "<<result.boundary.mesh.triangles.size()<<", tetrahedra: "<<result.mesh.cells.size()<<", volume: "<<std::setprecision(12)<<result.report.volume<<"\nReport: "<<prefix<<"_report.json\n";
    return result.complete()?0:2;
}
} // namespace fiber_surface
