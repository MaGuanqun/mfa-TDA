#pragma once
#include "marching_triangles.h"
#include "mesh_validation.h"
#include "root_io.h"
#include <iostream>

namespace marching_triangles {
// Common extraction/validation path for both explicit and MFA adapters.
inline void extract_isosurface(const FieldSetup<double>& setup,const std::string& root_file,
                              const std::string& singular_file,const std::string& output,bool allow_open) {
    const auto& options=setup.options;
    const auto seeds=read_roots(root_file);
    MarchingTriangles<double> mesher(setup.field,options);
    if(!singular_file.empty())mesher.set_degenerate_points(read_roots(singular_file));
    std::vector<std::vector<Point<double>>> vertices;
    std::vector<std::vector<Triangle>> triangles;
    std::cout<<"Loaded "<<seeds.size()<<" roots; seed edge length="<<options.step<<std::endl;
    if(!mesher.extract_all_sheets(seeds,vertices,triangles))throw std::runtime_error("No regular surface sheet found; check the seeds and isovalue");
    if(!mesher.diagnostic().empty())std::cerr<<mesher.diagnostic()<<std::endl;
    // Validate before publishing any of the meshes. Failure is observable
    // as a nonzero exit code, not a success message attached to an open OBJ.
    bool valid=true;
    for(size_t i=0;i<vertices.size();++i){
        const auto report=validate_mesh(vertices[i],triangles[i],options.step*1e-9);
        std::cout<<"sheet "<<i<<": vertices="<<vertices[i].size()<<" triangles="<<triangles[i].size()<<' '<<report.summary()<<std::endl;
        valid=valid&&report.valid()&&(allow_open||report.closed());
        double max_residual=0;
        for(const auto& p:vertices[i]) {
            const double residual=std::abs(setup.field.value(p)-options.iso_value);
            if(!std::isfinite(residual))max_residual=std::numeric_limits<double>::infinity();
            else max_residual=std::max(max_residual,residual);
            valid=valid&&(p.array()>=options.domain_min.array()).all()&&(p.array()<=options.domain_max.array()).all();
        }
        std::cout<<"sheet "<<i<<" max scalar residual="<<max_residual<<std::endl;
        valid=valid&&max_residual<=options.projection_tolerance*1.01;
        if(report.invalid_indices==0) {
            const auto quality=measure_mesh_quality(vertices[i],triangles[i],setup.field);
            std::cout<<"sheet "<<i<<" quality: "<<quality.summary()<<std::endl;
            valid=valid&&quality.reversed_faces==0;
            if(quality.triangles_below_5_degrees || quality.max_normal_error_degrees>45)
                std::cerr<<"sheet "<<i<<": poor triangle shapes or normal alignment remain; inspect the quality report and consider a finer resolution.\n";
        }
    }
    if(!valid)throw std::runtime_error("Mesh validation failed. No output written; reduce the step size or use --allow-open for an intentionally bounded surface");
    for(size_t i=0;i<vertices.size();++i){
        const std::string filename=output+std::to_string(i)+".obj";
        mesher.save_mesh_obj(filename,vertices[i],triangles[i]);
        std::cout<<"Wrote "<<filename<<std::endl;
    }
}
} // namespace marching_triangles
