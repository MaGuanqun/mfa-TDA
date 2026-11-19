#pragma once
#include "cli.h"
#include "volume.h"
namespace fiber_surface {
struct VolumeFiles {
    std::string prefix;
    VolumeStage stage;
    std::string boundary_input;
};

inline VolumeFiles volume_files(Arguments& args, const std::string& default_prefix) {
    VolumeFiles files{args.take("--output",default_prefix), VolumeStage::all, {}};
    const bool has_input=args.has("--boundary-input");
    const auto stage=args.take("--stage",has_input?"fill":"all");
    if (stage=="fill") {
        files.stage=VolumeStage::fill;
        files.boundary_input=args.take("--boundary-input",files.prefix+"_boundary.ply");
        if (files.boundary_input.empty()) throw std::invalid_argument("Expected a PLY path for --boundary-input");
    } else {
        if (has_input) throw std::invalid_argument("--boundary-input requires the fill stage");
        if (stage=="boundary") files.stage=VolumeStage::boundary;
        else if (stage!="all") throw std::invalid_argument("Expected --stage all, boundary, or fill");
    }
    return files;
}
inline void volume_options(Arguments& args,Field& field,RangeTriangle& triangle,Options& options,bool required=false) {
    const auto a=args.take("--q0"),b=args.take("--q1"),c=args.take("--q2");
    if(required || !a.empty() || !b.empty() || !c.empty()) {
        if(a.empty() || b.empty() || c.empty())throw std::invalid_argument("Supply all three range vertices: --q0, --q1, --q2");
        triangle=RangeTriangle(coordinates<2>(a),coordinates<2>(b),coordinates<2>(c));
    }
    if(args.has("--a") || args.has("--b"))throw std::invalid_argument("Use --q0 --q1 --q2 for a filled range triangle");
    common_options(args,field,options);
    for(const auto& s:triangle.sides)options.validate(field,s);
}
inline void volume_help() {
    std::cout<<"  --q0 u,v --q1 u,v --q2 u,v  counterclockwise range triangle\n"
        <<"  --min x,y,z --max x,y,z      spatial sub-box\n"
        <<"  --step H --fiber-step H --spatial-tolerance E --residual-tolerance E\n"
        <<"  --seed-grid N --max-triangles N (per side)\n"
        <<"  --stage all|boundary|fill    full extraction (default), export boundary, or fill from PLY\n"
        <<"  --boundary-input FILE        fill a saved boundary (default for fill: PREFIX_boundary.ply)\n"
        <<"  --output PREFIX              prefix for boundary meshes, tetrahedral VTK and JSON report\n"
        <<"Exit 0: requested stage passed; 1: input/runtime error; 2: unresolved partial result.\n";
}
} // namespace fiber_surface
