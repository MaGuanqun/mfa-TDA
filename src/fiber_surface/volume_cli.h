#pragma once
#include "cli.h"
#include "volume.h"
namespace fiber_surface {
inline void volume_options(Arguments& args,Field& field,RangeTriangle& triangle,Options& options,bool required=false) {
    const auto a=args.take("--q0"),b=args.take("--q1"),c=args.take("--q2");
    if(required || !a.empty() || !b.empty() || !c.empty()) {
        if(a.empty() || b.empty() || c.empty())throw std::invalid_argument("Supply all three range vertices: --q0, --q1, --q2");
        triangle=RangeTriangle(coordinates<2>(a),coordinates<2>(b),coordinates<2>(c));
    }
    if(args.has("--a") || args.has("--b"))throw std::invalid_argument("Use --q0 --q1 --q2 for a filled range triangle");
    Segment dummy=triangle.sides[0];common_options(args,field,dummy,options);
    for(const auto& s:triangle.sides)options.validate(field,s);
}
inline void volume_help() {
    std::cout<<"  --q0 u,v --q1 u,v --q2 u,v  counterclockwise range triangle\n"
        <<"  --min x,y,z --max x,y,z      spatial sub-box\n"
        <<"  --step H --fiber-step H --spatial-tolerance E --residual-tolerance E\n"
        <<"  --seed-grid N --max-triangles N (per side)\n"
        <<"  --output PREFIX  boundary OBJ/VTK/PLC, tetrahedral VTK, and JSON assessment\n"
        <<"Exit 0: validated volume; 1: invalid input/runtime error; 2: unresolved partial result.\n";
}
} // namespace fiber_surface
