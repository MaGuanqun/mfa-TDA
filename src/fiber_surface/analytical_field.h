#pragma once
#include "field.h"

namespace fiber_surface {
struct Example { Field field; Segment segment; };

inline Example analytical_field(const std::string& name) {
    Example e;
    if (name == "sphere" || name == "closed-sphere") {
        e.field.value = [](const Vec& x) { return Range(x.squaredNorm(), x[2]); };
        e.field.jacobian = [](const Vec& x) { Jacobian j; j.row(0)=2*x.transpose(); j.row(1)=Vec::UnitZ(); return j; };
        e.field.hessians = [](const Vec&) { return std::array<Hessian,2>{2*Hessian::Identity(), Hessian::Zero()}; };
        if (name == "closed-sphere") e.segment=Segment(Range(1,-2),Range(1,2));
    } else if (name == "strip" || name == "disconnected") {
        const bool disconnected = name == "disconnected";
        e.field.value = [disconnected](const Vec& x) { return Range(disconnected ? x[0]*x[0]+x[1] : x[0], x[1]); };
        e.field.jacobian = [disconnected](const Vec& x) {
            Jacobian j; j << (disconnected ? 2*x[0] : 1), (disconnected ? 1 : 0), 0, 0, 1, 0; return j;
        };
        e.field.hessians = [disconnected](const Vec&) {
            std::array<Hessian,2> h{Hessian::Zero(), Hessian::Zero()}; if (disconnected) h[0](0,0)=2; return h;
        };
        e.segment=Segment(Range(disconnected ? .25 : 0,0),Range(1,0));
    } else throw std::invalid_argument("Unknown analytical field: "+name+" (sphere, closed-sphere, strip, disconnected)");
    return e;
}
} // namespace fiber_surface
