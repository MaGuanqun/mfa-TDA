#pragma once
#include "analytical_field.h"
#include "volume.h"

namespace fiber_surface {
struct VolumeExample { Field field; RangeTriangle triangle; double exact_volume; };
inline VolumeExample volume_example(const std::string& name) {
    if (name=="hollow-cylinder" || name=="hollow-sphere") {
        const bool sphere=name=="hollow-sphere";
        Field f; f.box.min=Vec(-2,-2,sphere?-2:0);f.box.max=Vec(2,2,sphere?2:1);
        f.value=[sphere](const Vec& x){return Range(x[0]*x[0]+x[1]*x[1]+(sphere?x[2]*x[2]:0),0);};
        f.jacobian=[sphere](const Vec& x){Jacobian j=Jacobian::Zero();j.row(0)<<2*x[0],2*x[1],(sphere?2*x[2]:0);return j;};
        f.hessians=[sphere](const Vec&){Hessian h=2*Hessian::Identity();if(!sphere)h(2,2)=0;return std::array<Hessian,2>{h,Hessian::Zero()};};
        return {f,RangeTriangle(Range(.25,-1),Range(1.75,-1),Range(.25,1)),
            sphere?4*std::acos(-1.)/3*(1-.125):.75*std::acos(-1.)};
    }
    if (name=="prism" || name=="curved-prism") {
        const bool curved=name=="curved-prism";
        Field f;
        f.box.min=curved ? Vec(-2,-2,-1) : Vec(-1,-1,0);
        f.box.max=curved ? Vec(2,2,1) : Vec(2,2,1);
        f.value=[curved](const Vec& x) { const double t=curved ? .5*x[2]*x[2] : 0; return Range(x[0]+t,x[1]+t); };
        f.jacobian=[curved](const Vec& x) { Jacobian j; j<<1,0,(curved?x[2]:0),0,1,(curved?x[2]:0); return j; };
        f.hessians=[curved](const Vec&) { Hessian h=Hessian::Zero(); h(2,2)=curved?1:0; return std::array<Hessian,2>{h,h}; };
        return {f,RangeTriangle(),curved?1.:.5};
    }
    auto e=analytical_field(name);
    // Reuse the previous sphere and disconnected fields with filled selections.
    if (name=="sphere") {
        const double z=(std::sqrt(89.)-3)/8;
        const double volume=std::acos(-1.)*(1.25*z-.375*z*z-z*z*z/3+31./24-std::sqrt(2.)/3);
        return {e.field,RangeTriangle(Range(.5,-1),Range(2,-1),Range(.5,1)),volume};
    }
    // Integrate both x branches over z in [-2,2]: for x>=0, the
    // section is 0<=y<=min((x*x-.25)/2,1-x*x), .5<=x<=1.
    const double volume=name=="disconnected" ? 17./3-3*std::sqrt(3.) : 0;
    return {e.field,RangeTriangle(Range(.25,0),Range(1,0),Range(1,.25)),volume};
}
} // namespace fiber_surface
