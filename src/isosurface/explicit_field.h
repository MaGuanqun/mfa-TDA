#pragma once
#include "surface_field.h"
#include <iostream>
#include <string>
// The legacy analytic functions use these names unqualified. Keep the aliases
// local to their namespace instead of importing MFA/Torch just for Eigen types.
namespace closed_form_function {
template<class T> using VectorX = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using Eigen::VectorXd;
using Eigen::VectorXi;
using std::string;
}
#include "../closed_form_function.h"

namespace marching_triangles {
struct ExplicitField {
    SurfaceField<double> field;
    MeshingOptions<double> options;
    double base_step;
};
inline ExplicitField explicit_field(const std::string& name,double iso) {
    int type=0;
    if(name=="quartic_potential") type=1;
    else if(name=="quartic_potential_2") type=2;
    else if(name=="rotating_quartic_multiwell") type=3;
    else if(name=="ellipsoid") type=5;
    else if(name=="quartic_potential_3d") throw std::invalid_argument("quartic_potential_3d has a 4D domain; triangle extraction requires a 3D scalar field");
    else throw std::invalid_argument("Unknown explicit function: "+name);
    ExplicitField result;
    auto evaluate=[type](const Point<double>& p,const Eigen::VectorXi& deriv) {
        const Eigen::VectorXd point=p;
        Eigen::VectorXd out;
        switch(type) {
        case 1:closed_form_function::quartic_potential(point,out,deriv);break;
        case 2:closed_form_function::quartic_potential_2(point,out,deriv);break;
        case 3:closed_form_function::rotating_quartic_multiwell(point,out,deriv);break;
        case 5:closed_form_function::ellipsoid(point,out,deriv);break;
        }
        return out[0];
    };
    result.field.value=[evaluate](const Point<double>& p) {return evaluate(p,Eigen::VectorXi());};
    result.field.gradient=[evaluate](const Point<double>& p)->Point<double> {
        Point<double> g;
        for(int i=0;i<3;++i) {Eigen::VectorXi deriv=Eigen::VectorXi::Zero(3);deriv[i]=1;g[i]=evaluate(p,deriv);}
        return g;
    };
    result.options.iso_value=iso;
    result.options.domain_min=closed_form_function::domain_min(type);
    result.options.domain_max=closed_form_function::domain_max(type);
    if(type==5) {
        if(!std::isfinite(iso)||iso<=0) throw std::invalid_argument("ellipsoid requires an isovalue > 0 (its zero level is a single point)");
        result.options.domain_max=1.1*std::sqrt(iso)*Point<double>(1,2,3);
        result.options.domain_min=-result.options.domain_max;
    }
    // Use the final domain, including the ellipsoid's level-dependent padding.
    // The CLI resolution divisor further scales this initial edge length.
    result.base_step=(result.options.domain_max-result.options.domain_min).minCoeff()/20.0;
    result.options.step=result.base_step;
    result.options.min_step=result.base_step*.1;
    result.options.vertex_tolerance=result.base_step*1e-6;
    return result;
}
} // namespace marching_triangles
