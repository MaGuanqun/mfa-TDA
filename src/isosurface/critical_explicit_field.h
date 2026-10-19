#pragma once
#include "critical_field.h"
#include <string>

namespace critical_curves {
// b is a field offset, independent of the requested isovalue.
// Analytic answers are documented, never supplied as seeds to the algorithm.
inline Field explicit_field(const std::string& name, double b = 1.25) {
    if(!std::isfinite(b)) throw std::invalid_argument("The field offset must be finite");
    Field f;
    if(name == "crossing_line") {
        f.value = [b](const Point& p){return p.x()*p.y()+b;};
        f.gradient = [](const Point& p)->Point{return Point(p.y(),p.x(),0);};
        f.hessian = [](const Point&)->Hessian{Hessian h=Hessian::Zero();h(0,1)=h(1,0)=1;return h;};
    } else if(name == "crossing_circle") {
        f.value = [b](const Point& p){return p.z()*(p.x()*p.x()+p.y()*p.y()-1)+b;};
        f.gradient = [](const Point& p)->Point{return Point(2*p.x()*p.z(),2*p.y()*p.z(),p.x()*p.x()+p.y()*p.y()-1);};
        f.hessian = [](const Point& p)->Hessian{Hessian h;h << 2*p.z(),0,2*p.x(),0,2*p.z(),2*p.y(),2*p.x(),2*p.y(),0;return h;};
    } else if(name == "parallel_lines") {
        f.value = [b](const Point& p){return p.z()*(p.z()-p.y()*p.y()+.25)+b;};
        f.gradient = [](const Point& p)->Point{return Point(0,-2*p.y()*p.z(),2*p.z()-p.y()*p.y()+.25);};
        f.hessian = [](const Point& p)->Hessian{Hessian h=Hessian::Zero();h(1,1)=-2*p.z();h(1,2)=h(2,1)=-2*p.y();h(2,2)=2;return h;};
    } else if(name == "isolated_rank_two") {
        f.value = [b](const Point& p){return p.x()*p.x()-p.y()*p.y()+std::pow(p.z(),4)+b;};
        f.gradient = [](const Point& p)->Point{return Point(2*p.x(),-2*p.y(),4*std::pow(p.z(),3));};
        f.hessian = [](const Point& p)->Hessian{return Point(2,-2,12*p.z()*p.z()).asDiagonal();};
    } else if(name == "definite_line") {
        f.value = [b](const Point& p){return p.x()*p.x()+p.y()*p.y()+b;};
        f.gradient = [](const Point& p)->Point{return Point(2*p.x(),2*p.y(),0);};
        f.hessian = [](const Point&)->Hessian{return Point(2,2,0).asDiagonal();};
    } else if(name == "tilted_line") {
        // No critical points: a nonzero tangential derivative destroys the line.
        f.value = [b](const Point& p){return p.x()*p.y()+.01*p.z()+b;};
        f.gradient = [](const Point& p)->Point{return Point(p.y(),p.x(),.01);};
        f.hessian = [](const Point&)->Hessian{Hessian h=Hessian::Zero();h(0,1)=h(1,0)=1;return h;};
    } else throw std::invalid_argument("Unknown critical-curve example: "+name);
    return f;
}
} // namespace critical_curves
