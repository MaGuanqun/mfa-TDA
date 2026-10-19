#include "critical_explicit_field.h"
#include "critical_curves.h"
#include <iostream>

using namespace critical_curves;
static void require(bool ok,const std::string& message) {if(!ok) throw std::runtime_error(message);}

static TraceResult extract(const std::string& name, double level=1.25) {
    Options o; o.samples=5; o.step=.08;
    Problem p(explicit_field(name),o);
    const auto s=discover(p,level);
    auto r=trace(p,s.seeds,level);
    for(const auto& c:r.curves) for(const auto& x:c.points) {
        require(p.contains(x.position),"curve point outside domain");
        require(valid_sample(p,x.position,level),"curve sample failed residual/signature check");
    }
    return r;
}

int main() {
    try {
        const auto line=extract("crossing_line");
        require(line.curves.size()==1 && line.events.empty(),"expected one complete line");
        const auto& l=line.curves[0];
        require(!l.closed && l.begin=="boundary" && l.end=="boundary","line did not end on both faces");
        require(std::abs(l.points.front().position.z()-l.points.back().position.z())>3.999999,"line endpoints incomplete");
        for(const auto& p:l.points) require(p.position.head<2>().norm()<1e-6,"line position error");

        const auto circle=extract("crossing_circle");
        require(circle.curves.size()==1 && circle.curves[0].closed && circle.events.empty(),"expected one closed circle");
        double circumference=0;
        const auto& cp=circle.curves[0].points;
        for(size_t i=0;i<cp.size();++i) {
            require(std::abs(cp[i].position.z())<1e-6 && std::abs(cp[i].position.head<2>().norm()-1)<1e-6,"circle position error");
            circumference+=(cp[i].position-cp[(i+1)%cp.size()].position).norm();
        }
        require(std::abs(circumference-2*std::acos(-1.))<.01,"incorrect circle circumference");
        require(extract("parallel_lines").curves.size()==2,"parallel curves were lost or merged");
        require(extract("crossing_circle",1.4).curves.empty(),"wrong target level produced curves");
        require(extract("definite_line").curves.empty(),"definite critical manifold accepted as crossing");
        require(extract("tilted_line").curves.empty(),"nonzero tangential derivative ignored");

        Options options; options.samples=3; options.step=.08;
        Problem isolated(explicit_field("isolated_rank_two"),options);
        const std::vector<Point> origin{Point::Zero()};
        const auto candidate=discover(isolated,1.25,&origin);
        require(candidate.seeds.size()==1,"rank-two negative control was not a candidate");
        const auto rejected=trace(isolated,candidate.seeds,1.25);
        require(rejected.curves.empty() && !rejected.events.empty(),"isolated rank-two point reported as curve");

        // Rotate the closed curve, translate it, and change spatial/field units.
        // This exercises mixed derivatives and prevents coordinate-slice tracing.
        const Field base=explicit_field("crossing_circle",0);
        const Hessian rotation=Eigen::AngleAxisd(.67,Point(1,2,3).normalized()).toRotationMatrix();
        const Point center(10,-20,30);
        const double length=7,amplitude=1e-4,shift=.37;
        Field transformed;
        transformed.lower=center-Point::Constant(2*length);
        transformed.upper=center+Point::Constant(2*length);
        transformed.value=[=](const Point& x){return shift+amplitude*base.value(rotation.transpose()*(x-center)/length);};
        transformed.gradient=[=](const Point& x)->Point{return amplitude/length*rotation*base.gradient(rotation.transpose()*(x-center)/length);};
        transformed.hessian=[=](const Point& x)->Hessian{return amplitude/(length*length)*rotation*base.hessian(rotation.transpose()*(x-center)/length)*rotation.transpose();};
        options.step=.08*length;options.samples=5;
        Problem scaled(transformed,options);
        const auto scaled_search=discover(scaled,shift);
        const auto scaled_curves=trace(scaled,scaled_search.seeds,shift);
        require(scaled_curves.curves.size()==1 && scaled_curves.curves[0].closed,"rotated/scaled circle failed");
        for(const auto& p:scaled_curves.curves[0].points) {
            const Point x=rotation.transpose()*(p.position-center)/length;
            require(std::abs(x.z())<1e-6 && std::abs(x.head<2>().norm()-1)<1e-6,"rotated/scaled geometry error");
        }

        Options limited;limited.samples=3;limited.step=.08;limited.max_points=3;
        Problem incomplete(explicit_field("crossing_line"),limited);
        const auto incomplete_seeds=discover(incomplete,1.25,&origin);
        const auto incomplete_trace=trace(incomplete,incomplete_seeds.seeds,1.25);
        require(!incomplete_trace.events.empty(),"point-limit event was lost");
        std::cout << "Critical curves: line, circle, multiple components, selected value, negative controls, rotation, scaling, and limits passed\n";
        return 0;
    } catch(const std::exception& e) {std::cerr << "FAIL: " << e.what() << '\n';return 1;}
}
