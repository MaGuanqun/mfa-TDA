#pragma once
#include "output.h"
#include <sstream>

namespace fiber_surface {
class Arguments {
    std::map<std::string,std::string> values_;
public:
    bool has(const std::string& key) const { return values_.count(key)!=0; }
    Arguments(int argc,char** argv) {
        for (int i=1; i<argc; ++i) {
            std::string key=argv[i];
            if (key=="--help" || key=="-h" || key=="--inspect") { values_[key]="true"; continue; }
            if (key.rfind("--",0)!=0 || i+1==argc) throw std::invalid_argument("Expected --option value: "+key);
            if (values_.count(key)) throw std::invalid_argument("Repeated option: "+key);
            values_[key]=argv[++i];
        }
    }
    std::string take(const std::string& key, const std::string& fallback="") {
        auto it=values_.find(key); if (it==values_.end()) return fallback;
        const std::string value=it->second; values_.erase(it); return value;
    }
    double number(const std::string& key,double fallback) {
        const auto str=take(key); if (str.empty()) return fallback;
        size_t used=0; double x=std::stod(str,&used);
        if (used!=str.size() || !std::isfinite(x)) throw std::invalid_argument("Invalid number for "+key); return x;
    }
    int integer(const std::string& key,int fallback) {
        double x=number(key,fallback);
        if (x!=std::floor(x) || x<std::numeric_limits<int>::min() || x>std::numeric_limits<int>::max()) throw std::invalid_argument("Invalid integer for "+key);
        return int(x);
    }
    bool help() { const auto a=take("--help"),b=take("-h"); return !a.empty() || !b.empty(); }
    void finish() const { if (!values_.empty()) throw std::invalid_argument("Unknown option: "+values_.begin()->first); }
};
template<int N> Eigen::Matrix<double,N,1> coordinates(std::string str) {
    std::replace(str.begin(),str.end(),',',' '); std::istringstream in(str); Eigen::Matrix<double,N,1> x;
    for (int i=0; i<N; ++i) if (!(in>>x[i]) || !std::isfinite(x[i])) throw std::invalid_argument("Expected "+std::to_string(N)+" finite comma-separated coordinates");
    std::string extra; if (in>>extra) throw std::invalid_argument("Too many coordinates"); return x;
}
inline void common_options(Arguments& args, Field& field, Segment& segment, Options& options) {
    const auto a=args.take("--a"),b=args.take("--b");
    if (a.empty()!=b.empty()) throw std::invalid_argument("Supply both --a and --b");
    if (!a.empty()) segment=Segment(coordinates<2>(a),coordinates<2>(b));
    const auto min=args.take("--min"),max=args.take("--max");
    if (min.empty()!=max.empty()) throw std::invalid_argument("Supply both --min and --max");
    if (!min.empty()) {
        const Vec lo=coordinates<3>(min),hi=coordinates<3>(max);
        if (!field.box.contains(lo) || !field.box.contains(hi)) throw std::invalid_argument("Requested box is outside the model domain");
        field.box.min=lo; field.box.max=hi;
    }
    options.step=args.number("--step",options.step);
    options.min_step=args.number("--min-step",options.step/8);
    options.fiber_step=args.number("--fiber-step",options.step*.75);
    options.spatial_tolerance=args.number("--spatial-tolerance",options.step*.075);
    options.residual_tolerance=args.number("--residual-tolerance",options.residual_tolerance);
    options.derivative_tolerance=args.number("--derivative-tolerance",options.derivative_tolerance);
    options.min_angle_degrees=args.number("--min-angle",options.min_angle_degrees);
    options.max_normal_degrees=args.number("--max-normal-angle",options.max_normal_degrees);
    options.seed_grid=args.integer("--seed-grid",options.seed_grid);
    options.iterations=args.integer("--iterations",options.iterations);
    options.max_trace_steps=args.integer("--max-trace-steps",options.max_trace_steps);
    const int limit=args.integer("--max-triangles",int(options.max_triangles));
    if (limit<4) throw std::invalid_argument("--max-triangles must be at least 4"); options.max_triangles=size_t(limit);
    options.validate(field,segment);
}
inline void common_help() {
    std::cout<<"  --a u,v --b u,v        finite range segment\n"
        <<"  --min x,y,z --max x,y,z  spatial sub-box\n"
        <<"  --output PREFIX       OBJ, VTK, fiber VTK and JSON report\n"
        <<"  --step H --min-step H --fiber-step H (physical units)\n"
        <<"  --spatial-tolerance E --residual-tolerance E --derivative-tolerance E\n"
        <<"  --seed-grid N --iterations N --max-trace-steps N --max-triangles N\n"
        <<"Exit 0: complete discovered mesh; 1: input/runtime error; 2: unresolved/invalid partial result.\n";
}
} // namespace fiber_surface
