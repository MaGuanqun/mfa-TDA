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
inline Box coordnet_training_box(Arguments& args) {
    if (!args.has("--domain-min") || !args.has("--domain-max"))
        throw std::invalid_argument("Supply the full CoordNet training box with --domain-min and --domain-max");
    return {coordinates<3>(args.take("--domain-min")), coordinates<3>(args.take("--domain-max"))};
}
inline RangeCurve range_selection(Arguments& args, const Segment& fallback) {
    if (args.has("--polyline")) {
        if (args.has("--a") || args.has("--b")) throw std::invalid_argument("Choose either --polyline or --a and --b");
        const auto input=args.take("--polyline");
        if (input.empty() || input.back()==';') throw std::invalid_argument("Expected vertices in --polyline");
        std::istringstream values(input);
        std::vector<Range> points;
        std::string point;
        while (std::getline(values,point,';')) points.push_back(coordinates<2>(point));
        return RangeCurve(std::move(points));
    }
    if (args.has("--a")!=args.has("--b")) throw std::invalid_argument("Supply both --a and --b");
    if (!args.has("--a")) return fallback;
    const Range a=coordinates<2>(args.take("--a"));
    return Segment(a,coordinates<2>(args.take("--b")));
}
inline void common_options(Arguments& args, Field& field, Options& options) {
    const auto min=args.take("--min"),max=args.take("--max");
    if (min.empty()!=max.empty()) throw std::invalid_argument("Supply both --min and --max");
    if (!min.empty()) {
        const Vec lo=coordinates<3>(min),hi=coordinates<3>(max);
        if (!field.box.contains(lo) || !field.box.contains(hi)) throw std::invalid_argument("Requested box is outside the model domain");
        field.box.min=lo; field.box.max=hi;
    }
    const bool has_fiber_step=args.has("--fiber-step");
    const double fiber_step=args.number("--fiber-step",options.step);
    options.step=args.number("--step",fiber_step);
    if (has_fiber_step && fiber_step!=options.step)
        throw std::invalid_argument("--step and --fiber-step must agree. Use --step for both tracing and marching.");
    options.min_step=args.number("--min-step",options.step/8);
    options.spatial_tolerance=args.number("--spatial-tolerance",options.step*.075);
    if (args.has("--fiber-seed-radius"))
        options.fiber_seed_radius=args.number("--fiber-seed-radius",0);
    options.residual_tolerance=args.number("--residual-tolerance",options.residual_tolerance);
    options.corner_residual_tolerance=args.number("--fiber-field-tolerance",options.corner_residual_tolerance);
    if (args.has("--surface-solve-tolerance")) {
        options.surface_solve_tolerance=args.number("--surface-solve-tolerance",0);
        if (options.surface_solve_tolerance<=0) throw std::invalid_argument("--surface-solve-tolerance must be positive");
    }
    options.derivative_tolerance=args.number("--derivative-tolerance",options.derivative_tolerance);
    options.min_angle_degrees=args.number("--min-angle",options.min_angle_degrees);
    options.max_normal_degrees=args.number("--max-normal-angle",options.max_normal_degrees);
    options.seed_grid=args.integer("--seed-grid",options.seed_grid);
    options.iterations=args.integer("--iterations",options.iterations);
    options.max_trace_steps=args.integer("--max-trace-steps",options.max_trace_steps);
    const int limit=args.integer("--max-triangles",int(options.max_triangles));
    if (limit<4) throw std::invalid_argument("--max-triangles must be at least 4"); options.max_triangles=size_t(limit);
}
inline RangeCurve surface_options(Arguments& args, Field& field, const Segment& fallback, Options& options) {
    auto curve=range_selection(args,fallback);
    common_options(args,field,options);
    for (const auto& s:curve.segments) options.validate(field,s);
    return curve;
}
inline void common_help() {
    std::cout<<"  --a u,v --b u,v        finite range segment\n"
        <<"  --polyline 'u,v;u,v;...'  ordered range vertices; repeat the first to close\n"
        <<"  --min x,y,z --max x,y,z  spatial sub-box\n"
        <<"  --output PREFIX       OBJ, VTK, fiber VTK and JSON report\n"
        <<"  --step H                 shared fiber tracing and surface marching step (physical units)\n"
        <<"  --min-step H             minimum fiber retry step and marching height recovery threshold\n"
        <<"  --fiber-seed-radius R    endpoint seed duplicate-check radius (defaults to spatial tolerance)\n"
        <<"  --spatial-tolerance E --residual-tolerance E --derivative-tolerance E\n"
        <<"  --fiber-field-tolerance E    maximum range residual on a corner fiber (default 1e-9)\n"
        <<"  --surface-solve-tolerance E  surface root residual (default cap 1e-10 on triangle patches)\n"
        <<"  --seed-grid N         N seeds on the shortest axis; other counts = ceil(N * length / shortest)\n"
        <<"                        default 7 for all fields; 0 also selects 7\n"
        <<"                        uses the selected physical box; maximum 1000000 starts\n"
        <<"  --iterations N --max-trace-steps N --max-triangles N (per segment)\n"
        <<"Exit 0: complete discovered mesh; 1: input/runtime error; 2: unresolved/invalid partial result.\n";
}
} // namespace fiber_surface
