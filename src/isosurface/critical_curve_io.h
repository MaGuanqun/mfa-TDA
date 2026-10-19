#pragma once
#include "critical_curves.h"
#include "../opts.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace critical_curves {
inline void add_options(opts::Options& cli, Options& o) {
    cli >> opts::Option("samples",o.samples,"seed cells per axis (2..200; default 9)");
    cli >> opts::Option("seed-offset",o.seed_offset,"within-cell seed offset [0,1] (default 0.5)");
    cli >> opts::Option("iterations",o.iterations,"iterations per local solve (default 100)");
    cli >> opts::Option("gradient-tolerance",o.gradient_tolerance,"tolerance for L/A times gradient norm");
    cli >> opts::Option("value-tolerance",o.value_tolerance,"value residual / A tolerance; also candidate grouping width");
    cli >> opts::Option("hessian-absolute",o.hessian_absolute,"absolute rank threshold for L^2/A times Hessian");
    cli >> opts::Option("hessian-relative",o.hessian_relative,"relative Hessian rank threshold");
    cli >> opts::Option("field-scale",o.field_scale,"positive reference A; 0 estimates from a fixed grid");
    cli >> opts::Option("step",o.step,"maximum physical tracing step; 0 uses shortest box side / 80");
    cli >> opts::Option("min-step",o.min_step,"minimum physical step; 0 uses step / 1024");
    cli >> opts::Option("max-points",o.max_points,"point limit per directed trace");
    cli >> opts::Option("retries",o.retries,"retry limit per tracing step");
}

inline std::ofstream output_file(const std::string& filename) {
    const auto parent = std::filesystem::path(filename).parent_path();
    if(!parent.empty()) std::filesystem::create_directories(parent);
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(filename);
    out << std::setprecision(17);
    return out;
}

inline void write_report(std::ostream& out, const Problem& p, const SearchResult& s,
                         const std::string& input, std::optional<double> level) {
    out << "input=" << input << "\n";
    if(level) out << "isovalue=" << *level << "\n";
    out << "domain_min=" << p.field.lower.transpose() << "\ndomain_max=" << p.field.upper.transpose()
        << "\nlength_scale=" << p.length << "\nfield_scale=" << p.amplitude
        << "\nsamples_per_axis=" << p.options.samples << "\nseed_offset=" << p.options.seed_offset
        << "\niterations=" << p.options.iterations << "\nstep=" << p.options.step << "\nmin_step=" << p.options.min_step
        << "\nmax_points=" << p.options.max_points << "\nretries=" << p.options.retries
        << "\ngradient_tolerance=" << p.options.gradient_tolerance << "\nvalue_tolerance=" << p.options.value_tolerance
        << "\nhessian_absolute=" << p.options.hessian_absolute << "\nhessian_relative=" << p.options.hessian_relative
        << "\nstarts=" << s.starts << "\nfailed_solves=" << s.failed
        << "\nfull_rank_convergences=" << s.full_rank << "\ndefinite_rank_two_convergences=" << s.definite
        << "\ndegenerate_convergences=" << s.degenerate << "\nduplicate_seeds=" << s.duplicates
        << "\ncandidate_seeds=" << s.seeds.size() << "\nvalue_evaluations=" << p.values
        << "\ngradient_evaluations=" << p.gradients << "\nhessian_evaluations=" << p.hessians
        << "\ncoverage=finite_seed_search_not_a_completeness_certificate\n";
}

inline void print_search(const Problem& p, const SearchResult& s) {
    std::cout << "Scales: L=" << p.length << ", A=" << p.amplitude << "; starts=" << s.starts
              << ", candidate seeds=" << s.seeds.size() << ", failed solves=" << s.failed
              << ", full-rank=" << s.full_rank << ", definite rank-two=" << s.definite
              << ", other degenerate=" << s.degenerate << ", duplicate seeds=" << s.duplicates << '\n';
}

struct ValueGroup {
    double value, minimum, maximum;
    std::vector<Seed> seeds;
};

inline std::vector<ValueGroup> group_values(const Problem& p, const SearchResult& search) {
    auto sorted = search.seeds;
    std::stable_sort(sorted.begin(),sorted.end(),[](const Seed& a,const Seed& b){return a.value < b.value;});
    std::vector<ValueGroup> groups;
    for(const auto& seed: sorted) {
        // Bound the whole group's diameter; pairwise chaining can combine levels
        // that differ by many tolerances. Grouping does NOT imply connectivity.
        if(groups.empty() || seed.value-groups.back().minimum > p.amplitude*p.options.value_tolerance)
            groups.push_back({seed.value,seed.value,seed.value,{}});
        groups.back().seeds.push_back(seed); groups.back().maximum = seed.value;
    }
    for(auto& group: groups) group.value = group.seeds[group.seeds.size()/2].value;
    return groups;
}

inline void write_values(const std::string& prefix, const Problem& p, const SearchResult& search,
                         const std::string& input) {
    const auto groups = group_values(p,search);
    auto catalog = output_file(prefix+".values.csv");
    catalog << "candidate_id,value,minimum,maximum,seed_count,status\n";
    auto diagnostics = output_file(prefix+".seeds.csv");
    diagnostics << "candidate_id,x,y,z,value,scaled_gradient_norm,lambda0,lambda1,lambda2\n";
    for(size_t i=0;i<groups.size();++i) {
        const auto& group = groups[i];
        catalog << i << ',' << group.value << ',' << group.minimum << ',' << group.maximum << ','
                << group.seeds.size() << ",candidate_not_curve_validated\n";
        auto seeds = output_file(prefix+".value_"+std::to_string(i)+".seeds.txt");
        seeds << "# Candidate critical value " << group.value << "; positions require revalidation\n";
        for(const auto& seed: group.seeds) {
            seeds << seed.position.transpose() << '\n';
            diagnostics << i;
            for(int d=0;d<3;++d) diagnostics << ',' << seed.position[d];
            diagnostics << ',' << seed.value << ',' << seed.gradient_residual;
            for(int d=0;d<3;++d) diagnostics << ',' << seed.signature.eigenvalues[d];
            diagnostics << '\n';
        }
        std::cout << "Candidate " << i << ": value=" << std::setprecision(17) << group.value
                  << ", spread=" << group.maximum-group.minimum << ", seeds=" << group.seeds.size() << '\n';
    }
    auto report = output_file(prefix+".report.txt"); write_report(report,p,search,input,std::nullopt);
    report << "candidate_values=" << groups.size() << "\n";
    std::cout << "Wrote " << groups.size() << " candidate values to " << prefix << ".values.csv\n";
    if(groups.empty()) std::cout << "No crossing-curve candidates found at these search settings; this does not prove absence.\n";
}

inline void write_curves(const std::string& prefix, const Problem& p, const SearchResult& search,
                         const TraceResult& result, const std::string& input, double level) {
    auto obj = output_file(prefix+".obj");
    auto vtk = output_file(prefix+".vtk");
    auto points = output_file(prefix+".points.csv");
    auto curves = output_file(prefix+".curves.csv");
    auto events = output_file(prefix+".events.csv");
    auto report = output_file(prefix+".report.txt");
    obj << "# Critical-curve polylines; target value " << level << '\n';
    points << "curve_id,point_id,x,y,z,value,scaled_gradient_norm,lambda0,lambda1,lambda2\n";
    curves << "curve_id,point_count,closed,begin_status,end_status,min_value,max_value,max_scaled_gradient\n";
    events << "x,y,z,reason\n";
    size_t count = 0, line_entries = 0, closed = 0;
    for(const auto& c: result.curves) {count += c.points.size(); line_entries += 1+c.points.size()+(c.closed ? 1 : 0); closed += c.closed;}
    vtk << "# vtk DataFile Version 3.0\nCritical curves\nASCII\nDATASET POLYDATA\nPOINTS " << count << " double\n";
    size_t offset = 0;
    double worst_gradient = 0, worst_value = 0;
    for(size_t i=0;i<result.curves.size();++i) {
        const auto& c = result.curves[i];
        obj << "o curve_" << i << '\n';
        double lo = std::numeric_limits<double>::infinity(), hi = -lo, max_g = 0;
        for(size_t j=0;j<c.points.size();++j) {
            const auto& s = c.points[j];
            obj << "v " << s.position.transpose() << '\n'; vtk << s.position.transpose() << '\n';
            points << i << ',' << j;
            for(int d=0;d<3;++d) points << ',' << s.position[d];
            points << ',' << s.value << ',' << s.gradient_residual;
            for(int d=0;d<3;++d) points << ',' << s.eigenvalues[d];
            points << '\n';
            lo = std::min(lo,s.value); hi = std::max(hi,s.value); max_g = std::max(max_g,s.gradient_residual);
            worst_value = std::max(worst_value,std::abs(s.value-level));
        }
        obj << "l";
        for(size_t j=0;j<c.points.size();++j) obj << ' ' << offset+j+1;
        if(c.closed) obj << ' ' << offset+1;
        obj << '\n'; offset += c.points.size();
        curves << i << ',' << c.points.size() << ',' << c.closed << ',' << c.begin << ',' << c.end
               << ',' << lo << ',' << hi << ',' << max_g << '\n';
        worst_gradient = std::max(worst_gradient,max_g);
    }
    vtk << "LINES " << result.curves.size() << ' ' << line_entries << '\n'; offset = 0;
    for(const auto& c: result.curves) {
        vtk << c.points.size()+(c.closed ? 1 : 0);
        for(size_t j=0;j<c.points.size();++j) vtk << ' ' << offset+j;
        if(c.closed) vtk << ' ' << offset;
        vtk << '\n'; offset += c.points.size();
    }
    for(const auto& e: result.events) events << e.position.x() << ',' << e.position.y() << ',' << e.position.z() << ',' << e.reason << '\n';
    write_report(report,p,search,input,level);
    report << "curve_components=" << result.curves.size() << "\nclosed_components=" << closed
           << "\ncurve_points=" << count << "\ncovered_seeds=" << result.covered_seeds
           << "\nunresolved_events=" << result.events.size() << "\nmaximum_value_residual=" << worst_value
           << "\nmaximum_scaled_gradient=" << worst_gradient << '\n';
    std::cout << "Wrote " << result.curves.size() << " curve components (" << closed << " closed), "
              << count << " points, " << result.events.size() << " unresolved events to " << prefix << ".{obj,vtk,...}\n";
    if(result.curves.empty()) std::cout << "No supported critical curves traced at the selected value and tolerances.\n";
}
} // namespace critical_curves
