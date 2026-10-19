#pragma once

#include "critical_field.h"
#include <string>

namespace critical_curves {
struct Sample {
    Point position, tangent;
    double value, gradient_residual;
    Point eigenvalues;
};
struct Curve {
    std::vector<Sample> points;
    bool closed = false;
    std::string begin, end;
};
struct Event {
    Point position;
    std::string reason;
};
struct TraceResult {
    std::vector<Curve> curves;
    std::vector<Event> events;
    size_t covered_seeds = 0;
};

inline Sample sample(const Problem& problem, const Point& p, const Point& direction) {
    const auto s = classify(problem,p);
    Point t = s.eigenvectors.col(0);
    if(t.dot(direction) < 0) t = -t;
    return {p,t,problem.value(p),problem.gradient(p).norm(),s.eigenvalues};
}

inline bool valid_sample(const Problem& problem, const Point& p, double level) {
    if(!problem.contains(p)) return false;
    const double f = problem.value(p);
    const Point g = problem.gradient(p);
    return std::isfinite(f) && g.allFinite() &&
        std::abs(f-level) <= problem.amplitude*problem.options.value_tolerance &&
        g.norm() <= problem.options.gradient_tolerance && classify(problem,p).crossing;
}

// Keep Q fixed during each correction, so the plane is an equality constraint.
// All three gradient components participate, including the tangential one.
inline bool correct_plane(const Problem& problem, const Point& origin,
                          const Eigen::Matrix<double,3,2>& q, double radius, Point& p) {
    Eigen::Vector2d z = Eigen::Vector2d::Zero();
    p = origin;
    if(!problem.contains(p)) return false;
    for(int iteration=0;iteration<problem.options.iterations;++iteration) {
        const Point g = problem.gradient(p);
        if(!g.allFinite()) return false;
        if(g.norm() <= problem.options.gradient_tolerance) return true;
        const Eigen::Matrix<double,3,2> a = problem.hessian(p)*q;
        if(!a.allFinite()) return false;
        Eigen::JacobiSVD<Eigen::Matrix<double,3,2>> svd(a,Eigen::ComputeFullU | Eigen::ComputeFullV);
        if(svd.singularValues()[1] <= problem.options.hessian_absolute) return false;
        Eigen::Vector2d delta = problem.length*svd.solve(-g);
        bool advanced = false;
        for(int backtrack=0;backtrack<24;++backtrack) {
            const Eigen::Vector2d next_z = z+delta;
            const Point next = origin+q*next_z;
            if(next_z.norm() <= radius && problem.contains(next)) {
                const Point next_g = problem.gradient(next);
                if(next_g.allFinite() && next_g.squaredNorm() < g.squaredNorm()) {
                    z = next_z; p = next; advanced = true; break;
                }
            }
            delta *= .5;
        }
        if(!advanced) return false;
    }
    return problem.gradient(p).norm() <= problem.options.gradient_tolerance;
}

inline double segment_distance(const Point& p, const Point& a, const Point& b) {
    const Point d = b-a;
    if(d.squaredNorm() == 0) return (p-a).norm();
    const double u = std::clamp((p-a).dot(d)/d.squaredNorm(),0.,1.);
    return (p-(a+u*d)).norm();
}

inline bool covered(const Seed& seed, const std::vector<Curve>& curves, double distance) {
    for(const auto& curve: curves) {
        const size_t edges = curve.closed ? curve.points.size() : curve.points.size()-1;
        for(size_t i=0;i<edges;++i) {
            const auto& a = curve.points[i]; const auto& b = curve.points[(i+1)%curve.points.size()];
            if(std::abs(seed.signature.eigenvectors.col(0).dot(a.tangent)) > .95 &&
               segment_distance(seed.position,a.position,b.position) <= distance) return true;
        }
    }
    return false;
}

inline Curve trace_direction(const Problem& problem, const Seed& seed, double level, double sign) {
    Curve curve;
    const Point initial_direction = sign*seed.signature.eigenvectors.col(0);
    curve.points.push_back(sample(problem,seed.position,initial_direction));
    double step = problem.options.step, travelled = 0;
    for(int index=1;index<problem.options.max_points;++index) {
        const Sample current = curve.points.back();
        const Signature signature = classify(problem,current.position);
        const Eigen::Matrix<double,3,2> q = signature.eigenvectors.rightCols<2>();
        bool accepted = false;
        std::string failure = "correction_failed";
        for(int retry=0;retry<problem.options.retries;++retry) {
            double boundary_distance = std::numeric_limits<double>::infinity();
            int face = -1;
            for(int d=0;d<3;++d) {
                if(std::abs(current.tangent[d]) < 1e-14) continue;
                const double limit = current.tangent[d] > 0 ? problem.field.upper[d] : problem.field.lower[d];
                const double distance = (limit-current.position[d])/current.tangent[d];
                if(distance < boundary_distance) { boundary_distance = distance; face = d; }
            }
            // Finish transverse box contacts by correcting ON the selected face,
            // never by clamping a corrected curve point or evaluating outside.
            if(face >= 0 && boundary_distance <= step) {
                if(boundary_distance <= 1e-12*problem.length) { curve.end = "boundary"; return curve; }
                Point prediction = current.position+boundary_distance*current.tangent;
                prediction[face] = current.tangent[face] > 0 ? problem.field.upper[face] : problem.field.lower[face];
                Eigen::Matrix<double,3,2> face_basis;
                face_basis.col(0) = Point::Unit((face+1)%3);
                face_basis.col(1) = Point::Unit((face+2)%3);
                Point endpoint;
                if(correct_plane(problem,prediction,face_basis,step*.25,endpoint) &&
                   valid_sample(problem,endpoint,level)) {
                    const auto next = sample(problem,endpoint,current.tangent);
                    if(next.tangent.dot(current.tangent) >= .9 && (endpoint-current.position).norm() <= 1.25*step) {
                        curve.points.push_back(next); curve.end = "boundary"; return curve;
                    }
                }
            }
            const Point prediction = current.position+step*current.tangent;
            Point next_position;
            if(correct_plane(problem,prediction,q,step*.1,next_position)) {
                if(!valid_sample(problem,next_position,level)) failure = "value_or_hessian_check_failed";
                else {
                    const Sample next = sample(problem,next_position,current.tangent);
                    const double alignment = current.tangent.dot(next.tangent);
                    const double correction = (next_position-prediction).norm()/step;
                    if(alignment >= .9 && correction <= .1) {
                        travelled += (next_position-current.position).norm();
                        // Return to the seed with compatible orientation, after
                        // substantial progress. The OBJ/VTK closing edge reuses it.
                        if(curve.points.size() >= 8 && travelled > 4*problem.options.step &&
                           next.tangent.dot(initial_direction) > .95 &&
                           segment_distance(seed.position,current.position,next_position) <= .15*step &&
                           (next_position-seed.position).norm() < .75*step) {
                            if((current.position-seed.position).norm() > 1.25*step) curve.points.push_back(next);
                            curve.closed = true; curve.end = "closed"; return curve;
                        }
                        // Contacts with older parts of this trace remain unresolved;
                        // a small spatial distance does not authorize a branch merge.
                        for(size_t j=0;j+6<curve.points.size();++j)
                            if(segment_distance(next_position,curve.points[j].position,curve.points[j+1].position) < .05*step) {
                                curve.end = "self_contact"; return curve;
                            }
                        curve.points.push_back(next);
                        const double old_gap = current.eigenvalues.tail<2>().cwiseAbs().minCoeff();
                        const double new_gap = next.eigenvalues.tail<2>().cwiseAbs().minCoeff();
                        if(new_gap < .5*old_gap) step *= .5;
                        else if(alignment > .99 && correction < .025) step = std::min(problem.options.step,1.25*step);
                        accepted = true; break;
                    }
                    failure = "excessive_turn_or_correction";
                }
            }
            step *= .5;
            if(step < problem.options.min_step) { curve.end = failure+"_min_step"; return curve; }
        }
        if(!accepted) { curve.end = failure+"_retry_limit"; return curve; }
    }
    curve.end = "point_limit";
    return curve;
}

inline TraceResult trace(const Problem& problem, const std::vector<Seed>& seeds, double level) {
    if(!std::isfinite(level)) throw std::invalid_argument("A finite target isovalue is required");
    TraceResult result;
    for(const auto& seed: seeds) {
        if(!valid_sample(problem,seed.position,level)) {
            result.events.push_back({seed.position,"invalid_seed"}); continue;
        }
        // Tight coverage radius: tangent compatibility alone cannot distinguish
        // two close parallel curves. This is still a resolution-dependent test.
        if(covered(seed,result.curves,problem.options.step*.03)) { ++result.covered_seeds; continue; }
        Curve forward = trace_direction(problem,seed,level,1);
        if(forward.closed) { forward.begin = "closed"; result.curves.push_back(std::move(forward)); continue; }
        Curve backward = trace_direction(problem,seed,level,-1);
        if(backward.closed) { backward.begin = "closed"; result.curves.push_back(std::move(backward)); continue; }
        Curve combined; combined.begin = backward.end; combined.end = forward.end;
        for(auto i=backward.points.rbegin();i!=backward.points.rend();++i) {
            Sample s = *i; s.tangent = -s.tangent; combined.points.push_back(s);
        }
        combined.points.insert(combined.points.end(),forward.points.begin()+1,forward.points.end());
        double length = 0;
        for(size_t i=1;i<combined.points.size();++i) length += (combined.points[i].position-combined.points[i-1].position).norm();
        const bool boundary_to_boundary = combined.begin == "boundary" && combined.end == "boundary";
        if(combined.points.size() < 2 || (!boundary_to_boundary && length < 4*problem.options.step)) {
            result.events.push_back({seed.position,"unvalidated_short_trace:"+combined.begin+":"+combined.end}); continue;
        }
        if(combined.begin != "boundary") result.events.push_back({combined.points.front().position,combined.begin});
        if(combined.end != "boundary") result.events.push_back({combined.points.back().position,combined.end});
        result.curves.push_back(std::move(combined));
    }
    return result;
}
} // namespace critical_curves
