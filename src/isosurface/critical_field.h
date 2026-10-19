#pragma once

#include "surface_field.h"
#include <array>
#include <optional>
#include <vector>

namespace critical_curves {
using Point = Eigen::Vector3d;
using Hessian = Eigen::Matrix3d;

// Separate from SurfaceField: regular surface meshing does not need a Hessian.
struct Field {
    std::function<double(const Point&)> value;
    std::function<Point(const Point&)> gradient;
    std::function<Hessian(const Point&)> hessian;
    Point lower = Point::Constant(-2);
    Point upper = Point::Constant(2);
};

struct Options {
    int samples = 9;
    int iterations = 100;
    int max_points = 20000;             // per directed trace
    int retries = 20;
    double value_tolerance = 1e-9;      // dimensionless, after scaling by A
    double gradient_tolerance = 1e-9;   // after scaling by L/A
    double hessian_absolute = 1e-9;     // after scaling by L^2/A
    double hessian_relative = 1e-6;
    double field_scale = 0;             // 0: estimate once from a fixed 5^3 grid
    double step = 0;                    // physical distance; 0: shortest side / 80
    double min_step = 0;                // 0: step / 1024
    double seed_offset = .5;            // cell-centered seeds, in [0,1]

    void validate() const {
        if(samples < 2 || samples > 200 || iterations < 1 || max_points < 2 || retries < 1 ||
           !std::isfinite(value_tolerance) || value_tolerance <= 0 ||
           !std::isfinite(gradient_tolerance) || gradient_tolerance <= 0 ||
           !std::isfinite(hessian_absolute) || hessian_absolute <= 0 ||
           !std::isfinite(hessian_relative) || hessian_relative < 0 || hessian_relative >= 1 ||
           !std::isfinite(field_scale) || field_scale < 0 ||
           !std::isfinite(step) || step < 0 || !std::isfinite(min_step) || min_step < 0 ||
           (step > 0 && min_step > step) || !std::isfinite(seed_offset) || seed_offset < 0 || seed_offset > 1)
            throw std::invalid_argument("Invalid critical-curve sampling, tolerance, step, or iteration option");
    }
};

// Use an isotropic coordinate scale, preserving Euclidean tangent directions.
// All stored positions and step sizes remain in physical coordinates.
class Problem {
public:
    Field field;
    Options options;
    double length, amplitude;
    mutable size_t values = 0, gradients = 0, hessians = 0;

    Problem(Field f, Options o) : field(std::move(f)), options(o) {
        options.validate();
        if(!field.value || !field.gradient || !field.hessian || !field.lower.allFinite() ||
           !field.upper.allFinite() || !(field.upper.array() > field.lower.array()).all())
            throw std::invalid_argument("Critical curves require a finite box, value, gradient, and Hessian");
        const Point extent = field.upper-field.lower;
        length = extent.maxCoeff();
        amplitude = options.field_scale;
        if(amplitude == 0) {
            double lo = std::numeric_limits<double>::infinity(), hi = -lo, max_gradient = 0;
            for(int i=0;i<5;++i) for(int j=0;j<5;++j) for(int k=0;k<5;++k) {
                const Point p = field.lower+extent.cwiseProduct(Point(i,j,k)/4.);
                const double v = value(p);
                const Point g = field.gradient(p); ++gradients;
                if(!std::isfinite(v) || !g.allFinite()) throw std::runtime_error("Non-finite field during scale estimation");
                lo = std::min(lo,v); hi = std::max(hi,v); max_gradient = std::max(max_gradient,g.norm());
            }
            amplitude = std::max(hi-lo,length*max_gradient);
            if(amplitude == 0) amplitude = 1;
        }
        if(!std::isfinite(amplitude) || !std::isfinite(length) || amplitude <= 0)
            throw std::invalid_argument("Invalid field or coordinate scale");
        if(options.step == 0) options.step = extent.minCoeff()/80.;
        if(options.min_step == 0) options.min_step = options.step/1024.;
        options.validate();
    }
    bool contains(const Point& p) const {
        return p.allFinite() && (p.array() >= field.lower.array()).all() && (p.array() <= field.upper.array()).all();
    }
    double value(const Point& p) const { ++values; return field.value(p); }
    Point gradient(const Point& p) const { ++gradients; return field.gradient(p)*(length/amplitude); }
    Hessian hessian(const Point& p) const {
        ++hessians;
        const Hessian h = field.hessian(p)*(length*length/amplitude);
        return .5*(h+h.transpose());
    }
};

struct Signature {
    Point eigenvalues = Point::Zero(); // ordered by magnitude, not algebraically
    Hessian eigenvectors = Hessian::Identity();
    int rank = -1;
    bool crossing = false;
};

inline Signature classify(const Problem& problem, const Point& p) {
    Signature result;
    const Hessian h = problem.hessian(p);
    if(!h.allFinite()) return result;
    Eigen::SelfAdjointEigenSolver<Hessian> eig(h);
    if(eig.info() != Eigen::Success) return result;
    std::array<int,3> order{0,1,2};
    std::sort(order.begin(),order.end(),[&](int a,int b){return std::abs(eig.eigenvalues()[a]) < std::abs(eig.eigenvalues()[b]);});
    for(int i=0;i<3;++i) {
        result.eigenvalues[i] = eig.eigenvalues()[order[i]];
        result.eigenvectors.col(i) = eig.eigenvectors().col(order[i]);
    }
    const double threshold = problem.options.hessian_absolute +
        problem.options.hessian_relative*result.eigenvalues.cwiseAbs().maxCoeff();
    result.rank = (result.eigenvalues.cwiseAbs().array() > threshold).count();
    result.crossing = result.rank == 2 && result.eigenvalues[1]*result.eigenvalues[2] < 0;
    return result;
}

struct Seed {
    Point position;
    double value;
    double gradient_residual;
    Signature signature;
};

struct SearchResult {
    std::vector<Seed> seeds;
    size_t starts = 0, failed = 0, full_rank = 0, definite = 0, degenerate = 0, duplicates = 0;
};

// Rank-deficient Levenberg-Marquardt using an augmented SVD solve. Unlike the
// legacy critical_point/find_all_root.h Newton routine, singular H is allowed.
// The optional level adds F-c to the residual for target-level seed discovery.
inline bool refine_seed(const Problem& problem, Point& p, std::optional<double> level) {
    if(!problem.contains(p)) return false;
    double damping = 1e-6;
    for(int iteration=0;iteration<problem.options.iterations;++iteration) {
        const Point g = problem.gradient(p);
        const double f = level ? (problem.value(p)-*level)/problem.amplitude : 0;
        if(!g.allFinite() || !std::isfinite(f)) return false;
        if(g.norm() <= problem.options.gradient_tolerance && std::abs(f) <= problem.options.value_tolerance) return true;
        const Hessian h = problem.hessian(p);
        if(!h.allFinite()) return false;
        Eigen::Matrix<double,4,3> jacobian = Eigen::Matrix<double,4,3>::Zero();
        Eigen::Vector4d residual; residual << g,f;
        jacobian.topRows<3>() = h;
        if(level) jacobian.row(3) = g.transpose();
        const double scale = std::max(jacobian.squaredNorm(),1e-24);
        bool advanced = false;
        for(int trial=0;trial<12;++trial) {
            Eigen::Matrix<double,7,3> a;
            a.topRows<4>() = jacobian;
            a.bottomRows<3>() = std::sqrt(damping*scale)*Hessian::Identity();
            Eigen::Matrix<double,7,1> b = Eigen::Matrix<double,7,1>::Zero();
            b.head<4>() = -residual;
            Point delta = a.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
            if(!delta.allFinite()) return false;
            if(delta.norm() > .2) delta *= .2/delta.norm();
            for(int backtrack=0;backtrack<24;++backtrack) {
                const Point q = p+problem.length*delta;
                if(problem.contains(q)) {
                    const Point qg = problem.gradient(q);
                    const double qf = level ? (problem.value(q)-*level)/problem.amplitude : 0;
                    if(qg.allFinite() && std::isfinite(qf) && qg.squaredNorm()+qf*qf < residual.squaredNorm()) {
                        p = q; advanced = true; damping = std::max(1e-14,damping*.2); break;
                    }
                }
                delta *= .5;
            }
            if(advanced) break;
            damping *= 10;
        }
        if(!advanced) return false;
    }
    return problem.gradient(p).norm() <= problem.options.gradient_tolerance &&
        (!level || std::abs(problem.value(p)-*level) <= problem.amplitude*problem.options.value_tolerance);
}

inline SearchResult discover(const Problem& problem, std::optional<double> level = std::nullopt,
                             const std::vector<Point>* supplied = nullptr) {
    SearchResult result;
    const auto attempt = [&](Point p) {
        ++result.starts;
        if(!refine_seed(problem,p,level)) { ++result.failed; return; }
        const Signature s = classify(problem,p);
        if(!s.crossing) {
            if(s.rank == 3) ++result.full_rank;
            else if(s.rank == 2) ++result.definite;
            else ++result.degenerate;
            return;
        }
        const double value = problem.value(p);
        if(!std::isfinite(value)) { ++result.failed; return; }
        for(const auto& seed: result.seeds)
            if((seed.position-p).norm() <= problem.options.step*.05 &&
               std::abs(seed.signature.eigenvectors.col(0).dot(s.eigenvectors.col(0))) > .99) {
                ++result.duplicates; return;
            }
        result.seeds.push_back({p,value,problem.gradient(p).norm(),s});
    };
    if(supplied) { for(const auto& p: *supplied) attempt(p); }
    else {
        const int n = problem.options.samples;
        const Point extent = problem.field.upper-problem.field.lower;
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) for(int k=0;k<n;++k)
            attempt(problem.field.lower+extent.cwiseProduct((Point(i,j,k)+Point::Constant(problem.options.seed_offset))/n));
    }
    return result;
}
} // namespace critical_curves
