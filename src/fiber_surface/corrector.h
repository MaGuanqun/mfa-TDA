#pragma once
#include "field.h"

namespace fiber_surface {
class Corrector {
public:
    const Field& field;
    const Segment& segment;
    const Options& options;
    Corrector(const Field& f, const Segment& s, const Options& o) : field(f), segment(s), options(o) {}
    bool member(const Vec& x, int traced_endpoint=-1) const {
        if (!field.box.contains(x)) return false;
        const Range value=field.value(x);
        if (traced_endpoint==0 || traced_endpoint==1) {
            const Range& q=traced_endpoint==0 ? segment.a : segment.b;
            return (value-q).norm()<=options.corner_residual_tolerance;
        }
        return segment.contains_parameter(segment.parameter(value));
    }
    bool regular(const Jacobian& j) const {
        if (!j.allFinite()) return false;
        const Eigen::Matrix2d gram=j*j.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(gram);
        return eig.info()==Eigen::Success && eig.eigenvalues()[0] > options.derivative_tolerance*options.derivative_tolerance;
    }
    bool tangent(const Vec& x, Vec& t) const {
        if (!field.box.contains(x)) return false;
        const auto j=field.jacobian(x);
        if (!regular(j)) return false;
        t=Vec(j.row(0)).cross(Vec(j.row(1))).normalized(); return t.allFinite();
    }
    // Full minimum-norm Newton, optionally within a fixed box face/edge.
    // The caller handles boundary events; no model evaluation occurs outside the box.
    bool surface(Vec& x, double bound, unsigned faces=0, bool membership=true) const {
        const Vec origin=x; field.box.fix(x,faces);
        const double solve_tolerance=options.surface_solve_tolerance>0 ?
            std::min(options.residual_tolerance,options.surface_solve_tolerance) : options.residual_tolerance;
        for (int it=0; it<options.iterations; ++it) {
            if (!field.box.contains(x) || (x-origin).norm()>bound) return false;
            const Range f=field.value(x); const Vec g=segment.gradient(field.jacobian(x));
            if (!f.allFinite() || !g.allFinite() || g.norm()<=options.derivative_tolerance) return false;
            const double phi=segment.phi(f);
            Vec free=g; int free_coordinates=3;
            for (int d=0; d<3; ++d) if (faces & (3u<<(2*d))) { free[d]=0; --free_coordinates; }
            if (faces && free_coordinates && free.norm()<=options.derivative_tolerance) return false;
            if (std::abs(phi)<=solve_tolerance) return !membership || member(x);
            if (free.norm()<=options.derivative_tolerance) return false;
            x -= phi/free.squaredNorm()*free;
            field.box.fix(x,faces);
        }
        return false;
    }
    // Seed discovery uses a full minimum-norm Newton step.
    bool fiber_seed(Vec& x, const Range& q) const {
        if (!field.box.contains(x)) return false;
        Range r=field.value(x)-q;
        for (int it=0; it<options.iterations; ++it) {
            const Jacobian j=field.jacobian(x);
            if (!r.allFinite() || !regular(j)) return false;
            if (r.norm()<=options.corner_residual_tolerance) return true;
            const Vec trial=x-j.transpose()*(j*j.transpose()).ldlt().solve(r);
            if (!field.box.contains(trial)) return false;
            const Range next_r=field.value(trial)-q;
            if (!next_r.allFinite() || next_r.norm()>=r.norm()) return false;
            x=trial;
            r=next_r;
        }
        return false;
    }
    // F=q plus a fixed phase plane. At a contact, fixed box coordinates replace
    // the phase equation; QR also handles overdetermined edge/corner contacts.
    bool fiber(Vec& x, const Range& q, const Vec& phase, const Vec& anchor, double bound, unsigned faces=0) const {
        const Vec origin=x, fixed_anchor=anchor; field.box.fix(x,faces);
        for (int it=0; it<options.iterations; ++it) {
            if (!field.box.contains(x) || (x-origin).norm()>bound) return false;
            const Range r=field.value(x)-q; const Jacobian j=field.jacobian(x);
            if (!r.allFinite() || !regular(j)) return false;
            if (!faces) {
                const double p=phase.dot(x-fixed_anchor);
                if (r.norm()<=options.corner_residual_tolerance && std::abs(p)<=options.fiber_plane_tolerance()) return true;
                Hessian a; a.topRows<2>()=j; a.row(2)=phase.transpose();
                Eigen::FullPivLU<Hessian> lu(a); lu.setThreshold(options.derivative_tolerance);
                if (!lu.isInvertible()) return false;
                x += lu.solve(Vec(-r[0],-r[1],-p));
            } else {
                std::vector<int> free;
                for (int d=0; d<3; ++d) if (!(faces & (3u<<(2*d)))) free.push_back(d);
                if (free.empty()) return r.norm()<=options.corner_residual_tolerance;
                Eigen::MatrixXd a(2,free.size());
                for (size_t d=0; d<free.size(); ++d) a.col(d)=j.col(free[d]);
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a); qr.setThreshold(options.derivative_tolerance);
                if (qr.rank()!=int(free.size())) return false; // transversality on the face/edge
                if (r.norm()<=options.corner_residual_tolerance) return true;
                const Eigen::VectorXd delta=qr.solve(-r);
                if (!delta.allFinite()) return false;
                for (size_t d=0; d<free.size(); ++d) x[free[d]]+=delta[d];
                field.box.fix(x,faces);
            }
        }
        return false;
    }
};
} // namespace fiber_surface
