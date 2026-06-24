#pragma once
//
// Stable Feature Flow Fields (Weinkauf, Theisel, Van Gelder & Pang,
// IEEE TVCG 17(6), 2011) for critical-point tracking in 2D time-dependent
// vector fields (domain D = 3: 2D space + time, C = 2 components).
//
// This builds on the plain Feature Flow Field reimplementation in
// feature_flow_fields.h. The original files are kept intact so plain FFF can
// still be run; this is a separate, self-contained copy with the stabilization
// added on top.
//
// Background (plain FFF, Theisel & Seidel 2003):
//   - Input is a sampled vector field v on a regular space-time grid (here
//     v = grad_space(s), stored in a .vff file).
//   - v is reconstructed by MULTILINEAR interpolation; the in-cell ANALYTIC
//     derivatives are stacked into A(c,a) = d v_c / d x_a.
//   - The FFF is the cross product of the component gradients (eq. (5) of the
//     stable-FFF paper):
//         f[k] = (-1)^k * det( A_spatial with column k removed )
//     The feature line (path of v == 0) is a stream line of f, but its
//     neighborhood may be a source/saddle, so numerical drift is amplified.
//
// Stabilization (this file, paper sec. 4):
//   Add a correction field g that vanishes on the feature line but makes its
//   neighborhood an attracting sink, then integrate h = f + tau * g instead of
//   f. Since g == 0 on the feature line, h has the SAME feature lines as f
//   (no false positives), but self-corrects numerical error.
//
//   For v = (a, b) with gradients grad a, grad b (rows of A), f = grad a x grad b:
//     correction vector (paper eq. (6)):
//         G   = a * grad b - b * grad a          (G[k] = a*b_{x_k} - b*a_{x_k})
//         g   = (f x G) / ||f||
//       -> on a feature line, grad g has eigenvalues 0, -||f||, -||f|| (the
//          perpendicular plane is an isotropic sink; eigenvector f for 0).
//     adaptive strength (paper eqs. (23)/(25)), needs grad f:
//         s_f = div( f / ||f|| )
//         d_f = [ det(f, f_y, f_t); det(f_x, f, f_t); det(f_x, f_y, f) ]   (eq. 3)
//         p_f = (f . d_f) / ||f||^4
//         R   = sqrt( max(s_f^2 - 4 p_f, 0) )                 ( = Re sqrt(...) )
//         forward  branch: tau = (k + s_f)/2 + R/2 ,  h =  f + tau g
//         backward branch: tau = (k - s_f)/2 + R/2 ,  h = -f + tau g
//       with tau clamped to [0, tau_max]. k > 0 controls convergence strength.
//
//   grad f (i.e. second derivatives of v) is estimated by CENTRAL FINITE
//   DIFFERENCES of f -- the paper notes a rough estimate suffices since it only
//   affects the *strength* of the convergence, not the on-feature direction.
//
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <unordered_set>

// mfa/types.hpp defines the global-namespace VectorX/MatrixX alias templates
// that this header (and CP_Trace.h below) rely on; it also pulls in Eigen.
// Do NOT `using namespace Eigen;` here (see feature_flow_fields.h for why).
#include <mfa/types.hpp>

#include "CP_Trace.h"

namespace sfff
{

// ---------------------------------------------------------------------------
// Discrete vector field on a regular space-time grid.
//   - D domain axes, axis 0 varies fastest in memory, last axis = time.
//   - C components per node (= D-1 for a spatial-gradient field).
//   - data layout (AoS, components innermost):
//       offset = ( i0 + n0*( i1 + n1*( i2 + ... ) ) ) * C + c
// ---------------------------------------------------------------------------
template<typename T>
struct VectorFieldGrid
{
    int                 D = 0;      // number of domain axes (last = time)
    int                 C = 0;      // components per node (= D-1)
    std::vector<int>    n;          // grid resolution per axis (size D), axis 0 fastest
    VectorX<T>          dmin;       // domain minimum per axis (size D)
    VectorX<T>          dmax;       // domain maximum per axis (size D)
    VectorX<T>          spacing;    // (dmax-dmin)/(n-1) per axis (size D)
    std::vector<T>      data;       // C * prod(n) values

    void compute_spacing()
    {
        spacing.resize(D);
        for (int a = 0; a < D; ++a)
            spacing[a] = (n[a] > 1) ? (dmax[a] - dmin[a]) / static_cast<T>(n[a] - 1)
                                    : static_cast<T>(1);
    }

    // linear node index (x fastest, time slowest) for integer coords idx[0..D-1]
    inline std::size_t node_linear(const std::array<int, 4>& idx) const
    {
        std::size_t lin = 0;
        for (int a = D - 1; a >= 0; --a)
            lin = lin * static_cast<std::size_t>(n[a]) + static_cast<std::size_t>(idx[a]);
        return lin;
    }

    inline bool in_domain(const VectorX<T>& p, T tol = static_cast<T>(1e-9)) const
    {
        for (int a = 0; a < D; ++a)
        {
            if (p[a] < dmin[a] - tol * spacing[a]) return false;
            if (p[a] > dmax[a] + tol * spacing[a]) return false;
        }
        return true;
    }

    // Multilinear interpolation of v and/or the analytic in-cell derivatives.
    //   A (if non-null): C x D matrix, A(c,a) = d v_c / d x_a
    //   v (if non-null): C vector, interpolated value
    // Returns false if p is outside the domain.
    bool sample(const VectorX<T>& p, MatrixX<T>* A, VectorX<T>* v) const
    {
        if (!in_domain(p)) return false;

        std::array<int, 4> base{};
        std::array<T, 4>   frac{};
        for (int a = 0; a < D; ++a)
        {
            T local = (p[a] - dmin[a]) / spacing[a];
            int i = static_cast<int>(std::floor(local));
            if (i < 0)          i = 0;
            if (i > n[a] - 2)   i = std::max(0, n[a] - 2);
            base[a] = i;
            T f = local - static_cast<T>(i);
            if (f < T(0)) f = T(0);
            if (f > T(1)) f = T(1);
            frac[a] = f;
        }

        if (v) { v->setZero(C); }
        if (A) { A->setZero(C, D); }

        const int corners = 1 << D;
        for (int m = 0; m < corners; ++m)
        {
            // value weight for this corner
            T w = T(1);
            for (int a = 0; a < D; ++a)
                w *= ((m >> a) & 1) ? frac[a] : (T(1) - frac[a]);

            std::array<int, 4> idx{};
            for (int a = 0; a < D; ++a)
                idx[a] = base[a] + ((m >> a) & 1);

            const T* node = &data[node_linear(idx) * static_cast<std::size_t>(C)];

            if (v)
                for (int c = 0; c < C; ++c)
                    (*v)[c] += w * node[c];

            if (A)
            {
                for (int a = 0; a < D; ++a)
                {
                    // derivative weight wrt axis a: replace the (1-/+) factor of axis a
                    T dw = T(1);
                    for (int b = 0; b < D; ++b)
                    {
                        int bit = (m >> b) & 1;
                        if (b == a) dw *= bit ? T(1) : T(-1);
                        else        dw *= bit ? frac[b] : (T(1) - frac[b]);
                    }
                    dw /= spacing[a];
                    for (int c = 0; c < C; ++c)
                        (*A)(c, a) += dw * node[c];
                }
            }
        }
        return true;
    }

    // Plain feature flow field direction at p (un-normalized). Returns false if
    // outside. Uses only the first (D-1) spatial-gradient rows of A (see the
    // note in feature_flow_fields.h about non-square minors).
    bool fff_direction(const VectorX<T>& p, VectorX<T>& f) const
    {
        MatrixX<T> A;
        if (!sample(p, &A, nullptr)) return false;
        fff_from_A(A, f);
        return true;
    }

    // f from a sampled Jacobian A (shared by fff_direction and stable_fff_unit).
    void fff_from_A(const MatrixX<T>& A, VectorX<T>& f) const
    {
        const int Csp = D - 1;                 // spatial-gradient rows actually used
        f.resize(D);
        MatrixX<T> M(Csp, D - 1);
        for (int k = 0; k < D; ++k)
        {
            int col = 0;
            for (int j = 0; j < D; ++j)
                if (j != k) M.col(col++) = A.topRows(Csp).col(j);
            T det = M.determinant();
            f[k] = ((k & 1) ? T(-1) : T(1)) * det;
        }
    }

    // Unit-normalized plain feature flow direction (full space-time norm).
    bool fff_unit(const VectorX<T>& p, VectorX<T>& dir, T eps = static_cast<T>(1e-30)) const
    {
        if (!fff_direction(p, dir)) return false;
        T nrm = dir.norm();
        if (!(nrm > eps)) return false;
        dir /= nrm;
        return true;
    }

    // -----------------------------------------------------------------------
    // Stable feature flow field direction (Weinkauf et al. 2011), D == 3 case.
    //   forward == true  -> forward  branch  (h =  f + tau g)
    //   forward == false -> backward branch  (h = -f + tau g)
    // k > 0 sets convergence strength; tau is clamped to [0, tau_max]; grad f
    // is estimated with central finite differences using step fd_frac*spacing.
    // Falls back to the plain (normalized) +/- f direction when stabilization is
    // unavailable (D != 3, k <= 0, near a critical point of f, or finite
    // differences leaving the domain). Returns false only when p is outside the
    // domain or f is (near) zero.
    // -----------------------------------------------------------------------
    bool stable_fff_unit(const VectorX<T>& p, bool forward, T k, T tau_max,
                         T fd_frac, VectorX<T>& dir) const
    {
        MatrixX<T> A;
        VectorX<T> v;
        if (!sample(p, &A, &v)) return false;

        VectorX<T> f;
        fff_from_A(A, f);
        T fn = f.norm();
        if (!(fn > static_cast<T>(1e-30))) return false;

        // Stable correction is formulated for the 3D space-time case only.
        if (D != 3 || k <= T(0))
        {
            dir = (forward ? f : -f) / fn;
            return true;
        }

        // correction vector G = a*grad(b) - b*grad(a), with a = v[0], b = v[1].
        VectorX<T> G(3);
        for (int a = 0; a < 3; ++a)
            G[a] = v[0] * A(1, a) - v[1] * A(0, a);

        // g = (f x G) / ||f||  -- vanishes on the feature line, perpendicular sink.
        VectorX<T> g = cross3(f, G) / fn;

        // central-difference partials of the (un-normalized) f along each axis.
        VectorX<T> fx, fy, ft;
        bool ok = fd_partial(p, 0, fd_frac, fx)
               && fd_partial(p, 1, fd_frac, fy)
               && fd_partial(p, 2, fd_frac, ft);

        T tau = T(0);
        if (ok)
        {
            const T fn3 = fn * fn * fn;
            const T fn4 = fn * fn3;

            // s_f = div(f/||f||): sum_a [ f_{x_a}[a]/||f|| - f[a]*(f . f_{x_a})/||f||^3 ]
            const VectorX<T>* fa[3] = { &fx, &fy, &ft };
            T sf = T(0);
            for (int a = 0; a < 3; ++a)
            {
                T dot = f.dot(*fa[a]);
                sf += (*fa[a])[a] / fn - f[a] * dot / fn3;
            }

            // d_f and p_f = (f . d_f)/||f||^4
            VectorX<T> df(3);
            df[0] = det3(f, fy, ft);
            df[1] = det3(fx, f, ft);
            df[2] = det3(fx, fy, f);
            T pf = f.dot(df) / fn4;

            T disc = sf * sf - T(4) * pf;
            T R = (disc > T(0)) ? std::sqrt(disc) : T(0);
            tau = forward ? (k + sf) / T(2) + R / T(2)
                          : (k - sf) / T(2) + R / T(2);
            if (tau < T(0))      tau = T(0);
            if (tau > tau_max)   tau = tau_max;
        }

        VectorX<T> base = forward ? VectorX<T>(f) : VectorX<T>(-f);
        VectorX<T> h = base + tau * g;
        T hn = h.norm();
        if (!(hn > static_cast<T>(1e-30))) { dir = base / fn; return true; }
        dir = h / hn;
        return true;
    }

private:
    // 3-vector cross product.
    static inline VectorX<T> cross3(const VectorX<T>& a, const VectorX<T>& b)
    {
        VectorX<T> c(3);
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
        return c;
    }

    // determinant of the 3x3 matrix whose columns are c0, c1, c2.
    static inline T det3(const VectorX<T>& c0, const VectorX<T>& c1, const VectorX<T>& c2)
    {
        return c0[0] * (c1[1] * c2[2] - c1[2] * c2[1])
             - c1[0] * (c0[1] * c2[2] - c0[2] * c2[1])
             + c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]);
    }

    // Central finite-difference partial of f along `axis` (falls back to a
    // one-sided difference near the domain boundary). Returns false if neither
    // side nor the one-sided fallback can be evaluated.
    bool fd_partial(const VectorX<T>& p, int axis, T fd_frac, VectorX<T>& deriv) const
    {
        T h = fd_frac * spacing[axis];
        if (!(h > T(0))) return false;

        VectorX<T> pp = p, pm = p;
        pp[axis] += h;
        pm[axis] -= h;

        VectorX<T> fp, fm, f0;
        bool okp = fff_direction(pp, fp);
        bool okm = fff_direction(pm, fm);
        if (okp && okm) { deriv = (fp - fm) / (T(2) * h); return true; }
        if (!fff_direction(p, f0)) return false;
        if (okp) { deriv = (fp - f0) / h; return true; }
        if (okm) { deriv = (f0 - fm) / h; return true; }
        return false;
    }
};

// ---------------------------------------------------------------------------
// .vff binary format I/O (self-describing).
// Layout: "VFF1", uint32 dtype(0=f64,1=f32), uint32 D, uint32 C,
//         uint32 n[D], float64 dmin[D], float64 dmax[D], then data.
// ---------------------------------------------------------------------------
namespace vff_io
{
    template<typename T>
    bool load(const std::string& filename, VectorFieldGrid<T>& g)
    {
        std::ifstream in(filename, std::ios::binary);
        if (!in)
        {
            std::cerr << "vff_io::load: cannot open " << filename << std::endl;
            return false;
        }

        char magic[4];
        in.read(magic, 4);
        if (std::strncmp(magic, "VFF1", 4) != 0)
        {
            std::cerr << "vff_io::load: bad magic in " << filename << std::endl;
            return false;
        }

        std::uint32_t dtype = 0, D = 0, C = 0;
        in.read(reinterpret_cast<char*>(&dtype), sizeof(dtype));
        in.read(reinterpret_cast<char*>(&D), sizeof(D));
        in.read(reinterpret_cast<char*>(&C), sizeof(C));

        g.D = static_cast<int>(D);
        g.C = static_cast<int>(C);
        g.n.resize(g.D);
        std::vector<std::uint32_t> nn(g.D);
        in.read(reinterpret_cast<char*>(nn.data()), sizeof(std::uint32_t) * g.D);
        for (int a = 0; a < g.D; ++a) g.n[a] = static_cast<int>(nn[a]);

        std::vector<double> mn(g.D), mx(g.D);
        in.read(reinterpret_cast<char*>(mn.data()), sizeof(double) * g.D);
        in.read(reinterpret_cast<char*>(mx.data()), sizeof(double) * g.D);
        g.dmin.resize(g.D);
        g.dmax.resize(g.D);
        for (int a = 0; a < g.D; ++a) { g.dmin[a] = static_cast<T>(mn[a]); g.dmax[a] = static_cast<T>(mx[a]); }
        g.compute_spacing();

        std::size_t count = static_cast<std::size_t>(g.C);
        for (int a = 0; a < g.D; ++a) count *= static_cast<std::size_t>(g.n[a]);

        g.data.resize(count);
        if (dtype == 1) // float32
        {
            std::vector<float> tmp(count);
            in.read(reinterpret_cast<char*>(tmp.data()), sizeof(float) * count);
            for (std::size_t i = 0; i < count; ++i) g.data[i] = static_cast<T>(tmp[i]);
        }
        else // float64
        {
            std::vector<double> tmp(count);
            in.read(reinterpret_cast<char*>(tmp.data()), sizeof(double) * count);
            for (std::size_t i = 0; i < count; ++i) g.data[i] = static_cast<T>(tmp[i]);
        }

        if (!in)
        {
            std::cerr << "vff_io::load: unexpected EOF / read error in " << filename << std::endl;
            return false;
        }

        // Stable FFF for critical points needs the D-1 spatial-gradient
        // components (the first D-1 components, in axis order) and the field
        // value v itself. Accept either a spatial-gradient field (C == D-1) or a
        // full space-time gradient field (C == D, last component ignored).
        if (g.C == g.D)
        {
            std::cout << "vff_io::load: full gradient field (C=" << g.C
                      << "); using the first " << g.D - 1
                      << " spatial components for stable FFF." << std::endl;
        }
        else if (g.C < g.D - 1)
        {
            std::cerr << "vff_io::load: error, C(" << g.C << ") < D-1(" << g.D - 1
                      << ") -- not enough gradient components for FFF." << std::endl;
            return false;
        }
        return true;
    }

    template<typename T>
    bool save(const std::string& filename, const VectorFieldGrid<T>& g, bool as_float32 = false)
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            std::cerr << "vff_io::save: cannot open " << filename << std::endl;
            return false;
        }
        out.write("VFF1", 4);
        std::uint32_t dtype = as_float32 ? 1u : 0u;
        std::uint32_t D = static_cast<std::uint32_t>(g.D);
        std::uint32_t C = static_cast<std::uint32_t>(g.C);
        out.write(reinterpret_cast<const char*>(&dtype), sizeof(dtype));
        out.write(reinterpret_cast<const char*>(&D), sizeof(D));
        out.write(reinterpret_cast<const char*>(&C), sizeof(C));
        std::vector<std::uint32_t> nn(g.D);
        for (int a = 0; a < g.D; ++a) nn[a] = static_cast<std::uint32_t>(g.n[a]);
        out.write(reinterpret_cast<const char*>(nn.data()), sizeof(std::uint32_t) * g.D);
        std::vector<double> mn(g.D), mx(g.D);
        for (int a = 0; a < g.D; ++a) { mn[a] = static_cast<double>(g.dmin[a]); mx[a] = static_cast<double>(g.dmax[a]); }
        out.write(reinterpret_cast<const char*>(mn.data()), sizeof(double) * g.D);
        out.write(reinterpret_cast<const char*>(mx.data()), sizeof(double) * g.D);
        if (as_float32)
        {
            std::vector<float> tmp(g.data.size());
            for (std::size_t i = 0; i < g.data.size(); ++i) tmp[i] = static_cast<float>(g.data[i]);
            out.write(reinterpret_cast<const char*>(tmp.data()), sizeof(float) * tmp.size());
        }
        else
        {
            std::vector<double> tmp(g.data.size());
            for (std::size_t i = 0; i < g.data.size(); ++i) tmp[i] = static_cast<double>(g.data[i]);
            out.write(reinterpret_cast<const char*>(tmp.data()), sizeof(double) * tmp.size());
        }
        return static_cast<bool>(out);
    }
} // namespace vff_io

// ---------------------------------------------------------------------------
// RK4 stream-line integration of the STABLE feature flow field h = (+/-)f + tau g.
// Each branch (forward/backward) integrates with a positive arc-length step; the
// backward branch reverses the tangent via base = -f while keeping g attracting.
// ---------------------------------------------------------------------------
template<typename T>
struct Integrator
{
    const VectorFieldGrid<T>& g;
    T   step;          // fixed arc-length step (full space-time)
    int max_steps;     // maximum steps per direction
    T   k;             // convergence strength (paper's k > 0)
    T   tau_max;       // upper clamp for the adaptive tau
    T   fd_frac;       // finite-difference step as a fraction of grid spacing

    Integrator(const VectorFieldGrid<T>& grid, T step_size, int max_steps_,
               T k_strength, T tau_max_, T fd_frac_)
        : g(grid), step(step_size), max_steps(max_steps_),
          k(k_strength), tau_max(tau_max_), fd_frac(fd_frac_) {}

    // One RK4 step in the stable field for the given branch. Returns false to stop.
    bool rk4(const VectorX<T>& p, bool forward, VectorX<T>& p_next) const
    {
        VectorX<T> k1, k2, k3, k4, tmp;
        if (!g.stable_fff_unit(p, forward, k, tau_max, fd_frac, k1)) return false;
        tmp = p + (T(0.5) * step) * k1;
        if (!g.stable_fff_unit(tmp, forward, k, tau_max, fd_frac, k2)) return false;
        tmp = p + (T(0.5) * step) * k2;
        if (!g.stable_fff_unit(tmp, forward, k, tau_max, fd_frac, k3)) return false;
        tmp = p + step * k3;
        if (!g.stable_fff_unit(tmp, forward, k, tau_max, fd_frac, k4)) return false;

        VectorX<T> kk = (k1 + T(2) * k2 + T(2) * k3 + k4) / T(6);
        p_next = p + step * kk;
        if (!g.in_domain(p_next)) return false;
        return true;
    }

    // Integrate one branch from seed; appends new points (excludes seed) to out.
    void integrate(const VectorX<T>& seed, bool forward, std::vector<VectorX<T>>& out) const
    {
        VectorX<T> p = seed, pn;
        const T stall = static_cast<T>(1e-6) * step;
        for (int s = 0; s < max_steps; ++s)
        {
            if (!rk4(p, forward, pn)) break;
            if ((pn - p).norm() < stall) break;
            out.emplace_back(pn);
            p = pn;
        }
    }

    // Full stream line through seed = reverse(backward) + seed + forward.
    void trace(const VectorX<T>& seed, std::vector<VectorX<T>>& line) const
    {
        std::vector<VectorX<T>> back, fwd;
        integrate(seed, false, back);
        integrate(seed, true,  fwd);
        line.clear();
        line.reserve(back.size() + 1 + fwd.size());
        for (auto it = back.rbegin(); it != back.rend(); ++it) line.emplace_back(*it);
        line.emplace_back(seed);
        for (auto& q : fwd) line.emplace_back(q);
    }
};

// ---------------------------------------------------------------------------
// Spatial-temporal "covered" index (identical to the plain-FFF version).
// ---------------------------------------------------------------------------
template<typename T>
class CoveredIndex
{
public:
    CoveredIndex(int D, T spatial_eps, T temporal_eps)
        : D_(D), se_(spatial_eps), te_(temporal_eps) {}

    void insert(const VectorX<T>& p)
    {
        std::array<int, 4> cell;
        quantize(p, cell);
        occupied_.insert(key(cell));
    }

    void insert_line(const std::vector<VectorX<T>>& line)
    {
        for (const auto& p : line) insert(p);
    }

    bool covered(const VectorX<T>& p) const
    {
        std::array<int, 4> cell;
        quantize(p, cell);
        std::array<int, 4> nb = cell;
        return scan(cell, nb, 0);
    }

private:
    int D_;
    T   se_, te_;
    std::unordered_set<std::uint64_t> occupied_;

    inline void quantize(const VectorX<T>& p, std::array<int, 4>& cell) const
    {
        for (int a = 0; a < D_; ++a)
        {
            T eps = (a == D_ - 1) ? te_ : se_;
            cell[a] = static_cast<int>(std::lround(static_cast<double>(p[a] / eps)));
        }
    }

    inline std::uint64_t key(const std::array<int, 4>& cell) const
    {
        std::uint64_t h = 1469598103934665603ull; // FNV-style offset
        for (int a = 0; a < D_; ++a)
        {
            std::uint64_t v = static_cast<std::uint64_t>(static_cast<std::int64_t>(cell[a]));
            h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        }
        return h;
    }

    // recursive 3^D neighbourhood scan
    bool scan(const std::array<int, 4>& cell, std::array<int, 4>& nb, int axis) const
    {
        if (axis == D_)
            return occupied_.find(key(nb)) != occupied_.end();
        for (int d = -1; d <= 1; ++d)
        {
            nb[axis] = cell[axis] + d;
            if (scan(cell, nb, axis + 1)) return true;
        }
        return false;
    }
};

// ---------------------------------------------------------------------------
// Seed I/O: read critical points (one per row) from a CSV with a header line.
// Each row has D numeric columns (x, y, [z,] t) matching the grid domain axes.
// ---------------------------------------------------------------------------
template<typename T>
bool read_seeds_csv(const std::string& filename, int D, std::vector<VectorX<T>>& seeds)
{
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "read_seeds_csv: cannot open " << filename << std::endl;
        return false;
    }

    seeds.clear();
    std::string line;
    bool first = true;
    while (std::getline(in, line))
    {
        if (line.empty()) continue;

        std::vector<T> vals;
        std::stringstream ss(line);
        std::string tok;
        bool numeric = true;
        while (std::getline(ss, tok, ','))
        {
            try { vals.push_back(static_cast<T>(std::stod(tok))); }
            catch (...) { numeric = false; break; }
        }

        if (!numeric)
        {
            // header (or junk) line: only allowed as the very first line
            if (first) { first = false; continue; }
            continue;
        }
        first = false;

        if (static_cast<int>(vals.size()) < D) continue;
        VectorX<T> p(D);
        for (int a = 0; a < D; ++a) p[a] = vals[a];
        seeds.emplace_back(std::move(p));
    }
    return true;
}

// ---------------------------------------------------------------------------
// Track all critical points (same seeding algorithm as plain FFF):
//   sort seeds by time; for each seed, if already covered by a previously built
//   stream line, skip it; otherwise integrate its (stable) stream line in both
//   directions and register all its points as covered.
// ---------------------------------------------------------------------------
template<typename T>
void track_all(const VectorFieldGrid<T>& g,
               std::vector<VectorX<T>>&  seeds,
               T                         step,
               int                       max_steps,
               T                         spatial_eps,
               T                         temporal_eps,
               T                         k_strength,
               T                         tau_max,
               T                         fd_frac,
               std::vector<CP_Trace<T>>& traces)
{
    std::sort(seeds.begin(), seeds.end(),
              [](const VectorX<T>& a, const VectorX<T>& b) { return a[a.size() - 1] < b[b.size() - 1]; });

    Integrator<T>   integrator(g, step, max_steps, k_strength, tau_max, fd_frac);
    CoveredIndex<T> covered(g.D, spatial_eps, temporal_eps);

    traces.clear();
    traces.reserve(seeds.size());

    int skipped = 0;
    for (const auto& seed : seeds)
    {
        if (covered.covered(seed)) { ++skipped; continue; }

        std::vector<VectorX<T>> line;
        integrator.trace(seed, line);

        covered.insert_line(line);

        CP_Trace<T> tr;
        tr.traces = std::move(line);
        traces.emplace_back(std::move(tr));
    }

    std::cout << "seeds: " << seeds.size()
              << ", stream lines built: " << traces.size()
              << ", seeds skipped (already covered): " << skipped << std::endl;
}

} // namespace sfff
