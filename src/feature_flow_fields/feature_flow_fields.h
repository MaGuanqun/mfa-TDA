#pragma once
//
// Feature Flow Fields (Theisel & Seidel, VisSym 2003) for critical-point tracking.
//
// This is a *discrete* reimplementation, following the original paper:
//   - The input is a sampled vector field v on a regular space-time grid
//     (in this project v = grad_space(s), the spatial gradient of the scalar
//     field, stored in a .vff file -- see vff_io below).
//   - v is reconstructed at arbitrary points by MULTILINEAR interpolation
//     (C0 field). The derivatives of v needed for the feature flow field are
//     the ANALYTIC derivatives of that multilinear interpolant inside each
//     cell (not finite differences across nodes). Plugging these into eq. (3)
//     yields a (bi)quadratic, C^-1 feature flow field f, whose stream lines are
//     only G0 continuous (the "sharp corners" mentioned in the paper).
//   - Critical points are tracked by integrating stream lines of f (RK4),
//     forward and backward, from seed critical points extracted at the discrete
//     time steps. No Newton correction step is used (the paper has none).
//
// Feature flow field (generalized cross product / null vector of [grad v_1; ...]):
//   Stack the gradients of the components of v into a C x D matrix
//       A = [ grad v_1 ; ... ; grad v_C ]^T   (C = D-1 components, D domain axes,
//                                              last axis = time)
//   then  f[k] = (-1)^k * det( A with column k removed ),  k = 0..D-1.
//   For the 3D case (2D space + time) this reproduces eq. (3) exactly.
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

#include <Eigen/Dense>

// NOTE: do NOT `using namespace Eigen;` here. mfa (types.hpp) defines its own
// global-namespace VectorX/MatrixX alias templates, and Eigen 3.4 also defines
// Eigen::VectorX/MatrixX. Bringing Eigen's into the global namespace makes every
// unqualified VectorX/MatrixX ambiguous (even inside mfa's own headers). We rely
// on mfa's global typedefs instead, so mfa headers must be included first.

#include "CP_Trace.h"

namespace fff
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

    // Feature flow field direction at p (un-normalized). Returns false if outside.
    bool fff_direction(const VectorX<T>& p, VectorX<T>& f) const
    {
        MatrixX<T> A;
        if (!sample(p, &A, nullptr)) return false;

        f.resize(D);
        MatrixX<T> M(C, D - 1);
        for (int k = 0; k < D; ++k)
        {
            int col = 0;
            for (int j = 0; j < D; ++j)
                if (j != k) M.col(col++) = A.col(j);
            T det = M.determinant();
            f[k] = ((k & 1) ? T(-1) : T(1)) * det;
        }
        return true;
    }

    // Unit-normalized feature flow direction (full space-time norm).
    // Returns false if outside the domain or f is (near) zero (stagnation).
    bool fff_unit(const VectorX<T>& p, VectorX<T>& dir, T eps = static_cast<T>(1e-30)) const
    {
        if (!fff_direction(p, dir)) return false;
        T nrm = dir.norm();
        if (!(nrm > eps)) return false;
        dir /= nrm;
        return true;
    }
};

// ---------------------------------------------------------------------------
// .vff binary format I/O (self-describing). See feature_flow_fields.h header.
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

        if (g.C != g.D - 1)
            std::cerr << "vff_io::load: warning, C(" << g.C << ") != D-1(" << g.D - 1
                      << ") -- expected a spatial-gradient field." << std::endl;
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
// RK4 stream-line integration of the feature flow field (no correction step).
// ---------------------------------------------------------------------------
template<typename T>
struct Integrator
{
    const VectorFieldGrid<T>& g;
    T   step;          // fixed arc-length step (full space-time)
    int max_steps;     // maximum steps per direction

    Integrator(const VectorFieldGrid<T>& grid, T step_size, int max_steps_)
        : g(grid), step(step_size), max_steps(max_steps_) {}

    // One RK4 step. sign = +1 forward, -1 backward. Returns false to stop.
    bool rk4(const VectorX<T>& p, T sign, VectorX<T>& p_next) const
    {
        VectorX<T> k1, k2, k3, k4, tmp;
        if (!g.fff_unit(p, k1)) return false;
        tmp = p + sign * (T(0.5) * step) * k1;
        if (!g.fff_unit(tmp, k2)) return false;
        tmp = p + sign * (T(0.5) * step) * k2;
        if (!g.fff_unit(tmp, k3)) return false;
        tmp = p + sign * step * k3;
        if (!g.fff_unit(tmp, k4)) return false;

        VectorX<T> k = (k1 + T(2) * k2 + T(2) * k3 + k4) / T(6);
        p_next = p + sign * step * k;
        if (!g.in_domain(p_next)) return false;
        return true;
    }

    // Integrate one direction from seed; appends new points (excludes seed) to out.
    void integrate(const VectorX<T>& seed, T sign, std::vector<VectorX<T>>& out) const
    {
        VectorX<T> p = seed, pn;
        const T stall = static_cast<T>(1e-6) * step;
        for (int s = 0; s < max_steps; ++s)
        {
            if (!rk4(p, sign, pn)) break;
            if ((pn - p).norm() < stall) break;
            out.emplace_back(pn);
            p = pn;
        }
    }

    // Full stream line through seed = reverse(backward) + seed + forward.
    void trace(const VectorX<T>& seed, std::vector<VectorX<T>>& line) const
    {
        std::vector<VectorX<T>> back, fwd;
        integrate(seed, T(-1), back);
        integrate(seed, T(+1), fwd);
        line.clear();
        line.reserve(back.size() + 1 + fwd.size());
        for (auto it = back.rbegin(); it != back.rend(); ++it) line.emplace_back(*it);
        line.emplace_back(seed);
        for (auto& q : fwd) line.emplace_back(q);
    }
};

// ---------------------------------------------------------------------------
// Spatial-temporal "covered" index: quantize each point into a cell and test
// whether a candidate seed already lies on a previously built stream line.
// Spatial axes use spatial_eps, the time axis uses temporal_eps. Queries scan
// the 3^D neighbourhood so points within ~eps of an occupied cell count as
// covered.
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
// Track all critical points (the paper's seeding algorithm):
//   sort seeds by time; for each seed, if it is already covered by a previously
//   built stream line, skip it; otherwise integrate its stream line (both
//   directions) and register all its points as covered.
// ---------------------------------------------------------------------------
template<typename T>
void track_all(const VectorFieldGrid<T>& g,
               std::vector<VectorX<T>>&  seeds,
               T                         step,
               int                       max_steps,
               T                         spatial_eps,
               T                         temporal_eps,
               std::vector<CP_Trace<T>>& traces)
{
    std::sort(seeds.begin(), seeds.end(),
              [](const VectorX<T>& a, const VectorX<T>& b) { return a[a.size() - 1] < b[b.size() - 1]; });

    Integrator<T>   integrator(g, step, max_steps);
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

} // namespace fff
