//
// Spatial acceleration  d^2x/dt^2  of a critical point over a .ply mesh.
//
// For a scalar field f(x, y, t) the critical-point path x(t) satisfies
// grad_spatial f(x(t), t) = 0. Differentiating once gives the feature-flow
// velocity (the "x dot" of the figure)
//
//     dx/dt = - [H_s f]^{-1} (d/dt) grad_s f
//
// and differentiating a second time gives the spatial acceleration
//
//     d^2x_i/dt^2 = - sum_j (H_s^{-1})_{ij} (
//                       d^2/dt^2 (d f/dx_j)
//                     + 2 sum_k  (d^3 f / dx_j dx_k dt)  xdot_k
//                     +   sum_kl (d^3 f / dx_j dx_k dx_l) xdot_k xdot_l )
//
// where H_s is the 2x2 spatial Hessian [[f_xx, f_xy], [f_xy, f_yy]] (the 3rd
// domain axis is time) and xdot = dx/dt is the velocity above.
//
// Reads either an MFA model (a .mfa diy file) OR an INR model (a TorchScript
// .pt file + matching function name), plus a .ply mesh; it evaluates the
// velocity and acceleration at every mesh vertex and reports the min / max /
// average magnitudes.
//
// Inputs:
//   -f  .mfa  diy MFA input file               (MFA mode, the default)
//   -i  .pt   TorchScript INR model file       (INR mode; selected when set)
//   -n        INR function name (domain bounds), e.g. vortex_street_3d
//   -p  .ply  mesh whose vertices are sampled
//   -o        optional output CSV: x,y,t,xdot_x,xdot_y,acc_x,acc_y,|acc|
//
// MFA derivatives use the analytic mfa_extend::recover_mfa; INR derivatives use
// the exact nested-autograd INRModel::query_accel_derivs. Field is assumed to
// be 2 spatial dims + 1 time dim (domain dim 3).
//
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iomanip>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include <tbb/tbb.h>

#include "opts.h"
#include "block.hpp"
#include "mfa_extend.h"
#include "query_function.h"
#include "INRModel.h"

using namespace std;

// ---------------------------------------------------------------------------
// Minimal .ply reader (ascii + binary little/big endian).
//
// Only the vertex element is required; its first three float/double properties
// (conventionally x,y,z) are stored as a 3-vector. Any extra vertex properties
// and any other elements (edge, face, ...) are parsed and skipped so the reader
// stays robust to the various .ply files used in this project. The number of
// edges/faces, if present, is reported for convenience.
// ---------------------------------------------------------------------------
namespace ply_io
{
    enum class Format { ascii, binary_little, binary_big };

    enum class Type { i8, u8, i16, u16, i32, u32, f32, f64, unknown };

    static Type parse_type(const std::string& s)
    {
        if (s == "char"  || s == "int8")    return Type::i8;
        if (s == "uchar" || s == "uint8")   return Type::u8;
        if (s == "short" || s == "int16")   return Type::i16;
        if (s == "ushort"|| s == "uint16")  return Type::u16;
        if (s == "int"   || s == "int32")   return Type::i32;
        if (s == "uint"  || s == "uint32")  return Type::u32;
        if (s == "float" || s == "float32") return Type::f32;
        if (s == "double"|| s == "float64") return Type::f64;
        return Type::unknown;
    }

    static int type_size(Type t)
    {
        switch (t)
        {
            case Type::i8:  case Type::u8:  return 1;
            case Type::i16: case Type::u16: return 2;
            case Type::i32: case Type::u32: case Type::f32: return 4;
            case Type::f64: return 8;
            default: return 0;
        }
    }

    struct Property
    {
        std::string name;
        bool        is_list = false;
        Type        count_type = Type::unknown;   // list only
        Type        value_type = Type::unknown;
    };

    struct Element
    {
        std::string           name;
        size_t                count = 0;
        std::vector<Property> props;
    };

    static void byteswap(char* p, int n)
    {
        for (int i = 0; i < n / 2; ++i) std::swap(p[i], p[n - 1 - i]);
    }

    static double read_bin_scalar(std::istream& in, Type t, bool swap)
    {
        char buf[8];
        int  n = type_size(t);
        in.read(buf, n);
        if (swap) byteswap(buf, n);
        switch (t)
        {
            case Type::i8:  return static_cast<double>(*reinterpret_cast<int8_t*>(buf));
            case Type::u8:  return static_cast<double>(*reinterpret_cast<uint8_t*>(buf));
            case Type::i16: return static_cast<double>(*reinterpret_cast<int16_t*>(buf));
            case Type::u16: return static_cast<double>(*reinterpret_cast<uint16_t*>(buf));
            case Type::i32: return static_cast<double>(*reinterpret_cast<int32_t*>(buf));
            case Type::u32: return static_cast<double>(*reinterpret_cast<uint32_t*>(buf));
            case Type::f32: return static_cast<double>(*reinterpret_cast<float*>(buf));
            case Type::f64: return static_cast<double>(*reinterpret_cast<double*>(buf));
            default: return 0.0;
        }
    }

    // Read vertices (first 3 scalar props) from a .ply file.
    // Returns false on failure. n_edges / n_faces report counts of those elements.
    template<typename T>
    bool read(const std::string& filename, std::vector<VectorX<T>>& vertices,
              size_t& n_edges, size_t& n_faces)
    {
        n_edges = 0;
        n_faces = 0;
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (!in)
        {
            std::cerr << "ply_io::read: cannot open " << filename << std::endl;
            return false;
        }

        std::string line;
        std::getline(in, line);
        if (line.rfind("ply", 0) != 0)
        {
            std::cerr << "ply_io::read: " << filename << " is not a PLY file" << std::endl;
            return false;
        }

        Format format = Format::ascii;
        std::vector<Element> elements;

        // ---- header ----
        while (std::getline(in, line))
        {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            std::istringstream ss(line);
            std::string tok;
            ss >> tok;
            if (tok == "comment" || tok == "obj_info")
                continue;
            else if (tok == "format")
            {
                std::string f;
                ss >> f;
                if (f == "ascii")                       format = Format::ascii;
                else if (f == "binary_little_endian")   format = Format::binary_little;
                else if (f == "binary_big_endian")      format = Format::binary_big;
            }
            else if (tok == "element")
            {
                Element e;
                ss >> e.name >> e.count;
                elements.push_back(e);
            }
            else if (tok == "property")
            {
                if (elements.empty()) continue;
                Property p;
                std::string t1;
                ss >> t1;
                if (t1 == "list")
                {
                    std::string tcount, tvalue;
                    ss >> tcount >> tvalue >> p.name;
                    p.is_list    = true;
                    p.count_type = parse_type(tcount);
                    p.value_type = parse_type(tvalue);
                }
                else
                {
                    ss >> p.name;
                    p.value_type = parse_type(t1);
                }
                elements.back().props.push_back(p);
            }
            else if (tok == "end_header")
                break;
        }

        const bool swap = (format == Format::binary_big);   // host assumed little-endian

        // ---- body ----
        for (const Element& e : elements)
        {
            const bool is_vertex = (e.name == "vertex");
            if (e.name == "edge") n_edges = e.count;
            if (e.name == "face") n_faces = e.count;

            for (size_t i = 0; i < e.count; ++i)
            {
                if (format == Format::ascii)
                {
                    if (!std::getline(in, line)) { std::cerr << "ply_io::read: unexpected EOF\n"; return false; }
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    std::istringstream ss(line);
                    std::vector<double> scalars;
                    for (const Property& p : e.props)
                    {
                        if (p.is_list)
                        {
                            long cnt = 0; ss >> cnt;
                            for (long k = 0; k < cnt; ++k) { double v; ss >> v; }
                        }
                        else
                        {
                            double v; ss >> v;
                            scalars.push_back(v);
                        }
                    }
                    if (is_vertex)
                    {
                        VectorX<T> vtx(3);
                        for (int d = 0; d < 3; ++d)
                            vtx[d] = (d < (int)scalars.size()) ? static_cast<T>(scalars[d]) : T(0);
                        vertices.push_back(vtx);
                    }
                }
                else // binary
                {
                    std::vector<double> scalars;
                    for (const Property& p : e.props)
                    {
                        if (p.is_list)
                        {
                            double cnt = read_bin_scalar(in, p.count_type, swap);
                            for (long k = 0; k < (long)cnt; ++k)
                                read_bin_scalar(in, p.value_type, swap);
                        }
                        else
                        {
                            scalars.push_back(read_bin_scalar(in, p.value_type, swap));
                        }
                    }
                    if (is_vertex)
                    {
                        VectorX<T> vtx(3);
                        for (int d = 0; d < 3; ++d)
                            vtx[d] = (d < (int)scalars.size()) ? static_cast<T>(scalars[d]) : T(0);
                        vertices.push_back(vtx);
                    }
                    if (!in) { std::cerr << "ply_io::read: unexpected EOF in binary body\n"; return false; }
                }
            }
        }

        return true;
    }
}

// ---------------------------------------------------------------------------
// Acceleration computation.
//
// All derivatives are stored in physical-domain coordinates with domain axis
// order (x, y, t). The third-derivative vector `third` is ordered:
//   [fxxx, fxxy, fxyy, fyyy,   fxxt, fxyt, fyyt,   fxtt, fytt]
// ---------------------------------------------------------------------------

// Reciprocal condition number of a symmetric 2x2 matrix, rcond = |lambda_min| /
// |lambda_max| in [0, 1]. rcond == 0 means exactly singular; small rcond means
// near-degenerate (one eigenvalue ~ 0). Closed form for [[a,b],[b,c]].
template<typename T>
T sym2x2_rcond(const Eigen::Matrix<T, 2, 2>& H)
{
    const T a = H(0, 0), c = H(1, 1), b = T(0.5) * (H(0, 1) + H(1, 0));
    const T mean = T(0.5) * (a + c);
    const T diff = T(0.5) * (a - c);
    const T disc = std::sqrt(diff * diff + b * b);
    const T l1 = mean + disc;
    const T l2 = mean - disc;
    const T amax = std::max(std::abs(l1), std::abs(l2));
    const T amin = std::min(std::abs(l1), std::abs(l2));
    return (amax > T(0)) ? amin / amax : T(0);
}

// Returns true (and fills xdot, acc) when the spatial Hessian is well enough
// conditioned (reciprocal condition number > rcond_tol). The acceleration
// involves H^{-1} up to the cubic power, so ||x_tt|| ~ rcond^{-3}: a relative
// determinant / rcond test with a too-small tolerance lets near-degenerate
// (ill-conditioned) points through and produces enormous finite values. The
// threshold is therefore expressed directly as a reciprocal-condition-number
// cutoff (a critical point is "degenerate" exactly when an eigenvalue of the
// spatial Hessian vanishes, i.e. rcond -> 0).
template<typename T>
bool accel_from_derivs(const Eigen::Matrix<T, 2, 2>& H,
                       const Eigen::Matrix<T, 2, 1>& gt,
                       const Eigen::Matrix<T, 9, 1>& third,
                       T rcond_tol,
                       Eigen::Matrix<T, 2, 1>& xdot,
                       Eigen::Matrix<T, 2, 1>& acc)
{
    // Exclude (near-)degenerate critical points: x_tt = -H^{-1}(...) blows up as
    // the spatial Hessian becomes singular.
    if (sym2x2_rcond<T>(H) <= rcond_tol)
        return false;                                // degenerate spatial Hessian -> excluded

    const Eigen::Matrix<T, 2, 2> Hinv = H.inverse();

    // velocity:  xdot = -H^{-1} (d/dt grad_s f)
    xdot = -Hinv * gt;
    const T ux = xdot(0), uy = xdot(1);

    const T fxxx = third(0), fxxy = third(1), fxyy = third(2), fyyy = third(3);
    const T fxxt = third(4), fxyt = third(5), fyyt = third(6);
    const T fxtt = third(7), fytt = third(8);

    // R_j = d^2(df/dx_j)/dt^2
    //       + 2 sum_k d^3f/dx_j dx_k dt  xdot_k
    //       +   sum_kl d^3f/dx_j dx_k dx_l xdot_k xdot_l
    Eigen::Matrix<T, 2, 1> R;
    R(0) = fxtt
         + T(2) * (fxxt * ux + fxyt * uy)
         + (fxxx * ux * ux + T(2) * fxxy * ux * uy + fxyy * uy * uy);
    R(1) = fytt
         + T(2) * (fxyt * ux + fyyt * uy)
         + (fxxy * ux * ux + T(2) * fxyy * ux * uy + fyyy * uy * uy);

    // acceleration: x_tt = -H^{-1} R
    acc = -Hinv * R;
    return true;
}

// MFA: fill the spatial Hessian, time-derivative of the spatial gradient, and
// the spatial/time third derivatives from analytic mfa_extend::recover_mfa.
template<typename T>
bool mfa_accel_derivs(const Block<T>* b, const VectorX<T>& vertex,
                      Eigen::Matrix<T, 2, 2>& H,
                      Eigen::Matrix<T, 2, 1>& gt,
                      Eigen::Matrix<T, 9, 1>& third)
{
    const int D = static_cast<int>(b->core_mins.size());
    if (D < 3)
        return false;

    // clamp the query point into the block's domain
    VectorX<T> p(D);
    for (int d = 0; d < D; ++d)
    {
        T c = (d < vertex.size()) ? vertex[d] : T(0);
        if (c < b->core_mins[d]) c = b->core_mins[d];
        if (c > b->core_maxs[d]) c = b->core_maxs[d];
        p[d] = c;
    }

    VectorX<T> out(1);
    VectorXi   deriv = VectorXi::Zero(D);
    auto d3 = [&](int dx, int dy, int dt) -> T
    {
        deriv.setZero();
        deriv[0] = dx; deriv[1] = dy; deriv[2] = dt;   // axis order (x, y, t)
        mfa_extend::recover_mfa(b, p, out, deriv);
        return out[0];
    };

    const T fxx = d3(2, 0, 0), fxy = d3(1, 1, 0), fyy = d3(0, 2, 0);
    const T fxt = d3(1, 0, 1), fyt = d3(0, 1, 1);

    H(0, 0) = fxx; H(0, 1) = fxy; H(1, 0) = fxy; H(1, 1) = fyy;
    gt(0) = fxt;   gt(1) = fyt;

    third(0) = d3(3, 0, 0);   // fxxx
    third(1) = d3(2, 1, 0);   // fxxy
    third(2) = d3(1, 2, 0);   // fxyy
    third(3) = d3(0, 3, 0);   // fyyy
    third(4) = d3(2, 0, 1);   // fxxt
    third(5) = d3(1, 1, 1);   // fxyt
    third(6) = d3(0, 2, 1);   // fyyt
    third(7) = d3(1, 0, 2);   // fxtt
    third(8) = d3(0, 1, 2);   // fytt
    return true;
}

int main(int argc, char** argv)
{
    using T = double;

    diy::mpi::environment  env(argc, argv);
    diy::mpi::communicator world;

    string mfa_file   = "approx.mfa";   // MFA diy input file
    string inr_model  = "";             // TorchScript .pt; non-empty => INR mode
    string inr_name   = "";             // INR analytical function name (domain)
    string ply_file   = "mesh.ply";     // mesh whose vertices are sampled
    string out_file   = "";             // optional per-vertex CSV dump
    double rcond_tol  = 1e-3;           // degeneracy cutoff (reciprocal cond. number)
    bool   help       = false;

    opts::Options ops;
    ops >> opts::Option('f', "mfa",  mfa_file,  " input MFA .mfa diy file");
    ops >> opts::Option('i', "inr",  inr_model, " input INR TorchScript .pt model (selects INR mode)");
    ops >> opts::Option('n', "name", inr_name,  " INR function name (domain bounds), e.g. vortex_street_3d");
    ops >> opts::Option('p', "ply",  ply_file,  " input .ply mesh (vertices sampled)");
    ops >> opts::Option('o', "out",  out_file,  " optional output CSV: x,y,t,xdot_x,xdot_y,acc_x,acc_y,|acc|");
    ops >> opts::Option('t', "rcond", rcond_tol, " degeneracy cutoff: skip points whose spatial Hessian reciprocal condition number <= this (default 1e-3)");
    ops >> opts::Option('h', "help", help,      " show help");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    const bool use_inr = !inr_model.empty();

    // ---- read the .ply vertices ----
    std::vector<VectorX<T>> vertices;
    size_t n_edges = 0, n_faces = 0;
    if (!ply_io::read<T>(ply_file, vertices, n_edges, n_faces))
        return 1;
    std::cout << "read " << vertices.size() << " vertices, " << n_edges
              << " edges, " << n_faces << " faces from " << ply_file << std::endl;
    if (vertices.empty())
    {
        std::cerr << "no vertices -- nothing to evaluate" << std::endl;
        return 1;
    }

    const size_t N = vertices.size();
    std::vector<Eigen::Matrix<T, 2, 1>> xdot(N, Eigen::Matrix<T, 2, 1>::Zero());
    std::vector<Eigen::Matrix<T, 2, 1>> acc(N, Eigen::Matrix<T, 2, 1>::Zero());
    std::vector<char> valid(N, 0);

    auto t0 = std::chrono::high_resolution_clock::now();

    if (use_inr)
    {
        // -------- INR mode --------
        if (inr_name.empty())
        {
            std::cerr << "INR mode requires -n <function name> for the domain bounds" << std::endl;
            return 1;
        }
        INRModel<T> model(inr_name, inr_model);
        if (!model.isLoaded())
        {
            std::cerr << "failed to load INR model " << inr_model << std::endl;
            return 1;
        }
        const int D = static_cast<int>(model.domain_min.size());
        std::cout << "INR model '" << inr_name << "' domain dim " << D
                  << " min " << model.domain_min.transpose()
                  << " max " << model.domain_max.transpose() << std::endl;
        if (D < 3)
        {
            std::cerr << "spatial acceleration needs a (x,y,t) field (domain dim >= 3)" << std::endl;
            return 1;
        }

        Eigen::MatrixX<T> H;
        VectorX<T>        gt, third;
        for (size_t i = 0; i < N; ++i)
        {
            VectorX<T> p(D);
            for (int d = 0; d < D; ++d)
            {
                T c = (d < vertices[i].size()) ? vertices[i][d] : T(0);
                if (c < model.domain_min[d]) c = model.domain_min[d];
                if (c > model.domain_max[d]) c = model.domain_max[d];
                p[d] = c;
            }

            model.query_accel_derivs(p, H, gt, third);

            Eigen::Matrix<T, 2, 2> H2 = H.topLeftCorner(2, 2);
            Eigen::Matrix<T, 2, 1> g2 = gt.head(2);
            Eigen::Matrix<T, 9, 1> t9;
            for (int k = 0; k < 9; ++k) t9(k) = third(k);

            valid[i] = accel_from_derivs<T>(H2, g2, t9, T(rcond_tol), xdot[i], acc[i]) ? 1 : 0;
        }
    }
    else
    {
        // -------- MFA mode --------
        diy::FileStorage storage("./DIY.XXXXXX");
        diy::Master      master(world, -1, -1,
                                &Block<T>::create, &Block<T>::destroy,
                                &storage, &Block<T>::save, &Block<T>::load);
        diy::ContiguousAssigner assigner(world.size(), -1);

        diy::io::read_blocks(mfa_file.c_str(), world, assigner, master, &Block<T>::load);
        std::cout << master.size() << " blocks read from file " << mfa_file << std::endl;
        if (master.size() == 0)
        {
            std::cerr << "no blocks in " << mfa_file << std::endl;
            return 1;
        }

        master.foreach([&](Block<T>* b, const diy::Master::ProxyWithLink&)
        {
            std::cout << "MFA domain min " << b->core_mins.transpose()
                      << " max " << b->core_maxs.transpose() << std::endl;
            tbb::affinity_partitioner ap;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, N),
            [&](const tbb::blocked_range<size_t>& range)
            {
                Eigen::Matrix<T, 2, 2> H;
                Eigen::Matrix<T, 2, 1> gt;
                Eigen::Matrix<T, 9, 1> third;
                for (size_t i = range.begin(); i != range.end(); ++i)
                {
                    if (mfa_accel_derivs<T>(b, vertices[i], H, gt, third))
                        valid[i] = accel_from_derivs<T>(H, gt, third, T(rcond_tol), xdot[i], acc[i]) ? 1 : 0;
                    else
                        valid[i] = 0;
                }
            }, ap);
        });
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- report ----
    size_t n_valid = 0;
    T min_a = std::numeric_limits<T>::max();
    T max_a = std::numeric_limits<T>::lowest();
    T sum_a = T(0);
    T min_v = std::numeric_limits<T>::max();
    T max_v = std::numeric_limits<T>::lowest();
    T sum_v = T(0);
    for (size_t i = 0; i < N; ++i)
    {
        if (!valid[i]) continue;
        ++n_valid;
        const T a = acc[i].norm();
        const T v = xdot[i].norm();
        min_a = std::min(min_a, a); max_a = std::max(max_a, a); sum_a += a;
        min_v = std::min(min_v, v); max_v = std::max(max_v, v); sum_v += v;
    }

    std::cout << "spatial acceleration |d^2x/dt^2| over " << n_valid << " / " << N
              << " vertices (" << (N - n_valid)
              << " skipped: degenerate spatial Hessian, rcond <= " << rcond_tol << "):"
              << std::endl;
    if (n_valid > 0)
    {
        std::cout << std::setprecision(10);
        std::cout << "  velocity |dx/dt|     min = " << min_v
                  << "  max = " << max_v
                  << "  avg = " << (sum_v / static_cast<T>(n_valid)) << std::endl;
        std::cout << "  acceleration |x_tt|  min = " << min_a
                  << "  max = " << max_a
                  << "  avg = " << (sum_a / static_cast<T>(n_valid)) << std::endl;
    }
    std::cout << "evaluation time (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << std::endl;

    if (!out_file.empty())
    {
        std::ofstream out(out_file.c_str());
        if (out)
        {
            out << std::setprecision(12);
            out << "x,y,t,xdot_x,xdot_y,acc_x,acc_y,acc_norm,valid\n";
            for (size_t i = 0; i < N; ++i)
            {
                const T x = vertices[i].size() > 0 ? vertices[i][0] : T(0);
                const T y = vertices[i].size() > 1 ? vertices[i][1] : T(0);
                const T tt = vertices[i].size() > 2 ? vertices[i][2] : T(0);
                out << x << ',' << y << ',' << tt << ','
                    << xdot[i](0) << ',' << xdot[i](1) << ','
                    << acc[i](0) << ',' << acc[i](1) << ','
                    << (valid[i] ? acc[i].norm() : std::numeric_limits<T>::quiet_NaN()) << ','
                    << static_cast<int>(valid[i]) << '\n';
            }
            std::cout << "wrote per-vertex acceleration to " << out_file << std::endl;
        }
        else
            std::cerr << "failed to open " << out_file << " for writing" << std::endl;
    }

    return 0;
}
