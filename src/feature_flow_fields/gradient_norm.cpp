//
// Gradient-norm range over a .ply mesh.
//
// Reads either an MFA model (a .mfa diy file) OR an INR model (a TorchScript
// .pt file + the matching analytical-domain "function name"), plus a .ply file
// of vertices (and, optionally, edges/faces). It evaluates the gradient of the
// scalar field at every .ply vertex and reports the minimum, maximum and
// average gradient norm.
//
// The gradient is the full domain gradient (all domain axes, e.g. [fx, fy, ft]
// for a 2D-space + time field). Pass -s 1 to restrict the norm to the spatial
// axes only (all axes except the last).
//
// Inputs:
//   -f  .mfa  diy MFA input file               (MFA  mode, the default)
//   -i  .pt   TorchScript INR model file       (INR  mode; selected when set)
//   -n        INR function name (e.g. vortex_street_3d) -- needed for the
//             model's domain bounds when in INR mode
//   -p  .ply  mesh whose vertices are sampled
//   -s        spatial-only norm (0 = full gradient norm [default], 1 = drop the
//             last/time axis)
//   -o        optional output file: one float64 gradient norm per vertex
//
// Learns its MFA gradient evaluation from contour/*::compute_gradient and
// convert/write_gradient.cpp, and its INR evaluation from INRModel.h, following
// the same conventions used by critical_point_tracking and feature_flow_fields.
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

    // read one scalar of the given type from a binary stream into a double
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
            // strip trailing CR (files authored on Windows)
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

        const bool swap =
            (format == Format::binary_big);   // host assumed little-endian

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
// Gradient-norm helpers.
// ---------------------------------------------------------------------------

// MFA: physical-domain gradient via mfa_extend::recover_mfa (same scheme as
// contour/*::compute_gradient). Returns the (optionally spatial-only) norm.
template<typename T>
T mfa_grad_norm(const Block<T>* b, const VectorX<T>& vertex, bool spatial_only)
{
    const int D = static_cast<int>(b->core_mins.size());

    // build / clamp the query point into the block's domain
    VectorX<T> p(D);
    for (int d = 0; d < D; ++d)
    {
        T c = (d < vertex.size()) ? vertex[d] : T(0);
        if (c < b->core_mins[d]) c = b->core_mins[d];
        if (c > b->core_maxs[d]) c = b->core_maxs[d];
        p[d] = c;
    }

    VectorX<T> grad(D);
    VectorX<T> out(1);
    VectorXi   deriv = VectorXi::Zero(D);
    for (int d = 0; d < D; ++d)
    {
        deriv[d] = 1;
        mfa_extend::recover_mfa(b, p, out, deriv);
        grad[d] = out[0];
        deriv[d] = 0;
    }

    const int n = spatial_only ? std::max(1, D - 1) : D;
    return grad.head(n).norm();
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
    string out_file   = "";             // optional per-vertex norm dump (float64)
    int    spatial    = 0;              // 1 => norm over spatial axes only
    bool   help       = false;

    opts::Options ops;
    ops >> opts::Option('f', "mfa",      mfa_file,  " input MFA .mfa diy file");
    ops >> opts::Option('i', "inr",      inr_model, " input INR TorchScript .pt model (selects INR mode)");
    ops >> opts::Option('n', "name",     inr_name,  " INR function name (domain bounds), e.g. vortex_street_3d");
    ops >> opts::Option('p', "ply",      ply_file,  " input .ply mesh (vertices sampled for gradient norm)");
    ops >> opts::Option('s', "spatial",  spatial,   " 1: gradient norm over spatial axes only (drop last/time axis)");
    ops >> opts::Option('o', "out",      out_file,  " optional output file: one float64 gradient norm per vertex");
    ops >> opts::Option('h', "help",     help,      " show help");

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

    std::vector<T> norms(vertices.size(), T(0));

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

        // batched autograd gradient at all vertices (INR domain is x,y,t -> 3)
        std::vector<VectorX<T>> pts(vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            VectorX<T> p(D);
            for (int d = 0; d < D; ++d)
            {
                T c = (d < vertices[i].size()) ? vertices[i][d] : T(0);
                if (c < model.domain_min[d]) c = model.domain_min[d];
                if (c > model.domain_max[d]) c = model.domain_max[d];
                p[d] = c;
            }
            pts[i] = p;
        }

        std::vector<VectorX<T>> grads;
        model.eval_grad_batch_in_domain(pts, grads);
        for (size_t i = 0; i < grads.size(); ++i)
        {
            const int n = spatial ? std::max(1, (int)grads[i].size() - 1) : (int)grads[i].size();
            norms[i] = grads[i].head(n).norm();
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
            tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size()),
            [&](const tbb::blocked_range<size_t>& range)
            {
                for (size_t i = range.begin(); i != range.end(); ++i)
                    norms[i] = mfa_grad_norm<T>(b, vertices[i], spatial != 0);
            }, ap);
        });
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- report ----
    T min_norm = *std::min_element(norms.begin(), norms.end());
    T max_norm = *std::max_element(norms.begin(), norms.end());
    T sum_norm = std::accumulate(norms.begin(), norms.end(), T(0));
    T avg_norm = sum_norm / static_cast<T>(norms.size());

    std::cout << (spatial ? "spatial " : "full ") << "gradient norm over "
              << norms.size() << " vertices:" << std::endl;
    std::cout << "  min = " << min_norm << std::endl;
    std::cout << "  max = " << max_norm << std::endl;
    std::cout << "  avg = " << avg_norm << std::endl;
    std::cout << "evaluation time (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << std::endl;

    if (!out_file.empty())
    {
        std::ofstream out(out_file.c_str(), std::ios::binary);
        if (out)
        {
            out.write(reinterpret_cast<const char*>(norms.data()),
                      norms.size() * sizeof(T));
            std::cout << "wrote per-vertex norms to " << out_file << std::endl;
        }
        else
            std::cerr << "failed to open " << out_file << " for writing" << std::endl;
    }

    return 0;
}
