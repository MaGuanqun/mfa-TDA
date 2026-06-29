//
// Determinant of the Hessian at a set of input points.
//
// Reads a list of points from either a binary .dat file (the project's
// utility::writeMatrixVector format, as produced for degenerate/critical points
// -- see critical_point_tracking.cpp / Degenerate_case_tracing::read_degenerate_point)
// or a .csv file (one point per row), queries the scalar field from either an
// MFA model (.mfa diy file) or an INR model (TorchScript .pt + function name),
// builds the full domain Hessian at each point and writes its determinant.
//
// MFA Hessians use the analytic second derivatives of the MFA (same scheme as
// contour/*::compute_hessian). INR Hessians are central differences of the
// (autograd) gradient field returned by INRModel::eval_grad_batch_in_domain,
// which gives the full Hessian (including the time--time term).
//
// Inputs:
//   -f  .mfa  diy MFA input file                 (MFA mode, the default)
//   -i  .pt   TorchScript INR model file         (INR mode; selected when set)
//   -n        INR function name (domain bounds), e.g. vortex_street_3d
//   -p        input points file (.dat or .csv)
//   -c        input format override: auto (default) | dat | csv
//   -o        output .csv (point coords + det_hessian); if it ends with .dat,
//             a binary float64 array of determinants is written instead
//
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <numeric>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include <tbb/tbb.h>

#include "opts.h"
#include "block.hpp"
#include "mfa_extend.h"
#include "INRModel.h"
#include "../utility/utility_function.h"

using namespace std;

// ---------------------------------------------------------------------------
// Point readers.
// ---------------------------------------------------------------------------

// .csv : one point per row (comma separated). A non-numeric first line is
// treated as a header and skipped. Only the first D columns are used.
template<typename T>
bool read_points_csv(const std::string& filename, int D, std::vector<VectorX<T>>& points)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
        std::cerr << "read_points_csv: cannot open " << filename << std::endl;
        return false;
    }

    points.clear();
    std::string line;
    bool first = true;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
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
            if (first) { first = false; continue; }   // header line
            continue;
        }
        first = false;

        if (static_cast<int>(vals.size()) < D) continue;
        VectorX<T> p(D);
        for (int a = 0; a < D; ++a) p[a] = vals[a];
        points.emplace_back(std::move(p));
    }
    return true;
}

// .dat : utility::writeMatrixVector binary format (same as the degenerate /
// critical-point files). The first matrix holds one point per row; only the
// first D columns are used.
template<typename T>
bool read_points_dat(const std::string& filename, int D, std::vector<VectorX<T>>& points)
{
    std::ifstream probe(filename.c_str(), std::ios::binary);
    if (!probe)
    {
        std::cerr << "read_points_dat: cannot open " << filename << std::endl;
        return false;
    }
    probe.close();

    std::vector<MatrixX<T>> root;
    utility::loadMatrixVector(filename.c_str(), root);
    if (root.empty() || root[0].rows() == 0)
    {
        std::cerr << "read_points_dat: no points in " << filename << std::endl;
        return false;
    }
    if (root[0].cols() < D)
    {
        std::cerr << "read_points_dat: points have " << root[0].cols()
                  << " columns but " << D << " are required" << std::endl;
        return false;
    }

    points.resize(root[0].rows());
    for (int i = 0; i < root[0].rows(); ++i)
        points[i] = root[0].row(i).transpose().head(D);
    return true;
}

// ---------------------------------------------------------------------------
// Hessian determinant helpers.
// ---------------------------------------------------------------------------

// det(H) and det(H)/||H||^2 for one Hessian (||.|| is the Frobenius norm; the
// normalized value is 0 when ||H|| == 0 to avoid dividing by zero).
template<typename T>
inline void hessian_dets(const MatrixX<T>& H, T& det, T& det_norm)
{
    det = std::abs(H.determinant());
    const T n2 = H.squaredNorm();
    det_norm = (n2 > T(0)) ? det / n2 : T(0);
}

// MFA: determinant of the (spatial) Hessian from analytic second derivatives
// (mfa_extend::recover_mfa). hdim selects the leading sub-block: hdim = D-1 gives
// the spatial Hessian (drop the time axis), hdim = D the full Hessian.
// Returns det(H) in `det` and det(H)/||H||^2 in `det_norm`.
template<typename T>
void mfa_hessian_det(const Block<T>* b, const VectorX<T>& vin, int hdim,
                     T& det, T& det_norm)
{
    const int D = static_cast<int>(b->core_mins.size());

    VectorX<T> p(D);
    for (int d = 0; d < D; ++d)
    {
        T c = (d < vin.size()) ? vin[d] : T(0);
        if (c < b->core_mins[d]) c = b->core_mins[d];
        if (c > b->core_maxs[d]) c = b->core_maxs[d];
        p[d] = c;
    }

    MatrixX<T> H(hdim, hdim);
    VectorX<T> out(1);
    VectorXi   deriv = VectorXi::Zero(D);
    for (int i = 0; i < hdim; ++i)
        for (int j = i; j < hdim; ++j)
        {
            deriv.setZero();
            deriv[i] += 1;
            deriv[j] += 1;
            mfa_extend::recover_mfa(b, p, out, deriv);
            H(i, j) = out[0];
            H(j, i) = out[0];
        }

    hessian_dets<T>(H, det, det_norm);
}

// INR: determinant of the (spatial) Hessian using the exact nested-autograd
// Hessian (INRModel::query_hessian_autograd), in physical-domain coordinates.
// hdim selects the leading sub-block: hdim = D-1 -> spatial Hessian (e.g. 2x2),
// hdim = D -> full Hessian (including the time-time term).
// Fills dets[i] = det(H) and dets_norm[i] = det(H)/||H||^2.
template<typename T>
void inr_hessian_det(INRModel<T>& model, const std::vector<VectorX<T>>& points,
                     int hdim, std::vector<T>& dets, std::vector<T>& dets_norm)
{
    const int D = static_cast<int>(model.domain_min.size());

    dets.resize(points.size());
    dets_norm.resize(points.size());
    Eigen::MatrixX<T> H3;
    Eigen::MatrixX<T> H;
    for (size_t i = 0; i < points.size(); ++i)
    {
        VectorX<T> p(D);
        for (int d = 0; d < D; ++d)
        {
            T c = (d < points[i].size()) ? points[i][d] : T(0);
            if (c < model.domain_min[d]) c = model.domain_min[d];
            if (c > model.domain_max[d]) c = model.domain_max[d];
            p[d] = c;
        }

        model.query_hessian_autograd(p, H3);            // exact 3x3 Hessian (domain coords)
        H = H3.topLeftCorner(hdim, hdim);
        hessian_dets<T>(H, dets[i], dets_norm[i]);
    }
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
    using T = double;

    diy::mpi::environment  env(argc, argv);
    diy::mpi::communicator world;

    string mfa_file   = "approx.mfa";
    string inr_model  = "";
    string inr_name   = "";
    string point_file = "points.dat";
    string fmt        = "auto";          // auto | dat | csv
    string out_file   = "";
    int    spatial    = 1;               // 1: spatial Hessian (drop time axis); 0: full DxD
    bool   help       = false;

    opts::Options ops;
    ops >> opts::Option('f', "mfa",     mfa_file,   " input MFA .mfa diy file");
    ops >> opts::Option('i', "inr",     inr_model,  " input INR TorchScript .pt model (selects INR mode)");
    ops >> opts::Option('n', "name",    inr_name,   " INR function name (domain bounds), e.g. vortex_street_3d");
    ops >> opts::Option('p', "points",  point_file, " input points file (.dat or .csv)");
    ops >> opts::Option('c', "format",  fmt,        " input format: auto | dat | csv");
    ops >> opts::Option('s', "spatial", spatial,    " 1: spatial Hessian, e.g. 2x2 for 2D+time [default]; 0: full DxD Hessian");
    ops >> opts::Option('o', "out",     out_file,   " output .csv (coords + det); .dat => binary float64 dets");
    ops >> opts::Option('h', "help",    help,       " show help");

    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    const bool use_inr = !inr_model.empty();

    // resolve input format
    bool is_csv = false;
    if (fmt == "csv")      is_csv = true;
    else if (fmt == "dat") is_csv = false;
    else // auto: by extension
    {
        is_csv = (point_file.size() >= 4 &&
                  point_file.compare(point_file.size() - 4, 4, ".csv") == 0);
    }

    std::vector<VectorX<T>> points;
    std::vector<T>          dets;        // det(H)
    std::vector<T>          dets_norm;   // det(H)/||H||^2

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
        const int hdim = (spatial != 0) ? std::max(1, D - 1) : D;
        std::cout << "INR model '" << inr_name << "' domain dim " << D
                  << " min " << model.domain_min.transpose()
                  << " max " << model.domain_max.transpose() << std::endl;
        std::cout << "Hessian: " << hdim << "x" << hdim
                  << (spatial ? " (spatial)" : " (full)") << std::endl;

        bool ok = is_csv ? read_points_csv<T>(point_file, D, points)
                         : read_points_dat<T>(point_file, D, points);
        if (!ok) return 1;
        std::cout << "read " << points.size() << " points from " << point_file
                  << (is_csv ? " (csv)" : " (dat)") << std::endl;
        if (points.empty()) { std::cerr << "no points -- nothing to do" << std::endl; return 1; }

        inr_hessian_det<T>(model, points, hdim, dets, dets_norm);
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
        if (master.size() == 0) { std::cerr << "no blocks in " << mfa_file << std::endl; return 1; }

        master.foreach([&](Block<T>* b, const diy::Master::ProxyWithLink&)
        {
            const int D = static_cast<int>(b->core_mins.size());
            const int hdim = (spatial != 0) ? std::max(1, D - 1) : D;
            std::cout << "MFA domain dim " << D
                      << " min " << b->core_mins.transpose()
                      << " max " << b->core_maxs.transpose() << std::endl;
            std::cout << "Hessian: " << hdim << "x" << hdim
                      << (spatial ? " (spatial)" : " (full)") << std::endl;

            bool ok = is_csv ? read_points_csv<T>(point_file, D, points)
                             : read_points_dat<T>(point_file, D, points);
            if (!ok) return;
            std::cout << "read " << points.size() << " points from " << point_file
                      << (is_csv ? " (csv)" : " (dat)") << std::endl;
            if (points.empty()) return;

            dets.assign(points.size(), T(0));
            dets_norm.assign(points.size(), T(0));
            tbb::affinity_partitioner ap;
            tbb::parallel_for(tbb::blocked_range<size_t>(0, points.size()),
            [&](const tbb::blocked_range<size_t>& range)
            {
                for (size_t i = range.begin(); i != range.end(); ++i)
                    mfa_hessian_det<T>(b, points[i], hdim, dets[i], dets_norm[i]);
            }, ap);
        });

        if (points.empty() || dets.empty())
        {
            std::cerr << "no points evaluated" << std::endl;
            return 1;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- report (mean and max of det(H) and det(H)/||H||^2) ----
    const T n = static_cast<T>(dets.size());
    const T mean_det  = std::accumulate(dets.begin(),      dets.end(),      T(0)) / n;
    const T max_det   = *std::max_element(dets.begin(),      dets.end());
    const T mean_detn = std::accumulate(dets_norm.begin(), dets_norm.end(), T(0)) / n;
    const T max_detn  = *std::max_element(dets_norm.begin(), dets_norm.end());

    std::cout << std::setprecision(10);
    std::cout << "Hessian metrics over " << dets.size() << " points:" << std::endl;
    std::cout << "  det(H):           mean = " << mean_det  << "  max = " << max_det  << std::endl;
    std::cout << "  det(H)/||H||^2:   mean = " << mean_detn << "  max = " << max_detn << std::endl;
    std::cout << "evaluation time (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
              << std::endl;

    if(out_file.empty()){
        return 0;
    }

    // ---- write ----
    const bool out_dat = (out_file.size() >= 4 &&
                          out_file.compare(out_file.size() - 4, 4, ".dat") == 0);
    if (out_dat)
    {
        // binary float64: all det(H) values, followed by all det(H)/||H||^2 values
        std::ofstream out(out_file.c_str(), std::ios::binary);
        if (!out) { std::cerr << "cannot open " << out_file << " for writing" << std::endl; return 1; }
        out.write(reinterpret_cast<const char*>(dets.data()),      dets.size()      * sizeof(T));
        out.write(reinterpret_cast<const char*>(dets_norm.data()), dets_norm.size() * sizeof(T));
        std::cout << "wrote " << dets.size() << " det(H) + " << dets_norm.size()
                  << " det(H)/||H||^2 float64 values to " << out_file << std::endl;
    }
    else
    {
        std::ofstream out(out_file.c_str());
        if (!out) { std::cerr << "cannot open " << out_file << " for writing" << std::endl; return 1; }
        const int D = static_cast<int>(points[0].size());
        for (int d = 0; d < D; ++d) out << "x" << d << ",";
        out << "det_hessian,det_hessian_over_normsq\n";
        out << std::setprecision(17);
        for (size_t i = 0; i < points.size(); ++i)
        {
            for (int d = 0; d < D; ++d) out << points[i][d] << ",";
            out << dets[i] << "," << dets_norm[i] << "\n";
        }
        std::cout << "wrote per-point det(H) and det(H)/||H||^2 to " << out_file << std::endl;
    }

    return 0;
}
