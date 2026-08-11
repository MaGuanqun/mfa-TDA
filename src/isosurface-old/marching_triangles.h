#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <array>
#include <unordered_set>
#include <cstdint>

#include <Eigen/Dense>
#include <mfa/mfa.hpp>

#include "query_function.h"
#include "utility_function.h"
#include "block.hpp"

namespace marching_triangles
{

using Triangle = std::array<int, 3>;

inline Triangle make_triangle_key(int a, int b, int c) noexcept
{
    const int v[3] = {a, b, c};
    int min_i = 0;
    for (int i = 1; i < 3; ++i)
    {
        if (v[i] < v[min_i])
            min_i = i;
    }
    return {v[min_i], v[(min_i + 1) % 3], v[(min_i + 2) % 3]};
}

inline uint64_t pack_triangle(int a, int b, int c) noexcept
{
    return (static_cast<uint64_t>(a) << 42)
         | (static_cast<uint64_t>(b) << 21)
         | static_cast<uint64_t>(c);
}

inline uint64_t triangle_key_packed(int a, int b, int c)
{
    const Triangle tri = make_triangle_key(a, b, c);
    return pack_triangle(tri[0], tri[1], tri[2]);
}

/// Store triangle in boundary order (vk, v_mid, vk1); dedupe via make_triangle_key.
inline bool add_triangle_unique(std::vector<Triangle>& triangles,
                                std::unordered_set<uint64_t>& seen,
                                int a,
                                int b,
                                int c)
{
    const uint64_t key = triangle_key_packed(a, b, c);
    if (!seen.insert(key).second)
        return false;
    triangles.push_back({a, b, c});
    return true;
}

inline void seed_triangle_keys(const std::vector<Triangle>& triangles,
                               std::unordered_set<uint64_t>& seen)
{
    seen.clear();
    seen.reserve(triangles.size());
    for (const Triangle& tri : triangles)
        seen.insert(triangle_key_packed(tri[0], tri[1], tri[2]));
}

} // namespace marching_triangles

#include "crack_closing.h"

namespace marching_triangles
{

template<typename T>
class MarchingTriangles
{
public:
    MarchingTriangles(
        const VectorX<T>& domain_min_,
        const VectorX<T>& domain_max_,
        int function_type_,
        T iso_value_,
        T root_finding_epsilon_,
        int max_projection_itr_,
        T same_vertex_epsilon_,
        T hessian_rank_threshold_,
        T stop_curve_distance_,
        T d_min_,
        T d_step_,
        T projection_alpha_ = static_cast<T>(1.5),
        const Block<T>* b_ = nullptr,
        INRModel<T>* inr_model_ = nullptr,
        bool enable_crack_closing_ = false)
        : domain_min(domain_min_)
        , domain_max(domain_max_)
        , function_type(function_type_)
        , iso_value(iso_value_)
        , root_finding_epsilon(root_finding_epsilon_)
        , max_projection_itr(max_projection_itr_)
        , same_vertex_epsilon(same_vertex_epsilon_)
        , hessian_rank_threshold(hessian_rank_threshold_)
        , stop_curve_distance(stop_curve_distance_)
        , d_min(d_min_)
        , d_step(d_step_)
        , d_max(2.0 * d_step_)
        , projection_alpha(projection_alpha_)
        , b(b_)
        , inr_model(inr_model_)
        , enable_crack_closing(enable_crack_closing_)
    {}

    void set_degenerate_points(const std::vector<VectorX<T>>& pts) { degenerate_points = pts; }

    bool extract_all_sheets(const std::vector<VectorX<T>>& seed_roots,
    std::vector<std::vector<VectorX<T>>>& sheet_vertices,
    std::vector<std::vector<Triangle>>& sheet_triangles)
    {
        sheet_vertices.clear();
        sheet_triangles.clear();
        if (seed_roots.empty())
            return false;

        std::vector<bool> covered(seed_roots.size(), false);
        int sheet_id = 0;

        // while (true)
        // {
            int start_idx = 0;

            std::vector<VectorX<T>> verts;
            std::vector<Triangle> tris;
            std::vector<bool> local_covered = covered;

            if (!extract_one_sheet(seed_roots[start_idx], verts, tris, local_covered, seed_roots))
            {
                covered[start_idx] = true;
                // continue;
            }

            // mark_covered_seeds(verts, seed_roots, local_covered);
            // covered = local_covered;

            sheet_vertices.push_back(std::move(verts));
            sheet_triangles.push_back(std::move(tris));
            std::cout << "sheet " << sheet_id++ << " triangles " << sheet_triangles.back().size()
                      << " vertices " << sheet_vertices.back().size() << std::endl;
        // }

        return !sheet_vertices.empty();
    }




    static void save_mesh_obj(const std::string& filename,
                              const std::vector<VectorX<T>>& vertices,
                              const std::vector<Triangle>& triangles)
    {
        std::ofstream out(filename);
        if (!out.is_open())
        {
            std::cerr << "Cannot write " << filename << std::endl;
            return;
        }
        out << std::setprecision(15);
        for (const auto& v : vertices)
        {
            out << "v";
            for (int i = 0; i < v.size(); ++i)
                out << " " << v[i];
            out << "\n";
        }
        for (const auto& tri : triangles)
            out << "f " << (tri[0] + 1) << " " << (tri[1] + 1) << " " << (tri[2] + 1) << "\n";
        out.close();
    }

private:
    const VectorX<T> domain_min;
    const VectorX<T> domain_max;
    const int function_type;
    const T iso_value;
    const T root_finding_epsilon;
    const int max_projection_itr;
    const T same_vertex_epsilon;
    const T hessian_rank_threshold;
    const T stop_curve_distance;
    const T d_min;
    const T d_step;
    const T d_max;
    const T projection_alpha;
    const Block<T>* b;
    INRModel<T>* inr_model;
    const bool enable_crack_closing;

    std::vector<VectorX<T>> degenerate_points;

    // Per-vertex incident triangle lists (cleared each sheet). Speeds up Delaunay checks.
    mutable std::vector<std::vector<int>> vertex_triangle_incident_;
    mutable std::vector<int> tri_visit_stamp_;
    mutable std::vector<int> tri_visit_candidates_;
    mutable int tri_visit_generation_ = 1;

    void clear_triangle_index() const
    {
        vertex_triangle_incident_.clear();
        tri_visit_stamp_.clear();
        tri_visit_candidates_.clear();
        tri_visit_generation_ = 1;
    }

    /// One or more closed vertex rings; each loop[i] -> loop[(i+1)%n] is a directed frontier edge.
    struct BoundaryLoops
    {
        std::vector<std::vector<int>> loops;

        void clear() { loops.clear(); }

        void prune_small_loops()
        {
            std::vector<std::vector<int>> kept;
            kept.reserve(loops.size());
            for (auto& loop : loops)
            {
                if (loop.size() >= 3)
                    kept.push_back(std::move(loop));
            }
            loops.swap(kept);
        }

        bool is_vertex_on_loop(int vid) const
        {
            for (const auto& loop : loops)
            {
                for (int v : loop)
                {
                    if (v == vid)
                        return true;
                }
            }
            return false;
        }

        static bool find_directed_edge(const std::vector<int>& loop, int from, int to, size_t& edge_i)
        {
            if (loop.size() < 3)
                return false;
            for (size_t i = 0; i < loop.size(); ++i)
            {
                const size_t j = (i + 1) % loop.size();
                if (loop[i] == from && loop[j] == to)
                {
                    edge_i = i;
                    return true;
                }
            }
            return false;
        }

        size_t edge_count() const
        {
            size_t count = 0;
            for (const auto& loop : loops)
            {
                if (loop.size() >= 3)
                    count += loop.size();
            }
            return count;
        }

        bool edge_at(size_t global_edge, size_t& loop_i, size_t& edge_i) const
        {
            size_t offset = 0;
            for (size_t li = 0; li < loops.size(); ++li)
            {
                const size_t n = loops[li].size() >= 3 ? loops[li].size() : 0;
                if (global_edge < offset + n)
                {
                    loop_i = li;
                    edge_i = global_edge - offset;
                    return true;
                }
                offset += n;
            }
            return false;
        }

        size_t global_edge_index(size_t loop_i, size_t edge_i) const
        {
            size_t g = 0;
            for (size_t li = 0; li < loop_i; ++li)
            {
                if (loops[li].size() >= 3)
                    g += loops[li].size();
            }
            return g + edge_i;
        }

        bool get_edge_vertices(size_t loop_i, size_t edge_i, int& vk, int& vk1) const
        {
            if (loop_i >= loops.size())
                return false;
            const std::vector<int>& loop = loops[loop_i];
            if (loop.size() < 3 || edge_i >= loop.size())
                return false;
            vk = loop[edge_i];
            vk1 = loop[(edge_i + 1) % loop.size()];
            return true;
        }

        int vertex_before(size_t loop_i, size_t edge_i) const
        {
            if (loop_i >= loops.size())
                return -1;
            const std::vector<int>& loop = loops[loop_i];
            if (loop.size() < 3 || edge_i >= loop.size())
                return -1;
            return loop[(edge_i + loop.size() - 1) % loop.size()];
        }

        int vertex_after(size_t loop_i, size_t edge_i) const
        {
            if (loop_i >= loops.size())
                return -1;
            const std::vector<int>& loop = loops[loop_i];
            if (loop.size() < 3 || edge_i >= loop.size())
                return -1;
            return loop[(edge_i + 2) % loop.size()];
        }

        bool can_consume_edge_with_new_vertex(size_t loop_i, int vk, int vk1) const
        {
            if (loop_i >= loops.size())
                return false;
            size_t ei = 0;
            return find_directed_edge(loops[loop_i], vk, vk1, ei);
        }

        bool can_consume_edge_predecessor_wedge(size_t loop_i, int vk, int vk1, int vm1) const
        {
            if (loop_i >= loops.size())
                return false;
            const std::vector<int>& loop = loops[loop_i];
            size_t ei = 0;
            if (!find_directed_edge(loop, vk, vk1, ei))
                return false;
            return loop.size() >= 4
                && loop[(ei + loop.size() - 1) % loop.size()] == vm1;
        }

        bool can_consume_edge_successor_wedge(size_t loop_i, int vk, int vk1, int vp2) const
        {
            if (loop_i >= loops.size())
                return false;
            const std::vector<int>& loop = loops[loop_i];
            size_t ei = 0;
            if (!find_directed_edge(loop, vk, vk1, ei))
                return false;
            const size_t n = loop.size();
            return loop[(ei + 1) % n] == vk1 && loop[(ei + 2) % n] == vp2;
        }

        /// Used edge vk->vk1 is removed; new edges vk->vp and vp->vk1 are added.
        bool consume_edge_with_new_vertex(size_t loop_i, int vk, int vk1, int vp)
        {
            if (loop_i >= loops.size())
                return false;
            std::vector<int>& loop = loops[loop_i];
            size_t ei = 0;
            if (!find_directed_edge(loop, vk, vk1, ei))
                return false;
            loop.insert(loop.begin() + static_cast<std::vector<int>::difference_type>(ei + 1), vp);
            return true;
        }

        /// Wedge with predecessor vm1 on loop ... vm1, vk, vk1: remove used edge vk->vk1 and vertex vk.
        bool consume_edge_predecessor_wedge(size_t loop_i, int vk, int vk1, int vm1)
        {
            if (loop_i >= loops.size())
                return false;
            std::vector<int>& loop = loops[loop_i];
            size_t ei = 0;
            if (!find_directed_edge(loop, vk, vk1, ei))
                return false;
            if (loop[(ei + loop.size() - 1) % loop.size()] != vm1)
                return false;
            if (loop.size() < 4)
                return false;
            loop.erase(loop.begin() + static_cast<std::vector<int>::difference_type>(ei));
            return true;
        }

        /// Wedge with successor vp2: remove edge vk->vk1; loop becomes ... vk1, vk, vp2 ...
        bool consume_edge_successor_wedge(size_t loop_i, int vk, int vk1, int vp2)
        {
            if (loop_i >= loops.size())
                return false;
            std::vector<int>& loop = loops[loop_i];
            size_t ei = 0;
            if (!find_directed_edge(loop, vk, vk1, ei))
                return false;
            const size_t n = loop.size();
            if (loop[(ei + 1) % n] != vk1 || loop[(ei + 2) % n] != vp2)
                return false;

            std::vector<int> next;
            next.reserve(n);
            for (size_t i = 0; i < n; ++i)
            {
                if (i != ei)
                    next.push_back(loop[i]);
            }
            for (size_t i = 0; i < next.size(); ++i)
            {
                if (next[i] == vk1)
                {
                    next.insert(next.begin() + static_cast<std::vector<int>::difference_type>(i + 1), vk);
                    break;
                }
            }
            if (next.size() < 3)
                return false;
            loop.swap(next);
            return true;
        }

        std::vector<int> largest_loop() const
        {
            size_t best = 0;
            size_t best_i = 0;
            for (size_t i = 0; i < loops.size(); ++i)
            {
                if (loops[i].size() > best)
                {
                    best = loops[i].size();
                    best_i = i;
                }
            }
            if (best >= 3)
                return loops[best_i];
            return {};
        }
    };

    /// True if some triangle already uses va->vb as its frontier (tri[0]->tri[2]).
    bool directed_frontier_edge_taken(int va,
                                      int vb,
                                      const std::vector<Triangle>& triangles) const
    {
        if (va < 0 || vb < 0)
            return false;
        for (const int v : std::array<int, 2>{va, vb})
        {
            if (v >= static_cast<int>(vertex_triangle_incident_.size()))
                continue;
            for (const int ti : vertex_triangle_incident_[static_cast<size_t>(v)])
            {
                const Triangle& tri = triangles[static_cast<size_t>(ti)];
                if (tri[0] == va && tri[2] == vb)
                    return true;
            }
        }
        return false;
    }

    /// Triangle stored as (vk, v_mid, vk1). Reject duplicate frontier edge or opposite orientation.
    bool boundary_triangle_orientation_ok(int vk,
                                        int v_mid,
                                        int vk1,
                                        const std::vector<VectorX<T>>& vertices,
                                        const std::vector<Triangle>& triangles) const
    {
        const Vector3<T> n_cand = triangle_normal(vertices[vk], vertices[v_mid], vertices[vk1]);
        if (n_cand.squaredNorm() <= std::numeric_limits<T>::epsilon())
            return false;

        if (directed_frontier_edge_taken(vk, vk1, triangles))
            return false;

        collect_incident_triangles(vk, vk1, v_mid);
        for (const int ti : tri_visit_candidates_)
        {
            const Triangle& tri = triangles[static_cast<size_t>(ti)];
            int shared = 0;
            for (const int v : tri)
            {
                if (v == vk || v == v_mid || v == vk1)
                    ++shared;
            }
            if (shared < 2)
                continue;

            const Vector3<T> n_tri = triangle_normal(
                vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
            if (n_tri.squaredNorm() <= std::numeric_limits<T>::epsilon())
                continue;
            if (n_tri.dot(n_cand) <= 0)
                return false;
        }
        return true;
    }

    /// Store triangle in boundary order (a,b,c) = (vk, v_mid, vk1).
    bool register_triangle(int a,
                           int b,
                           int c,
                           const std::vector<VectorX<T>>& vertices,
                           std::vector<Triangle>& triangles) const
    {
        if (!boundary_triangle_orientation_ok(a, b, c, vertices, triangles))
            return false;
        if (triangle_exists_adjacent(a, b, c, triangles))
            return false;

        const int ti = static_cast<int>(triangles.size());
        triangles.push_back({a, b, c});
        for (const int v : triangles[static_cast<size_t>(ti)])
        {
            if (v >= static_cast<int>(vertex_triangle_incident_.size()))
                vertex_triangle_incident_.resize(static_cast<size_t>(v) + 1);
            vertex_triangle_incident_[static_cast<size_t>(v)].push_back(ti);
        }
        if (ti >= static_cast<int>(tri_visit_stamp_.size()))
            tri_visit_stamp_.resize(static_cast<size_t>(ti) + 1, 0);
        return true;
    }

    bool triangle_exists_adjacent(int a, int b, int c, const std::vector<Triangle>& triangles) const
    {
        const uint64_t key = triangle_key_packed(a, b, c);
        collect_incident_triangles(a, b, c);
        for (const int ti : tri_visit_candidates_)
        {
            const Triangle& tri = triangles[static_cast<size_t>(ti)];
            if (triangle_key_packed(tri[0], tri[1], tri[2]) == key)
                return true;
        }
        return false;
    }

    static void compact_duplicate_triangles(std::vector<Triangle>& triangles)
    {
        std::unordered_set<uint64_t> seen;
        seen.reserve(triangles.size());
        std::vector<Triangle> unique;
        unique.reserve(triangles.size());
        for (const Triangle& tri : triangles)
        {
            const uint64_t key = triangle_key_packed(tri[0], tri[1], tri[2]);
            if (seen.insert(key).second)
                unique.push_back(tri);
        }
        triangles.swap(unique);
    }

    void collect_incident_triangles(int v0, int v1, int vp) const
    {
        tri_visit_candidates_.clear();
        if (++tri_visit_generation_ == 0)
        {
            tri_visit_generation_ = 1;
            std::fill(tri_visit_stamp_.begin(), tri_visit_stamp_.end(), 0);
        }

        const int verts[3] = {v0, v1, vp};
        for (int v : verts)
        {
            if (v < 0 || v >= static_cast<int>(vertex_triangle_incident_.size()))
                continue;
            for (const int ti : vertex_triangle_incident_[static_cast<size_t>(v)])
            {
                if (ti >= static_cast<int>(tri_visit_stamp_.size()))
                    tri_visit_stamp_.resize(static_cast<size_t>(ti) + 1, 0);
                if (tri_visit_stamp_[static_cast<size_t>(ti)] == tri_visit_generation_)
                    continue;
                tri_visit_stamp_[static_cast<size_t>(ti)] = tri_visit_generation_;
                tri_visit_candidates_.push_back(ti);
            }
        }
    }

    void query_value(const VectorX<T>& p, T& value) const
    {
        VectorX<T> out(1);
        query_function::query_function(p, out, function_type, b, VectorXi(), inr_model);
        value = out[0] - iso_value;
    }

    void compute_gradient(const VectorX<T>& p, VectorX<T>& grad) const
    {
        grad.resize(p.size());
        VectorX<T> out(1);
        VectorXi deriv = VectorXi::Zero(p.size());
        for (int i = 0; i < p.size(); ++i)
        {
            deriv[i] = 1;
            query_function::query_function(p, out, function_type, b, deriv);
            grad[i] = out[0];
            deriv[i] = 0;
        }
    }

    void compute_hessian(const VectorX<T>& p, MatrixX<T>& H) const
    {
        const int n = p.size();
        H.resize(n, n);
        H.setZero();
        VectorX<T> out(1);
        VectorXi deriv = VectorXi::Zero(n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = i; j < n; ++j)
            {
                deriv.setZero();
                deriv[i] += 1;
                deriv[j] += 1;
                query_function::query_function(p, out, function_type, b, deriv);
                H(i, j) = out[0];
                H(j, i) = out[0];
                deriv[i] = 0;
                deriv[j] = 0;
            }
        }
    }

    bool is_full_rank_point(const VectorX<T>& p) const
    {
        MatrixX<T> H;
        compute_hessian(p, H);
        Eigen::SelfAdjointEigenSolver<MatrixX<T>> es(H);
        if (es.info() != Eigen::Success)
            return false;
        const T min_ev = es.eigenvalues().minCoeff();
        return std::abs(min_ev) > hessian_rank_threshold;
    }

    bool project_to_surface(VectorX<T>& p) const
    {
        if (!utility::In_Domain(p, domain_min, domain_max))
            return false;

        T fval = 0;
        query_value(p, fval);
        if (std::abs(fval) < root_finding_epsilon)
            return true;

        for (int itr = 0; itr < max_projection_itr; ++itr)
        {
            VectorX<T> grad;
            compute_gradient(p, grad);
            const T g2 = grad.squaredNorm();
            if (g2 < std::numeric_limits<T>::epsilon())
                return false;

            p -=  fval * grad / g2;

            if (!utility::In_Domain(p, domain_min, domain_max))
                return false;

            query_value(p, fval);
            if (std::abs(fval) < root_finding_epsilon)
                return true;
        }
        return std::abs(fval) < root_finding_epsilon;
    }

    int find_or_add_vertex(std::vector<VectorX<T>>& vertices, const VectorX<T>& p) const
    {
        const T eps2 = same_vertex_epsilon * same_vertex_epsilon;
        for (int i = 0; i < static_cast<int>(vertices.size()); ++i)
        {
            if ((vertices[i] - p).squaredNorm() < eps2)
                return i;
        }
        vertices.push_back(p);
        return static_cast<int>(vertices.size()) - 1;
    }

    bool near_stop_set(const VectorX<T>& p) const
    {
        const T d2 = stop_curve_distance * stop_curve_distance;
        for (const auto& q : degenerate_points)
        {
            if ((p - q).squaredNorm() < d2)
                return true;
        }
        return false;
    }

    // 3D cross product (domain is always 3D; VectorX kept for query_function API).
    static VectorX<T> cross3(const VectorX<T>& u, const VectorX<T>& v)
    {
        VectorX<T> result(3);
        result[0] = u[1] * v[2] - u[2] * v[1];
        result[1] = u[2] * v[0] - u[0] * v[2];
        result[2] = u[0] * v[1] - u[1] * v[0];
        return result;
    }


    bool circumsphere(const VectorX<T>& a, const VectorX<T>& b, const VectorX<T>& c, VectorX<T>& center, T& radius) const
    {
        const VectorX<T> v1 = b - a;
        const VectorX<T> v2 = c - a;
        const VectorX<T> n = cross3(v1, v2);
        const T denom = 2 * n.squaredNorm();
        if (denom < std::numeric_limits<T>::epsilon())
            return false;

        const VectorX<T> c3 = a
            + (v1.squaredNorm() * cross3(v2, n) - v2.squaredNorm() * cross3(v1, n)) / denom;
        center = c3;
        radius = (c3 - a).norm();
        return true;
    }

    Vector3<T> triangle_normal(const VectorX<T>& a, const VectorX<T>& b, const VectorX<T>& c) const
    {
        VectorX<T> x1=b - a;
        VectorX<T> x2=c - a;
        return cross3(x1, x2);
    }

    bool delaunay_constraint_ok(const std::vector<VectorX<T>>& vertices,
                                const std::vector<Triangle>& triangles,
                                int v0, int v1, int vp) const
    {
        if (vertices[v0].size() != 3)
            return false;

        VectorX<T> center;
        T radius = 0;
        if (!circumsphere(vertices[v0], vertices[v1], vertices[vp], center, radius))
            return false;

        const Vector3<T> n_new = triangle_normal(vertices[v0], vertices[v1], vertices[vp]);
        if (n_new.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;

        const T r2 = radius * radius;
        collect_incident_triangles(v0, v1, vp);

        for (const int ti : tri_visit_candidates_)
        {
            const Triangle& tri = triangles[static_cast<size_t>(ti)];
            for (int j = 0; j < 3; ++j)
            {
                const int vid = tri[j];
                if (vid == v0 || vid == v1 || vid == vp)
                    continue;
                if ((vertices[vid] - center).squaredNorm() < r2)
                {
                    const Vector3<T> n_old = triangle_normal(
                        vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]);
                    if (n_old.dot(n_new) <= 0)
                        continue;
                    return false;
                }
            }
        }
        return true;
    }

    T edge_length(const VectorX<T>& a, const VectorX<T>& b) const
    {
        return (a - b).norm();
    }

    T blended_step(T d) const
    {
        if (d < d_min)
            d = static_cast<T>(0.75) * d + static_cast<T>(0.25) * d_min;
        return std::min(d, d_max);
    }

    bool triangle_edges_within_limit(const std::vector<VectorX<T>>& vertices,
                                     int a, int b, int c) const
    {
        const T max_e2 = d_max * d_max;
        return (vertices[a] - vertices[b]).squaredNorm() <= max_e2
            && (vertices[b] - vertices[c]).squaredNorm() <= max_e2
            && (vertices[c] - vertices[a]).squaredNorm() <= max_e2;
    }

    bool tangent_frame_at(const VectorX<T>& p, VectorX<T>& t1, VectorX<T>& t2) const
    {
        VectorX<T> grad;
        compute_gradient(p, grad);
        if (grad.norm() < std::numeric_limits<T>::epsilon())
            return false;

        VectorX<T> axis = VectorX<T>::Zero(grad.size());
        axis[0] = 1;
        t1 = cross3(grad, axis);
        if (t1.squaredNorm() < std::numeric_limits<T>::epsilon())
        {
            axis[0] = 0;
            axis[1] = 1;
            t1 = cross3(grad, axis);
        }
        if (t1.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;
        t1.normalize();
        t2 = cross3(grad, t1);
        if (t2.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;
        t2.normalize();
        return true;
    }

    bool try_wedge_predecessor_triangle(std::vector<VectorX<T>>& vertices,
                                        std::vector<Triangle>& triangles,
                                        BoundaryLoops& boundary,
                                        size_t loop_i,
                                        size_t edge_i,
                                        int vk,
                                        int vk1,
                                        int vm1,
                                        size_t& edge_cursor) const
    {
        if (vm1 == vk || vm1 == vk1)
            return false;
        if (boundary.vertex_before(loop_i, edge_i) != vm1)
            return false;
        if (!triangle_edges_within_limit(vertices, vk, vm1, vk1))
            return false;
        if (!delaunay_constraint_ok(vertices, triangles, vk, vm1, vk1))
            return false;
        if (!boundary.can_consume_edge_predecessor_wedge(loop_i, vk, vk1, vm1))
            return false;
        if (!boundary_triangle_orientation_ok(vk, vm1, vk1, vertices, triangles))
            return false;
        if (triangle_exists_adjacent(vk, vm1, vk1, triangles))
            return false;

        const std::vector<std::vector<int>> backup = boundary.loops;
        if (!boundary.consume_edge_predecessor_wedge(loop_i, vk, vk1, vm1))
        {
            boundary.loops = backup;
            return false;
        }
        if (!register_triangle(vk, vm1, vk1, vertices, triangles))
        {
            boundary.loops = backup;
            return false;
        }

        const size_t n_edges = boundary.edge_count();
        edge_cursor = n_edges > 0 ? boundary.global_edge_index(loop_i, edge_i) % n_edges : 0;
        return true;
    }

    bool try_wedge_successor_triangle(std::vector<VectorX<T>>& vertices,
                                      std::vector<Triangle>& triangles,
                                      BoundaryLoops& boundary,
                                      size_t loop_i,
                                      size_t edge_i,
                                      int vk,
                                      int vk1,
                                      int vp2,
                                      size_t& edge_cursor) const
    {
        if (vp2 == vk || vp2 == vk1)
            return false;
        if (boundary.vertex_after(loop_i, edge_i) != vp2)
            return false;
        if (!triangle_edges_within_limit(vertices, vk1, vk, vp2))
            return false;
        if (!delaunay_constraint_ok(vertices, triangles, vk1, vk, vp2))
            return false;
        if (!boundary.can_consume_edge_successor_wedge(loop_i, vk, vk1, vp2))
            return false;
        if (!boundary_triangle_orientation_ok(vk1, vk, vp2, vertices, triangles))
            return false;
        if (triangle_exists_adjacent(vk1, vk, vp2, triangles))
            return false;

        const std::vector<std::vector<int>> backup = boundary.loops;
        if (!boundary.consume_edge_successor_wedge(loop_i, vk, vk1, vp2))
        {
            boundary.loops = backup;
            return false;
        }
        if (!register_triangle(vk1, vk, vp2, vertices, triangles))
        {
            boundary.loops = backup;
            return false;
        }

        const size_t n_edges = boundary.edge_count();
        edge_cursor = n_edges > 0 ? boundary.global_edge_index(loop_i, edge_i) % n_edges : 0;
        return true;
    }

    /// Paper-style seed: equilateral triangle on the surface, closed boundary ring of 3 edges.
    bool build_initial_seed_triangle(const VectorX<T>& seed,
                                     std::vector<VectorX<T>>& vertices,
                                     std::vector<Triangle>& triangles,
                                     BoundaryLoops& boundary) const
    {
        VectorX<T> p0 = seed;
        if (!project_to_surface(p0))
            return false;

        VectorX<T> t1, t2;
        if (!tangent_frame_at(p0, t1, t2))
            return false;

        const T d0 = blended_step(d_step);
        const T half = static_cast<T>(-0.5);
        const T h = static_cast<T>(std::sqrt(3) / 2);
        VectorX<T> p1 = p0 + d0 * t1;
        VectorX<T> p2 = p0 + d0 * (half * t1 + h * t2);
        if (!project_to_surface(p1) || !project_to_surface(p2))
            return false;

        const int i0 = find_or_add_vertex(vertices, p0);
        const int i1 = find_or_add_vertex(vertices, p1);
        const int i2 = find_or_add_vertex(vertices, p2);
        if (!triangle_edges_within_limit(vertices, i0, i1, i2))
            return false;

        if (!register_triangle(i0, i1, i2, vertices, triangles))
            return false;

        boundary.clear();
        boundary.loops.push_back({i0, i1, i2});
        return true;
    }

    bool try_add_new_vertex_triangle(std::vector<VectorX<T>>& vertices,
                                     std::vector<Triangle>& triangles,
                                     BoundaryLoops& boundary,
                                     size_t loop_i,
                                     size_t edge_i,
                                     int vk,
                                     int vk1,
                                     int vp,
                                     size_t& edge_cursor) const
    {
        if (vp == vk || vp == vk1)
            return false;
        if (boundary.is_vertex_on_loop(vp))
            return false;
        if (!triangle_edges_within_limit(vertices, vk, vk1, vp))
            return false;
        if (!delaunay_constraint_ok(vertices, triangles, vk, vp, vk1))
            return false;
        if (!boundary.can_consume_edge_with_new_vertex(loop_i, vk, vk1))
            return false;
        if (!boundary_triangle_orientation_ok(vk, vp, vk1, vertices, triangles))
            return false;
        if (triangle_exists_adjacent(vk, vp, vk1, triangles))
            return false;

        const std::vector<std::vector<int>> backup = boundary.loops;
        if (!boundary.consume_edge_with_new_vertex(loop_i, vk, vk1, vp))
        {
            boundary.loops = backup;
            return false;
        }
        if (!register_triangle(vk, vp, vk1, vertices, triangles))
        {
            boundary.loops = backup;
            return false;
        }

        const size_t n_edges = boundary.edge_count();
        edge_cursor = n_edges > 0 ? (boundary.global_edge_index(loop_i, edge_i) + 1) % n_edges : 0;
        return true;
    }

    bool try_grow_triangle(std::vector<VectorX<T>>& vertices,
                           std::vector<Triangle>& triangles,
                           BoundaryLoops& boundary,
                           size_t& edge_cursor) const
    {
        const size_t n_edges = boundary.edge_count();
        if (n_edges == 0)
            return false;

        size_t loop_i = 0;
        size_t edge_i = 0;
        if (!boundary.edge_at(edge_cursor % n_edges, loop_i, edge_i))
            return false;

        int vk = 0;
        int vk1 = 0;
        if (!boundary.get_edge_vertices(loop_i, edge_i, vk, vk1))
            return false;
        if (vk < 0 || vk1 < 0
            || vk >= static_cast<int>(vertices.size())
            || vk1 >= static_cast<int>(vertices.size()))
            return false;

        const VectorX<T>& xk = vertices[vk];
        const VectorX<T>& xk1 = vertices[vk1];
        VectorX<T> xs = (xk + xk1) * static_cast<T>(0.5);
        if (!project_to_surface(xs))
            return false;

        VectorX<T> grad;
        compute_gradient(xs, grad);
        VectorX<T> e = xk1 - xk;
        VectorX<T> t = cross3(e, grad);
        if (t.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;
        t.normalize();

        const T e_cur = edge_length(xk, xk1);
        const int vm1 = boundary.vertex_before(loop_i, edge_i);
        const int vp2 = boundary.vertex_after(loop_i, edge_i);
        T e_prev = e_cur;
        T e_next = e_cur;
        if (vm1 >= 0 && vm1 < static_cast<int>(vertices.size()))
            e_prev = edge_length(vertices[vm1], xk);
        if (vp2 >= 0 && vp2 < static_cast<int>(vertices.size()))
            e_next = edge_length(xk1, vertices[vp2]);

        const T d = blended_step(static_cast<T>(std::sqrt(3) / 2) * (e_prev + e_cur + e_next) / 3);

        VectorX<T> xp = xs + d * t;
        if (project_to_surface(xp))
        {
            const int vp = find_or_add_vertex(vertices, xp);
            if (!boundary.is_vertex_on_loop(vp)
                && try_add_new_vertex_triangle(
                       vertices, triangles, boundary, loop_i, edge_i, vk, vk1, vp, edge_cursor))
                return true;
        }

        if (vm1 >= 0
            && try_wedge_predecessor_triangle(
                   vertices, triangles, boundary, loop_i, edge_i, vk, vk1, vm1, edge_cursor))
            return true;

        if (vp2 >= 0
            && try_wedge_successor_triangle(
                   vertices, triangles, boundary, loop_i, edge_i, vk, vk1, vp2, edge_cursor))
            return true;

        return false;
    }

    bool grow_sheet(std::vector<VectorX<T>>& vertices,
                    std::vector<Triangle>& triangles,
                    BoundaryLoops& boundary,
                    int max_triangles) const
    {
        triangles.reserve(static_cast<size_t>(max_triangles));
        int stagnation = 0;
        size_t edge_cursor = 0;

        while (static_cast<int>(triangles.size()) < max_triangles && boundary.edge_count() > 0)
        {
            const size_t n_edges = boundary.edge_count();

            if (try_grow_triangle(vertices, triangles, boundary, edge_cursor))
            {
                stagnation = 0;
            }
            else
            {
                ++edge_cursor;
                if (edge_cursor >= n_edges)
                {
                    edge_cursor = 0;
                    ++stagnation;
                    if (stagnation > static_cast<int>(n_edges))
                        break;
                }
            }
        }
        return !triangles.empty();
    }

    void mark_covered_seeds(const std::vector<VectorX<T>>& mesh_vertices,
                            const std::vector<VectorX<T>>& seeds,
                            std::vector<bool>& covered) const
    {
        const T eps2 = same_vertex_epsilon * same_vertex_epsilon;
        for (int i = 0; i < static_cast<int>(seeds.size()); ++i)
        {
            if (covered[i])
                continue;
            for (const auto& v : mesh_vertices)
            {
                if ((seeds[i] - v).squaredNorm() < eps2)
                {
                    covered[i] = true;
                    break;
                }
            }
        }
    }

    bool extract_one_sheet(const VectorX<T>& seed,
                             std::vector<VectorX<T>>& vertices,
                             std::vector<Triangle>& triangles,
                             std::vector<bool>& covered,
                             const std::vector<VectorX<T>>& all_seeds)
    {
        vertices.clear();
        triangles.clear();
        clear_triangle_index();
        BoundaryLoops boundary;

        if (!build_initial_seed_triangle(seed, vertices, triangles, boundary))
            return false;

        const int max_triangles = 1000000;
        if (!grow_sheet(vertices, triangles, boundary, max_triangles))
            return false;

        if (enable_crack_closing)
        {
            const std::vector<int> boundary_verts = boundary.largest_loop();
            if (boundary_verts.size() >= 3)
            {
                crack_closing::CrackCloser<T> closer(same_vertex_epsilon, d_max, triangles);
                const size_t tri_before = triangles.size();
                closer.close_cracks(vertices, triangles, boundary_verts, true);
                std::cout << "crack closing: +" << (triangles.size() - tri_before) << " triangles" << std::endl;
            }
        }

        compact_duplicate_triangles(triangles);

        mark_covered_seeds(vertices, all_seeds, covered);
        return true;
    }
};

} // namespace marching_triangles
