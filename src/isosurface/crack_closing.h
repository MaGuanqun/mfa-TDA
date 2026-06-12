#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <cstdint>

#include <Eigen/Dense>

namespace crack_closing
{

using Triangle = marching_triangles::Triangle;

/// Ordered vertex ring; consecutive vertices form boundary edges (closed if is_closed).
template<typename T>
struct Contour
{
    std::vector<int> verts;
    bool is_closed = true;
};

template<typename T>
class CrackCloser
{
public:
    static constexpr int kSmallContourMaxEdges = 6;

    CrackCloser(T same_vertex_epsilon_, T max_edge_length_, const std::vector<Triangle>& existing_triangles)
        : same_vertex_epsilon(same_vertex_epsilon_)
        , max_edge_length(max_edge_length_)
        , max_edge_length_sq(max_edge_length_ * max_edge_length_)
    {
        marching_triangles::seed_triangle_keys(existing_triangles, triangle_keys_);
    }

    /// Original Akkouche crack-fixing: contour splitting + contour meshing on boundary Le.
    void close_cracks(const std::vector<Eigen::VectorX<T>>& vertices,
                      std::vector<Triangle>& triangles,
                      const std::vector<int>& boundary_verts,
                      bool boundary_is_closed) const
    {
        if (boundary_verts.size() < 3)
            return;

        Contour<T> contour;
        contour.verts = boundary_verts;
        contour.is_closed = boundary_is_closed;

        if (!contour.is_closed)
        {
            const T d = (vertices[contour.verts.front()] - vertices[contour.verts.back()]).norm();
            if (d < same_vertex_epsilon * 10 && contour.verts.size() >= 3)
                contour.is_closed = true;
        }

        std::vector<Contour<T>> leaves;
        leaves.push_back(contour);
        split_all_contours(vertices, triangles, leaves);
        contour_meshing(vertices, triangles, leaves);
    }

private:
    T same_vertex_epsilon;
    T max_edge_length;
    T max_edge_length_sq;
    mutable std::unordered_set<uint64_t> triangle_keys_;

    bool triangle_edges_within_limit(const std::vector<Eigen::VectorX<T>>& vertices,
                                     int a, int b, int c) const
    {
        return (vertices[a] - vertices[b]).squaredNorm() <= max_edge_length_sq
            && (vertices[b] - vertices[c]).squaredNorm() <= max_edge_length_sq
            && (vertices[c] - vertices[a]).squaredNorm() <= max_edge_length_sq;
    }

    int fan_apex_index(const std::vector<Eigen::VectorX<T>>& vertices, const Contour<T>& contour) const
    {
        Eigen::Matrix<T, 3, 1> centroid = Eigen::Matrix<T, 3, 1>::Zero();
        for (int vid : contour.verts)
            centroid += to3(vertices[vid]);
        centroid /= static_cast<T>(contour.verts.size());

        int best = contour.verts[0];
        T best_d2 = std::numeric_limits<T>::max();
        for (int vid : contour.verts)
        {
            const T d2 = (to3(vertices[vid]) - centroid).squaredNorm();
            if (d2 < best_d2)
            {
                best_d2 = d2;
                best = vid;
            }
        }
        return best;
    }

    int num_edges(const Contour<T>& c) const
    {
        const int m = static_cast<int>(c.verts.size());
        return c.is_closed ? m : std::max(0, m - 1);
    }

    int vert_at(const Contour<T>& c, int k) const
    {
        const int m = static_cast<int>(c.verts.size());
        if (c.is_closed)
            return c.verts[(k % m + m) % m];
        return c.verts[k];
    }

    int next_index(const Contour<T>& c, int k) const
    {
        return c.is_closed ? k + 1 : k + 1;
    }

    Eigen::Matrix<T, 3, 1> to3(const Eigen::VectorX<T>& p) const
    {
        Eigen::Matrix<T, 3, 1> q = Eigen::Matrix<T, 3, 1>::Zero();
        for (int i = 0; i < std::min(3, static_cast<int>(p.size())); ++i)
            q[i] = p[i];
        return q;
    }

    T dist_point_segment(const Eigen::Matrix<T, 3, 1>& p,
                         const Eigen::Matrix<T, 3, 1>& a,
                         const Eigen::Matrix<T, 3, 1>& b) const
    {
        const Eigen::Matrix<T, 3, 1> ab = b - a;
        const T ab2 = ab.squaredNorm();
        if (ab2 < std::numeric_limits<T>::epsilon())
            return (p - a).norm();
        T t = (p - a).dot(ab) / ab2;
        t = std::max(static_cast<T>(0), std::min(static_cast<T>(1), t));
        return (p - (a + t * ab)).norm();
    }

    /// Prefer the growth triangle (third vertex not on the crack contour).
    bool find_adjacent_triangle(const std::vector<Triangle>& triangles,
                                int a,
                                int b,
                                const std::vector<int>& contour_vertex_set,
                                int& tri_index,
                                bool& edge_is_ab) const
    {
        int fallback_idx = -1;
        bool fallback_ab = true;

        for (int t = 0; t < static_cast<int>(triangles.size()); ++t)
        {
            const auto& tri = triangles[t];
            for (int i = 0; i < 3; ++i)
            {
                const int v0 = tri[i];
                const int v1 = tri[(i + 1) % 3];
                const int v2 = tri[(i + 2) % 3];
                bool is_ab = false;
                bool on_edge = false;
                if (v0 == a && v1 == b)
                {
                    is_ab = true;
                    on_edge = true;
                }
                else if (v0 == b && v1 == a)
                {
                    is_ab = false;
                    on_edge = true;
                }
                if (!on_edge)
                    continue;

                const bool interior_w = contour_vertex_set.empty() ||
                    (v2 >= static_cast<int>(contour_vertex_set.size())) ||
                    (contour_vertex_set[static_cast<size_t>(v2)] == 0);
                if (interior_w)
                {
                    tri_index = t;
                    edge_is_ab = is_ab;
                    return true;
                }
                if (fallback_idx < 0)
                {
                    fallback_idx = t;
                    fallback_ab = is_ab;
                }
            }
        }
        if (fallback_idx >= 0)
        {
            tri_index = fallback_idx;
            edge_is_ab = fallback_ab;
            return true;
        }
        return false;
    }

    static std::vector<int> build_contour_vertex_set(const Contour<T>& contour)
    {
        int max_vid = -1;
        for (int vid : contour.verts)
            max_vid = std::max(max_vid, vid);
        if (max_vid < 0)
            return {};
        std::vector<int> on_contour(static_cast<size_t>(max_vid) + 1, 0);
        for (int vid : contour.verts)
            on_contour[static_cast<size_t>(vid)] = 1;
        return on_contour;
    }

    Eigen::Matrix<T, 3, 1> triangle_normal_at(const std::vector<Eigen::VectorX<T>>& vertices,
                                              const Triangle& tri) const
    {
        return (to3(vertices[tri[1]]) - to3(vertices[tri[0]]))
            .cross(to3(vertices[tri[2]]) - to3(vertices[tri[0]]));
    }

    bool faces_open_edge(const std::vector<Eigen::VectorX<T>>& vertices,
                         const std::vector<Triangle>& triangles,
                         const Contour<T>& contour,
                         const std::vector<int>& contour_vertex_set,
                         int edge_k,
                         int candidate_vert) const
    {
        const int vk = vert_at(contour, edge_k);
        const int vk1 = vert_at(contour, next_index(contour, edge_k));
        const Eigen::Matrix<T, 3, 1> a = to3(vertices[vk]);
        const Eigen::Matrix<T, 3, 1> b = to3(vertices[vk1]);
        const Eigen::Matrix<T, 3, 1> mid = (a + b) * static_cast<T>(0.5);
        const Eigen::Matrix<T, 3, 1> to_xs = to3(vertices[candidate_vert]) - mid;

        int tri_idx = -1;
        bool ab = true;
        if (!find_adjacent_triangle(triangles, vk, vk1, contour_vertex_set, tri_idx, ab))
            return to_xs.norm() > std::numeric_limits<T>::epsilon();

        const Eigen::Matrix<T, 3, 1> n = triangle_normal_at(vertices, triangles[tri_idx]);
        Eigen::Matrix<T, 3, 1> e = b - a;
        if (!ab)
            e = -e;

        Eigen::Matrix<T, 3, 1> open_dir = e.cross(n);
        if (open_dir.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;
        open_dir.normalize();
        return open_dir.dot(to_xs) > 0;
    }

    bool is_forbidden_split_vertex(const Contour<T>& contour, int edge_k, int vert_index) const
    {
        const int m = static_cast<int>(contour.verts.size());
        if (vert_index == edge_k || vert_index == next_index(contour, edge_k))
            return true;

        if (contour.is_closed)
        {
            const int km1 = (edge_k - 1 + m) % m;
            const int kp2 = (edge_k + 2) % m;
            return vert_index == km1 || vert_index == kp2;
        }

        if (edge_k > 0 && vert_index == edge_k - 1)
            return true;
        if (edge_k + 2 < m && vert_index == edge_k + 2)
            return true;
        return false;
    }

    int find_splitting_vertex(const std::vector<Eigen::VectorX<T>>& vertices,
                              const std::vector<Triangle>& triangles,
                              const Contour<T>& contour,
                              const std::vector<int>& contour_vertex_set,
                              int edge_k) const
    {
        const int m = static_cast<int>(contour.verts.size());
        const int vk = vert_at(contour, edge_k);
        const int vk1 = vert_at(contour, next_index(contour, edge_k));
        const Eigen::Matrix<T, 3, 1> a = to3(vertices[vk]);
        const Eigen::Matrix<T, 3, 1> b = to3(vertices[vk1]);

        int best = -1;
        T best_dist = std::numeric_limits<T>::max();

        for (int j = 0; j < m; ++j)
        {
            if (is_forbidden_split_vertex(contour, edge_k, j))
                continue;

            const int xs = contour.verts[j];
            if (xs == vk || xs == vk1)
                continue;

            if (!faces_open_edge(vertices, triangles, contour, contour_vertex_set, edge_k, xs))
                continue;

            const T d = dist_point_segment(to3(vertices[xs]), a, b);
            if (d > max_edge_length)
                continue;
            if (!triangle_edges_within_limit(vertices, vk, vk1, xs))
                continue;
            if (d < best_dist)
            {
                best_dist = d;
                best = j;
            }
        }
        return best;
    }

    void extract_subcontour(const Contour<T>& contour, int start, int end, Contour<T>& out) const
    {
        out.verts.clear();
        const int m = static_cast<int>(contour.verts.size());
        if (m == 0)
            return;

        int i = start;
        out.verts.push_back(contour.verts[i]);
        while (i != end)
        {
            if (!contour.is_closed && i + 1 >= m)
                break;
            i = contour.is_closed ? (i + 1) % m : i + 1;
            out.verts.push_back(contour.verts[i]);
            if (!contour.is_closed && i >= m - 1)
                break;
        }
        out.is_closed = contour.is_closed && out.verts.size() >= 3;
    }

    bool split_contour(const Contour<T>& contour, int edge_k, int split_index,
                       Contour<T>& sub_a, Contour<T>& sub_b) const
    {
        const int m = static_cast<int>(contour.verts.size());
        const int kp1 = contour.is_closed ? (edge_k + 1) % m : edge_k + 1;

        extract_subcontour(contour, kp1, split_index, sub_a);
        extract_subcontour(contour, split_index, edge_k, sub_b);

        return sub_a.verts.size() >= 3 && sub_b.verts.size() >= 3;
    }

    bool triangle_orientation_ok(const std::vector<Eigen::VectorX<T>>& vertices,
                                 const std::vector<Triangle>& triangles,
                                 const std::vector<int>& contour_vertex_set,
                                 int v0,
                                 int v1,
                                 int v2) const
    {
        const Eigen::Matrix<T, 3, 1> n_new =
            (to3(vertices[v1]) - to3(vertices[v0])).cross(to3(vertices[v2]) - to3(vertices[v0]));
        if (n_new.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;

        int tri_idx = -1;
        bool ab = true;
        if (find_adjacent_triangle(triangles, v0, v1, contour_vertex_set, tri_idx, ab))
        {
            const Eigen::Matrix<T, 3, 1> n_old = triangle_normal_at(vertices, triangles[tri_idx]);
            return n_new.dot(n_old) > 0;
        }
        return true;
    }

    void fan_triangulate(const std::vector<Eigen::VectorX<T>>& vertices,
                         std::vector<Triangle>& triangles,
                         const Contour<T>& contour,
                         const std::vector<int>& contour_vertex_set) const
    {
        const int m = static_cast<int>(contour.verts.size());
        if (m < 3)
            return;

        const int apex = fan_apex_index(vertices, contour);
        for (int i = 0; i < m; ++i)
        {
            if (contour.verts[i] == apex)
                continue;
            const int v1 = contour.verts[i];
            const int v2 = contour.verts[(i + 1) % m];
            if (v1 == apex || v2 == apex)
                continue;
            if (!triangle_edges_within_limit(vertices, apex, v1, v2))
                continue;
            if (!triangle_orientation_ok(vertices, triangles, contour_vertex_set, apex, v1, v2))
                continue;
            marching_triangles::add_triangle_unique(triangles, triangle_keys_, apex, v1, v2);
        }
    }

    /// Recursive divide-and-conquer split (paper §3.2.2). Leaf contours (|L_e| < 6 or unsplittable) go to out_leaves.
    void contour_splitting(const std::vector<Eigen::VectorX<T>>& vertices,
                           std::vector<Triangle>& triangles,
                           const Contour<T>& contour,
                           std::vector<Contour<T>>& out_leaves,
                           int edge_start = 0) const
    {
        const int m = static_cast<int>(contour.verts.size());
        const int n_edges = num_edges(contour);
        if (m < 3 || n_edges == 0)
            return;

        if (n_edges < kSmallContourMaxEdges)
        {
            out_leaves.push_back(contour);
            return;
        }

        const std::vector<int> contour_vertex_set = build_contour_vertex_set(contour);

        for (int t = 0; t < n_edges; ++t)
        {
            const int k = contour.is_closed ? (edge_start + t) % n_edges : edge_start + t;
            if (!contour.is_closed && k >= n_edges)
                break;

            const int split_idx = find_splitting_vertex(vertices, triangles, contour, contour_vertex_set, k);
            if (split_idx < 0)
                continue;

            const int vk = vert_at(contour, k);
            const int vk1 = vert_at(contour, next_index(contour, k));
            const int vs = contour.verts[split_idx];

            if (!triangle_orientation_ok(vertices, triangles, contour_vertex_set, vk, vk1, vs))
                continue;

            Contour<T> sub_a, sub_b;
            if (!split_contour(contour, k, split_idx, sub_a, sub_b))
                continue;

            marching_triangles::add_triangle_unique(triangles, triangle_keys_, vk, vk1, vs);

            contour_splitting(vertices, triangles, sub_a, out_leaves, 0);
            contour_splitting(vertices, triangles, sub_b, out_leaves, 0);
            return;
        }

        out_leaves.push_back(contour);
    }

    void split_all_contours(const std::vector<Eigen::VectorX<T>>& vertices,
                            std::vector<Triangle>& triangles,
                            std::vector<Contour<T>>& contours) const
    {
        std::vector<Contour<T>> leaves;
        for (const auto& c : contours)
            contour_splitting(vertices, triangles, c, leaves, 0);
        contours = std::move(leaves);
    }

    T contour_pair_distance(const std::vector<Eigen::VectorX<T>>& vertices,
                            const Contour<T>& c1,
                            const Contour<T>& c2,
                            int& best_i,
                            int& best_j) const
    {
        T best = std::numeric_limits<T>::max();
        best_i = best_j = -1;
        for (int i = 0; i < static_cast<int>(c1.verts.size()); ++i)
        {
            for (int j = 0; j < static_cast<int>(c2.verts.size()); ++j)
            {
                const T d = (vertices[c1.verts[i]] - vertices[c2.verts[j]]).norm();
                if (d < best)
                {
                    best = d;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        return best;
    }

    Contour<T> merge_contours(const Contour<T>& c1, int i, const Contour<T>& c2, int j) const
    {
        Contour<T> merged;
        merged.is_closed = true;
        const int n1 = static_cast<int>(c1.verts.size());
        const int n2 = static_cast<int>(c2.verts.size());

        for (int t = 0; t < n1; ++t)
            merged.verts.push_back(c1.verts[(i + t) % n1]);
        if (n2 > 2)
        {
            for (int t = 1; t < n2 - 1; ++t)
                merged.verts.push_back(c2.verts[(j + t) % n2]);
        }
        else if (n2 == 2)
        {
            merged.verts.push_back(c2.verts[(j + 1) % n2]);
        }

        return merged;
    }

    bool try_merge_pair(const std::vector<Eigen::VectorX<T>>& vertices,
                        std::vector<Triangle>& triangles,
                        const Contour<T>& c1,
                        const Contour<T>& c2,
                        int i,
                        int j,
                        Contour<T>& merged_out) const
    {
        const int n1 = static_cast<int>(c1.verts.size());
        const int n2 = static_cast<int>(c2.verts.size());
        if (n1 < 2 || n2 < 2)
            return false;

        std::vector<int> merge_vertex_set = build_contour_vertex_set(c1);
        for (int vid : c2.verts)
        {
            if (vid >= static_cast<int>(merge_vertex_set.size()))
                merge_vertex_set.resize(static_cast<size_t>(vid) + 1, 0);
            merge_vertex_set[static_cast<size_t>(vid)] = 1;
        }

        const int im1 = (i - 1 + n1) % n1;
        const int jp1 = (j + 1) % n2;
        const int a = c1.verts[im1];
        const int b = c1.verts[i];
        const int c = c2.verts[j];

        if (!triangle_orientation_ok(vertices, triangles, merge_vertex_set, a, b, c))
        {
            const int a2 = c1.verts[i];
            const int b2 = c2.verts[j];
            const int c2v = c2.verts[jp1];
            if (!triangle_orientation_ok(vertices, triangles, merge_vertex_set, a2, b2, c2v))
                return false;
            marching_triangles::add_triangle_unique(triangles, triangle_keys_, a2, b2, c2v);
            merged_out = merge_contours(c2, j, c1, i);
            return true;
        }

        marching_triangles::add_triangle_unique(triangles, triangle_keys_, a, b, c);
        merged_out = merge_contours(c1, i, c2, j);
        return true;
    }

    /// Mesh small holes, then merge remaining contour pairs with fixing triangles.
    void contour_meshing(const std::vector<Eigen::VectorX<T>>& vertices,
                         std::vector<Triangle>& triangles,
                         std::vector<Contour<T>>& contours) const
    {
        // Phase 1: mesh simple holes (six or fewer boundary edges, paper §3.2.2).
        bool progress = true;
        while (progress)
        {
            progress = false;
            for (auto it = contours.begin(); it != contours.end();)
            {
                if (num_edges(*it) <= kSmallContourMaxEdges)
                {
                    if (it->is_closed)
                    {
                        const std::vector<int> vset = build_contour_vertex_set(*it);
                        fan_triangulate(vertices, triangles, *it, vset);
                    }
                    it = contours.erase(it);
                    progress = true;
                }
                else
                {
                    ++it;
                }
            }
        }

        // Phase 2: merge facing contour pairs with a fixing triangle (tubular sections).
        while (contours.size() > 1)
        {
            T best_d = std::numeric_limits<T>::max();
            int best_a = -1, best_b = -1, best_i = -1, best_j = -1;

            for (int a = 0; a < static_cast<int>(contours.size()); ++a)
            {
                for (int b = a + 1; b < static_cast<int>(contours.size()); ++b)
                {
                    int i = -1, j = -1;
                    const T d = contour_pair_distance(vertices, contours[a], contours[b], i, j);
                    if (d < best_d)
                    {
                        best_d = d;
                        best_a = a;
                        best_b = b;
                        best_i = i;
                        best_j = j;
                    }
                }
            }

            if (best_a < 0)
                break;

            Contour<T> merged;
            if (!try_merge_pair(vertices, triangles, contours[best_a], contours[best_b],
                                best_i, best_j, merged))
                break;

            const int hi = std::max(best_a, best_b);
            const int lo = std::min(best_a, best_b);
            contours.erase(contours.begin() + hi);
            contours.erase(contours.begin() + lo);
            contours.push_back(std::move(merged));
        }

        // Phase 3: if one contour remains, close it as a simple hole (paper §3.2.2).
        for (const auto& c : contours)
        {
            if (c.is_closed)
            {
                const std::vector<int> vset = build_contour_vertex_set(c);
                fan_triangulate(vertices, triangles, c, vset);
            }
        }
    }
};

} // namespace crack_closing
