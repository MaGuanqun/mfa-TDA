#pragma once

#include <vector>
#include <deque>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iomanip>

#include <Eigen/Dense>
#include <mfa/mfa.hpp>

#include "query_function.h"
#include "utility_function.h"
#include "block.hpp"

namespace marching_triangles
{

template<typename T>
struct Triangle
{
    int v0, v1, v2;
};

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
        bool enable_crack_closing_ = true)
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
        , projection_alpha(projection_alpha_)
        , b(b_)
        , inr_model(inr_model_)
        , enable_crack_closing(enable_crack_closing_)
    {}

    void set_degenerate_points(const std::vector<VectorX<T>>& pts) { degenerate_points = pts; }

    bool extract_all_sheets(const std::vector<VectorX<T>>& seed_roots,
    std::vector<std::vector<VectorX<T>>>& sheet_vertices,
    std::vector<std::vector<Triangle<T>>>& sheet_triangles)
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
            std::vector<Triangle<T>> tris;
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
                              const std::vector<Triangle<T>>& triangles)
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
            out << "f " << (tri.v0 + 1) << " " << (tri.v1 + 1) << " " << (tri.v2 + 1) << "\n";
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
    const T projection_alpha;
    const Block<T>* b;
    INRModel<T>* inr_model;
    const bool enable_crack_closing;

    std::vector<VectorX<T>> degenerate_points;

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

    // 3D cross product (points are always 3D; VectorX is used for API compatibility only).

    static VectorX<T> cross3(const VectorX<T>& u, const Vector3<T>& v)
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
                                const std::vector<Triangle<T>>& triangles,
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

        // here need to be fastened
        for (const auto& tri : triangles)
        {
            const int ids[3] = {tri.v0, tri.v1, tri.v2};
            for (int vid : ids)
            {
                if (vid == v0 || vid == v1 || vid == vp)
                    continue;
                if ((vertices[vid] - center).squaredNorm() < r2)
                {
                    const Vector3<T> n_old = triangle_normal(vertices[tri.v0], vertices[tri.v1], vertices[tri.v2]);
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
        if (d >= d_min)
            return d;
        return static_cast<T>(0.75) * d + static_cast<T>(0.25) * d_min;
    }

    bool build_initial_edge(const VectorX<T>& seed,
                            std::vector<VectorX<T>>& vertices,
                            std::deque<int>& boundary) const
    {
        VectorX<T> p0 = seed;
        if (!project_to_surface(p0))
            return false;

        VectorX<T> grad;
        compute_gradient(p0, grad);
        if (grad.norm() < std::numeric_limits<T>::epsilon())
            return false;

        VectorX<T> axis_direction = VectorX<T>::Zero(grad.size());
        axis_direction[0] = 1;

        VectorX<T> tangent = cross3(grad, axis_direction);
        if (tangent.squaredNorm() < std::numeric_limits<T>::epsilon())
        {
            axis_direction[0] = 0;
            axis_direction[1] = 1;
            tangent = cross3(grad, axis_direction);
        }
        tangent.normalize();

        const T d0 = blended_step(d_step);
        VectorX<T> p1 = p0 + d0 * tangent;
        if (!project_to_surface(p1))
            return false;

        const int i0 = find_or_add_vertex(vertices, p0);
        const int i1 = find_or_add_vertex(vertices, p1);
        boundary.clear();
        boundary.push_back(i0);
        boundary.push_back(i1);
        return true;
    }

    bool try_grow_triangle(std::vector<VectorX<T>>& vertices,
                           std::vector<Triangle<T>>& triangles,
                           std::deque<int>& boundary,
                           size_t& edge_index) const
    {
        if (boundary.size() < 2)
            return false;

        const size_t n = boundary.size();
        const size_t n_edges = n - 1;
        size_t k = edge_index % n_edges;
        const int vk = boundary[k];
        const int vk1 = boundary[k + 1];

        const VectorX<T>& xk = vertices[vk];
        const VectorX<T>& xk1 = vertices[vk1];
        const VectorX<T> xm = (xk + xk1) * static_cast<T>(0.5);

        VectorX<T> xs = xm;
        if (!project_to_surface(xs))
            return false;

        VectorX<T> grad;
        compute_gradient(xs, grad);
        VectorX<T> temp=xk1 - xk;
        VectorX<T> t = cross3(temp, grad);
        if (t.squaredNorm() < std::numeric_limits<T>::epsilon())
            return false;
        t.normalize();

        const T e_cur = edge_length(xk, xk1);
        T e_prev = e_cur;
        T e_next = e_cur;
        if (k > 0)
            e_prev = edge_length(vertices[boundary[k - 1]], xk);
        if (k + 2 < n)
            e_next = edge_length(xk1, vertices[boundary[k + 2]]);

        T d_bar = (e_prev + e_cur + e_next) / 3;
        T d = blended_step(static_cast<T>(std::sqrt(3) / 2) * d_bar);

        VectorX<T> xp = xs + d * t;
        if (!project_to_surface(xp))
            return false;

        // if (near_stop_set(xp) || near_stop_set(xs))
        //     return false;

        const int vp = find_or_add_vertex(vertices, xp);
        if (delaunay_constraint_ok(vertices, triangles, vk, vk1, vp))
        {
            triangles.push_back({vk, vk1, vp});
            boundary.insert(boundary.begin() + static_cast<std::deque<int>::difference_type>(k + 1), vp);
            edge_index = (k + 1) % (boundary.size() - 1);
            return true;
        }

        // Paper §3.2 step 3.2: try triangles using only existing boundary vertices.
        if (k > 0)
        {
            const int vm1 = boundary[k - 1];
            if (delaunay_constraint_ok(vertices, triangles, vm1, vk, vk1))
            {
                triangles.push_back({vm1, vk, vk1});
                boundary.erase(boundary.begin() + static_cast<std::deque<int>::difference_type>(k));
                edge_index = k % std::max<size_t>(1, boundary.size() - 1);
                return true;
            }
        }
        if (k + 2 < n)
        {
            const int vp2 = boundary[k + 2];
            if (delaunay_constraint_ok(vertices, triangles, vk, vk1, vp2))
            {
                triangles.push_back({vk, vk1, vp2});
                boundary.erase(boundary.begin() + static_cast<std::deque<int>::difference_type>(k + 1));
                edge_index = k % std::max<size_t>(1, boundary.size() - 1);
                return true;
            }
        }

        return false;
    }

    bool grow_sheet(std::vector<VectorX<T>>& vertices,
                    std::vector<Triangle<T>>& triangles,
                    std::deque<int>& boundary,
                    int max_triangles) const
    {
        int stagnation = 0;
        size_t edge_cursor = 0;

        while (static_cast<int>(triangles.size()) < max_triangles && boundary.size() >= 2)
        {
            bool stop_boundary = false;
            // for (int vid : boundary)
            // {
            //     if (near_stop_set(vertices[vid]))
            //     {
            //         stop_boundary = true;
            //         break;
            //     }
            // }
            // if (stop_boundary)
            //     break;

            const size_t n_edges = boundary.size() - 1;
            if (n_edges == 0)
                break;

            if (try_grow_triangle(vertices, triangles, boundary, edge_cursor))
            {
                stagnation = 0;
                const size_t new_n_edges = boundary.size() - 1;
                if (new_n_edges > 0)
                    edge_cursor = (edge_cursor + 1) % new_n_edges;
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

            if(triangles.size() == 2000)
            {
                std::cout<<"triangles.size() "<<triangles.size()<<std::endl;
                return true;
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
                             std::vector<Triangle<T>>& triangles,
                             std::vector<bool>& covered,
                             const std::vector<VectorX<T>>& all_seeds)
    {
        vertices.clear();
        triangles.clear();
        std::deque<int> boundary;

        if (!build_initial_edge(seed, vertices, boundary))
            return false;

        const int max_triangles = 1000000;
        if (!grow_sheet(vertices, triangles, boundary, max_triangles))
            return false;

        std::cout<<"enable_crack_closing "<<enable_crack_closing<<std::endl;

        if (enable_crack_closing && boundary.size() >= 3)
        {
            std::vector<int> boundary_verts(boundary.begin(), boundary.end());
            const T gap = (vertices[boundary_verts.front()] - vertices[boundary_verts.back()]).norm();
            const bool boundary_closed = gap < same_vertex_epsilon * 10;

            crack_closing::CrackCloser<T> closer(same_vertex_epsilon);
            const size_t tri_before = triangles.size();
            closer.close_cracks(vertices, triangles, boundary_verts, boundary_closed);
            std::cout << "crack closing: +" << (triangles.size() - tri_before) << " triangles" << std::endl;
        }

        mark_covered_seeds(vertices, all_seeds, covered);
        return true;
    }
};

} // namespace marching_triangles
