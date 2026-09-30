#pragma once

#include "mesh_validation.h"
#include <numeric>

namespace marching_triangles {

// Scale-invariant mean-ratio quality: 1 for an equilateral triangle, 0 at
// degeneracy. A nonzero area is a validity test, not a useful quality threshold.
template<class T>
T triangle_quality(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
    const T sum = (b-a).squaredNorm() + (c-b).squaredNorm() + (a-c).squaredNorm();
    if (!(sum > 0)) return 0;
    return T(2)*std::sqrt(T(3))*(b-a).cross(c-a).norm()/sum;
}

template<class T>
T triangle_normal_alignment(const SurfaceField<T>& field, const Point<T>& a,
                            const Point<T>& b, const Point<T>& c) {
    const Point<T> cross = (b-a).cross(c-a);
    const T length = cross.stableNorm();
    if (!(length > 0)) return -1;
    const Point<T> normal = cross/length;
    T alignment = 1;
    for (const Point<T>& p : std::array<Point<T>,4>{a,b,c,(a+b+c)/3}) {
        const Point<T> gradient = field.gradient(p);
        const T norm = gradient.stableNorm();
        if (!gradient.allFinite() || !(norm > 0)) return -1;
        alignment = std::min(alignment, normal.dot(gradient/norm));
    }
    return alignment;
}

struct MeshQuality {
    double min_angle_degrees = 180;
    double min_quality = 1;
    double mean_quality = 0;
    double max_normal_error_degrees = 0;
    size_t triangles_below_5_degrees = 0;
    size_t reversed_faces = 0;

    std::string summary() const {
        std::ostringstream out;
        out << "min_angle_deg=" << min_angle_degrees
            << " triangles_below_5deg=" << triangles_below_5_degrees
            << " min_quality=" << min_quality << " mean_quality=" << mean_quality
            << " max_normal_error_deg=" << max_normal_error_degrees
            << " reversed_faces=" << reversed_faces;
        return out.str();
    }
};

template<class T>
MeshQuality measure_mesh_quality(const std::vector<Point<T>>& vertices,
                                 const std::vector<Triangle>& triangles,
                                 const SurfaceField<T>& field) {
    MeshQuality result;
    const double to_degrees = 180/std::acos(-1.0);
    for (const auto& t : triangles) {
        const Point<T>& a = vertices[t[0]];
        const Point<T>& b = vertices[t[1]];
        const Point<T>& c = vertices[t[2]];
        const double q = triangle_quality(a,b,c);
        result.min_quality = std::min(result.min_quality,q);
        result.mean_quality += q;
        double angle = 180;
        for (int i=0; i<3; ++i) {
            const Point<T> u = vertices[t[(i+1)%3]]-vertices[t[i]];
            const Point<T> v = vertices[t[(i+2)%3]]-vertices[t[i]];
            angle = std::min(angle,std::atan2(double(u.cross(v).norm()),double(u.dot(v)))*to_degrees);
        }
        result.min_angle_degrees = std::min(result.min_angle_degrees,angle);
        if (angle < 5) ++result.triangles_below_5_degrees;
        const T cosine = triangle_normal_alignment(field,a,b,c);
        if (cosine <= 0) ++result.reversed_faces;
        result.max_normal_error_degrees = std::max(result.max_normal_error_degrees,
            std::acos(std::clamp(double(cosine),-1.0,1.0))*to_degrees);
    }
    if (!triangles.empty()) result.mean_quality /= triangles.size();
    return result;
}

// Optional refinement after a sheet is closed. It changes diagonals and moves
// existing vertices along the implicit surface; it neither caps holes nor
// changes connectivity between components. This is an addition to the paper's
// two meshing stages, whose crack triangulation alone has no quality bound.
template<class T>
class MeshQualityOptimizer {
public:
    MeshQualityOptimizer(const SurfaceField<T>& field, const MeshingOptions<T>& options,
                         const std::vector<Point<T>>& stop_points)
        : field_(field), options_(options), projector_(field,options), stop_points_(stop_points) {}

    void improve(TriangleMesh<T>& mesh) {
        if (options_.quality_iterations == 0 || !mesh.boundary.empty() || mesh.triangles.empty()) return;
        if (!validate_mesh(mesh.vertices,mesh.triangles,epsilon(),false).valid()) return;
        incident_.assign(mesh.vertices.size(),{});
        boxes_.resize(mesh.triangles.size());
        for (size_t i=0; i<mesh.triangles.size(); ++i) {
            for (int v : mesh.triangles[i]) incident_[v].insert(static_cast<int>(i));
            boxes_[i] = bounds(mesh,mesh.triangles[i]);
        }
        for (int pass=0; pass<options_.quality_iterations; ++pass) {
            const size_t flips = flip_edges(mesh);
            const size_t moves = relocate_vertices(mesh);
            if (flips == 0 && moves == 0) break;
        }
        flip_edges(mesh);
    }

private:
    struct Bounds { Point<T> low, high; };
    const SurfaceField<T>& field_;
    const MeshingOptions<T>& options_;
    SurfaceProjector<T> projector_;
    const std::vector<Point<T>>& stop_points_;
    std::vector<std::set<int>> incident_;
    std::vector<Bounds> boxes_;

    T epsilon() const { return options_.step*T(1e-9); }
    Bounds bounds(const TriangleMesh<T>& mesh,const Triangle& t) const {
        Bounds box{mesh.vertices[t[0]],mesh.vertices[t[0]]};
        for (int v : t) {
            box.low = box.low.cwiseMin(mesh.vertices[v]);
            box.high = box.high.cwiseMax(mesh.vertices[v]);
        }
        return box;
    }
    bool overlaps(const Bounds& a,const Bounds& b) const {
        return !(a.low.array()>b.high.array()+epsilon()).any() &&
               !(b.low.array()>a.high.array()+epsilon()).any();
    }
    T score(const TriangleMesh<T>& mesh,const Triangle& t) const {
        const auto& a=mesh.vertices[t[0]];
        const auto& b=mesh.vertices[t[1]];
        const auto& c=mesh.vertices[t[2]];
        return std::min(triangle_quality(a,b,c),triangle_normal_alignment(field_,a,b,c));
    }
    bool surface_ok(const TriangleMesh<T>& mesh,const Triangle& t) const {
        const Point<T>& a=mesh.vertices[t[0]];
        const Point<T>& b=mesh.vertices[t[1]];
        const Point<T>& c=mesh.vertices[t[2]];
        const T h=std::max({(a-b).norm(),(b-c).norm(),(c-a).norm()});
        if (h>2*options_.step || score(mesh,t)<=T(1e-8)) return false;
        if (std::min({(a-b).norm(),(b-c).norm(),(c-a).norm()})<=options_.vertex_tolerance) return false;
        for (Point<T> p : std::array<Point<T>,4>{(a+b)/2,(b+c)/2,(c+a)/2,(a+b+c)/3}) {
            for (const auto& stop : stop_points_)
                if ((p-stop).norm()<options_.stop_distance) return false;
            if (!projector_.project(p,h*T(.35))) return false;
        }
        return true;
    }
    bool patch_ok(const TriangleMesh<T>& mesh,const std::vector<Triangle>& patch,
                  const std::vector<int>& replaced) const {
        std::vector<Bounds> patch_boxes;
        for (const auto& t : patch) {
            if (!surface_ok(mesh,t)) return false;
            patch_boxes.push_back(bounds(mesh,t));
        }
        for (size_t i=0; i<patch.size(); ++i)
            for (size_t j=i+1; j<patch.size(); ++j)
                if (triangles_conflict(mesh.vertices,patch[i],patch[j],epsilon())) return false;
        Bounds all=patch_boxes.front();
        for (const auto& box : patch_boxes) {
            all.low=all.low.cwiseMin(box.low);
            all.high=all.high.cwiseMax(box.high);
        }
        for (size_t i=0; i<mesh.triangles.size(); ++i) {
            if (!overlaps(all,boxes_[i]) ||
                std::find(replaced.begin(),replaced.end(),static_cast<int>(i))!=replaced.end()) continue;
            for (size_t j=0; j<patch.size(); ++j)
                if (overlaps(patch_boxes[j],boxes_[i]) &&
                    triangles_conflict(mesh.vertices,patch[j],mesh.triangles[i],epsilon())) return false;
        }
        return true;
    }
    size_t flip_edges(TriangleMesh<T>& mesh) {
        std::vector<Edge> edges;
        for (const auto& e : mesh.edges) edges.push_back(e.first);
        size_t changed=0;
        for (const auto& edge : edges) {
            auto found=mesh.edges.find(edge);
            if (found==mesh.edges.end() || found->second.second<0) continue;
            const auto use=found->second;
            const int a=use.from,b=use.to,f=use.face,g=use.second;
            int c=-1,d=-1;
            for (int v : mesh.triangles[f]) if (v!=a && v!=b) c=v;
            for (int v : mesh.triangles[g]) if (v!=a && v!=b) d=v;
            if (c<0 || d<0 || c==d || mesh.edges.count(edge_key(c,d))) continue;
            const Triangle first{c,d,b},second{d,c,a};
            const T before=std::min(score(mesh,mesh.triangles[f]),score(mesh,mesh.triangles[g]));
            const T after=std::min(score(mesh,first),score(mesh,second));
            if (after<=before+T(1e-4) || !patch_ok(mesh,{first,second},{f,g})) continue;
            for (int v : mesh.triangles[f]) incident_[v].erase(f);
            for (int v : mesh.triangles[g]) incident_[v].erase(g);
            mesh.replace_pair(f,g,first,second);
            for (int v : first) incident_[v].insert(f);
            for (int v : second) incident_[v].insert(g);
            boxes_[f]=bounds(mesh,first);
            boxes_[g]=bounds(mesh,second);
            ++changed;
        }
        return changed;
    }
    size_t relocate_vertices(TriangleMesh<T>& mesh) {
        std::vector<std::pair<T,int>> order;
        for (size_t v=0; v<mesh.vertices.size(); ++v) {
            T worst=1;
            for (int f : incident_[v]) worst=std::min(worst,score(mesh,mesh.triangles[f]));
            order.emplace_back(worst,static_cast<int>(v));
        }
        std::sort(order.begin(),order.end());
        size_t changed=0;
        for (const auto& item : order) {
            const int v=item.second;
            const Point<T> original=mesh.vertices[v];
            std::set<int> neighbors;
            std::vector<int> replaced(incident_[v].begin(),incident_[v].end());
            std::vector<Triangle> patch;
            T before=1;
            for (int f : replaced) {
                patch.push_back(mesh.triangles[f]);
                before=std::min(before,score(mesh,mesh.triangles[f]));
                for (int n : mesh.triangles[f]) if (n!=v) neighbors.insert(n);
            }
            if (neighbors.empty()) continue;
            Point<T> average=Point<T>::Zero();
            T length=0;
            for (int n : neighbors) {
                average+=mesh.vertices[n];
                length+=(mesh.vertices[n]-original).norm();
            }
            average/=T(neighbors.size());
            length/=T(neighbors.size());
            const Point<T> normal=field_.gradient(original).normalized();
            Point<T> direction=average-original;
            direction-=normal*direction.dot(normal);
            if (direction.norm()>T(.35)*length) direction*=T(.35)*length/direction.norm();
            for (T fraction : {T(1),T(.5),T(.25),T(.125)}) {
                Point<T> candidate=original+fraction*direction;
                if (!projector_.project(candidate,length*T(.35))) continue;
                if ((candidate-original).norm()>length*T(.5)) continue;
                bool stopped=false;
                for (const auto& stop : stop_points_)
                    if ((candidate-stop).norm()<options_.stop_distance) stopped=true;
                if (stopped) continue;
                mesh.vertices[v]=candidate;
                T after=1;
                for (const auto& t : patch) after=std::min(after,score(mesh,t));
                if (after>before+T(1e-4) && patch_ok(mesh,patch,replaced)) {
                    for (int f : replaced) boxes_[f]=bounds(mesh,mesh.triangles[f]);
                    ++changed;
                    break;
                }
                mesh.vertices[v]=original;
            }
        }
        return changed;
    }
};
} // namespace marching_triangles
