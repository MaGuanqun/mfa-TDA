#pragma once

#include "mesh_geometry.h"
#include <functional>

namespace crack_closing {
// Section 3.2.2. All insertion goes through the mesher's geometry/topology
// validator. A failed insertion leaves the contour intact. No new vertices.
template<class T> class CrackCloser {
public:
    using Mesh = marching_triangles::TriangleMesh<T>;
    using Triangle = marching_triangles::Triangle;
    using Insert = std::function<bool(const Triangle&)>;
    using EligibleEdge = std::function<bool(int,int)>;

    static void close_cracks(Mesh& mesh, const Insert& insert, size_t max_triangles,
                             const EligibleEdge& eligible={}, T max_edge=std::numeric_limits<T>::infinity()) {
        while (!mesh.boundary.empty() && mesh.triangles.size() < max_triangles) {
            const auto loops=mesh.boundary_loops();
            std::vector<Contour> contours;
            for (const auto& vertices:loops) contours.emplace_back(mesh,vertices,eligible,max_edge);
            bool progress=false;
            // First split long contours at a facing, non-neighboring vertex.
            for (const auto& loop:contours) {
                if(loop.ids.size()<=6) continue;
                if (try_contour(mesh,loop,insert,true,max_edge)) {progress=true;break;}
            }
            if(progress) continue;
            // Mesh simple holes by validated ears; a blind fan may cross a
            // concave contour or introduce a third face on an existing edge.
            for(const auto& loop:contours) {
                if(try_contour(mesh,loop,insert,false,max_edge)) {progress=true;break;}
            }
            if(progress) continue;
            // Facing contours on tubular sections need joining, not capping.
            for(size_t i=0;i<loops.size()&&!progress;++i)
                for(size_t j=0;j<loops.size()&&!progress;++j) {
                    if(i==j) continue;
                    for(size_t k:contours[i].active) {
                        const int a=loops[i][k],b=loops[i][(k+1)%loops[i].size()];
                        for(int c:nearest(mesh,a,b,contours[j].candidates(mesh,a,max_edge)))
                            if(insert({a,c,b})) {progress=true;break;}
                        if(progress) break;
                    }
                }
            if(!progress) break;
        }
    }
private:
    struct Contour {
        const std::vector<int>& ids;
        std::vector<size_t> active, by_x;
        Contour(const Mesh& mesh,const std::vector<int>& loop,const EligibleEdge& eligible,T max_edge)
            : ids(loop) {
            for (size_t k=0;k<ids.size();++k) {
                if (!eligible || eligible(ids[k],ids[(k+1)%ids.size()])) active.push_back(k);
                if (std::isfinite(max_edge)) by_x.push_back(k);
            }
            std::sort(by_x.begin(),by_x.end(),[&](size_t i,size_t j) {
                return mesh.vertices[ids[i]].x()<mesh.vertices[ids[j]].x();
            });
        }
        std::vector<int> candidates(const Mesh& mesh,int a,T max_edge,
                                    size_t edge=std::numeric_limits<size_t>::max()) const {
            std::vector<int> result;
            auto add=[&](size_t j) {
                const size_t n=ids.size();
                if (edge<n && (j==edge || j==(edge+1)%n || j==(edge+n-1)%n || j==(edge+2)%n)) return;
                result.push_back(ids[j]);
            };
            if (by_x.empty()) {
                for (size_t j=0;j<ids.size();++j) add(j);
            } else {
                // The caller rejects longer edges. Round the slab outward,
                // then retain its original distance and vertex-ID ordering.
                const T inf=std::numeric_limits<T>::infinity();
                const T radius=std::nextafter(max_edge,inf), x=mesh.vertices[a].x();
                const T low=std::nextafter(x-radius,-inf),high=std::nextafter(x+radius,inf);
                auto it=std::lower_bound(by_x.begin(),by_x.end(),low,[&](size_t j,T bound) {
                    return mesh.vertices[ids[j]].x()<bound;
                });
                for (;it!=by_x.end() && mesh.vertices[ids[*it]].x()<=high;++it)
                    if ((mesh.vertices[ids[*it]]-mesh.vertices[a]).norm()<=radius) add(*it);
            }
            return result;
        }
    };
    static std::vector<int> nearest(const Mesh& mesh,int a,int b,const std::vector<int>& vertices) {
        std::vector<std::pair<T,int>> candidates;
        for(int c:vertices) if(c!=a&&c!=b) {
            const auto q=marching_triangles::closest_on_segment(mesh.vertices[c],mesh.vertices[a],mesh.vertices[b]);
            candidates.emplace_back((mesh.vertices[c]-q).squaredNorm(),c);
        }
        std::sort(candidates.begin(),candidates.end());
        std::vector<int> result;
        for(const auto& c:candidates) result.push_back(c.second);
        return result;
    }
    static bool try_contour(Mesh& mesh,const Contour& contour,const Insert& insert,bool split,T max_edge) {
        const auto& loop=contour.ids;
        const size_t n=loop.size();
        if(n<3) return false;
        for(size_t k:contour.active) {
            const int a=loop[k],b=loop[(k+1)%n];
            std::vector<int> candidates;
            if(split) candidates=contour.candidates(mesh,a,max_edge,k);
            else candidates={loop[(k+n-1)%n],loop[(k+2)%n]};
            for(int c:nearest(mesh,a,b,candidates))
                if(insert({a,c,b})) return true;
        }
        return false;
    }
};
} // namespace crack_closing
