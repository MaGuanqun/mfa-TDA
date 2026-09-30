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

    static void close_cracks(Mesh& mesh, const Insert& insert, size_t max_triangles) {
        while (!mesh.boundary.empty() && mesh.triangles.size() < max_triangles) {
            const auto loops=mesh.boundary_loops();
            bool progress=false;
            // First split long contours at a facing, non-neighboring vertex.
            for (const auto& loop:loops) {
                if(loop.size()<=6) continue;
                if (try_contour(mesh,loop,insert,true)) {progress=true;break;}
            }
            if(progress) continue;
            // Mesh simple holes by validated ears; a blind fan may cross a
            // concave contour or introduce a third face on an existing edge.
            for(const auto& loop:loops) {
                if(try_contour(mesh,loop,insert,false)) {progress=true;break;}
            }
            if(progress) continue;
            // Facing contours on tubular sections need joining, not capping.
            for(size_t i=0;i<loops.size()&&!progress;++i)
                for(size_t j=0;j<loops.size()&&!progress;++j) {
                    if(i==j) continue;
                    for(size_t k=0;k<loops[i].size()&&!progress;++k) {
                        const int a=loops[i][k],b=loops[i][(k+1)%loops[i].size()];
                        for(int c:nearest(mesh,a,b,loops[j]))
                            if(insert({a,c,b})) {progress=true;break;}
                    }
                }
            if(!progress) break;
        }
    }
private:
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
    static bool try_contour(Mesh& mesh,const std::vector<int>& loop,const Insert& insert,bool split) {
        const size_t n=loop.size();
        if(n<3) return false;
        for(size_t k=0;k<n;++k) {
            const int a=loop[k],b=loop[(k+1)%n];
            std::vector<int> candidates;
            if(split) {
                for(size_t j=0;j<n;++j)
                    if(j!=k && j!=(k+1)%n && j!=(k+n-1)%n && j!=(k+2)%n) candidates.push_back(loop[j]);
            } else {candidates={loop[(k+n-1)%n],loop[(k+2)%n]};}
            for(int c:nearest(mesh,a,b,candidates))
                if(insert({a,c,b})) return true;
        }
        return false;
    }
};
} // namespace crack_closing
