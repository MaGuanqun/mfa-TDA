#pragma once
#include "mesh_geometry.h"
#include <sstream>

namespace marching_triangles {
struct MeshValidation {
    size_t boundary_edges=0,nonmanifold_edges=0,inconsistent_edges=0;
    size_t duplicate_faces=0,degenerate_faces=0,invalid_indices=0,unused_vertices=0;
    size_t nonmanifold_vertices=0,intersections=0;
    long euler_characteristic=0;
    bool valid() const {return !(nonmanifold_edges||inconsistent_edges||duplicate_faces||degenerate_faces||invalid_indices||unused_vertices||nonmanifold_vertices||intersections);}
    bool closed() const {return valid()&&boundary_edges==0;}
    std::string summary() const {
        std::ostringstream out;
        out<<"boundary_edges="<<boundary_edges<<" nonmanifold_edges="<<nonmanifold_edges
           <<" inconsistent_edges="<<inconsistent_edges<<" duplicate_faces="<<duplicate_faces
           <<" degenerate_faces="<<degenerate_faces<<" invalid_indices="<<invalid_indices
           <<" unused_vertices="<<unused_vertices<<" nonmanifold_vertices="<<nonmanifold_vertices
           <<" intersections="<<intersections<<" Euler="<<euler_characteristic;
        return out.str();
    }
};
// Reconstruct adjacency from the exported arrays, independently of growth state.
template<class T> MeshValidation validate_mesh(const std::vector<Point<T>>& v,const std::vector<Triangle>& t,
                                               T eps,bool check_intersections=true) {
    MeshValidation result;
    std::map<Edge,std::vector<int>> edges;
    std::set<Triangle> faces;
    std::vector<std::vector<Edge>> links(v.size());
    for(const auto& f:t) {
        bool valid=true;
        for(int i:f)if(i<0||static_cast<size_t>(i)>=v.size())valid=false;
        if(!valid){++result.invalid_indices;continue;}
        if(!faces.insert(make_triangle_key(f[0],f[1],f[2])).second)++result.duplicate_faces;
        const Point<T> a=v[f[1]]-v[f[0]],b=v[f[2]]-v[f[0]];
        if(!v[f[0]].allFinite()||!a.allFinite()||!b.allFinite()||a.cross(b).norm()<=eps*std::max(a.norm(),b.norm()))++result.degenerate_faces;
        for(int j=0;j<3;++j) {
            int p=f[j],q=f[(j+1)%3];
            edges[edge_key(p,q)].push_back(p<q?1:-1);
            links[p].emplace_back(q,f[(j+2)%3]);
        }
    }
    for(const auto& e:edges) {
        if(e.second.size()==1)++result.boundary_edges;
        if(e.second.size()>2)++result.nonmanifold_edges;
        if(e.second.size()==2&&e.second[0]==e.second[1])++result.inconsistent_edges;
    }
    for(const auto& link:links) {
        if(link.empty()){++result.unused_vertices;continue;}
        std::map<int,std::vector<int>> graph;
        for(auto e:link){graph[e.first].push_back(e.second);graph[e.second].push_back(e.first);}
        bool bad=false;int ends=0;
        for(const auto& x:graph){if(x.second.size()==1)++ends;else if(x.second.size()!=2)bad=true;}
        if(ends!=0&&ends!=2)bad=true;
        std::set<int> seen;std::vector<int> stack{graph.begin()->first};
        while(!stack.empty()){int i=stack.back();stack.pop_back();if(!seen.insert(i).second)continue;for(int j:graph[i])stack.push_back(j);}
        if(seen.size()!=graph.size())bad=true;
        if(bad)++result.nonmanifold_vertices;
    }
    if(check_intersections&&!result.invalid_indices&&!result.degenerate_faces)
        for(size_t i=0;i<t.size();++i)for(size_t j=i+1;j<t.size();++j)
            if(triangles_conflict(v,t[i],t[j],eps))++result.intersections;
    result.euler_characteristic=static_cast<long>(v.size())-static_cast<long>(edges.size())+static_cast<long>(t.size());
    return result;
}
} // namespace marching_triangles
