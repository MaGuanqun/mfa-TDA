#pragma once

#include "surface_field.h"
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>

namespace marching_triangles {
using Triangle = std::array<int, 3>;
using Edge = std::pair<int, int>;
inline Edge edge_key(int a, int b) { return {std::min(a,b), std::max(a,b)}; }
inline Triangle make_triangle_key(int a, int b, int c) {
    Triangle t{a,b,c};
    std::sort(t.begin(),t.end());
    return t;
}

template<class T> Point<T> closest_on_segment(const Point<T>& p, const Point<T>& a, const Point<T>& b) {
    const Point<T> e = b-a;
    const T n = e.squaredNorm();
    if (n == 0) return a;
    return a + std::clamp((p-a).dot(e)/n,T(0),T(1))*e;
}

template<class T> Point<T> closest_on_triangle(const Point<T>& p, const Point<T>& a,
                                              const Point<T>& b, const Point<T>& c) {
    const Point<T> ab=b-a, ac=c-a, ap=p-a;
    const T d1=ab.dot(ap), d2=ac.dot(ap);
    if (d1<=0 && d2<=0) return a;
    const Point<T> bp=p-b;
    const T d3=ab.dot(bp), d4=ac.dot(bp);
    if (d3>=0 && d4<=d3) return b;
    const T vc=d1*d4-d3*d2;
    if (vc<=0 && d1>=0 && d3<=0) return a+(d1/(d1-d3))*ab;
    const Point<T> cp=p-c;
    const T d5=ab.dot(cp), d6=ac.dot(cp);
    if (d6>=0 && d5<=d6) return c;
    const T vb=d5*d2-d1*d6;
    if (vb<=0 && d2>=0 && d6<=0) return a+(d2/(d2-d6))*ac;
    const T va=d3*d6-d5*d4;
    if (va<=0 && d4-d3>=0 && d5-d6>=0) return b+((d4-d3)/((d4-d3)+(d5-d6)))*(c-b);
    const T denom=va+vb+vc;
    if (denom == 0) return closest_on_segment(p,a,b);
    return a+(vb/denom)*ab+(vc/denom)*ac;
}

template<class T> bool circumcircle(const Point<T>& a,const Point<T>& b,const Point<T>& c,Point<T>& center,T& r2) {
    const Point<T> u=b-a,v=c-a,n=u.cross(v);
    const T n2=n.squaredNorm();
    if (n2 <= T(1e-24)*u.squaredNorm()*v.squaredNorm()) return false;
    center=a+(u.squaredNorm()*v.cross(n)+v.squaredNorm()*n.cross(u))/(2*n2);
    r2=(center-a).squaredNorm();
    return center.allFinite() && std::isfinite(r2);
}

// True for any contact outside the shared topological vertex/edge. Includes
// coplanar overlaps and T-junctions, which segment/plane tests alone miss.
template<class T> bool triangles_conflict(const std::vector<Point<T>>& v, const Triangle& a,
                                          const Triangle& b, T eps) {
    Point<T> alo=v[a[0]],ahi=alo,blo=v[b[0]],bhi=blo;
    std::vector<int> shared;
    for (int i:a) {
        alo=alo.cwiseMin(v[i]); ahi=ahi.cwiseMax(v[i]);
        if (std::find(b.begin(),b.end(),i)!=b.end()) shared.push_back(i);
    }
    for (int i:b) {blo=blo.cwiseMin(v[i]);bhi=bhi.cwiseMax(v[i]);}
    if ((alo.array()>bhi.array()+eps).any() || (blo.array()>ahi.array()+eps).any()) return false;
    if (shared.size()==3) return true;
    Point<T> na=(v[a[1]]-v[a[0]]).cross(v[a[2]]-v[a[0]]).normalized();
    Point<T> nb=(v[b[1]]-v[b[0]]).cross(v[b[2]]-v[b[0]]).normalized();
    auto allowed=[&](const Point<T>& p) {
        if (shared.size()==2) return (p-closest_on_segment(p,v[shared[0]],v[shared[1]])).norm()<=eps;
        return shared.size()==1 && (p-v[shared[0]]).norm()<=eps;
    };
    const bool coplanar=na.cross(nb).norm()<T(1e-8) && std::abs(na.dot(v[b[0]]-v[a[0]]))<=eps;
    if (shared.size()==2 && !coplanar) return false;
    auto inside=[&](const Point<T>& p,const Triangle& t) {
        return (p-closest_on_triangle(p,v[t[0]],v[t[1]],v[t[2]])).norm()<=eps;
    };
    if (coplanar) {
        for(int i:a) if(!allowed(v[i]) && inside(v[i],b)) return true;
        for(int i:b) if(!allowed(v[i]) && inside(v[i],a)) return true;
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
            const Point<T> p=v[a[i]], q=v[b[j]], u=v[a[(i+1)%3]]-p, w=v[b[(j+1)%3]]-q;
            const T det=na.dot(u.cross(w));
            if(std::abs(det)>eps*(u.norm()+w.norm())) {
                const T s=na.dot((q-p).cross(w))/det, t=na.dot((q-p).cross(u))/det;
                if(s>=0 && s<=1 && t>=0 && t<=1 && !allowed(p+s*u)) return true;
            } else {
                for(const Point<T>& x:std::array<Point<T>,2>{p,p+u})
                    if((x-closest_on_segment(x,q,Point<T>(q+w))).norm()<=eps && !allowed(x)) return true;
                for(const Point<T>& x:std::array<Point<T>,2>{q,q+w})
                    if((x-closest_on_segment(x,p,Point<T>(p+u))).norm()<=eps && !allowed(x)) return true;
            }
        }
        // Identical coplanar interiors with all vertices on shared edges.
        const Point<T> ca=(v[a[0]]+v[a[1]]+v[a[2]])/3;
        const Point<T> cb=(v[b[0]]+v[b[1]]+v[b[2]])/3;
        return (!allowed(ca)&&inside(ca,b)) || (!allowed(cb)&&inside(cb,a));
    }
    auto crosses=[&](const Triangle& s,const Triangle& t,const Point<T>& n) {
        for(int i=0;i<3;++i) {
            const Point<T> p=v[s[i]], q=v[s[(i+1)%3]];
            const T dp=n.dot(p-v[t[0]]), dq=n.dot(q-v[t[0]]);
            if(std::abs(dp)<=eps && inside(p,t) && !allowed(p)) return true;
            if((dp>eps && dq>eps)||(dp < -eps && dq < -eps)) continue;
            if(std::abs(dp-dq)>std::numeric_limits<T>::epsilon()*(q-p).norm()) {
                const T u=dp/(dp-dq);
                if(u>=0 && u<=1) {
                    const Point<T> x=p+u*(q-p);
                    if(inside(x,t)&&!allowed(x)) return true;
                }
            }
        }
        return false;
    };
    return crosses(a,b,nb)||crosses(b,a,na);
}

// Chord triangles on a curved sheet can cover the same patch without literally
// intersecting in R^3. Test local same-facing faces in the candidate's plane as
// well. The distance bound excludes remote parts of the surface.
template<class T> bool surface_overlap(const std::vector<Point<T>>& v,const Triangle& a,
                                       const Triangle& b,T eps,const Point<T>& direction=Point<T>::Zero()) {
    const Point<T> n=direction.squaredNorm()>0?Point<T>(direction.normalized()):Point<T>((v[a[1]]-v[a[0]]).cross(v[a[2]]-v[a[0]]).normalized());
    const Point<T> nb=(v[b[1]]-v[b[0]]).cross(v[b[2]]-v[b[0]]).normalized();
    if(n.dot(nb)<=0) return false;
    const T h=std::max({(v[a[1]]-v[a[0]]).norm(),(v[a[2]]-v[a[0]]).norm(),(v[a[2]]-v[a[1]]).norm()});
    Point<T> lo=v[a[0]],hi=lo,blo=v[b[0]],bhi=blo;
    for(int i:a){lo=lo.cwiseMin(v[i]);hi=hi.cwiseMax(v[i]);}
    for(int i:b){blo=blo.cwiseMin(v[i]);bhi=bhi.cwiseMax(v[i]);}
    if((lo.array()>bhi.array()+h*T(.5)).any()||(blo.array()>hi.array()+h*T(.5)).any()) return false;
    T distance=std::numeric_limits<T>::max();
    for(int i:b) distance=std::min(distance,std::abs(n.dot(v[i]-v[a[0]])));
    if(distance>h*T(.5)) return false;
    std::vector<Point<T>> projected;
    for(int i:a) projected.push_back(v[i]-n*n.dot(v[i]-v[a[0]]));
    Triangle other;
    for(int i=0;i<3;++i) {
        auto it=std::find(a.begin(),a.end(),b[i]);
        if(it!=a.end())other[i]=static_cast<int>(it-a.begin());
        else {
            other[i]=static_cast<int>(projected.size());
            projected.push_back(v[b[i]]-n*n.dot(v[b[i]]-v[a[0]]));
        }
    }
    return triangles_conflict(projected,Triangle{0,1,2},other,eps);
}

template<class T> struct TriangleMesh {
    struct EdgeUse { int from, to, face; int second = -1; };
    std::vector<Point<T>> vertices;
    std::vector<Triangle> triangles;
    std::map<Edge,EdgeUse> edges;
    // Directed boundary halfedges have the same direction as their incident face.
    std::map<Edge,int> boundary;
    std::set<Triangle> faces;

    bool topology_ok(const Triangle& t) const {
        if(t[0]==t[1]||t[1]==t[2]||t[2]==t[0] || faces.count(make_triangle_key(t[0],t[1],t[2]))) return false;
        for(int i=0;i<3;++i) {
            const int a=t[i],b=t[(i+1)%3];
            auto it=edges.find(edge_key(a,b));
            if(it!=edges.end() && (it->second.second>=0 || it->second.from==a)) return false;
        }
        return true;
    }
    void add(const Triangle& t) {
        const int face=static_cast<int>(triangles.size());
        triangles.push_back(t);
        register_face(t,face);
    }
    // Replace a checked two-face patch without changing face IDs. The caller
    // must first verify the new diagonal, winding, and geometric intersections.
    void replace_pair(int first,int second,const Triangle& a,const Triangle& b) {
        unregister_face(triangles[first],first);
        unregister_face(triangles[second],second);
        triangles[first]=a;
        triangles[second]=b;
        register_face(a,first);
        register_face(b,second);
    }
private:
    void register_face(const Triangle& t,int face) {
        faces.insert(make_triangle_key(t[0],t[1],t[2]));
        for(int i=0;i<3;++i) {
            const int a=t[i],b=t[(i+1)%3];
            auto it=edges.find(edge_key(a,b));
            if(it==edges.end()) {edges.emplace(edge_key(a,b),EdgeUse{a,b,face,-1});boundary[{a,b}]=face;}
            else {it->second.second=face;boundary.erase({b,a});}
        }
    }
    void unregister_face(const Triangle& t,int face) {
        faces.erase(make_triangle_key(t[0],t[1],t[2]));
        for(int i=0;i<3;++i) {
            const Edge key=edge_key(t[i],t[(i+1)%3]);
            auto it=edges.find(key);
            auto& use=it->second;
            if(use.second<0) {
                boundary.erase({use.from,use.to});
                edges.erase(it);
            } else {
                if(use.face==face) {
                    std::swap(use.from,use.to);
                    use.face=use.second;
                }
                use.second=-1;
                boundary[{use.from,use.to}]=use.face;
            }
        }
    }
public:
    // Follow halfedges through incident faces, not a guessed single successor
    // per vertex: contour splits and joins can visit the same vertex twice.
    std::vector<std::vector<int>> boundary_loops() const {
        std::set<Edge> remaining;
        for(const auto& e:boundary) remaining.insert(e.first);
        std::vector<std::vector<int>> loops;
        while(!remaining.empty()) {
            const Edge start=*remaining.begin(); Edge current=start;
            std::vector<int> loop;
            do {
                if(!remaining.erase(current)) throw std::runtime_error("Broken boundary halfedge cycle");
                loop.push_back(current.first);
                int face=boundary.at(current), pivot=current.second, previous=current.first;
                for(size_t count=0;;++count) {
                    if(count>triangles.size()) throw std::runtime_error("Broken vertex link");
                    const Triangle& t=triangles[face]; int next=-1;
                    for(int k=0;k<3;++k) if(t[k]==pivot) next=t[(k+1)%3];
                    if(next<0||next==previous) throw std::runtime_error("Inconsistent face orientation");
                    if(boundary.count({pivot,next})) {current={pivot,next};break;}
                    const auto& use=edges.at(edge_key(pivot,next));
                    previous=next;face=use.face==face?use.second:use.face;
                    if(face<0) throw std::runtime_error("Missing boundary halfedge");
                }
            } while(current!=start);
            // A splitting triangle can touch an existing boundary vertex in
            // a second fan. Its halfedge walk then visits that vertex twice.
            // Split this weak contour there, as in paper Fig. 8, instead of
            // treating it as one polygon and adding spurious handles.
            std::vector<int> path;
            loop.push_back(loop.front());
            for(int vertex:loop) {
                auto repeated=std::find(path.begin(),path.end(),vertex);
                if(repeated==path.end()) path.push_back(vertex);
                else {
                    if(path.end()-repeated<3) throw std::runtime_error("Degenerate boundary contour");
                    loops.emplace_back(repeated,path.end());
                    path.erase(repeated+1,path.end());
                }
            }
        }
        return loops;
    }
};
} // namespace marching_triangles
