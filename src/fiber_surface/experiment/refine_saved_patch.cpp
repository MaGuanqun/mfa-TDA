// Subdivide a proposed local repair while preserving all of its exposed edges.
#include "../pipeline.h"
#include "../volume.h"
#include "../corrector.h"
#include <iostream>

using namespace fiber_surface;

namespace {
struct InternalEdge {
    int label = 0;
    int neighbor = 0;
    int midpoint = -1;
};

void refine(VolumeBoundary& patch, const Field& field, const RangeTriangle& range, Options options) {
    options.surface_solve_tolerance = std::min(1e-10, options.residual_tolerance);
    auto& mesh = patch.mesh;
    std::map<Edge, InternalEdge> edges;
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const int label = patch.patches[i];
        if (label < 1 || label > 3) throw std::invalid_argument("Expected range-side faces");
        const auto& t = mesh.triangles[i];
        for (int j = 0; j < 3; ++j) {
            auto& edge = edges[edge_key(t[j], t[(j+1)%3])];
            if (!edge.label) edge.label = label;
            else if (!edge.neighbor) edge.neighbor = label;
            else throw std::invalid_argument("Nonmanifold repair edge");
        }
    }
    for (auto& [ids, edge] : edges) {
        if (!edge.neighbor) continue;
        const Vec a = mesh.vertices[ids.first], b = mesh.vertices[ids.second];
        const Vec start = volume_midpoint(a,b);
        Vec p = start;
        Corrector solver(field, range.sides[edge.label-1], options);
        if (edge.label == edge.neighbor) {
            if (!solver.surface(p, options.spatial_tolerance) ||
                !range.contains(field.value(p), options.residual_tolerance)) continue;
        } else {
            const int corner = (edge.label%3+1 == edge.neighbor ? edge.neighbor : edge.label)-1;
            const Range& q = range.corners[corner];
            if ((field.value(a)-q).norm() > Options::fiber_field_tolerance ||
                (field.value(b)-q).norm() > Options::fiber_field_tolerance)
                throw std::invalid_argument("Expected a fitted shared corner edge");
            if (!solver.fiber(p, q, (b-a).normalized(), start, options.spatial_tolerance)) {
                p = start;
                if (!solver.fiber_seed(p,q) || (p-start).norm() > options.spatial_tolerance) continue;
            }
        }
        if (!p.allFinite() || (p-a).norm() <= 1e-8*options.step || (p-b).norm() <= 1e-8*options.step) continue;
        edge.midpoint = static_cast<int>(mesh.vertices.size());
        mesh.vertices.push_back(p);
    }

    std::vector<Triangle> faces;
    std::vector<int> labels;
    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        const auto& t = mesh.triangles[i];
        std::array<int,3> m;
        for (int j = 0; j < 3; ++j) m[j] = edges.at(edge_key(t[j],t[(j+1)%3])).midpoint;
        auto add = [&](int a, int b, int c) {
            faces.push_back({a,b,c});
            labels.push_back(patch.patches[i]);
        };
        const int missing = std::count(m.begin(),m.end(),-1);
        if (missing == 3) add(t[0],t[1],t[2]);
        else if (missing == 0) {
            for (int j = 0; j < 3; ++j) add(t[j],m[j],m[(j+2)%3]);
            add(m[0],m[1],m[2]);
        } else if (missing == 2) {
            int j = int(std::find_if(m.begin(),m.end(),[](int v){return v>=0;}) - m.begin());
            add(t[j],m[j],t[(j+2)%3]);
            add(m[j],t[(j+1)%3],t[(j+2)%3]);
        } else {
            int j = int(std::find(m.begin(),m.end(),-1) - m.begin());
            int a = t[j], b = t[(j+1)%3], c = t[(j+2)%3];
            int bc = m[(j+1)%3], ca = m[(j+2)%3];
            add(c,ca,bc);
            if ((mesh.vertices[a]-mesh.vertices[bc]).squaredNorm() <
                (mesh.vertices[b]-mesh.vertices[ca]).squaredNorm()) {
                add(a,b,bc); add(a,bc,ca);
            } else {
                add(a,b,ca); add(b,bc,ca);
            }
        }
    }
    mesh.triangles.swap(faces);
    patch.patches.swap(labels);
}
}

int main(int argc,char** argv) {
    try {
        if (argc!=4) throw std::invalid_argument("Use: fiber_refine_saved_patch STAGES INPUT.ply OUTPUT.ply");
        const auto config = read_config(argv[1]);
        if (config.curve.size()!=3 || !config.curve.closed()) throw std::invalid_argument("Expected a closed range triangle");
        const auto field = load_field(config,true);
        RangeTriangle range(config.curve.point(0),config.curve.point(1),config.curve.point(2));
        auto patch = read_boundary_ply(argv[2]);
        refine(patch,field.field,range,config.options);
        write_boundary_ply(argv[3],patch,false);
    } catch(const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
