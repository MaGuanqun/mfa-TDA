#pragma once
#include "output.h"

namespace fiber_surface {
// Difference-based affine combinations preserve a shared box coordinate
// exactly and avoid cancellation when the geometry is far from the origin.
inline Vec volume_midpoint(const Vec& a,const Vec& b) {return a+(b-a)/2;}
inline Vec volume_centroid(const Vec& a,const Vec& b,const Vec& c) {return a+((b-a)+(c-a))/3;}
inline std::array<Vec,7> volume_triangle_samples(const Vec& a,const Vec& b,const Vec& c) {
    return {a,b,c,volume_midpoint(a,b),volume_midpoint(b,c),volume_midpoint(c,a),volume_centroid(a,b,c)};
}
struct RangeTriangle {
    std::array<Range,3> corners;
    std::array<Segment,3> sides;
    RangeTriangle(Range a=Range(0,0), Range b=Range(1,0), Range c=Range(0,1));
    bool contains(const Range& value, double tolerance=0) const;
};

struct VolumeBoundary {
    Mesh mesh;
    // 1..3: range edges; 4..9: xmin,xmax,ymin,ymax,zmin,zmax.
    std::vector<int> patches;
    std::vector<Event> unresolved;
    std::array<size_t,3> corner_components{};
    size_t box_contact_triangles=0;
};
struct BoundaryAssessment {
    marching_triangles::MeshValidation topology;
    std::vector<int> shells;
    std::vector<double> signed_volumes;
    std::vector<Vec> holes;
    size_t constraint_failures=0, orientation_failures=0, coverage_misses=0;
    double side_residual=0, corner_residual=0, inequality_violation=0, spatial_error=0;
    bool valid=false;
};
using Tetrahedron=std::array<int,4>;
struct VolumeMesh {
    std::vector<Vec> vertices;
    std::vector<Tetrahedron> cells;
    std::vector<Triangle> boundary_faces;
    std::vector<int> parent_facets;
};
struct VolumeAssessment {
    marching_triangles::MeshValidation output_topology;
    size_t invalid_cells=0, invalid_faces=0, conformity_failures=0, intersections=0;
    double volume=0, boundary_volume=0, min_volume=0, max_volume=0;
    double side_residual=0, inequality_violation=0, spatial_error=0;
    bool valid=false;
};
struct VolumeResult {
    VolumeBoundary boundary;
    BoundaryAssessment boundary_report;
    VolumeMesh mesh;
    VolumeAssessment report;
    std::string failed_stage;
    std::array<double,4> seconds{};
    bool complete() const { return failed_stage.empty() && boundary_report.valid && report.valid; }
};

VolumeBoundary assemble_volume_boundary(const Field&, const RangeTriangle&, const Options&);
BoundaryAssessment validate_volume_boundary(VolumeBoundary&, const Field&, const RangeTriangle&, const Options&);
VolumeMesh fill_volume(const VolumeBoundary&, const BoundaryAssessment&);
VolumeAssessment assess_volume(const VolumeMesh&, const VolumeBoundary&, const BoundaryAssessment&,
                               const Field&, const RangeTriangle&, const Options&);
VolumeResult extract_volume(const Field&, const RangeTriangle&, const Options&);
void write_volume(const std::string&, const VolumeResult&, const Field&, const RangeTriangle&, const Options&);
int run_volume(const Field&, const RangeTriangle&, const Options&, const std::string&);
// Winding membership of a closed shell, used for containment and coverage.
bool inside_shell(const Vec&, const Mesh&, const std::vector<int>& faces);
bool correct_volume_sample(Vec&, int patch, const Field&, const RangeTriangle&, const Options&);
} // namespace fiber_surface
