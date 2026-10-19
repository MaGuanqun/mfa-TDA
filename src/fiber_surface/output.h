#pragma once
#include "field.h"
#include "../isosurface/mesh_validation.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

namespace fiber_surface {
struct Samples {
    size_t count=0, violations=0;
    double sum=0, maximum=0, min_s=std::numeric_limits<double>::infinity(), max_s=-std::numeric_limits<double>::infinity();
    void add(const Vec& x, const Field& f, const Segment& s, const Options& o) {
        const Range value=f.value(x); const double residual=std::abs(s.phi(value)), parameter=s.parameter(value);
        ++count; sum+=residual; maximum=std::max(maximum,residual);
        min_s=std::min(min_s,parameter); max_s=std::max(max_s,parameter);
        if (!value.allFinite() || parameter < -o.residual_tolerance/s.length || parameter > 1+o.residual_tolerance/s.length) ++violations;
    }
    double mean() const { return count ? sum/count : 0; }
};
struct Report {
    marching_triangles::MeshValidation topology;
    Samples vertices, midpoints, centroids;
    double endpoint_residual=0, area=0, min_angle=180, median_angle=0, p05_angle=0;
    size_t unexpected_boundary=0, incorrect_feature_incidence=0, components=0;
    std::vector<int> face_components;
    bool complete(const Result& r) const {
        return !r.mesh.triangles.empty() && topology.valid() && r.unresolved.empty() &&
            !unexpected_boundary && !incorrect_feature_incidence && !vertices.violations && !midpoints.violations && !centroids.violations;
    }
};

Report evaluate(const Result& r, const Field& field, const Segment& segment, const Options& options);
void write_outputs(const std::string& prefix, const Result& r, const Report& report,
                   const Field& field, const Segment& segment, const Options& options);
int evaluate_and_write(Result result, const Field& field, const Segment& segment,
                       const Options& options, const std::string& prefix);
int run_and_write(const Field& field, const Segment& segment, const Options& options, const std::string& prefix);
} // namespace fiber_surface
