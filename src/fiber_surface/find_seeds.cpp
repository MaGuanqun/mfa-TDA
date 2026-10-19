#include "corrector.h"
#include <chrono>

namespace fiber_surface {
Seeds find_seeds(const Field& field, const Segment& segment, const Options& options) {
    options.validate(field,segment);
    const auto start=std::chrono::steady_clock::now();
    const auto grid=discovery_grid(field.box,options.seed_grid);
    Corrector corrector(field,segment,options);
    Seeds seeds;
    seeds.starting_points=grid.size();
    for (int endpoint=0; endpoint<2; ++endpoint)
        for (Vec seed:grid)
            if (corrector.fiber_seed(seed,endpoint ? segment.b : segment.a)) seeds.endpoints[endpoint].push_back(seed);
    for (Vec seed:grid)
        if (corrector.surface(seed,(field.box.max-field.box.min).norm())) seeds.surface.push_back(seed);
    // Preserve discovery order. Coverage/component identity is established by
    // tracing and marching, rather than by proximity-only seed deduplication.
    seeds.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return seeds;
}
} // namespace fiber_surface
