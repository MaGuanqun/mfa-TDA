#include "corrector.h"
#include <chrono>
#include <tbb/parallel_for.h>

namespace fiber_surface {
Seeds find_seeds(const Field& field, const Segment& segment, const Options& options) {
    options.validate(field,segment);
    const auto start=std::chrono::steady_clock::now();
    const auto shape=discovery_shape(field,options);
    const auto grid=discovery_grid(field.box,shape);
    Corrector corrector(field,segment,options);
    Seeds seeds;
    seeds.starting_points=grid.size();
    seeds.grid_shape=shape;
    const double bound=(field.box.max-field.box.min).norm();
    // Each worker owns its output slot. Compact in grid order after the join;
    // tracing/marching depend on this order when identifying covered components.
    for (int target=0; target<3; ++target) {
        std::vector<Vec> corrected(grid.size());
        std::vector<unsigned char> accepted(grid.size(),0);
        tbb::parallel_for(size_t(0),grid.size(),[&](size_t i) {
            Vec seed=grid[i];
            accepted[i]=target==2 ? corrector.surface(seed,bound) : corrector.fiber_seed(seed,target ? segment.b : segment.a);
            if (accepted[i]) corrected[i]=seed;
        });
        auto& output=target==2 ? seeds.surface : seeds.endpoints[target];
        for (size_t i=0; i<grid.size(); ++i) if (accepted[i]) output.push_back(corrected[i]);
    }
    // Preserve discovery order. Coverage/component identity is established by
    // tracing and marching, rather than by proximity-only seed deduplication.
    seeds.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return seeds;
}
} // namespace fiber_surface
