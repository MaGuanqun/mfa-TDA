#include "cli.h"
#include "extract.h"
#include "pipeline.h"
#include "volume.h"

int main(int argc, char** argv) {
    using namespace fiber_surface;
    try {
        Arguments args(argc, argv);
        if (args.help()) {
            std::cout << "fiber_volume_from_patches --run-dir DIRECTORY --output PREFIX\n"
                      << "Assemble, validate and fill a saved counterclockwise three-segment loop.\n";
            return 0;
        }
        const fs::path run=args.take("--run-dir");
        const auto prefix=args.take("--output");
        args.finish();
        if (run.empty() || prefix.empty()) throw std::invalid_argument("Supply --run-dir and --output");
        const auto config=read_config(run);
        const auto& curve=config.curve;
        if (curve.size()!=3 || !curve.closed()) throw std::invalid_argument("Expected a closed three-segment loop");
        const RangeTriangle range(curve.point(0),curve.point(1),curve.point(2));
        const auto loaded=load_field(config,true);
        const auto opts=curve_options(curve,config.options);
        std::vector<Result> patches;
        for (int i=0;i<curve.size();++i)
            patches.push_back(read_checkpoint(patch_directory(run,curve,i),config,Stage::stitching,i).result);
        auto result=volume_from_patches(loaded.field,range,opts,join_curve_patches(curve,patches));
        write_volume(prefix,result,loaded.field,range,opts);
        std::cout << "Volume complete: " << result.complete() << ", boundary triangles: "
                  << result.boundary.mesh.triangles.size() << ", tetrahedra: " << result.mesh.cells.size() << '\n';
        return result.complete()?0:2;
    } catch (const std::exception& e) {
        std::cerr << "fiber_volume_from_patches: " << e.what() << '\n';
        return 1;
    }
}
