#include "cli.h"
#include "extract.h"
#include "pipeline.h"
#include "volume.h"

int main(int argc, char** argv) {
    using namespace fiber_surface;
    try {
        Arguments args(argc, argv);
        if (args.help()) {
            std::cout << "fiber_volume_from_patches --run-dir DIRECTORY --output PREFIX [--stage all|surface]\n"
                      << "Assemble, validate and fill a saved counterclockwise three-segment loop.\n"
                      << "The surface stage exports just the assembled PLY for an independent surface audit.\n";
            return 0;
        }
        const fs::path run=args.take("--run-dir");
        const auto prefix=args.take("--output");
        const auto stage=args.take("--stage","all");
        args.finish();
        if(stage!="all" && stage!="surface") throw std::invalid_argument("Expected --stage all or surface");
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
        auto joined=join_curve_patches(curve,patches);
        patches.clear();
        if(stage=="surface") {
            const auto boundary=assemble_volume_boundary(loaded.field,range,opts,std::move(joined));
            if(!fs::path(prefix).parent_path().empty()) fs::create_directories(fs::path(prefix).parent_path());
            write_boundary_ply(prefix+"_boundary.ply",boundary,false);
            std::cout<<"Assembled "<<boundary.mesh.triangles.size()<<" triangles, "
                     <<boundary.unresolved.size()<<" unresolved events.\n";
            for(const auto& event:boundary.unresolved)
                std::cout<<event.position.transpose()<<"  "<<event.reason<<'\n';
            return 0;
        }
        auto result=volume_from_patches(loaded.field,range,opts,std::move(joined));
        write_volume(prefix,result,loaded.field,range,opts);
        std::cout << "Volume complete: " << result.complete() << ", boundary triangles: "
                  << result.boundary.mesh.triangles.size() << ", tetrahedra: " << result.mesh.cells.size() << '\n';
        return result.complete()?0:2;
    } catch (const std::exception& e) {
        std::cerr << "fiber_volume_from_patches: " << e.what() << '\n';
        return 1;
    }
}
