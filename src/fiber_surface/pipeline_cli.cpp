#include "pipeline.h"
#include "analytical_field.h"
#include "cli.h"
#include "extract.h"
#include <chrono>

namespace fiber_surface {
namespace {
RunConfig initial_config(Arguments& args) {
    RunConfig c;
    c.model=args.take("--model","analytical");
    if (c.model=="analytical") {
        c.example=args.take("--example","sphere"); c.segment=analytical_field(c.example).segment;
    } else if (c.model=="mfa") {
        c.example.clear(); c.first=args.take("--first"); c.second=args.take("--second",c.first);
        if (c.first.empty()) throw std::invalid_argument("--first is required for MFA");
        if (!args.has("--a") || !args.has("--b")) throw std::invalid_argument("MFA seeds require --a and --b range endpoints");
        c.first=fs::absolute(c.first).lexically_normal().string(); c.second=fs::absolute(c.second).lexically_normal().string();
        c.variable1=args.integer("--variable1",0); c.variable2=args.integer("--variable2",c.first==c.second?1:0);
        c.block1=args.integer("--block1",-1); c.block2=args.integer("--block2",c.block1);
        c.first_digest=file_digest(c.first); c.second_digest=file_digest(c.second);
    } else throw std::invalid_argument("--model must be analytical or mfa");
    return c;
}
void evaluation_manifest(const fs::path& run,const RunConfig& config) {
    const auto target=run/"evaluation.conf";
    std::ofstream out(target);
    out<<"FIBER_SURFACE_EVALUATION 1\nconfig "<<config_signature(config)
        <<"\ninput "<<file_digest(state_path(run,Stage::stitching))<<"\nEND\n";
    out.close(); if (!out) throw std::runtime_error("Cannot write "+target.string());
}
}

int stage_main(Stage stage,int argc,char** argv) {
    try {
        Arguments args(argc,argv);
        if (args.help()) {
            std::cout<<"Fiber-surface stage "<<int(stage)<<" --run-dir DIRECTORY\n";
            std::cout<<"All outputs are saved directly in DIRECTORY.\n";
            if (stage==Stage::seeds) {
                std::cout<<"  --model analytical --example sphere|closed-sphere|strip|disconnected\n"
                    <<"  --model mfa --first FILE [--second FILE] --a u,v --b u,v\n"
                    <<"  --variable1 N --variable2 N --block1 N --block2 N\n"
                    <<"  All resolved parameters are saved in run.conf.\n";
                common_help();
            } else std::cout<<"Reads run.conf and preceding stage files. To change parameters, restart at find_seeds.\n";
            return 0;
        }
        const std::string directory=args.take("--run-dir");
        if (directory.empty()) throw std::invalid_argument("--run-dir is required");
        const fs::path run=fs::absolute(directory).lexically_normal();
        if (stage==Stage::seeds) {
            auto config=initial_config(args);
            auto loaded=load_field(config,false);
            if (config.model=="mfa") config.options.step=(loaded.field.box.max-loaded.field.box.min).minCoeff()/20;
            common_options(args,loaded.field,config.segment,config.options); args.finish();
            config.box=loaded.field.box;
            if (config.options.max_triangles>10000000) throw std::invalid_argument("Staged runs support at most 10000000 triangles");
            const auto seeds=find_seeds(loaded.field,config.segment,config.options);
            write_config(run,config); write_seeds(run,config,seeds);
            std::cout<<"Surface candidates="<<seeds.surface.size()<<" endpoint candidates="<<seeds.endpoints[0].size()<<','<<seeds.endpoints[1].size()
                <<"\nSaved "<<run<<'\n';
            return 0;
        }
        args.finish();
        const auto config=read_config(run);
        auto loaded=load_field(config,true);
        if (stage==Stage::fibers) {
            const auto seeds=read_seeds(run,config);
            validate_seeds(seeds,loaded.field,config.segment,config.options);
            Checkpoint state;
            state.result=trace_endpoint_fibers(loaded.field,config.segment,config.options,seeds);
            write_checkpoint(run,config,stage,state); write_stage_visuals(run,stage,state);
            write_stage_report(run,stage,state,state.result.seconds-seeds.seconds);
            std::cout<<"Traced "<<state.result.fibers.size()<<" endpoint components; unresolved="<<state.result.unresolved.size()<<'\n';
            return state.result.unresolved.empty()?0:2;
        }
        if (stage==Stage::marching) {
            const auto seeds=read_seeds(run,config);
            validate_seeds(seeds,loaded.field,config.segment,config.options);
            auto fibers=read_checkpoint(run,config,Stage::fibers);
            const double prior=fibers.result.seconds;
            auto state=Extractor(loaded.field,config.segment,config.options).march(seeds,std::move(fibers.result));
            write_checkpoint(run,config,stage,state); write_stage_visuals(run,stage,state);
            write_stage_report(run,stage,state,state.result.seconds-prior);
            std::cout<<"Marching saved "<<state.result.mesh.triangles.size()<<" triangles; run stitch_surface next\n";
            // Exposed fronts are expected at this stage, not final failures.
            return 0;
        }
        if (stage==Stage::stitching) {
            auto input=read_checkpoint(run,config,Stage::marching); const double prior=input.result.seconds;
            auto state=Extractor(loaded.field,config.segment,config.options).finish(std::move(input));
            write_checkpoint(run,config,stage,state); write_stage_visuals(run,stage,state);
            write_stage_report(run,stage,state,state.result.seconds-prior);
            std::cout<<"Stitching saved "<<state.result.mesh.triangles.size()<<" triangles; unresolved="<<state.result.unresolved.size()<<'\n';
            return state.result.unresolved.empty()?0:2;
        }
        auto input=read_checkpoint(run,config,Stage::stitching);
        fs::create_directories(run);
        const int status=evaluate_and_write(std::move(input.result),loaded.field,config.segment,config.options,(run/"surface").string());
        // Keep the final report at a predictable name in addition to the legacy
        // prefix convention used by the preview and one-shot executables.
        fs::copy_file(run/"surface_report.json",run/"report.json",fs::copy_options::overwrite_existing);
        evaluation_manifest(run,config);
        return status;
    } catch (const std::exception& e) {
        std::cerr<<"fiber-surface stage "<<int(stage)<<": "<<e.what()<<'\n'; return 1;
    }
}
} // namespace fiber_surface
