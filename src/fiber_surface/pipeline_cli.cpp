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
        c.example=args.take("--example","sphere"); c.curve=analytical_field(c.example).segment;
    } else if (c.model=="mfa") {
        c.example.clear(); c.first=args.take("--first"); c.second=args.take("--second",c.first);
        if (c.first.empty()) throw std::invalid_argument("--first is required for MFA");
        if (!args.has("--polyline") && (!args.has("--a") || !args.has("--b"))) throw std::invalid_argument("MFA seeds require --polyline or --a and --b range endpoints");
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
        <<"\ninput";
    for (int i=0; i<config.curve.size(); ++i)
        out<<' '<<file_digest(state_path(patch_directory(run,config.curve,i),Stage::stitching));
    out<<"\nEND\n";
    out.close(); if (!out) throw std::runtime_error("Cannot write "+target.string());
}
}

int stage_main(Stage stage,int argc,char** argv) {
    try {
        Arguments args(argc,argv);
        if (args.help()) {
            std::cout<<"Fiber-surface stage "<<int(stage)<<" --run-dir DIRECTORY\n";
            std::cout<<"Multi-segment runs save per-segment checkpoints under DIRECTORY.\n";
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
            config.curve=surface_options(args,loaded.field,config.curve.segments.front(),config.options);
            args.finish();
            config.box=loaded.field.box;
            if (config.options.max_triangles>10000000) throw std::invalid_argument("Staged runs support at most 10000000 triangles per segment");
            const auto seeds=find_curve_seeds(loaded.field,config.curve,curve_options(config.curve,config.options));
            write_config(run,config);
            for (int i=0; i<config.curve.size(); ++i) {
                write_seeds(patch_directory(run,config.curve,i),config,seeds[i]);
                std::cout<<"Segment "<<i<<" surface candidates="<<seeds[i].surface.size()
                    <<" endpoint candidates="<<seeds[i].endpoints[0].size()<<','<<seeds[i].endpoints[1].size()<<'\n';
            }
            return 0;
        }
        args.finish();
        const auto config=read_config(run);
        auto loaded=load_field(config,true);
        const auto& curve=config.curve;
        const auto options=curve_options(curve,config.options);
        std::vector<Seeds> seeds;
        if (stage==Stage::fibers || stage==Stage::marching) {
            for (int i=0; i<curve.size(); ++i) {
                seeds.push_back(read_seeds(patch_directory(run,curve,i),config));
                validate_seeds(seeds.back(),loaded.field,curve.segments[i],options);
            }
        }
        std::vector<Result> patches;
        if (stage==Stage::fibers) patches=trace_curve_fibers(loaded.field,curve,options,seeds);
        int status=0;
        for (int i=0; i<curve.size(); ++i) {
            const auto dir=patch_directory(run,curve,i);
            if (stage==Stage::evaluation) {
                patches.push_back(read_checkpoint(dir,config,Stage::stitching,i).result);
                continue;
            }
            Checkpoint state;
            double prior=0;
            if (stage==Stage::fibers) {
                prior=seeds[i].seconds;
                state.result=std::move(patches[i]);
            } else {
                const auto previous=stage==Stage::marching ? Stage::fibers : Stage::marching;
                auto input=read_checkpoint(dir,config,previous,i);
                prior=input.result.seconds;
                Extractor extractor(loaded.field,curve.segments[i],options);
                state=stage==Stage::marching ? extractor.march(seeds[i],std::move(input.result))
                                             : extractor.finish(std::move(input),curve.size()>1);
            }
            write_checkpoint(dir,config,stage,state,i);
            write_stage_visuals(dir,stage,state);
            write_stage_report(dir,stage,state,state.result.seconds-prior);
            std::cout<<"Segment "<<i<<" triangles="<<state.result.mesh.triangles.size()
                <<" fibers="<<state.result.fibers.size()<<" unresolved="<<state.result.unresolved.size()<<'\n';
            if (stage!=Stage::marching && !state.result.unresolved.empty()) status=2;
        }
        if (stage!=Stage::evaluation) return status;
        auto surface=curve.size()==1 ? std::move(patches.front()) : join_curve_patches(curve,patches);
        status=evaluate_and_write(std::move(surface),loaded.field,curve,options,(run/"surface").string());
        fs::copy_file(run/"surface_report.json",run/"report.json",fs::copy_options::overwrite_existing);
        evaluation_manifest(run,config);
        return status;
    } catch (const std::exception& e) {
        std::cerr<<"fiber-surface stage "<<int(stage)<<": "<<e.what()<<'\n'; return 1;
    }
}
} // namespace fiber_surface
