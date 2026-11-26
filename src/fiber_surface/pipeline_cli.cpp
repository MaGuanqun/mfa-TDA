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
    } else if (c.model=="mfa" || c.model=="coordnet") {
        c.example.clear(); c.first=args.take("--first");
        if (c.first.empty()) throw std::invalid_argument("--first is required for "+c.model);
        if (!args.has("--polyline") && (!args.has("--a") || !args.has("--b")))
            throw std::invalid_argument("Supply --polyline or --a and --b range endpoints");
        c.second=args.take("--second",c.model=="mfa" ? c.first : "");
        if (c.second.empty()) throw std::invalid_argument("--second is required for "+c.model);
        c.first=fs::absolute(c.first).lexically_normal().string(); c.second=fs::absolute(c.second).lexically_normal().string();
        if (c.model=="coordnet") {
            c.training_box=coordnet_training_box(args);
        } else {
            c.variable1=args.integer("--variable1",0); c.variable2=args.integer("--variable2",c.first==c.second?1:0);
            c.block1=args.integer("--block1",-1); c.block2=args.integer("--block2",c.block1);
        }
        c.first_digest=file_digest(c.first); c.second_digest=file_digest(c.second);
    } else throw std::invalid_argument("--model must be analytical, mfa or coordnet");
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

// Allow separate Slurm jobs to discover sides of a prepared curve. The second
// run may change tracing and marching distances, while using the same roots.
void discover_segment(const fs::path& run, int side, const std::string& also_run) {
    const auto config=read_config(run);
    if (side<0 || side>=config.curve.size()) throw std::invalid_argument("Segment index outside the saved range curve");
    const auto options=curve_options(config.curve,config.options);
    std::optional<RunConfig> other;
    if (!also_run.empty()) {
        other=read_config(also_run);
        auto comparable=*other;
        // Permit only these downstream controls to differ. The residual also
        // participates in Newton, so compare its effective seed threshold.
        const auto target=curve_options(other->curve,other->options);
        auto seed_tolerance=[](const Options& o) {
            return o.surface_solve_tolerance>0 ? std::min(o.surface_solve_tolerance,o.residual_tolerance) : o.residual_tolerance;
        };
        auto& o=comparable.options;
        o.step=config.options.step;
        o.min_step=config.options.min_step;
        o.spatial_tolerance=config.options.spatial_tolerance;
        o.fiber_seed_radius=config.options.fiber_seed_radius;
        o.residual_tolerance=config.options.residual_tolerance;
        if (seed_tolerance(options)!=seed_tolerance(target) || config_signature(comparable)!=config_signature(config))
            throw std::invalid_argument("Seed discovery settings differ between the two run directories");
    }
    auto loaded=load_field(config,true);
    const bool last=!config.curve.closed() && side+1==config.curve.size();
    const auto seeds=find_seeds(loaded.field,config.curve.segments[side],options,last);
    write_seeds(patch_directory(run,config.curve,side),config,seeds);
    if (other) write_seeds(patch_directory(also_run,other->curve,side),*other,seeds);
    std::cout<<"Segment "<<side<<" surface candidates="<<seeds.surface.size()
             <<" endpoint candidates="<<seeds.endpoints[0].size()<<','<<seeds.endpoints[1].size()<<'\n';
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
                    <<"  --model coordnet --first U.pt --second V.pt --domain-min x,y,t --domain-max x,y,t\n"
                    <<"    Use float64 scalar exports and range endpoints in model-output units.\n"
                    <<"  All resolved parameters are saved in run.conf.\n";
                std::cout<<"  --max-seed-starts N  seed lattice allocation limit (default 1000000)\n";
                std::cout<<"  --prepare-only  save run.conf for parallel discovery, then run this executable\n"
                         <<"    with --run-dir DIRECTORY --segment I [--also-run-dir COMPATIBLE_DIRECTORY]\n";
                common_help();
            } else {
                std::cout<<"Reads run.conf and preceding stage files. To change parameters, restart at find_seeds.\n";
                if (stage==Stage::marching || stage==Stage::stitching)
                    std::cout<<"  --segment I   process only zero-based range segment I\n";
            }
            return 0;
        }
        const std::string directory=args.take("--run-dir");
        if (directory.empty()) throw std::invalid_argument("--run-dir is required");
        const fs::path run=fs::absolute(directory).lexically_normal();
        if (stage==Stage::seeds) {
            if (args.has("--segment")) {
                const int side=args.integer("--segment",0);
                const auto other=args.take("--also-run-dir");
                args.finish();
                discover_segment(run,side,other);
                return 0;
            }
            const bool prepare=!args.take("--prepare-only").empty();
            auto config=initial_config(args);
            config.options.max_seed_starts=args.integer("--max-seed-starts",config.options.max_seed_starts);
            if (config.options.max_seed_starts<1 || config.options.max_seed_starts>10000000)
                throw std::invalid_argument("Staged runs support a seed allocation limit from 1 to 10000000");
            auto loaded=load_field(config,false);
            if (config.model!="analytical") config.options.step=(loaded.field.box.max-loaded.field.box.min).minCoeff()/20;
            config.curve=surface_options(args,loaded.field,config.curve.segments.front(),config.options);
            args.finish();
            config.box=loaded.field.box;
            if (config.options.max_triangles>10000000) throw std::invalid_argument("Staged runs support at most 10000000 triangles per segment");
            if (prepare) {
                write_config(run,config);
                return 0;
            }
            const auto seeds=find_curve_seeds(loaded.field,config.curve,curve_options(config.curve,config.options));
            write_config(run,config);
            for (int i=0; i<config.curve.size(); ++i) {
                write_seeds(patch_directory(run,config.curve,i),config,seeds[i]);
                std::cout<<"Segment "<<i<<" surface candidates="<<seeds[i].surface.size()
                    <<" endpoint candidates="<<seeds[i].endpoints[0].size()<<','<<seeds[i].endpoints[1].size()<<'\n';
            }
            return 0;
        }
        const bool single = (stage==Stage::marching || stage==Stage::stitching) && args.has("--segment");
        const int side = single ? args.integer("--segment",0) : 0;
        args.finish();
        const auto config=read_config(run);
        if (side<0 || side>=config.curve.size()) throw std::invalid_argument("Segment index outside the saved range curve");
        auto loaded=load_field(config,true);
        const auto& curve=config.curve;
        const auto options=curve_options(curve,config.options);
        std::vector<Seeds> seeds(curve.size());
        if (stage==Stage::fibers || stage==Stage::marching) {
            for (int i=0; i<curve.size(); ++i) {
                if (single && i!=side) continue;
                seeds[i]=read_seeds(patch_directory(run,curve,i),config);
                validate_seeds(seeds[i],loaded.field,curve.segments[i],options);
            }
        }
        std::vector<Result> patches;
        if (stage==Stage::fibers) patches=trace_curve_fibers(loaded.field,curve,options,seeds);
        int status=0;
        for (int i=0; i<curve.size(); ++i) {
            if (single && i!=side) continue;
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
