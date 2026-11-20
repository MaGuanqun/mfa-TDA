#include "pipeline.h"
#include "analytical_field.h"
#include "coordnet_field.h"
#include "mfa_field.h"

namespace fiber_surface {
namespace {
struct MpiRuntime { diy::mpi::environment environment; diy::mpi::communicator world; };
}
LoadedField load_field(const RunConfig& config, bool restrict_box) {
    LoadedField loaded;
    if (config.model=="mfa" || config.model=="coordnet") {
        if (file_digest(config.first)!=config.first_digest || file_digest(config.second)!=config.second_digest) {
            const std::string kind=config.model=="mfa" ? "MFA" : "CoordNet";
            throw std::runtime_error(kind+" input changed since seed discovery; restart at find_seeds");
        }
    }
    if (config.model=="analytical") loaded.field=analytical_field(config.example).field;
    else if (config.model=="coordnet")
        loaded.field=coordnet_field(config.first,config.second,config.training_box);
    else if (config.model=="mfa") {
        auto runtime=std::make_shared<MpiRuntime>();
        auto first=std::make_shared<marching_triangles::MfaModel>(config.first,runtime->world,config.block1);
        auto second=config.first==config.second && config.block1==config.block2 ? first :
            std::make_shared<marching_triangles::MfaModel>(config.second,runtime->world,config.block2);
        loaded.runtime=runtime;
        loaded.field=mfa_field(first,second,config.variable1,config.variable2);
    } else throw std::invalid_argument("--model must be analytical, mfa or coordnet");
    if (restrict_box) {
        if (!loaded.field.box.contains(config.box.min) || !loaded.field.box.contains(config.box.max))
            throw std::runtime_error("Saved spatial box is outside the model domain");
        loaded.field.box=config.box;
        for (const auto& s:config.curve.segments) config.options.validate(loaded.field,s);
    }
    return loaded;
}
} // namespace fiber_surface
