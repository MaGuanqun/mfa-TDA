#include "pipeline.h"
#include "analytical_field.h"
#include "mfa_field.h"

namespace fiber_surface {
namespace {
struct MpiRuntime { diy::mpi::environment environment; diy::mpi::communicator world; };
}
LoadedField load_field(const RunConfig& config, bool restrict_box) {
    LoadedField loaded;
    if (config.model=="analytical") loaded.field=analytical_field(config.example).field;
    else if (config.model=="mfa") {
        if (config.first.empty() || config.second.empty()) throw std::invalid_argument("Both MFA paths are required");
        if (file_digest(config.first)!=config.first_digest || file_digest(config.second)!=config.second_digest)
            throw std::runtime_error("MFA input changed since seed discovery; restart at find_seeds");
        auto runtime=std::make_shared<MpiRuntime>();
        auto first=std::make_shared<marching_triangles::MfaModel>(config.first,runtime->world,config.block1);
        auto second=config.first==config.second && config.block1==config.block2 ? first :
            std::make_shared<marching_triangles::MfaModel>(config.second,runtime->world,config.block2);
        loaded.runtime=runtime;
        loaded.field=mfa_field(first,second,config.variable1,config.variable2);
    } else throw std::invalid_argument("--model must be analytical or mfa");
    if (restrict_box) {
        if (!loaded.field.box.contains(config.box.min) || !loaded.field.box.contains(config.box.max))
            throw std::runtime_error("Saved spatial box is outside the model domain");
        loaded.field.box=config.box;
        config.options.validate(loaded.field,config.segment);
    }
    return loaded;
}
} // namespace fiber_surface
