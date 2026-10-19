#pragma once
#include "field.h"
#include <filesystem>
#include <memory>

namespace fiber_surface {
namespace fs = std::filesystem;
enum class Stage { seeds=1, fibers=2, marching=3, stitching=4, evaluation=5 };

struct RunConfig {
    std::string model="analytical", example="sphere", first, second;
    std::string first_digest, second_digest;
    int variable1=0, variable2=0, block1=-1, block2=-1;
    Box box;
    Segment segment;
    Options options;
};
struct LoadedField {
    // Destroy model evaluators before finalizing MPI. Analytical runs need no MPI.
    std::shared_ptr<void> runtime;
    Field field;
};

LoadedField load_field(const RunConfig& config, bool restrict_box);
std::string file_digest(const fs::path& path);
std::string config_signature(const RunConfig& config);
void write_config(const fs::path& run, const RunConfig& config);
RunConfig read_config(const fs::path& run);
void write_seeds(const fs::path& run, const RunConfig& config, const Seeds& seeds);
Seeds read_seeds(const fs::path& run, const RunConfig& config);
void validate_seeds(const Seeds& seeds, const Field& field, const Segment& segment, const Options& options);
fs::path state_path(const fs::path& run, Stage stage);
void write_checkpoint(const fs::path& run, const RunConfig& config, Stage stage, const Checkpoint& state);
Checkpoint read_checkpoint(const fs::path& run, const RunConfig& config, Stage stage);
void write_stage_visuals(const fs::path& run, Stage stage, const Checkpoint& state);
void write_stage_report(const fs::path& run, Stage stage, const Checkpoint& state, double stage_seconds);
int stage_main(Stage stage, int argc, char** argv);
} // namespace fiber_surface
