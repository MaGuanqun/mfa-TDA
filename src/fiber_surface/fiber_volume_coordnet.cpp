#include "coordnet_field.h"
#include "volume_cli.h"

int main(int argc, char** argv) {
    using namespace fiber_surface;
    try {
        Arguments args(argc, argv);
        if (args.help()) {
            std::cout << "fiber_volume_coordnet --first U.pt --second V.pt --domain-min x,y,t --domain-max x,y,t\n"
                      << "  Supply the full physical training box and two float64 scalar TorchScript models.\n"
                      << "  Range vertices use model-output units. --inspect samples the selected box.\n";
            volume_help();
            return 0;
        }
        const auto first = args.take("--first"), second = args.take("--second");
        if (first.empty() || second.empty())
            throw std::invalid_argument("Supply both CoordNet scalar models with --first and --second");
        auto field = coordnet_field(first, second, coordnet_training_box(args));
        return run_model_volume(args, field, "fiber_volume_coordnet");
    } catch (const std::exception& e) {
        std::cerr << "fiber_volume_coordnet: " << e.what() << '\n';
        return 1;
    }
}
