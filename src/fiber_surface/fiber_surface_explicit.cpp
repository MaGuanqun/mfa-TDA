#include "analytical_field.h"
#include "cli.h"

int main(int argc,char** argv) {
    using namespace fiber_surface;
    try {
        Arguments args(argc,argv);
        if (args.help()) { std::cout<<"fiber_surface_explicit --example sphere|closed-sphere|strip|disconnected (default sphere)\n"; common_help(); return 0; }
        auto example=analytical_field(args.take("--example","sphere"));
        const auto output=args.take("--output","fiber_surface"); Options options;
        common_options(args,example.field,example.segment,options); args.finish();
        return run_and_write(example.field,example.segment,options,output);
    } catch (const std::exception& e) { std::cerr<<"fiber_surface_explicit: "<<e.what()<<'\n'; return 1; }
}
