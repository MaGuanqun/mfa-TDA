#include "volume_cli.h"
#include "volume_examples.h"
int main(int argc,char** argv) {
    using namespace fiber_surface;
    try {
        Arguments args(argc,argv);
        if(args.help()){std::cout<<"fiber_volume_explicit --example prism|curved-prism|sphere|disconnected|hollow-sphere|hollow-cylinder (default prism)\n";volume_help();return 0;}
        const auto name=args.take("--example","prism");auto e=volume_example(name);
        Options options; if(name=="curved-prism" || name=="sphere")options.residual_tolerance=.005;
        const auto prefix=args.take("--output","fiber_volume");
        volume_options(args,e.field,e.triangle,options);args.finish();
        return run_volume(e.field,e.triangle,options,prefix);
    }catch(const std::exception& e){std::cerr<<"fiber_volume_explicit: "<<e.what()<<'\n';return 1;}
}
