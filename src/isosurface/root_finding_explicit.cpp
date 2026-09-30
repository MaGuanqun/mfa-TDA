#include "explicit_field.h"
#include "find_initial_root.h"
#include "root_io.h"
#include "../opts.h"

int main(int argc,char** argv) {
    using namespace marching_triangles;
    try {
        std::string name="ellipsoid",filename="root.dat",shrink="0-1-0-1-0-1";
        double resolution=1,epsilon=1e-8,iso=std::numeric_limits<double>::quiet_NaN();
        int iterations=80,samples=20;bool help=false;
        opts::Options cli;
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"input_function_name",name,"3D closed-form function name");
        cli >> opts::Option('s',"root_file",filename,"exact output path (no automatic suffix)");
        cli >> opts::Option('v',"function_value",iso,"isovalue (default: ellipsoid=1, other functions=0)");
        cli >> opts::Option('g',"spatial_step_size",resolution,"resolution divisor: root separation = shortest domain dimension / (20 * g)");
        cli >> opts::Option('x',"root_finding_epsilon",epsilon,"projection residual tolerance");
        cli >> opts::Option('m',"max_itr",iterations,"projection iteration limit");
        cli >> opts::Option('k',"shrink_range",shrink,"legacy option; only full range supported");
        cli >> opts::Option("samples",samples,"initial samples per axis (2 to 500)");
        if(!cli.parse(argc,argv)){std::cerr<<cli;return 1;}
        if(help){std::cout<<cli;return 0;}
        if(!std::isfinite(resolution)||resolution<=0)throw std::invalid_argument("Resolution must be positive");
        if(shrink!="0-1-0-1-0-1")throw std::invalid_argument("Non-default shrink ranges are unsupported");
        if(std::isnan(iso))iso=name=="ellipsoid"?1:0;
        auto setup=explicit_field(name,iso);
        setup.options.projection_tolerance=epsilon;setup.options.projection_iterations=iterations;
        const auto roots=find_initial_roots(setup.field,setup.options,samples,setup.base_step/resolution);
        if(roots.empty())throw std::runtime_error("No regular roots found; check the isovalue or increase --samples");
        write_roots(filename,roots);
        std::cout<<"Wrote "<<roots.size()<<" roots to "<<filename<<'\n';
        return 0;
    } catch(const std::exception& e){std::cerr<<"root_finding_explicit: "<<e.what()<<'\n';return 1;}
}
