#include "mfa_field.h"
#include "find_initial_root.h"
#include "root_io.h"
#include "../opts.h"
#include <algorithm>

int main(int argc,char** argv) {
    using namespace marching_triangles;
    diy::mpi::environment env(argc,argv);
    diy::mpi::communicator world;
    try {
        std::string name="model.mfa",filename="root.dat",shrink="0-1-0-1-0-1";
        double resolution=1,epsilon=1e-8,iso=.5;
        int iterations=80,samples=20;bool help=false;
        int variable=0,block_index=-1;
        opts::Options cli;
        cli >> opts::Option("variable",variable,"zero-based scalar science-variable index");
        cli >> opts::Option("block",block_index,"zero-based block index; required for multiblock files");
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"infile",name,"3D MFA model (.mfa DIY file)");
        cli >> opts::Option('s',"root_file",filename,"exact output path (no automatic suffix)");
        cli >> opts::Option('v',"function_value",iso,"isovalue (default: 0.5)");
        cli >> opts::Option('g',"spatial_step_size",resolution,"resolution divisor: root separation = shortest domain dimension / (20 * g)");
        cli >> opts::Option('x',"root_finding_epsilon",epsilon,"projection residual tolerance");
        cli >> opts::Option('m',"max_itr",iterations,"projection iteration limit");
        cli >> opts::Option('k',"shrink_range",shrink,"legacy option; only full range supported");
        cli >> opts::Option("samples",samples,"initial samples per axis (2 to 500)");
        if(!cli.parse(argc,argv)){std::cerr<<cli;return 1;}
        if(help){std::cout<<cli;return 0;}
        if(!std::isfinite(resolution)||resolution<=0)throw std::invalid_argument("Resolution must be positive");
        if(shrink!="0-1-0-1-0-1")throw std::invalid_argument("Non-default shrink ranges are unsupported");
        auto model=std::make_shared<MfaModel>(name,world,block_index);
        auto setup=mfa_field(model,iso,variable);
        setup.options.projection_tolerance=epsilon;setup.options.projection_iterations=iterations;
        auto roots=find_initial_roots(setup.field,setup.options,samples,setup.base_step/resolution);
        if(roots.empty())throw std::runtime_error("No regular roots found; check the isovalue or increase --samples");
        // Start well inside the domain: a boundary seed can create a tiny
        // partial sheet before another seed reaches the same component.
        const auto clearance=[&](const Point<double>& p) {
            return (p-setup.options.domain_min).cwiseMin(setup.options.domain_max-p).minCoeff();
        };
        std::stable_sort(roots.begin(),roots.end(),[&](const auto& a,const auto& b){return clearance(a)>clearance(b);});
        write_roots(filename,roots);
        std::cout<<"Wrote "<<roots.size()<<" roots to "<<filename<<'\n';
        return 0;
    } catch(const std::exception& e){std::cerr<<"root_finding_mfa: "<<e.what()<<'\n';return 1;}
}
