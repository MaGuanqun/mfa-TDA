#include "critical_mfa_field.h"
#include "critical_curve_io.h"

int main(int argc,char** argv) {
    using namespace critical_curves;
    diy::mpi::environment env(argc,argv);diy::mpi::communicator world;
    try {
        Options options;
        std::string filename="model.mfa",prefix="critical_values";
        int variable=0,block=-1;bool help=false;
        opts::Options cli;add_options(cli,options);
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"infile",filename,"3D MFA model");
        cli >> opts::Option('o',"output",prefix,"output prefix for candidate values, seed files, and diagnostics");
        cli >> opts::Option("variable",variable,"zero-based scalar variable");
        cli >> opts::Option("block",block,"zero-based block; required for multiblock files");
        if(!cli.parse(argc,argv)){std::cerr<<cli;return 1;}
        if(help){std::cout<<cli;return 0;}
        auto model=std::make_shared<marching_triangles::MfaModel>(filename,world,block);
        Problem problem(mfa_field(model,variable),options);
        const auto search=discover(problem);
        print_search(problem,search);
        write_values(prefix,problem,search,filename+";variable="+std::to_string(variable)+";block="+std::to_string(block));
        return 0;
    }catch(const std::exception& e){std::cerr<<"critical_values_mfa: "<<e.what()<<'\n';return 1;}
}
