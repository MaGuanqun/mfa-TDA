#include "critical_mfa_field.h"
#include "critical_curve_io.h"
#include "root_io.h"

int main(int argc,char** argv) {
    using namespace critical_curves;
    diy::mpi::environment env(argc,argv);diy::mpi::communicator world;
    try {
        Options options;
        std::string filename="model.mfa",prefix="critical_curves",seed_file;
        double level=std::numeric_limits<double>::quiet_NaN();
        int variable=0,block=-1;bool help=false;
        opts::Options cli;add_options(cli,options);
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"infile",filename,"3D MFA model");
        cli >> opts::Option('v',"iso",level,"required target isovalue from the candidate catalog or user");
        cli >> opts::Option('s',"seeds",seed_file,"optional value_ID.seeds.txt file; revalidated on this field and level");
        cli >> opts::Option('o',"output",prefix,"output prefix for OBJ, VTK, CSV, and report files");
        cli >> opts::Option("variable",variable,"zero-based scalar variable");
        cli >> opts::Option("block",block,"zero-based block; required for multiblock files");
        if(!cli.parse(argc,argv)){std::cerr<<cli;return 1;}
        if(help){std::cout<<cli;return 0;}
        if(!std::isfinite(level))throw std::invalid_argument("Supply a finite target value using --iso (-v)");
        auto model=std::make_shared<marching_triangles::MfaModel>(filename,world,block);
        Problem problem(mfa_field(model,variable),options);
        std::vector<Point> supplied;
        if(!seed_file.empty())supplied=marching_triangles::read_roots(seed_file);
        const auto search=discover(problem,level,seed_file.empty()?nullptr:&supplied);
        print_search(problem,search);
        const auto curves=trace(problem,search.seeds,level);
        write_curves(prefix,problem,search,curves,filename+";variable="+std::to_string(variable)+";block="+std::to_string(block)+";seeds="+seed_file,level);
        return 0;
    }catch(const std::exception& e){std::cerr<<"critical_curves_mfa: "<<e.what()<<'\n';return 1;}
}
