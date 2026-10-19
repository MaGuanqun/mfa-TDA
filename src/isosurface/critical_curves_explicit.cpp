#include "critical_explicit_field.h"
#include "critical_curve_io.h"
#include "root_io.h"

int main(int argc,char** argv) {
    using namespace critical_curves;
    try {
        Options options;
        std::string name="crossing_circle",prefix="critical_circle",seed_file;
        double offset=1.25,level=std::numeric_limits<double>::quiet_NaN();
        bool help=false;
        opts::Options cli; add_options(cli,options);
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"function",name,"crossing_line, crossing_circle, parallel_lines, isolated_rank_two, definite_line, tilted_line");
        cli >> opts::Option("offset",offset,"additive field offset b; default 1.25");
        cli >> opts::Option('v',"iso",level,"target value; default is the known critical value b");
        cli >> opts::Option('o',"output",prefix,"output prefix for OBJ, VTK, CSV, and report files");
        cli >> opts::Option('s',"seeds",seed_file,"optional xyz seed file; default searches a grid at the selected value");
        if(!cli.parse(argc,argv)) {std::cerr << cli;return 1;}
        if(help) {std::cout << cli;return 0;}
        if(std::isnan(level)) level=offset;
        if(!std::isfinite(level)) throw std::invalid_argument("The isovalue must be finite");
        Problem problem(explicit_field(name,offset),options);
        std::vector<Point> supplied;
        if(!seed_file.empty()) supplied=marching_triangles::read_roots(seed_file);
        const auto search=discover(problem,level,seed_file.empty()?nullptr:&supplied);
        print_search(problem,search);
        const auto curves=trace(problem,search.seeds,level);
        write_curves(prefix,problem,search,curves,name+";offset="+std::to_string(offset)+";seeds="+seed_file,level);
        return 0;
    } catch(const std::exception& e) {std::cerr << "critical_curves_explicit: " << e.what() << '\n';return 1;}
}
