#include "mfa_field.h"
#include "extract_isosurface.h"
#include "../opts.h"

int main(int argc,char** argv) {
    using namespace marching_triangles;
    diy::mpi::environment env(argc,argv);
    diy::mpi::communicator world;
    try {
        std::string name="model.mfa",root_file="root.dat",singular_file,output="isosurface_sheet";
        std::string shrink="0-1-0-1-0-1";
        double resolution=1,iso=.5,epsilon=1e-8;
        double min_ratio=.1,stop_distance=.05,legacy_hessian=1e-10;
        int projection_iterations=80,triangle_limit=200000,quality_iterations=8;
        bool help=false,legacy_close=false,no_close=false,allow_open=false;
        int variable=0,block_index=-1;
        opts::Options cli;
        cli >> opts::Option("variable",variable,"zero-based scalar science-variable index");
        cli >> opts::Option("block",block_index,"zero-based block index; required for multiblock files");
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"infile",name,"3D MFA model (.mfa DIY file)");
        cli >> opts::Option('s',"root_file",root_file,"binary root matrix or text xyz file");
        cli >> opts::Option('o',"output_mesh_prefix",output,"output OBJ prefix per connected sheet");
        cli >> opts::Option('v',"function_value",iso,"isovalue (default: 0.5)");
        cli >> opts::Option('g',"spatial_step_size",resolution,"resolution divisor: seed edge length = shortest domain dimension / (20 * g)");
        cli >> opts::Option('x',"root_finding_epsilon",epsilon,"surface projection residual tolerance");
        cli >> opts::Option('m',"max_projection_itr",projection_iterations,"projection iteration limit");
        cli >> opts::Option('n',"d_min_ratio",min_ratio,"minimum adaptive step / seed edge length");
        cli >> opts::Option('d',"singular_point_file",singular_file,"optional stop-set points");
        cli >> opts::Option('c',"stop_curve_distance",stop_distance,"stop distance from supplied singular points");
        cli >> opts::Option('H',"hessian_rank_threshold",legacy_hessian,"legacy option; regularity uses the gradient");
        cli >> opts::Option('k',"shrink_range",shrink,"legacy option; only full range 0-1-0-1-0-1 supported");
        cli >> opts::Option('r',"enable_crack_closing",legacy_close,"compatibility flag; crack closing is enabled by default");
        cli >> opts::Option("no-crack-closing",no_close,"debug: stop after growth");
        cli >> opts::Option("allow-open",allow_open,"allow a valid mesh with boundary (e.g. a clipped surface)");
        cli >> opts::Option("max-triangles",triangle_limit,"maximum triangles per sheet");
        cli >> opts::Option("quality-iterations",quality_iterations,"closed-surface quality improvement passes (default 8; 0 disables; max 100)");
        if(!cli.parse(argc,argv)){std::cerr<<cli;return 1;}
        if(help){std::cout<<cli;return 0;}
        if(shrink!="0-1-0-1-0-1")throw std::invalid_argument("Non-default shrink ranges are unsupported; supply explicit bounds through the library API");
        if(!std::isfinite(resolution)||resolution<=0||triangle_limit<4)throw std::invalid_argument("Resolution must be positive and max-triangles at least 4");
        auto model=std::make_shared<MfaModel>(name,world,block_index);
        auto setup=mfa_field(model,iso,variable);
        auto& options=setup.options;
        options.step=setup.base_step/resolution;options.min_step=options.step*min_ratio;
        options.vertex_tolerance=options.step*1e-6;options.projection_tolerance=epsilon;
        options.projection_iterations=projection_iterations;options.stop_distance=stop_distance;
        options.max_triangles=static_cast<size_t>(triangle_limit);options.close_cracks=!no_close;
        options.quality_iterations=quality_iterations;
        options.validate();
        extract_isosurface(setup,root_file,singular_file,output,allow_open);
        return 0;
    } catch(const std::exception& e){std::cerr<<"isosurface_mfa: "<<e.what()<<'\n';return 1;}
}
