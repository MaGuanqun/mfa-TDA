#include "explicit_field.h"
#include "marching_triangles.h"
#include "mesh_validation.h"
#include "root_io.h"
#include "../opts.h"

int main(int argc,char** argv) {
    using namespace marching_triangles;
    try {
        std::string name="ellipsoid",root_file="root.dat",singular_file,output="isosurface_sheet";
        std::string shrink="0-1-0-1-0-1";
        double resolution=1,iso=std::numeric_limits<double>::quiet_NaN(),epsilon=1e-8;
        double min_ratio=.1,stop_distance=.05,legacy_hessian=1e-10;
        int projection_iterations=80,triangle_limit=200000,quality_iterations=8;
        bool help=false,legacy_close=false,no_close=false,allow_open=false;
        opts::Options cli;
        cli >> opts::Option('h',"help",help,"show help");
        cli >> opts::Option('f',"input_function_name",name,"3D closed-form function name");
        cli >> opts::Option('s',"root_file",root_file,"binary root matrix or text xyz file");
        cli >> opts::Option('o',"output_mesh_prefix",output,"output OBJ prefix per connected sheet");
        cli >> opts::Option('v',"function_value",iso,"isovalue (default: ellipsoid=1, other functions=0)");
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
        if(std::isnan(iso))iso=name=="ellipsoid"?1:0;
        auto setup=explicit_field(name,iso);
        auto& options=setup.options;
        options.step=setup.base_step/resolution;options.min_step=options.step*min_ratio;
        options.vertex_tolerance=options.step*1e-6;options.projection_tolerance=epsilon;
        options.projection_iterations=projection_iterations;options.stop_distance=stop_distance;
        options.max_triangles=static_cast<size_t>(triangle_limit);options.close_cracks=!no_close;
        options.quality_iterations=quality_iterations;
        options.validate();
        const auto seeds=read_roots(root_file);
        MarchingTriangles<double> mesher(setup.field,options);
        if(!singular_file.empty())mesher.set_degenerate_points(read_roots(singular_file));
        std::vector<std::vector<Point<double>>> vertices;
        std::vector<std::vector<Triangle>> triangles;
        std::cout<<"Loaded "<<seeds.size()<<" roots; seed edge length="<<options.step<<'\n';
        if(!mesher.extract_all_sheets(seeds,vertices,triangles))throw std::runtime_error("No regular surface sheet found; check the seeds and isovalue");
        if(!mesher.diagnostic().empty())std::cerr<<mesher.diagnostic()<<'\n';
        // Validate before publishing any of the meshes. Failure is observable
        // as a nonzero exit code, not a success message attached to an open OBJ.
        bool valid=true;
        for(size_t i=0;i<vertices.size();++i){
            const auto report=validate_mesh(vertices[i],triangles[i],options.step*1e-9);
            std::cout<<"sheet "<<i<<": vertices="<<vertices[i].size()<<" triangles="<<triangles[i].size()<<' '<<report.summary()<<'\n';
            valid=valid&&report.valid()&&(allow_open||report.closed());
            if(report.invalid_indices==0) {
                const auto quality=measure_mesh_quality(vertices[i],triangles[i],setup.field);
                std::cout<<"sheet "<<i<<" quality: "<<quality.summary()<<'\n';
                valid=valid&&quality.reversed_faces==0;
                if(quality.triangles_below_5_degrees || quality.max_normal_error_degrees>45)
                    std::cerr<<"sheet "<<i<<": poor triangle shapes or normal alignment remain; inspect the quality report and consider a finer resolution.\n";
            }
        }
        if(!valid)throw std::runtime_error("Mesh validation failed. No output written; reduce the step size or use --allow-open for an intentionally bounded surface");
        for(size_t i=0;i<vertices.size();++i){
            const std::string filename=output+std::to_string(i)+".obj";
            mesher.save_mesh_obj(filename,vertices[i],triangles[i]);
            std::cout<<"Wrote "<<filename<<'\n';
        }
        return 0;
    } catch(const std::exception& e){std::cerr<<"isosurface_explicit: "<<e.what()<<'\n';return 1;}
}
