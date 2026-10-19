#include "mfa_field.h"
#include "cli.h"

int main(int argc,char** argv) {
    using namespace fiber_surface;
    diy::mpi::environment environment(argc,argv); diy::mpi::communicator world;
    try {
        Arguments args(argc,argv);
        if (args.help()) {
            std::cout<<"fiber_surface_mfa --first FILE [--second FILE] --variable1 0 --variable2 0\n"
                <<"  Omit --second to select two scalar variables from --first.\n"
                <<"  --block1 N --block2 N select blocks explicitly in multiblock files.\n"
                <<"  --inspect prints sampled range and physical derivative checks, without extraction.\n";
            common_help(); return 0;
        }
        const auto first=args.take("--first"),second=args.take("--second",first);
        if (first.empty()) throw std::invalid_argument("--first MFA file is required");
        if (!args.has("--inspect") && (!args.has("--a") || !args.has("--b")))
            throw std::invalid_argument("MFA extraction requires --a and --b range endpoints; use --inspect to sample the range first");
        const int variable1=args.integer("--variable1",0),variable2=args.integer("--variable2",first==second?1:0);
        const int block1=args.integer("--block1",-1),block2=args.integer("--block2",block1);
        auto model1=std::make_shared<marching_triangles::MfaModel>(first,world,block1);
        auto model2=first==second && block1==block2 ? model1 : std::make_shared<marching_triangles::MfaModel>(second,world,block2);
        Field field=mfa_field(model1,model2,variable1,variable2); Segment segment; Options options;
        options.step=(field.box.max-field.box.min).minCoeff()/20;
        const auto output=args.take("--output","fiber_surface_mfa"); const bool inspect=!args.take("--inspect").empty();
        common_options(args,field,segment,options); args.finish();
        if (inspect) {
            Range lo=Range::Constant(std::numeric_limits<double>::infinity()),hi=-lo;
            for (const Vec& x:discovery_grid(field.box,options.seed_grid)) { const Range f=field.value(x); lo=lo.cwiseMin(f); hi=hi.cwiseMax(f); }
            std::cout<<std::setprecision(17)<<"Sampled range min: "<<lo.transpose()<<" max: "<<hi.transpose()<<'\n';
            const Vec center=(field.box.min+field.box.max)/2;
            std::cout<<"Center xyz: "<<center.transpose()<<" F: "<<field.value(center).transpose()<<"\nJ:\n"<<field.jacobian(center)<<'\n';
            double error=0;
            for (int d=0; d<3; ++d) {
                const double h=(field.box.max[d]-field.box.min[d])*1e-5; Vec a=center,b=center; a[d]-=h; b[d]+=h;
                const Range fd=(field.value(b)-field.value(a))/(2*h);
                error=std::max(error,(fd-field.jacobian(center).col(d)).norm()/std::max(1.,fd.norm()));
            }
            std::cout<<"Center finite-difference gradient relative error: "<<error<<'\n'; return error<1e-4 ? 0 : 2;
        }
        return run_and_write(field,segment,options,output);
    } catch (const std::exception& e) { std::cerr<<"fiber_surface_mfa: "<<e.what()<<'\n'; return 1; }
}
