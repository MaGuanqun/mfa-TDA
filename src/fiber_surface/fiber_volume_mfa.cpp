#include "volume_cli.h"
#include "mfa_field.h"
int main(int argc,char** argv) {
    using namespace fiber_surface;
    diy::mpi::environment environment(argc,argv);diy::mpi::communicator world;
    try {
        Arguments args(argc,argv);
        if(args.help()) {
            std::cout<<"fiber_volume_mfa --first FILE [--second FILE] --variable1 0 --variable2 0\n"
                <<"  Omit --second to use two scalar variables from one file.\n"
                <<"  --block1 N --block2 N select blocks; --inspect samples the selected box.\n";
            volume_help();return 0;
        }
        const auto first=args.take("--first"),second=args.take("--second",first);
        if(first.empty())throw std::invalid_argument("--first MFA file is required");
        const int v1=args.integer("--variable1",0),v2=args.integer("--variable2",first==second?1:0);
        const int b1=args.integer("--block1",-1),b2=args.integer("--block2",b1);
        auto m1=std::make_shared<marching_triangles::MfaModel>(first,world,b1);
        auto m2=first==second && b1==b2?m1:std::make_shared<marching_triangles::MfaModel>(second,world,b2);
        auto field=mfa_field(m1,m2,v1,v2);Options options;options.step=(field.box.max-field.box.min).minCoeff()/20;
        const auto prefix=args.take("--output","fiber_volume_mfa");const bool inspect=!args.take("--inspect").empty();RangeTriangle triangle;
        volume_options(args,field,triangle,options,!inspect);args.finish();
        if(inspect) {
            Range lo=Range::Constant(std::numeric_limits<double>::infinity()),hi=-lo;
            for(const Vec& x:discovery_grid(field.box,discovery_shape(field,options))) {const Range y=field.value(x);lo=lo.cwiseMin(y);hi=hi.cwiseMax(y);}
            std::cout<<std::setprecision(17)<<"Sampled range min: "<<lo.transpose()<<" max: "<<hi.transpose()<<'\n';return 0;
        }
        return run_volume(field,triangle,options,prefix);
    }catch(const std::exception& e){std::cerr<<"fiber_volume_mfa: "<<e.what()<<'\n';return 1;}
}
