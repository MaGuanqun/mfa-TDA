#include "mfa_field.h"
#include "marching_triangles.h"
#include "mesh_validation.h"
#include <filesystem>

using namespace marching_triangles;
using P=Point<double>;
static void require(bool pass,const std::string& message) {if(!pass)throw std::runtime_error(message);}

static mfa::MFA_Data<double>* polynomial(int degree,int first,int components,const P& origin,const P& extent) {
    auto* data=new mfa::MFA_Data<double>(Eigen::VectorXi::Constant(3,degree),size_t(1),first,first+components-1);
    auto& t=data->tmesh.tensor_prods[0];const int n=degree+1,count=n*n*n;
    t.nctrl_pts=Eigen::VectorXi::Constant(3,n);t.ctrl_pts.resize(count,components);t.weights=Eigen::VectorXd::Ones(count);
    t.knot_mins.assign(3,0);t.knot_maxs.assign(3,2*n-1);t.level=0;
    data->tmesh.all_knots.resize(3);data->tmesh.all_knot_levels.resize(3);
    for(int d=0;d<3;++d) {
        data->tmesh.all_knots[d].assign(n,0);data->tmesh.all_knots[d].insert(data->tmesh.all_knots[d].end(),n,1);
        data->tmesh.all_knot_levels[d].assign(2*n,0);
    }
    for(int i=0;i<count;++i) {
        int index=i;P u;double value=0;
        for(int d=0;d<3;++d){int j=index%n;index/=n;u[d]=double(j)/degree;value+=j==1?-1.:1.;}
        if(components==3)t.ctrl_pts.row(i)=(origin+extent.cwiseProduct(u)).transpose();
        else t.ctrl_pts(i,0)=value; // Bernstein coefficients of sum_d (2*u_d-1)^2.
    }
    data->tmesh.tensor_knot_idxs(t);return data;
}
static void check_gradients(const FieldSetup<double>& setup) {
    const P range=setup.options.domain_max-setup.options.domain_min;
    double worst=0;
    for(int i=1;i<=9;++i) {
        const P p=setup.options.domain_min+range.cwiseProduct(P(.07*i,.11+.07*i,.93-.07*i));
        const P grad=setup.field.gradient(p);
        for(int d=0;d<3;++d) {
            const double h=range[d]*1e-6;P lo=p,hi=p;lo[d]-=h;hi[d]+=h;
            const double fd=(setup.field.value(hi)-setup.field.value(lo))/(2*h);
            worst=std::max(worst,std::abs(fd-grad[d])/std::max(1.,std::abs(fd)));
        }
    }
    require(worst<1e-5,"MFA gradient disagrees with physical-coordinate finite differences");
    std::cout<<"MFA finite-difference relative error: "<<worst<<std::endl;
}
int main(int argc,char** argv) {
    diy::mpi::environment env(argc,argv);diy::mpi::communicator world;
    try {
        if(argc>1) {
            auto setup=mfa_field(std::make_shared<MfaModel>(argv[1],world),.5);
            check_gradients(setup);return 0;
        }
        const auto dir=std::filesystem::current_path()/"mfa_test_files";std::filesystem::create_directories(dir);
        diy::Master writer(world,1,-1,&Block<double>::create,&Block<double>::destroy,nullptr,&Block<double>::save,&Block<double>::load);
        auto* b=new Block<double>();b->dom_dim=3;b->pt_dim=4;b->save_datasets=0;
        const P origin(2,-5,7),extent(6,8,10);
        b->bounds_mins=b->core_mins=origin;b->bounds_maxs=b->core_maxs=origin+extent;
        std::vector<mfa::MFA_Data<double>*> vars{polynomial(2,3,1,origin,extent)};
        b->mfa=new mfa::MFA<double>(3,0,polynomial(1,0,3,origin,extent),vars);
        writer.add(0,b,new diy::RegularLink<diy::Bounds<double>>);
        const auto filename=(dir/"current.mfa").string();
        diy::io::write_blocks(filename,world,writer,&Block<double>::save);
        auto setup=mfa_field(std::make_shared<MfaModel>(filename,world),.5);
        for(int i=0;i<=10;++i) {
            const P u(i/10.,.2,.7),p=origin+extent.cwiseProduct(u);
            require(std::abs(setup.field.value(p)-(2*u-P::Ones()).squaredNorm())<1e-12,"Wrong MFA scalar value");
            require((setup.field.gradient(p)-(8*u-P::Constant(4)).cwiseQuotient(extent)).norm()<1e-12,"Wrong physical MFA gradient");
        }
        check_gradients(setup);
        require(std::isnan(setup.field.value(origin-P::Ones())),"MFA evaluated outside its domain");
        setup.options.step=.8;setup.options.min_step=.4;
        MarchingTriangles<double> mesher(setup.field,setup.options);
        std::vector<std::vector<P>> vertices;std::vector<std::vector<Triangle>> triangles;
        const P seed=origin+extent*.5+P(extent[0]*std::sqrt(.5)/2,0,0);
        require(mesher.extract_all_sheets({seed},vertices,triangles)&&vertices.size()==1,"MFA ellipsoid extraction failed");
        const auto report=validate_mesh(vertices[0],triangles[0],1e-9);
        require(report.closed()&&report.euler_characteristic==2,"MFA ellipsoid topology: "+report.summary());
        for(const auto& p:vertices[0])require(std::abs(setup.field.value(p)-.5)<=1.01e-8,"MFA mesh vertex off surface");
        bool rejected=false;
        try {mfa_field(std::make_shared<MfaModel>(filename,world),.5,1);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"Out-of-range scalar variable was accepted");
        std::cout<<"Current MFA loading, scalar values, gradients, domain guards, and ellipsoid extraction passed"<<std::endl;
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<std::endl;return 1;}
}
