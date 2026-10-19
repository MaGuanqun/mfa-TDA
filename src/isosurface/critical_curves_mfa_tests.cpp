#include "critical_mfa_field.h"
#include "critical_explicit_field.h"
#include "critical_curve_io.h"

using namespace critical_curves;
static void require(bool ok,const std::string& message) {if(!ok) throw std::runtime_error(message);}

// An exact cubic Bernstein representation of z*(x*x+y*y-1)+0.37,
// with unequal physical axis lengths. No fitting error masks adapter errors.
static mfa::MFA_Data<double>* polynomial(bool geometry) {
    const Point origin(-2,-3,-1),extent(4,6,2);
    const int degree=geometry?1:3,n=degree+1;
    auto* data=new mfa::MFA_Data<double>(Eigen::VectorXi::Constant(3,degree),size_t(1),geometry?0:3,geometry?2:3);
    auto& t=data->tmesh.tensor_prods[0];
    t.nctrl_pts=Eigen::VectorXi::Constant(3,n);
    t.ctrl_pts.resize(n*n*n,geometry?3:1);t.weights=Eigen::VectorXd::Ones(n*n*n);
    t.knot_mins.assign(3,0);t.knot_maxs.assign(3,2*n-1);t.level=0;
    data->tmesh.all_knots.resize(3);data->tmesh.all_knot_levels.resize(3);
    for(int d=0;d<3;++d) {
        data->tmesh.all_knots[d].assign(n,0);data->tmesh.all_knots[d].insert(data->tmesh.all_knots[d].end(),n,1);
        data->tmesh.all_knot_levels[d].assign(2*n,0);
    }
    for(int index=0;index<n*n*n;++index) {
        int remainder=index;Point linear,squared;
        for(int d=0;d<3;++d) {
            const int j=remainder%n;remainder/=n;
            linear[d]=origin[d]+extent[d]*j/degree;
            if(!geometry) squared[d]=origin[d]*origin[d]+2*origin[d]*extent[d]*j/degree+
                extent[d]*extent[d]*j*(j-1)/(degree*(degree-1));
        }
        if(geometry)t.ctrl_pts.row(index)=linear.transpose();
        else t.ctrl_pts(index,0)=linear.z()*(squared.x()+squared.y()-1)+.37;
    }
    data->tmesh.tensor_knot_idxs(t);
    return data;
}

int main(int argc,char** argv) {
    diy::mpi::environment env(argc,argv);diy::mpi::communicator world;
    try {
        const std::filesystem::path filename=argc>1?argv[1]:"critical_curve_test_files/circle.mfa";
        if(!filename.parent_path().empty())std::filesystem::create_directories(filename.parent_path());
        diy::Master writer(world,1,-1,&Block<double>::create,&Block<double>::destroy,nullptr,&Block<double>::save,&Block<double>::load);
        auto* block=new Block<double>();block->dom_dim=3;block->pt_dim=4;block->save_datasets=0;
        block->core_mins=block->bounds_mins=Point(-2,-3,-1);block->core_maxs=block->bounds_maxs=Point(2,3,1);
        std::vector<mfa::MFA_Data<double>*> vars{polynomial(false)};
        block->mfa=new mfa::MFA<double>(3,0,polynomial(true),vars);
        writer.add(0,block,new diy::RegularLink<diy::Bounds<double>>);
        diy::io::write_blocks(filename.string(),world,writer,&Block<double>::save);
        auto model=std::make_shared<marching_triangles::MfaModel>(filename.string(),world);
        const Field field=mfa_field(model);
        const Field exact=explicit_field("crossing_circle",.37);
        for(int i=0;i<7;++i)for(int j=0;j<7;++j)for(int k=0;k<7;++k) {
            const Point x=field.lower+(field.upper-field.lower).cwiseProduct(Point(i,j,k)/6.);
            require(std::abs(field.value(x)-exact.value(x))<1e-11,"MFA value mismatch");
            require((field.gradient(x)-exact.gradient(x)).norm()<1e-10,"MFA gradient mismatch");
            require((field.hessian(x)-exact.hessian(x)).norm()<1e-10,"MFA Hessian mismatch (including mixed terms and physical scaling)");
        }
        // Check analytic MFA Hessians against finite differences away from knots.
        const Point check(.31,-.47,.23);
        for(int d=0;d<3;++d) {
            const Point delta=1e-5*Point::Unit(d);
            const Point finite_difference=(field.gradient(check+delta)-field.gradient(check-delta))/(2e-5);
            require((finite_difference-field.hessian(check).col(d)).norm()<1e-8,"Hessian finite-difference mismatch");
        }
        Options options;options.samples=5;options.step=.06;
        Problem problem(field,options);
        const auto search=discover(problem); // no input isovalue
        const auto groups=group_values(problem,search);
        require(groups.size()==1 && std::abs(groups[0].value-.37)<1e-9,"MFA critical value discovery failed");
        const auto curves=trace(problem,groups[0].seeds,groups[0].value);
        require(curves.curves.size()==1 && curves.curves[0].closed && curves.events.empty(),"MFA closed curve tracing failed");
        for(const auto& p:curves.curves[0].points)
            require(std::abs(p.position.z())<1e-6 && std::abs(p.position.head<2>().norm()-1)<1e-6,"MFA curve geometry error");

        // Adding tangential variation to a line destroys its critical set.
        // This tests discovery's full-gradient residual independently of tracing.
        Problem tilted(explicit_field("tilted_line"),options);
        require(discover(tilted).seeds.empty(),"Discovery accepted nonzero tangential gradient");
        std::cout<<"MFA exact derivatives, mixed Hessians, physical scaling, automatic value discovery and circle tracing passed\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
