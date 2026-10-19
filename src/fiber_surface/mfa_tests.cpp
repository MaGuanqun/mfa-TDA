#include "mfa_field.h"
#include "output.h"
#include "extract.h"
#include <filesystem>

using namespace fiber_surface;
static void require(bool ok,const std::string& message) { if (!ok) throw std::runtime_error(message); }

// Exact Bernstein coefficients, serialized by the project's current DIY writer.
// Unequal extents detect missing physical-coordinate derivative scaling.
static mfa::MFA_Data<double>* polynomial(int kind,const Vec& origin,const Vec& extent) {
    const int degree=kind==0?1:2, n=degree+1, components=kind==0?3:1, first=kind==0?0:kind+2;
    auto* data=new mfa::MFA_Data<double>(Eigen::VectorXi::Constant(3,degree),size_t(1),first,first+components-1);
    auto& t=data->tmesh.tensor_prods[0];
    t.nctrl_pts=Eigen::VectorXi::Constant(3,n); t.ctrl_pts.resize(n*n*n,components); t.weights=Eigen::VectorXd::Ones(n*n*n);
    t.knot_mins.assign(3,0); t.knot_maxs.assign(3,2*n-1); t.level=0;
    data->tmesh.all_knots.resize(3); data->tmesh.all_knot_levels.resize(3);
    for (int d=0; d<3; ++d) {
        data->tmesh.all_knots[d].assign(n,0); data->tmesh.all_knots[d].insert(data->tmesh.all_knots[d].end(),n,1);
        data->tmesh.all_knot_levels[d].assign(2*n,0);
    }
    for (int i=0; i<n*n*n; ++i) {
        int index=i; Vec u; double square=0;
        for (int d=0; d<3; ++d) {
            const int j=index%n; index/=n; u[d]=double(j)/degree;
            square+=origin[d]*origin[d]+2*origin[d]*extent[d]*u[d]+extent[d]*extent[d]*j*(j-1)/2.;
        }
        if (kind==0) t.ctrl_pts.row(i)=(origin+extent.cwiseProduct(u)).transpose();
        else t.ctrl_pts(i,0)=kind==1 ? square : origin[2]+extent[2]*u[2];
    }
    data->tmesh.tensor_knot_idxs(t); return data;
}
int main(int argc,char** argv) {
    diy::mpi::environment environment(argc,argv); diy::mpi::communicator world;
    try {
        const auto dir=std::filesystem::current_path()/"fiber_surface_mfa_fixture";
        std::filesystem::create_directories(dir);
        diy::Master writer(world,1,-1,&Block<double>::create,&Block<double>::destroy,nullptr,&Block<double>::save,&Block<double>::load);
        auto* b=new Block<double>(); b->dom_dim=3; b->pt_dim=5; b->save_datasets=0;
        const Vec origin(-2,-3,-4),extent(4,6,8);
        b->bounds_mins=b->core_mins=origin; b->bounds_maxs=b->core_maxs=origin+extent;
        std::vector<mfa::MFA_Data<double>*> vars{polynomial(1,origin,extent),polynomial(2,origin,extent)};
        b->mfa=new mfa::MFA<double>(3,0,polynomial(0,origin,extent),vars);
        writer.add(0,b,new diy::RegularLink<diy::Bounds<double>>);
        const auto file=(dir/"sphere.mfa").string(); diy::io::write_blocks(file,world,writer,&Block<double>::save);
        auto model=std::make_shared<marching_triangles::MfaModel>(file,world);
        Field field=mfa_field(model,model,0,1);
        for (const auto& x:discovery_grid(field.box,4)) {
            require((field.value(x)-Range(x.squaredNorm(),x[2])).norm()<1e-10,"MFA values");
            Jacobian j; j.row(0)=2*x.transpose(); j.row(1)=Vec::UnitZ();
            require((field.jacobian(x)-j).norm()<1e-10,"Physical MFA gradients");
            const auto h=field.hessians(x);
            require((h[0]-2*Hessian::Identity()).norm()<1e-10 && h[1].norm()<1e-10,"Physical MFA Hessians");
        }
        bool rejected=false; try { mfa_field(model,model,0,2); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"Invalid MFA variable accepted");
        // Restrict to the same physical box as the analytical example.
        field.box.min=Vec::Constant(-2); field.box.max=Vec::Constant(2);
        Segment segment; Options options;
        auto result=Extractor(field,segment,options).run(); auto report=evaluate(result,field,segment,options);
        write_outputs((dir/"sphere").string(),result,report,field,segment,options);
        require(report.complete(result) && report.components==1 && report.topology.euler_characteristic==1,"MFA truncated sphere extraction");
        require(result.fibers.size()==1 && result.fibers.front().closed,"MFA circle boundary");
        std::cout<<"MFA values, physical gradients/Hessians, loading, and bounded surface extraction passed\n"; return 0;
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
