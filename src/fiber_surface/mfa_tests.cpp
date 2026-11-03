#include "mfa_field.h"
#include "analytical_field.h"
#include "output.h"
#include "extract.h"
#include <filesystem>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

using namespace fiber_surface;
static void require(bool ok,const std::string& message) { if (!ok) throw std::runtime_error(message); }

// Exact polynomial B-spline coefficients, serialized by the current DIY writer.
// Unequal extents detect missing physical-coordinate derivative scaling.
static mfa::MFA_Data<double>* polynomial(int kind,const Vec& origin,const Vec& extent,
                                      Eigen::Vector3i degrees=Eigen::Vector3i::Constant(2),SeedShape spans={1,1,1},int first=-1) {
    if (kind==0) degrees.setOnes();
    const Eigen::Vector3i controls=degrees+Eigen::Vector3i(spans[0],spans[1],spans[2]);
    const int components=kind==0?3:1, total=controls.prod();
    if (first<0) first=kind==0?0:kind+2;
    auto* data=new mfa::MFA_Data<double>(degrees,size_t(1),first,first+components-1);
    auto& t=data->tmesh.tensor_prods[0];
    t.nctrl_pts=controls; t.ctrl_pts.resize(total,components); t.weights=Eigen::VectorXd::Ones(total);
    t.knot_mins.assign(3,0); t.knot_maxs.resize(3); t.level=0;
    data->tmesh.all_knots.resize(3); data->tmesh.all_knot_levels.resize(3);
    for (int d=0; d<3; ++d) {
        auto& knots=data->tmesh.all_knots[d];
        knots.assign(degrees[d]+1,0);
        for (int i=1; i<spans[d]; ++i) knots.push_back(double(i)/spans[d]);
        knots.insert(knots.end(),degrees[d]+1,1);
        data->tmesh.all_knot_levels[d].assign(knots.size(),0);
        t.knot_maxs[d]=int(knots.size())-1;
    }
    for (int i=0; i<total; ++i) {
        int index=i; Vec u; double square=0;
        for (int d=0; d<3; ++d) {
            const int j=index%controls[d]; index/=controls[d];
            double sum=0,squares=0;
            for (int k=1; k<=degrees[d]; ++k) { const double knot=data->tmesh.all_knots[d][j+k]; sum+=knot; squares+=knot*knot; }
            u[d]=sum/degrees[d];
            const double quadratic=degrees[d]>1 ? (sum*sum-squares)/(degrees[d]*(degrees[d]-1)) : 0;
            square+=origin[d]*origin[d]+2*origin[d]*extent[d]*u[d]+extent[d]*extent[d]*quadratic;
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
        auto* b=new Block<double>(); b->dom_dim=3; b->pt_dim=7; b->save_datasets=0;
        const Vec origin(-2,-3,-4),extent(4,6,8);
        b->bounds_mins=b->core_mins=origin; b->bounds_maxs=b->core_maxs=origin+extent;
        std::vector<mfa::MFA_Data<double>*> vars{polynomial(1,origin,extent),polynomial(2,origin,extent),
            polynomial(1,origin,extent,Eigen::Vector3i(3,4,5),{2,3,4},5),
            polynomial(2,origin,extent,Eigen::Vector3i(4,3,3),{3,2,5},6)};
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
        Options sampling;
        require(discovery_shape(field,sampling)==SeedShape{7,11,14},"Default MFA counts scale with physical lengths");
        auto multi=mfa_field(model,model,2,3);
        require(discovery_shape(multi,sampling)==discovery_shape(field,sampling),"MFA counts are independent of degree and span count");
        sampling.seed_grid=0;
        require(discovery_shape(multi,sampling)==SeedShape{7,11,14},"Zero selects default proportional MFA counts");
        multi.box.max=origin+Vec(4,2,4);
        require(discovery_shape(multi,sampling)==SeedShape{14,7,14},"The selected box determines the shortest axis");
        sampling.seed_grid=4;
        require(discovery_shape(multi,sampling)==SeedShape{8,4,8},"Explicit shortest-axis MFA count");
        auto analytical=analytical_field("sphere").field; analytical.box=multi.box;
        require(discovery_grid(multi.box,discovery_shape(multi,sampling))==
                discovery_grid(analytical.box,discovery_shape(analytical,sampling)),"MFA and analytical fields use identical starting grids");
        const auto samples=discovery_grid(field.box,SeedShape{13,14,15});
        tbb::global_control control(tbb::global_control::max_allowed_parallelism,4);
        tbb::parallel_for(size_t(0),samples.size(),[&](size_t i) {
            const Vec& x=samples[i];
            for (const Field* candidate:{&field,&multi}) {
                require((candidate->value(x)-Range(x.squaredNorm(),x[2])).norm()<1e-9,"Parallel MFA values");
                Jacobian j; j.row(0)=2*x.transpose(); j.row(1)=Vec::UnitZ();
                require((candidate->jacobian(x)-j).norm()<1e-9,"Parallel MFA gradients");
                const auto h=candidate->hessians(x);
                require((h[0]-2*Hessian::Identity()).norm()<1e-9 && h[1].norm()<1e-9,"Parallel MFA Hessians");
            }
        });
        bool rejected=false; try { mfa_field(model,model,0,4); } catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"Invalid MFA variable accepted");
        // Restrict to the same physical box as the analytical example.
        field.box.min=Vec::Constant(-2); field.box.max=Vec::Constant(2);
        Segment segment; Options options;
        Seeds serial;
        {
            tbb::global_control single(tbb::global_control::max_allowed_parallelism,1);
            serial=find_seeds(field,segment,options);
        }
        const auto parallel=find_seeds(field,segment,options);
        require(serial.surface==parallel.surface && serial.endpoints==parallel.endpoints,"MFA parallel seed ordering/values");
        auto result=Extractor(field,segment,options).run(); auto report=evaluate(result,field,segment,options);
        write_outputs((dir/"sphere").string(),result,report,field,segment,options);
        require(report.complete(result) && report.components==1 && report.topology.euler_characteristic==1,"MFA truncated sphere extraction");
        require(result.fibers.size()==1 && result.fibers.front().closed,"MFA circle boundary");
        std::cout<<"MFA values, physical gradients/Hessians, loading, and bounded surface extraction passed\n"; return 0;
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
