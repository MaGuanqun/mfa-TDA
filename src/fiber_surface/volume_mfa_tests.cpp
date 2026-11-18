#include "mfa_field.h"
#include "volume_examples.h"
#include <filesystem>
using namespace fiber_surface;
static void require(bool ok,const std::string& why){if(!ok)throw std::runtime_error(why);}
// Exact tensor-product Bernstein coefficients, written with the current MFA
// serializer. This tests the adapter on the same field, without fitting error.
static mfa::MFA_Data<double>* polynomial(const Box& box,int component,bool curved) {
    const bool geometry=component<0;const int degree=geometry?1:2,n=degree+1,count=n*n*n;
    auto* data=new mfa::MFA_Data<double>(Eigen::Vector3i::Constant(degree),size_t(1),geometry?0:3+component,geometry?2:3+component);
    auto& t=data->tmesh.tensor_prods[0];t.nctrl_pts=Eigen::Vector3i::Constant(n);
    t.ctrl_pts.resize(count,geometry?3:1);t.weights=Eigen::VectorXd::Ones(count);
    t.knot_mins.assign(3,0);t.knot_maxs.assign(3,2*n-1);t.level=0;
    data->tmesh.all_knots.resize(3);data->tmesh.all_knot_levels.resize(3);
    for(int d=0;d<3;++d){data->tmesh.all_knots[d].assign(n,0);data->tmesh.all_knots[d].insert(data->tmesh.all_knots[d].end(),n,1);data->tmesh.all_knot_levels[d].assign(2*n,0);}
    const Vec extent=box.max-box.min;
    for(int i=0;i<count;++i) {
        int index=i;Eigen::Vector3i k;for(int d=0;d<3;++d){k[d]=index%n;index/=n;}
        const Vec u=k.cast<double>()/degree,x=box.min+extent.cwiseProduct(u);
        if(geometry)t.ctrl_pts.row(i)=x.transpose();
        else {
            const double z2=box.min[2]*box.min[2]+2*box.min[2]*extent[2]*u[2]+extent[2]*extent[2]*(k[2]==2?1.:0.);
            t.ctrl_pts(i,0)=x[component]+(curved?.5*z2:0.);
        }
    }
    data->tmesh.tensor_knot_idxs(t);return data;
}
int main(int argc,char** argv) {
    diy::mpi::environment environment(argc,argv);diy::mpi::communicator world;
    try {
        for(const std::string name:{"prism","curved-prism"}) {
            const auto e=volume_example(name);const bool curved=name=="curved-prism";
            const auto dir=std::filesystem::current_path()/"fiber_volume_mfa_fixture";std::filesystem::create_directories(dir);
            const auto file=(dir/(name+".mfa")).string();
            {
                diy::Master writer(world,1,-1,&Block<double>::create,&Block<double>::destroy,nullptr,&Block<double>::save,&Block<double>::load);
                auto* b=new Block<double>();b->dom_dim=3;b->pt_dim=5;b->save_datasets=0;
                b->bounds_mins=b->core_mins=e.field.box.min;b->bounds_maxs=b->core_maxs=e.field.box.max;
                std::vector<mfa::MFA_Data<double>*> vars{polynomial(e.field.box,0,curved),polynomial(e.field.box,1,curved)};
                b->mfa=new mfa::MFA<double>(3,0,polynomial(e.field.box,-1,false),vars);
                writer.add(0,b,new diy::RegularLink<diy::Bounds<double>>);diy::io::write_blocks(file,world,writer,&Block<double>::save);
            }
            auto model=std::make_shared<marching_triangles::MfaModel>(file,world);const auto f=mfa_field(model,model,0,1);
            for(const Vec& x:discovery_grid(f.box,5)) {
                require((f.value(x)-e.field.value(x)).norm()<1e-11,"Polynomial MFA value mismatch");
                require((f.jacobian(x)-e.field.jacobian(x)).norm()<1e-11,"Polynomial MFA Jacobian mismatch");
                const auto h=f.hessians(x),reference=e.field.hessians(x);
                require((h[0]-reference[0]).norm()<1e-11 && (h[1]-reference[1]).norm()<1e-11,"Polynomial MFA Hessian mismatch");
            }
            Options o;if(curved)o.residual_tolerance=.005;
            auto result=extract_volume(f,e.triangle,o);write_volume((dir/name).string(),result,f,e.triangle,o);
            require(result.complete(),"Incomplete MFA "+name);
            require(std::abs(result.report.volume-e.exact_volume)<(curved?.002:1e-9),"MFA volume differs from exact reference");
            std::cout<<"MFA "<<name<<" volume "<<result.report.volume<<'\n';
        }
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
