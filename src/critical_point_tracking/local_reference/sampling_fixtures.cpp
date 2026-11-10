// Exact linear fields in both model formats detect axis order, normalization,
// local-coordinate resets, and batch/file-layout mistakes in the sampler.
#include "sampling_grid.h"
#include "../../isosurface/mfa_model.h"
#include <torch/script.h>
#include <filesystem>

using namespace local_reference;
mfa::MFA_Data<double>* linear(const Dataset& d, bool geometry) {
    const Eigen::Vector3i degree=Eigen::Vector3i::Ones();
    auto* data=new mfa::MFA_Data<double>(degree,size_t(1),geometry?0:3,geometry?2:3);
    auto& tensor=data->tmesh.tensor_prods[0];
    tensor.nctrl_pts=Eigen::Vector3i::Constant(2);
    tensor.ctrl_pts.resize(8,geometry?3:1);
    tensor.weights=Eigen::VectorXd::Ones(8);
    tensor.knot_mins.assign(3,0); tensor.knot_maxs.assign(3,3); tensor.level=0;
    data->tmesh.all_knots.assign(3,{0,0,1,1});
    data->tmesh.all_knot_levels.assign(3,{0,0,0,0});
    for(int i=0;i<8;++i) {
        const double ux=i%2, uy=(i/2)%2, ut=i/4;
        if(geometry) {
            tensor.ctrl_pts(i,0)=d.min[0]+ux*(d.max[0]-d.min[0]);
            tensor.ctrl_pts(i,1)=d.min[1]+uy*(d.max[1]-d.min[1]);
            tensor.ctrl_pts(i,2)=d.min[2]+ut*(d.max[2]-d.min[2]);
        } else tensor.ctrl_pts(i,0)=5*(2*ux-1)+3*(2*uy-1)+2*(2*ut-1);
    }
    data->tmesh.tensor_knot_idxs(tensor);
    return data;
}

int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::invalid_argument("Expected fixture directory");
        const std::filesystem::path dir(argv[1]);
        std::filesystem::create_directories(dir);
        torch::jit::Module module("LinearField");
        module.define(R"JIT(
def forward(self, p: Tensor) -> Tensor:
    return (2.0 * p[:, 0] + 3.0 * p[:, 1] + 5.0 * p[:, 2]).unsqueeze(1)
)JIT");
        module.save((dir/"linear.pt").string());
        diy::mpi::environment environment(argc,argv);
        diy::mpi::communicator world;
        for(const auto& name:{"boussinesq_3d","vortex_street_3d"}) {
            const auto d=dataset(name);
            diy::Master writer(world,1,-1,&Block<double>::create,&Block<double>::destroy,
                               nullptr,&Block<double>::save,&Block<double>::load);
            auto* b=new Block<double>();
            b->dom_dim=3;b->pt_dim=4;b->save_datasets=0;
            Eigen::Vector3d lo,hi;
            for(int a=0;a<3;++a) {lo[a]=d.min[a];hi[a]=d.max[a];}
            b->bounds_mins=b->core_mins=lo;b->bounds_maxs=b->core_maxs=hi;
            std::vector<mfa::MFA_Data<double>*> variables{linear(d,false)};
            b->mfa=new mfa::MFA<double>(3,0,linear(d,true),variables);
            writer.add(0,b,new diy::RegularLink<diy::Bounds<double>>);
            diy::io::write_blocks((dir/(std::string(name)+".mfa")).string(),world,writer,&Block<double>::save);
        }
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
