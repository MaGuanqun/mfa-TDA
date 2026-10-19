#pragma once
#include "surface_field.h"
#include "mfa_model.h"

namespace marching_triangles {
// Serial evaluator: reusable MFA decoder buffers avoid allocating a decoder for
// every value/derivative query. A field owns its model through shared_ptr.
class MfaEvaluator {
public:
    MfaEvaluator(std::shared_ptr<MfaModel> model,int variable)
        : model_(std::move(model)),data_(checked_data(model_->block(),variable)),
          decoder_(data_,0),scratch_(decoder_),param_(3),value_(1),gradient_(3) {
        scratch_.ResizeDers(1);
        const auto& b=model_->block();
        const auto& geom=b.mfa->geom();
        if(geom.dim()!=3 || geom.tmesh.tensor_prods.size()!=1)
            throw std::invalid_argument("MFA isosurfacing requires a 3D Cartesian geometry model");
        const auto& tensor=geom.tmesh.tensor_prods[0];
        if(!tensor.weights.allFinite() || (tensor.weights.array()-1).abs().maxCoeff()>1e-12)
            throw std::invalid_argument("Weighted MFA geometry is not supported");
        Eigen::VectorXd zero=Eigen::VectorXd::Zero(3),one=Eigen::VectorXd::Ones(3),lo(3),hi(3);
        b.mfa->DecodeGeom(zero,lo);b.mfa->DecodeGeom(one,hi);
        origin_=lo;extent_=hi-lo;
        if(!origin_.allFinite() || !extent_.allFinite() || (extent_.array()<=0).any())
            throw std::invalid_argument("Invalid MFA geometry bounds");
        const double tolerance=1e-9*std::max(1.,extent_.maxCoeff());
        // An affine Cartesian map has control points at the Greville abscissae.
        // Check all controls instead of assuming parametric and physical axes agree.
        for(Eigen::Index i=0;i<tensor.ctrl_pts.rows();++i) {
            Eigen::Index index=i;Point<double> expected;
            for(int d=0;d<3;++d) {
                const int j=index%tensor.nctrl_pts[d];index/=tensor.nctrl_pts[d];
                if(geom.p[d]<1)throw std::invalid_argument("MFA geometry degree must be positive");
                double u=0;for(int k=1;k<=geom.p[d];++k)u+=geom.tmesh.all_knots[d][j+k];
                expected[d]=origin_[d]+extent_[d]*u/geom.p[d];
            }
            if((tensor.ctrl_pts.row(i).transpose()-expected).norm()>tolerance)
                throw std::invalid_argument("Non-Cartesian MFA geometry requires an inverse geometry map and is not supported");
        }
        if(b.core_mins.size()!=3 || b.core_maxs.size()!=3 || !b.core_mins.allFinite() || !b.core_maxs.allFinite() ||
           (b.core_maxs.array()<=b.core_mins.array()).any() ||
           (b.core_mins.array()<origin_.array()-tolerance).any() ||
           (b.core_maxs.array()>(origin_+extent_).array()+tolerance).any())
            throw std::invalid_argument("MFA core bounds are outside the geometry domain");
        std::cout<<"MFA domain: "<<b.core_mins.transpose()<<" -> "<<b.core_maxs.transpose()
                 <<"; variable="<<variable<<" degree="<<data_.p.transpose()
                 <<" controls="<<data_.tmesh.tensor_prods[0].nctrl_pts.transpose()<<std::endl;
    }
    double value(const Point<double>& p) {
        if(!parameter(p))return std::numeric_limits<double>::quiet_NaN();
        decoder_.FastVolPt(param_,value_,scratch_,data_.tmesh.tensor_prods[0]);
        return value_[0];
    }
    Point<double> gradient(const Point<double>& p) {
        if(!parameter(p))return Point<double>::Constant(std::numeric_limits<double>::quiet_NaN());
        decoder_.FastGrad(param_,scratch_,data_.tmesh.tensor_prods[0],gradient_);
        // MFA differentiates with respect to u in [0,1]. Apply dx/du scaling.
        return gradient_.cwiseQuotient(extent_);
    }
    // Critical-curve continuation needs a continuous Hessian. Checking knot
    // multiplicity is conservative: low-continuity splines whose coefficients
    // happen to cancel the discontinuity are not specially recognized.
    void enable_hessian() {
        if(hessian_scratch_) return;
        const auto& b=model_->block();
        for(int d=0;d<3;++d) {
            const auto& knots=data_.tmesh.all_knots[d];
            const double lo=(b.core_mins[d]-origin_[d])/extent_[d];
            const double hi=(b.core_maxs[d]-origin_[d])/extent_[d];
            const double endpoint_epsilon=64*std::numeric_limits<double>::epsilon();
            for(size_t i=0;i<knots.size();) {
                size_t j=i+1;while(j<knots.size() && knots[j]==knots[i])++j;
                if(knots[i]>std::max(0.,lo)+endpoint_epsilon && knots[i]<std::min(1.,hi)-endpoint_epsilon &&
                   data_.p[d]-static_cast<int>(j-i)<2)
                    throw std::invalid_argument("Critical curves require C2 MFA splines: axis "+std::to_string(d)+
                        ", interior knot "+std::to_string(knots[i])+", degree "+std::to_string(data_.p[d])+
                        ", multiplicity "+std::to_string(j-i));
                i=j;
            }
        }
        hessian_scratch_=std::make_unique<mfa::DecodeInfo<double>>(data_,Eigen::VectorXi::Constant(3,2));
    }
    Eigen::Matrix3d hessian(const Point<double>& p) {
        enable_hessian();
        if(!parameter(p))return Eigen::Matrix3d::Constant(std::numeric_limits<double>::quiet_NaN());
        Eigen::Matrix3d h;
        for(int i=0;i<3;++i)for(int j=i;j<3;++j) {
            Eigen::VectorXi deriv=Eigen::VectorXi::Zero(3);++deriv[i];++deriv[j];
            double entry=0;
            if((deriv.array()<=data_.p.array()).all()) {
                decoder_.VolPt(param_,value_,*hessian_scratch_,data_.tmesh.tensor_prods[0],deriv);
                entry=value_[0]/(extent_[i]*extent_[j]);
            }
            h(i,j)=h(j,i)=entry;
        }
        return h;
    }
    const Block<double>& block() const {return model_->block();}
private:
    static const mfa::MFA_Data<double>& checked_data(const Block<double>& b,int variable) {
        if(!b.mfa || b.dom_dim!=3 || variable<0 || variable>=b.mfa->nvars())
            throw std::invalid_argument("Invalid 3D MFA model or --variable index");
        const auto& data=b.mfa->var(variable);
        if(data.dim()!=1 || data.tmesh.tensor_prods.size()!=1)
            throw std::invalid_argument("MFA isosurfacing requires a scalar variable with one tensor product");
        const auto& t=data.tmesh.tensor_prods[0];
        if(!t.ctrl_pts.size() || !t.ctrl_pts.allFinite() || !t.weights.size() || !t.weights.allFinite() ||
           (t.weights.array()-1).abs().maxCoeff()>1e-12 || (data.p.array()<1).any())
            throw std::invalid_argument("MFA isosurfacing requires finite, unweighted splines of positive degree");
        return data;
    }
    bool parameter(const Point<double>& p) {
        if(!p.allFinite())return false;
        param_=(p-origin_).cwiseQuotient(extent_);
        const double eps=64*std::numeric_limits<double>::epsilon();
        if((param_.array() < -eps).any() || (param_.array() > 1+eps).any())return false;
        param_=param_.cwiseMax(0).cwiseMin(1);return true;
    }
    std::shared_ptr<MfaModel> model_;
    const mfa::MFA_Data<double>& data_;
    mfa::Decoder<double> decoder_;
    mfa::FastDecodeInfo<double> scratch_;
    std::unique_ptr<mfa::DecodeInfo<double>> hessian_scratch_;
    Point<double> origin_,extent_;
    Eigen::VectorXd param_,value_,gradient_;
};

inline FieldSetup<double> mfa_field(std::shared_ptr<MfaModel> model,double iso,int variable=0) {
    if(!std::isfinite(iso))throw std::invalid_argument("The isovalue must be finite");
    auto evaluator=std::make_shared<MfaEvaluator>(std::move(model),variable);
    FieldSetup<double> setup;
    setup.field.value=[evaluator](const Point<double>& p){return evaluator->value(p);};
    setup.field.gradient=[evaluator](const Point<double>& p){return evaluator->gradient(p);};
    setup.options.domain_min=evaluator->block().core_mins;
    setup.options.domain_max=evaluator->block().core_maxs;
    setup.options.iso_value=iso;
    setup.base_step=(setup.options.domain_max-setup.options.domain_min).minCoeff()/20.;
    setup.options.step=setup.base_step;setup.options.min_step=setup.base_step*.1;
    setup.options.vertex_tolerance=setup.base_step*1e-6;
    return setup;
}
} // namespace marching_triangles
