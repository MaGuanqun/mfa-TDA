#pragma once
#include "field.h"
#include "../isosurface/mfa_field.h"

namespace fiber_surface {
inline Field mfa_field(std::shared_ptr<marching_triangles::MfaModel> first,
                       std::shared_ptr<marching_triangles::MfaModel> second, int variable1=0, int variable2=0) {
    auto f=std::make_shared<marching_triangles::MfaEvaluator>(std::move(first),variable1);
    auto g=std::make_shared<marching_triangles::MfaEvaluator>(std::move(second),variable2);
    f->enable_hessian(); g->enable_hessian();
    Field field;
    field.box.min=f->block().core_mins.cwiseMax(g->block().core_mins);
    field.box.max=f->block().core_maxs.cwiseMin(g->block().core_maxs);
    if ((field.box.max.array()<=field.box.min.array()).any())
        throw std::invalid_argument("The two MFA fields have no common physical box");
    field.value=[f,g](const Vec& x) { return Range(f->value(x),g->value(x)); };
    field.jacobian=[f,g](const Vec& x) { Jacobian j; j.row(0)=f->gradient(x).transpose(); j.row(1)=g->gradient(x).transpose(); return j; };
    field.hessians=[f,g](const Vec& x) { return std::array<Hessian,2>{f->hessian(x),g->hessian(x)}; };
    return field;
}
} // namespace fiber_surface
