#pragma once
#include "field.h"
#include "../isosurface/mfa_field.h"
#include <tbb/enumerable_thread_specific.h>

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
    using Evaluator=marching_triangles::MfaEvaluator;
    using Pair=std::pair<std::shared_ptr<Evaluator>,std::shared_ptr<Evaluator>>;
    auto workers=std::make_shared<tbb::enumerable_thread_specific<Pair>>([f,g] {
        return Pair{std::make_shared<Evaluator>(*f),std::make_shared<Evaluator>(*g)};
    });
    field.value=[workers](const Vec& x) { auto& e=workers->local(); return Range(e.first->value(x),e.second->value(x)); };
    field.jacobian=[workers](const Vec& x) { auto& e=workers->local(); Jacobian j; j.row(0)=e.first->gradient(x).transpose(); j.row(1)=e.second->gradient(x).transpose(); return j; };
    field.hessians=[workers](const Vec& x) { auto& e=workers->local(); return std::array<Hessian,2>{e.first->hessian(x),e.second->hessian(x)}; };
    return field;
}
} // namespace fiber_surface
