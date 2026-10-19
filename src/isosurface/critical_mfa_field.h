#pragma once
#include "critical_field.h"
#include "mfa_field.h"

namespace critical_curves {
inline Field mfa_field(std::shared_ptr<marching_triangles::MfaModel> model,int variable=0) {
    auto evaluator=std::make_shared<marching_triangles::MfaEvaluator>(std::move(model),variable);
    evaluator->enable_hessian();
    Field field;
    field.lower=evaluator->block().core_mins;field.upper=evaluator->block().core_maxs;
    field.value=[evaluator](const Point& p){return evaluator->value(p);};
    field.gradient=[evaluator](const Point& p){return evaluator->gradient(p);};
    field.hessian=[evaluator](const Point& p){return evaluator->hessian(p);};
    return field;
}
} // namespace critical_curves
