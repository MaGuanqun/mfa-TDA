#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>
#include <mfa/mfa.hpp>
#include "opts.h"
#include "block.hpp"
#include "../utility/utility_function.h"
#include "mfa_extend.h"
#include "closed_form_function.h"
#include "INRModel.h"

namespace query_function
{

    // we extend the query to entire space. If a point is outside the domain, use the info of the closest span.
    template<typename T>
    void query_function(
            const VectorX<T>&   p,                  // parameters of point to decode
            VectorX<T>&         out,const int function_type=0, const Block<T>* b=nullptr, const VectorXi&     deriv = VectorXi(), INRModel<T>* inr_model = nullptr)
    {
        switch (function_type)
        {
        case -1:
            // inr_model->query(p, out, deriv);
            std::cerr<<"INR model querying not implemented in this function"<<std::endl;
            break;
        case 0:
            mfa_extend::recover_mfa(b, p, out, deriv);
            break;
        case 1:
            closed_form_function::quartic_potential(p, out, deriv);
            break;
        case 2:
            closed_form_function::quartic_potential_2(p, out, deriv);
            break;
        case 3:
            closed_form_function::rotating_quartic_multiwell(p, out, deriv);
            break;
        case 4:
            closed_form_function::quartic_potential_3d(p, out, deriv);
            break;
        default:
            std::cerr<<"invalid function type"<<std::endl;
            exit(0);
            break;
        }
    }

}