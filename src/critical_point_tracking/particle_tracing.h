#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

#include "utility_function.h"
#include "RK4.h"


namespace particle_tracing{
 template<typename T>
    bool tracing_one_direction(T time_step, T spatial_step_size, const Block<T>* b, VectorX<T>& initial, std::vector<VectorX<T>>& result, bool upper_search,
    T hessian_det_epsilon, T gradient_epsilon, int max_itr, T d_max_square, const VectorX<T>&         block_min,
    const VectorX<T>&         block_max)
    {
        result.clear();
        VectorX<T> p_new = VectorX<T>::Zero(initial.size());
        VectorX<T> p_old = initial;

        bool valid = RK4::RK4_correction(b,initial, p_new,time_step,spatial_step_size,hessian_det_epsilon,gradient_epsilon,upper_search,max_itr, d_max_square);

        result.emplace_back(initial);

        // std::cout<<"start rkf 45 correction "<<valid<<std::endl;

        if(!valid)
        {
            return true;
        }

        int step=1;

        while (utility::In_Domain(p_new, block_min,block_max))
        {            
            result.emplace_back(p_new);
            p_old = p_new;

            valid = RK4::RK4_correction(b,p_old, p_new,time_step,spatial_step_size,hessian_det_epsilon,gradient_epsilon,upper_search,max_itr, d_max_square);

            if(!valid)
            {
                break;
            }

            step++;
        }
            
        return true;
    }
}