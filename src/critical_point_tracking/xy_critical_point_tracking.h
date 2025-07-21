#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

#include "particle_tracing.h"
#include "xy_critical_point_finding.h"
#include "CP_Trace.h"
#include "degenerate_case_tracing.h"

namespace xy_cp_tracking{

    template<typename T>
    bool tracing_single_cpt(T time_step, T spatial_step_size, const Block<T>* b, VectorX<T>& initial, std::vector<VectorX<T>>& result, int correction_max_itr, 
    T hessian_det_epsilon, T gradient_epsilon, T d_max_square)
    {
        result.clear();
        
        particle_tracing::tracing_one_direction(time_step,spatial_step_size,b,initial,result,false,hessian_det_epsilon,gradient_epsilon,correction_max_itr,d_max_square,
            b->core_mins, b->core_maxs);
        std::reverse(result.begin(),result.end());


        std::vector<VectorX<T>> temp_result;
        particle_tracing::tracing_one_direction(time_step,spatial_step_size,b,initial,temp_result,true,hessian_det_epsilon,gradient_epsilon,correction_max_itr,d_max_square, b->core_mins, b->core_maxs);
        if(temp_result.size()>1)
        {
            result.pop_back();
            result.insert(result.end(),temp_result.begin(),temp_result.end());
        }

        return true;
    }




    template<typename T>
    bool check_duplication(CP_Trace<T>& trace_in_span, std::vector<VectorX<T>>& result, T threshold_square)
    {
        
        // std::cout<<"check duplication "<<threshold_square<<std::endl;
        if(result.empty())
        {
            return true;
        }
        for(auto i=0;i<trace_in_span.traces.size();++i)
        {

            if((result[0]-trace_in_span.traces[i][0]).squaredNorm()<threshold_square)
            {
                if((result.back()-trace_in_span.traces[i].back()).squaredNorm()<threshold_square)
                {                            
                    return true;
                }
            }                                                
            
        }
            
        
        return false;
    }



    template<typename T>
    void find_trace(T time_step, T spatial_step_size,int max_step, const Block<T>* b, std::vector<VectorX<T>>& initial,
    std::vector<CP_Trace<T>>& traces, 
    T hessian_det_epsilon, T gradient_epsilon,int correction_max_itr)
    {
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);
        auto domain_range = b->core_maxs - b->core_mins;


        int distance_stop_itr = 5;
        T span_size= domain_range.cwiseQuotient(span_num.cast<T>()).head(span_num.size()-1).squaredNorm();
        T d_max_square = distance_stop_itr*distance_stop_itr*span_size;

        tbb::affinity_partitioner ap;
  

        tbb::parallel_for(tbb::blocked_range<size_t>(0,initial.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i = r.begin(); i != r.end(); ++i)
            {  

                tracing_single_cpt(time_step,spatial_step_size,b,initial[i],traces[i].traces,correction_max_itr,hessian_det_epsilon,gradient_epsilon,d_max_square);           
            }
        },ap);



    }



}