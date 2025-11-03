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
#include "find_boundary_roots.h"
#include "CP_Trace.h"
#include "degenerate_case_tracing.h"


template<typename T>
class Boundary_critical_point_tracking{

private:

    Block<T>* b;
    const int function_type;
    const VectorX<T> core_mins;
    const VectorX<T> core_maxs;
    int correction_max_itr;
    T gradient_epsilon;
    T time_step;
    T spatial_step_size;
    INRModel<T>* inr_model=nullptr;
    
    
    bool tracing_single_cpt(VectorX<T>& initial, std::vector<VectorX<T>>& result,
    T d_max_square)
    {
        result.clear();

        std::vector<VectorX<T>> temp_result;
        particle_tracing::tracing_one_direction(time_step,spatial_step_size,initial,temp_result,true,gradient_epsilon,correction_max_itr,d_max_square, core_mins, core_maxs, core_mins, core_maxs,function_type,b,inr_model);
        // if(temp_result.size()>1)
        // {
        //     result.pop_back();
            result.insert(result.end(),temp_result.begin(),temp_result.end());
        // }

        return true;
    }




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



public:
    Boundary_critical_point_tracking(const VectorX<T> domain_min, const VectorX<T> domain_max, T gradient_epsi, T time_step_, T spatial_step_size_, int func_type=0, int correction_max_it=50, Block<T>* block=nullptr, INRModel<T>* inr_model_=nullptr): core_mins(domain_min), core_maxs(domain_max),  b(block), function_type(func_type),correction_max_itr(correction_max_it), gradient_epsilon(gradient_epsi), time_step(time_step_), spatial_step_size(spatial_step_size_), inr_model(inr_model_) {}

    ~Boundary_critical_point_tracking(){}



    void find_trace(std::vector<VectorX<T>>& initial,
    std::vector<CP_Trace<T>>& traces)
    {
        // auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        // VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        // VectorXi number_in_every_domain; //span
        // utility::obtain_number_in_every_domain(span_num,number_in_every_domain);
        // auto domain_range = b->core_maxs - b->core_mins;


        // int distance_stop_itr = 1;
        // T span_size= domain_range.cwiseQuotient(span_num.cast<T>()).head(span_num.size()-1).squaredNorm();
        T d_max_square = 25*spatial_step_size*spatial_step_size;

        tbb::affinity_partitioner ap;
  

        // VectorX<T> test_point(3);
        // test_point<<-1.9325, -2.0, 0.45603;

        tbb::parallel_for(tbb::blocked_range<size_t>(0,initial.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i = r.begin(); i != r.end(); ++i)
            {  
                // if((initial[i]-test_point).norm()>0.0001)
                // {
                //     continue;
                // }

                // std::cout<<"start tracing on certain point"<<i<<std::endl;

                // std::cout<<"tracing point "<<i<<" "<<initial[i].transpose()<<std::endl;

                tracing_single_cpt(initial[i],traces[i].traces,d_max_square);           
            }
        },ap);



    }



};