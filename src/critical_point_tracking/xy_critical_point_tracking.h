#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

#include "utility_function.h"
#include "mfa_extend.h"
#include "CP_Trace.h"

namespace xy_cp_tracking{
    template<typename T>
    //compute dx/dt, dy/dt
    bool compute_gradient(const Block<T>* b, VectorX<T>& p, VectorX<T>& gradient, T hessian_det_epsilon)
    {
        int domain_dim = 2;
        VectorX<T> dev_f(domain_dim);
        MatrixX<T> hessian(domain_dim,domain_dim);

        VectorX<T> f_vector(1);
        VectorX<T> dev_f_vector(1);    

        // std::cout<<local_domain_range.transpose()<<std::endl;

        VectorXi deriv(3);

        for(int i=0;i<domain_dim;i++)
        {
            for(int j=i;j<domain_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                mfa_extend::recover_mfa(b,p,dev_f_vector,deriv);
                hessian(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
                hessian(i,j) = dev_f_vector[0];
            }
        }

        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]=1;
            deriv[2]=1;
            mfa_extend::recover_mfa(b,p,f_vector,deriv);
            dev_f[i] = f_vector[0];
        }

        T determinant = hessian.data()[3]*hessian.data()[0]-hessian.data()[1]*hessian.data()[2];
        if(std::abs(determinant) < hessian_det_epsilon)
        {
            return false;                
        }

        MatrixX<T> inv(domain_dim,domain_dim);
        inv.data()[0]=hessian.data()[3];
        inv.data()[1]=-hessian.data()[1];
        inv.data()[2]=-hessian.data()[2];
        inv.data()[3]=hessian.data()[0];           
        inv/=determinant;

        gradient = -inv*dev_f;
        return true;
    }

    template<typename T>
    bool compute_direction(const Block<T>* b, VectorX<T>& p, VectorX<T>& direction, T hessian_det_epsilon, T gradient_epsilon)
    {
        VectorX<T> gradient(2);
        if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
        {
            return false;
        }

        if(!compute_gradient(b,p,gradient,hessian_det_epsilon))
        {
            return false;
        }
        // T norm = gradient.norm();
        // if(norm<gradient_epsilon)
        // {
        //     return false;
        // }
        // gradient /= norm;
        direction.head(2) = gradient;
        direction[2] = 1.0;

        return true;
    }

    template<typename T>
    bool RKF45(const Block<T>* b, VectorX<T>& p,VectorX<T>& result, T step_size, T hessian_det_epsilon, T gradient_epsilon, bool upper_search)
    {
        VectorX<T> k1(3);
        if(!compute_direction(b,p, k1,hessian_det_epsilon,gradient_epsilon))
        {
            return false;
        }
        VectorX<T> k2(3);
        VectorX<T> p2;
        if(upper_search)
        {
            p2 = p+0.5*step_size*k1;
        }
        else
        {
            p2 = p-0.5*step_size*k1;
        }
        if(!compute_direction(b,p2, k2,hessian_det_epsilon,gradient_epsilon))
        {
            return false;
        }
        VectorX<T> k3(3);
        VectorX<T> p3;
        if(upper_search)
        {
            p3 = p+0.5*step_size*k2;
        }
        else
        {
            p3 = p-0.5*step_size*k2;
        }

        if(!compute_direction(b,p3, k3,hessian_det_epsilon,gradient_epsilon))
        {
            return false;
        }
        VectorX<T> k4(3);
        VectorX<T> p4;
        if(upper_search)
        {
            p4 = p+step_size*k3;
        }
        else
        {
            p4 = p-step_size*k3;
        }

        if(!compute_direction(b,p4, k4,hessian_det_epsilon,gradient_epsilon))
        {
            return false;
        }
        VectorX<T> k = (k1+2*k2+2*k3+k4)/6.0;
        k[2] = 1.0;
        if(upper_search)
        {
            result = p + step_size*k;
        }
        else
        {
            result = p - step_size*k;
        }
        if(!utility::In_Domain(result,b->core_mins,b->core_maxs))
        {
            return false;
        }
        return true;

    }

    template<typename T>
    bool tracing_one_direction(T step_size, const Block<T>* b, VectorX<T>& initial, std::vector<VectorX<T>>& result, bool upper_search, std::vector<std::vector<T>>& span_range, T threshold_correction,
    T hessian_det_epsilon, T gradient_epsilon)
    {
        result.clear();
        VectorX<T> p_new = VectorX<T>::Zero(initial.size());
        VectorX<T> p_old = initial;

        bool valid = RKF45(b,initial, p_new,step_size,hessian_det_epsilon,gradient_epsilon,upper_search);

        result.emplace_back(initial);

        bool work_on_tracing = true;

        if(!valid)
        {
            work_on_tracing = false;
            return true;
        }

        int step=1;


        while (utility::InBlock(span_range,p_new))
        {            
            result.emplace_back(p_new);
            p_old = p_new;

            valid = RKF45(b,p_old, p_new,step_size,hessian_det_epsilon,gradient_epsilon,upper_search);

            if(!valid)
            {
                work_on_tracing = false;
                break;
            }

            step++;
        }
            
        return true;
    }

    template<typename T>
    bool tracing_single_cpt(T step_size, const Block<T>* b, VectorX<T>& initial, std::vector<VectorX<T>>& result, std::vector<std::vector<T>>& span_range, T threshold_correction, int correction_max_itr, 
    T hessian_det_epsilon, T gradient_epsilon)
    {
        result.clear();
        
        tracing_one_direction(step_size,b,initial,result,false,span_range,threshold_correction,hessian_det_epsilon,gradient_epsilon);
        std::reverse(result.begin(),result.end());

        std::vector<VectorX<T>> temp_result;
        tracing_one_direction(step_size,b,initial,temp_result,true,span_range,threshold_correction,hessian_det_epsilon,gradient_epsilon);
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
    void tracing_cpt(T step_size, const Block<T>* b, std::vector<VectorX<T>>& initial, CP_Trace<T>& trace, std::vector<std::vector<T>>& span_range, T threshold_correction, int correction_max_itr,
    T hessian_det_epsilon, T gradient_epsilon, T trace_threshold_square)
    {
         std::vector<VectorX<T>> result;

        for(int i =0;i<initial.size();i++)
        {
            tracing_single_cpt(step_size,b,initial[i],result,span_range,threshold_correction,correction_max_itr,hessian_det_epsilon,gradient_epsilon);

            // check duplicaiton
            if(!check_duplication(trace, result, trace_threshold_square))
            {
                trace.traces.emplace_back(result);
            }


           
        }

    }


    template<typename T>
    void find_trace(T step_size,int max_step, const Block<T>* b, std::vector<std::vector<VectorX<T>>>& initial,
    std::vector<size_t>& span_index, std::vector<CP_Trace<T>>& traces, 
    T hessian_det_epsilon, T gradient_epsilon,T threshold_correction,int correction_max_itr,
    T trace_threshold_square)
    {
        auto& tc = b->mfa->var(0).tmesh.tensor_prods[0];
        VectorXi span_num = tc.nctrl_pts-b->mfa->var(0).p;

        VectorXi number_in_every_domain; //span
        utility::obtain_number_in_every_domain(span_num,number_in_every_domain);
        auto domain_range = b->core_maxs - b->core_mins;
        tbb::affinity_partitioner ap;
        tbb::parallel_for(tbb::blocked_range<size_t>(0,initial.size()),[&](const tbb::blocked_range<size_t>& r)
        {
            for(auto i = r.begin(); i != r.end(); ++i)
            {  
                VectorXi span_domain_index(b->mfa->var(0).p.size());
                //decode the span index in every dimension
                utility::obtainDomainIndex(span_index[i],span_domain_index,number_in_every_domain);

                std::vector<std::vector<T>> span_range(number_in_every_domain.size());
                for(int j=0;j<span_domain_index.size();++j)
                {    
                    span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[j][span_domain_index[j]+b->mfa->var(0).p[j]]*domain_range[j]+b->core_mins[j]);
                    span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[j][span_domain_index[j]+1+b->mfa->var(0).p[j]]*domain_range[j]+b->core_mins[j]);
                }

                tracing_cpt(step_size,b,initial[i],traces[i],span_range,threshold_correction,correction_max_itr,hessian_det_epsilon,gradient_epsilon,trace_threshold_square);                
            }
        },ap);
    }

}