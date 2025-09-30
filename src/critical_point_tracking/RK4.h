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
#include "tracking_derivatives.h"

namespace RK4
{
    template<typename T>
    //compute dx/dt, dy/dt
    bool compute_gradient(VectorX<T>& p, VectorX<T>& gradient, const int function_type=0, const Block<T>* b=nullptr)
    {
        int domain_dim = p.size()-1;
        VectorXi deriv(p.size());
        VectorX<T> dev_f(domain_dim);
        MatrixX<T> hessian(domain_dim,domain_dim);
        VectorX<T> f_vector(1);

        tracking_derivatives::compute_Hessian(p, hessian, domain_dim, function_type, b);

        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]=1;
            deriv[2]=1;
            query_function::query_function(p, f_vector, function_type, b, deriv);
            dev_f[i] = f_vector[0];
        }

        Eigen::ColPivHouseholderQR<MatrixX<T>> qr(hessian);
        // if(qr.rank() < domain_dim)
        // {
        //     std::cout<<hessian<<std::endl;
        //     std::cout<<"Hessian is not full rank"<<std::endl;
        //     return false;
        // }
        gradient =-1.0*qr.solve(dev_f);

        return true;
    }

    template<typename T>
    bool compute_direction(VectorX<T>& p, VectorX<T>& direction, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> gradient;
        // std::cout<<"p in compute direction "<<  p.transpose() <<std::endl;
        if(!utility::In_Domain(p,core_mins,core_maxs))
        {

            return false;
        }

        if(!compute_gradient(p,gradient,function_type,b))
        {
            return false;
        }

        direction.head(p.size()-1) = gradient;
        direction[p.size()-1] = 1.0;

        return true;
    }

    template<typename T>
    bool RK4(VectorX<T>& p,VectorX<T>& result, T step_size, bool upper_search, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> k1(p.size());
        if(!compute_direction(p, k1,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        VectorX<T> k2(p.size());
        VectorX<T> p2;
        if(upper_search)
        {
            p2 = p+0.5*step_size*k1;
        }
        else
        {
            p2 = p-0.5*step_size*k1;
        }
        if(!compute_direction(p2, k2,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        VectorX<T> k3(p.size());
        VectorX<T> p3;
        if(upper_search)
        {
            p3 = p+0.5*step_size*k2;
        }
        else
        {
            p3 = p-0.5*step_size*k2;
        }

        if(!compute_direction(p3, k3,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        VectorX<T> k4(p.size());
        VectorX<T> p4;
        if(upper_search)
        {
            p4 = p+step_size*k3;
        }
        else
        {
            p4 = p-step_size*k3;
        }
        if(!compute_direction(p4, k4,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        VectorX<T> k = (k1+2*k2+2*k3+k4)/6.0;
        k[p.size()-1] = 1.0;
        if(upper_search)
        {
            result = p + step_size*k;
        }
        else
        {
            result = p - step_size*k;
        }
        if(!utility::In_Domain(result,core_mins,core_maxs))
        {
            return false;
        }

        return true;

    }

         // newton method with single initial_point
    template<typename T>
    bool correction_newton(VectorX<T>& input_point, int max_itr, T d_max_square, T root_finding_epsilon, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        int itr_num=0;

        MatrixX<T> dev_f;
        VectorX<T> f;
        tracking_derivatives::compute_gradient(input_point, f,function_type,b);

        if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon)
        {
            return false;
        }

        VectorX<T> p = input_point;

        T temp_rec;

        while(itr_num<max_itr)
        {
            tracking_derivatives::compute_Hessian(p, dev_f, input_point.size()-1, function_type, b);

            Eigen::ColPivHouseholderQR<MatrixX<T>> qr(dev_f);

            if(qr.rank() < dev_f.cols())
            {
                return false;
            }

            p.head(p.size()-1)-= qr.solve(f);                
            
            
            
            if(!utility::In_Domain(p,core_mins,core_maxs))
            {
                return false;
            }
            

            if((input_point-p).squaredNorm()>d_max_square)
            {
                return false;
            }

            tracking_derivatives::compute_gradient(p, f,function_type,b);


 
            if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon){             
                input_point = p;
                return true;
            }
            
            itr_num++;
        }

        return false;

    }

    template<typename T>
    void normalize_m(VectorX<T>& m)
    {
        T norm=m.head(m.size()-1).norm();
        m = m/norm;
    }

    template<typename T>
    // RKF45 with fixed spatial step size
    bool RK4_normalized_step(VectorX<T>& p,VectorX<T>& result, T step_size, bool upper_search,const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> k1(p.size());
        if(!compute_direction(p, k1,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        normalize_m(k1);
        
        VectorX<T> k2(p.size());
        VectorX<T> p2;
        if(upper_search)
        {
            p2 = p+0.5*step_size*k1;
        }
        else
        {
            p2 = p-0.5*step_size*k1;
        }
        if(!compute_direction(p2, k2,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        normalize_m(k2);

        VectorX<T> k3(p.size());
        VectorX<T> p3;
        if(upper_search)
        {
            p3 = p+0.5*step_size*k2;
        }
        else
        {
            p3 = p-0.5*step_size*k2;
        }
        if(!compute_direction(p3, k3,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        normalize_m(k3);

        VectorX<T> k4(p.size());
        VectorX<T> p4;
        if(upper_search)
        {
            p4 = p+step_size*k3;
        }
        else
        {
            p4 = p-step_size*k3;
        }
        if(!compute_direction(p4, k4,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        normalize_m(k4);

        VectorX<T> k = (k1+2*k2+2*k3+k4)/6.0;
        
        if(upper_search)
        {
            result = p + step_size*k;
        }
        else
        {
            result = p - step_size*k;
        }
        if(!utility::In_Domain(result,core_mins,core_maxs))
        {
            return false;
        }

        return true;

    }

    template<typename T>
    bool RK4_choose_direction(VectorX<T>& p,VectorX<T>& result, T time_step, T sptial_step_size,  T gradient_epsilon, bool upper_search,int max_itr, T d_max_square, bool first_fixed_time,const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        

        if(first_fixed_time)
        {
            if(RK4(p,result,time_step,upper_search,core_mins,core_maxs,function_type,b))
            {
            
                if((result.head(result.size()-1)-p.head(p.size()-1)).squaredNorm()>sptial_step_size*sptial_step_size)
                {
                if(RK4_normalized_step(p,result,sptial_step_size,upper_search,core_mins,core_maxs,function_type,b))
                {
                        correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                        return true;
                }
                }
                else
                {
                    correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                    return true;
                }
            }
            else
            {
                if(RK4_normalized_step(p,result,sptial_step_size,upper_search,core_mins,core_maxs,function_type,b))
                {
                    correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                    return true;
                }
            }
        }
        else
        {
            if(RK4_normalized_step(p,result,sptial_step_size,upper_search,core_mins,core_maxs,function_type,b))
            {
            
                if(std::abs(result[result.size()-1]-p[p.size()-1])>time_step)
                {
                    if(RK4(p,result,time_step,upper_search,core_mins,core_maxs,function_type,b))
                    {
                            correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                            return true;
                    }
                }
                else
                {
                    correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                    return true;
                }
            }
            else
            {
                if(RK4(p,result,time_step,upper_search,core_mins,core_maxs,function_type,b))
                {
                    correction_newton(result, max_itr, d_max_square, gradient_epsilon, core_mins,core_maxs,function_type,b);
                    return true;
                }
            }

        }
        return false;
    }

    // determine we should fix time step or spatial step size
    template<typename T>
    bool determine_fixed_space_time(VectorX<T>& p, T time_step, T sptial_step_size,  bool& fixed_time,const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> m(p.size());
        if(!compute_direction(p, m,core_mins,core_maxs,function_type,b))
        {
            return false;
        }
        if(m.head(m.size()-1).squaredNorm() < sptial_step_size*sptial_step_size/(time_step*time_step))
        {
            fixed_time = true;
        }
        else
        {
            fixed_time = false;
        }

        return true;
    }

    template<typename T>
    bool RK4_correction(VectorX<T>& p,VectorX<T>& result, T time_step, T sptial_step_size,  T gradient_epsilon, bool upper_search,int max_itr, T d_max_square, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        bool fixed_time;
        if(!determine_fixed_space_time(p,time_step,sptial_step_size,fixed_time,core_mins,core_maxs,function_type,b))
        {
            return false;
        }

        // std::cout<<"determined fixed time "<<fixed_time<<std::endl;

        if(RK4_choose_direction(p,result,time_step,sptial_step_size,gradient_epsilon,upper_search,max_itr,d_max_square,fixed_time,core_mins,core_maxs,function_type,b))
        {
            return true;
        }

        return false;
    } 



}