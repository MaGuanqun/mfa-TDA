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
#include "mfa_extend.h"

namespace RK4
{
    template<typename T>
    //compute dx/dt, dy/dt
    bool compute_gradient(const Block<T>* b, VectorX<T>& p, VectorX<T>& gradient, T hessian_det_epsilon)
    {
        int domain_dim = b->dom_dim-1;
        VectorXi deriv(b->dom_dim);
        VectorX<T> dev_f(domain_dim);
        MatrixX<T> hessian(domain_dim,domain_dim);
        VectorX<T> f_vector(1);

        find_boundary_roots::compute_Hessian(b, p, hessian, b->dom_dim-1);

        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]=1;
            deriv[2]=1;
            mfa_extend::recover_mfa(b,p,f_vector,deriv);
            dev_f[i] = f_vector[0];
        }

        T determinant;
        if(domain_dim==2){
            determinant = hessian.data()[3]*hessian.data()[0]-hessian.data()[1]*hessian.data()[2];
        }
        else{
            determinant = hessian.determinant();
        }
        if(std::abs(determinant) < hessian_det_epsilon)
        {
            return false;                
        }

        if(domain_dim==2){
            MatrixX<T> inv(domain_dim,domain_dim);
            inv.data()[0]=hessian.data()[3];
            inv.data()[1]=-hessian.data()[1];
            inv.data()[2]=-hessian.data()[2];
            inv.data()[3]=hessian.data()[0];           
            inv/=determinant;
            gradient = -inv*dev_f;
        }
        else{
            gradient =-1.0*hessian.colPivHouseholderQr().solve(dev_f);
        }

        return true;
    }

    template<typename T>
    bool compute_direction(const Block<T>* b, VectorX<T>& p, VectorX<T>& direction, T hessian_det_epsilon)
    {
        VectorX<T> gradient;
        if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
        {
            return false;
        }

        if(!compute_gradient(b,p,gradient,hessian_det_epsilon))
        {
            return false;
        }

        direction.head(b->dom_dim-1) = gradient;
        direction[b->dom_dim-1] = 1.0;

        return true;
    }

    template<typename T>
    bool RK4(const Block<T>* b, VectorX<T>& p,VectorX<T>& result, T step_size, T hessian_det_epsilon, bool upper_search)
    {
        VectorX<T> k1(b->dom_dim);
        if(!compute_direction(b,p, k1,hessian_det_epsilon))
        {
            return false;
        }
        VectorX<T> k2(b->dom_dim);
        VectorX<T> p2;
        if(upper_search)
        {
            p2 = p+0.5*step_size*k1;
        }
        else
        {
            p2 = p-0.5*step_size*k1;
        }
        if(!compute_direction(b,p2, k2,hessian_det_epsilon))
        {
            return false;
        }
        VectorX<T> k3(b->dom_dim);
        VectorX<T> p3;
        if(upper_search)
        {
            p3 = p+0.5*step_size*k2;
        }
        else
        {
            p3 = p-0.5*step_size*k2;
        }
        if(!compute_direction(b,p3, k3,hessian_det_epsilon))
        {
            return false;
        }
        VectorX<T> k4(b->dom_dim);
        VectorX<T> p4;
        if(upper_search)
        {
            p4 = p+step_size*k3;
        }
        else
        {
            p4 = p-step_size*k3;
        }
        if(!compute_direction(b,p4, k4,hessian_det_epsilon))
        {
            return false;
        }
        VectorX<T> k = (k1+2*k2+2*k3+k4)/6.0;
        k[b->dom_dim-1] = 1.0;
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

         // newton method with single initial_point
    template<typename T>
    bool correction_newton(const Block<T>* b, VectorX<T>& input_point, int max_itr, T d_max_square, T root_finding_epsilon, T hessian_det_epsilon)
    {
        int itr_num=0;

        MatrixX<T> dev_f;
        VectorX<T> f;
        find_boundary_roots::compute_gradient(b, input_point, f);

        if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon)
        {
            return false;
        }

        VectorX<T> p = input_point;

        T temp_rec;

        while(itr_num<max_itr)
        {
            find_boundary_roots::compute_Hessian(b, p, dev_f, b->dom_dim-1);

            T determinant = dev_f.determinant();

            if(std::abs(determinant) < hessian_det_epsilon)
            {
                return false;                
            }
            else
            {  
                VectorX<T> tem;
                if(p.size()==3)
                {               
                        MatrixX<T> inv(2,2);
                        inv.data()[0]=dev_f.data()[3];
                        inv.data()[1]=-dev_f.data()[1];
                        inv.data()[2]=-dev_f.data()[2];
                        inv.data()[3]=dev_f.data()[0];           
                        inv/=determinant;
                        tem= inv*f;
                }
                else
                {
                    tem = dev_f.colPivHouseholderQr().solve(f);
                }


                p.head(b->dom_dim-1)-= tem;                
            }
            
            
            if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
            {
                return false;
            }
            

            if((input_point-p).squaredNorm()>d_max_square)
            {
                return false;
            }

            find_boundary_roots::compute_gradient(b, p, f);


 
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
    bool RK4_normalized_step(const Block<T>* b, VectorX<T>& p,VectorX<T>& result, T step_size, T hessian_det_epsilon, bool upper_search)
    {
        VectorX<T> k1(b->dom_dim);
        if(!compute_direction(b,p, k1,hessian_det_epsilon))
        {
            return false;
        }
        normalize_m(k1);
        
        VectorX<T> k2(b->dom_dim);
        VectorX<T> p2;
        if(upper_search)
        {
            p2 = p+0.5*step_size*k1;
        }
        else
        {
            p2 = p-0.5*step_size*k1;
        }
        if(!compute_direction(b,p2, k2,hessian_det_epsilon))
        {
            return false;
        }
        normalize_m(k2);

        VectorX<T> k3(b->dom_dim);
        VectorX<T> p3;
        if(upper_search)
        {
            p3 = p+0.5*step_size*k2;
        }
        else
        {
            p3 = p-0.5*step_size*k2;
        }
        if(!compute_direction(b,p3, k3,hessian_det_epsilon))
        {
            return false;
        }
        normalize_m(k3);

        VectorX<T> k4(b->dom_dim);
        VectorX<T> p4;
        if(upper_search)
        {
            p4 = p+step_size*k3;
        }
        else
        {
            p4 = p-step_size*k3;
        }
        if(!compute_direction(b,p4, k4,hessian_det_epsilon))
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
        if(!utility::In_Domain(result,b->core_mins,b->core_maxs))
        {
            return false;
        }

        return true;

    }

    template<typename T>
    bool RK4_correction(const Block<T>* b, VectorX<T>& p,VectorX<T>& result, T time_step, T sptial_step_size, T hessian_det_epsilon, T gradient_epsilon, bool upper_search,int max_itr, T d_max_square, bool first_fixed_time=true)
    {
        
        if(first_fixed_time)
        {
            if(RK4(b,p,result,time_step,hessian_det_epsilon,upper_search))
            {
            
                if((result.head(result.size()-1)-p.head(p.size()-1)).squaredNorm()>sptial_step_size*sptial_step_size)
                {
                if(RK4_normalized_step(b,p,result,sptial_step_size,hessian_det_epsilon,upper_search))
                {
                        // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                        return true;
                }
                }
                else
                {
                    // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                    return true;
                }
            }
            else
            {
                if(RK4_normalized_step(b,p,result,sptial_step_size,hessian_det_epsilon,upper_search))
                {
                    // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                    return true;
                }
            }
        }
        else
        {
            if(RK4_normalized_step(b,p,result,sptial_step_size,hessian_det_epsilon,upper_search))
            {
            
                if(std::abs(result[result.size()-1]-p[p.size()-1])>time_step)
                {
                    if(RK4(b,p,result,time_step,hessian_det_epsilon,upper_search))
                    {
                            // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                            return true;
                    }
                }
                else
                {
                    // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                    return true;
                }
            }
            else
            {
                if(RK4(b,p,result,time_step,hessian_det_epsilon,upper_search))
                {
                    // correction_newton(b, result, max_itr, d_max_square, gradient_epsilon, hessian_det_epsilon);
                    return true;
                }
            }

        }
        return false;
    }





}