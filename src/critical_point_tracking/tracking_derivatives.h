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
#include "query_function.h"

namespace tracking_derivatives
{
    
    int remap(int idx, int remove_idx) {
        return idx - (idx > remove_idx);
    }


    template<typename T>
    void compute_Hessian(VectorX<T>& p, MatrixX<T>& dev_f, int removed_dom, const int function_type=0, const Block<T>* b=nullptr)
    {

        int domain_dim = p.size()-1;
        dev_f.resize(domain_dim,domain_dim);
        
        int ori_domain_dim = p.size();

        VectorX<T> dev_f_vector(1);    

        // std::cout<<local_domain_range.transpose()<<std::endl;

        VectorXi deriv(p.size());

        if(removed_dom==domain_dim)
        {
            for(int i=0;i<domain_dim;i++)
            {
                for(int j=i;j<domain_dim;j++)
                {
                    deriv.setZero();
                    deriv[i]+=1;
                    deriv[j]+=1;

                    query_function::query_function(p, dev_f_vector, function_type, b, deriv);

                    
                    dev_f(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
                    dev_f(i,j) = dev_f_vector[0];
                }
            }
        }
        else
        {
            for(int i=0;i<domain_dim;i++)
            {
                for(int j=i;j<ori_domain_dim;j++)
                {
                    if(j==i && j==removed_dom)
                    {
                        continue;
                    }
                    deriv.setZero();
                    deriv[i]+=1;
                    deriv[j]+=1;

                    query_function::query_function(p, dev_f_vector, function_type, b, deriv);
                   
                    if(j==domain_dim)
                    {
                        dev_f(i,domain_dim-1) = dev_f_vector[0];
                        
                    }
                    else if(j==removed_dom)
                    {
                        dev_f(removed_dom,i) = dev_f_vector[0];
                    }
                    else
                    {
                        dev_f(i,remap(j,removed_dom)) = dev_f_vector[0];
                        dev_f(j,remap(i,removed_dom)) = dev_f_vector[0];
                    }                    
                }
            }
        }

    }

        //gradient exclude last dimension
    template<typename T>
    void compute_gradient(VectorX<T>& p, VectorX<T>& f,const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> f_vector(1);
        int domain_dim = p.size()-1;
        VectorXi deriv(p.size());
        f.resize(domain_dim);
        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            query_function::query_function(p, f_vector, function_type, b, deriv);
            
            f[i] = f_vector[0];
        }
    }


}