// for a 3D dataset, find critical point on every xy plane.
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


namespace find_2d_roots
{

    // compute derivative in a single span. here the third dimension is fixed. So only compute derivative in x and y directions
    template<typename T>
    void compute_f_dev_f(const Block<T>* b, VectorX<T>& p, VectorX<T>& f, MatrixX<T>& dev_f)
    {
        int domain_dim = b->dom_dim-1;
        dev_f.resize(domain_dim,domain_dim);
        f.resize(domain_dim);

        VectorX<T> f_vector(1);
        VectorX<T> dev_f_vector(1);    

        // std::cout<<local_domain_range.transpose()<<std::endl;

        VectorXi deriv(b->dom_dim);

        for(int i=0;i<domain_dim;i++)
        {
            for(int j=i;j<domain_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                mfa_extend::recover_mfa(b, p,dev_f_vector,deriv);
                dev_f(j,i) = dev_f_vector[0];// / (local_domain_range[i]*local_domain_range[j]);
                dev_f(i,j) = dev_f_vector[0];
            }
        }

        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            mfa_extend::recover_mfa(b, p,f_vector, deriv);
            f[i] = f_vector[0];
        }

    }
    // newton method with single initial_point
    template<typename T>
    bool newton(const Block<T>* b,VectorX<T>& result, VectorX<T>& p, int max_itr, std::vector<std::vector<T>>& span_range,
                    T d_max_square, VectorX<T>& center,
                    T root_finding_epsilon, T hessian_det_epsilon)
    {
        int itr_num=0;


        MatrixX<T> dev_f;
        VectorX<T> f;
        compute_f_dev_f(b,p,f,dev_f);

        if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon)
        {
            result = p;
            return true;
        }

        result=p;

        T temp_rec;

        VectorX<T> pre_point = p;
        VectorX<T> intersection = p;

        while(itr_num<max_itr)
        {

            // if(dev_f.squaredNorm() <ROOT_FINDING_EPSILON*ROOT_FINDING_EPSILON)
            // {
            //     return false;
            // }



            T determinant = dev_f.data()[3]*dev_f.data()[0]-dev_f.data()[1]*dev_f.data()[2];
            if(std::abs(determinant) < hessian_det_epsilon)
            {
                return false;                
            }
            else
            {  
                MatrixX<T> inv(2,2);
                inv.data()[0]=dev_f.data()[3];
                inv.data()[1]=-dev_f.data()[1];
                inv.data()[2]=-dev_f.data()[2];
                inv.data()[3]=dev_f.data()[0];           
                inv/=determinant;
                p.head(2) -= inv * f;
            }
            

            if((p.head(2)-center).squaredNorm()>d_max_square)
            {
                return false;
            }

            if(!utility::InBlock(span_range,p))
            {
                return false;
            }


            compute_f_dev_f(b,p,f,dev_f);


            if(itr_num>0){
                if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon){             
                    
                    result = p;
                    return true;
                }
            }

            // result = p;
            itr_num++;
        }

        return false;

    }

    template<typename T>
    bool newRoot(VectorX<T>& z, std::vector<VectorX<T>>& root_so_far, T threshold)
    {
        for(int i=0;i<root_so_far.size();++i)
        {
            if((z-root_so_far[i]).squaredNorm()<threshold*threshold)
            {
                return false;
            }
        }
        return true;
    }


    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
        T root_finding_grad_epsilon, T same_root_epsilon,
        int t_initial_number, T hessian_det_epsilon) { 

        root.clear();
        
        VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=100;

        int distance_stop_itr = 5;

        
        auto domain_range = b->core_maxs-b->core_mins;
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        std::vector<std::vector<T>> span_range(2);


        VectorX<T> center(2);
        for(int i=0;i<span_range.size();++i)
        {    
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);

            center[i]=(span_range[i][0]+span_range[i][1])*0.5;
        }   



        // compute distance to terminate iteration
        T d_max_square=0;
        for(auto i=span_range.begin();i!=span_range.end();++i)
        {  
            d_max_square+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
        }

        d_max_square*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2

        std::vector<std::vector<T>>initial_point;


        utility::compute_initial_points(initial_point,b->mfa->var(0).p,span_range);
        // only for first two dimensions

        std::vector<T> t_span_range;
        t_span_range.emplace_back(b->mfa->var(0).tmesh.all_knots[2][span_index[2]]*domain_range[2]+b->core_mins[2]);
        t_span_range.emplace_back(b->mfa->var(0).tmesh.all_knots[2][span_index[2]+1]*domain_range[2]+b->core_mins[2]);
        span_range.emplace_back(t_span_range);

        std::vector<T> initial_t;
        utility::compute_initial_points(initial_t,t_initial_number,span_range[2]);


        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point[i].size();
        }

        int num_initial_point = num_initial_point_every_domain.prod();

        VectorX<T> next_root; 

        // std::cout<<"initial_point "<<num_initial_point<<std::endl;

        // for(int i=0;i<initial_point.size();++i)
        // {
        //     std::cout<<initial_point[i]<<" ";
        // }
        // std::cout<<std::endl;
        // std::cout << "Press Enter to continue...";
        // std::cin.get(); // Waits for the user to press Enter

        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size()+1);
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);




        for(int g=0;g<initial_t.size();g++)
        {  
            std::vector<VectorX<T>> root_in_original_domain;

            for(int i=0;i<num_initial_point;++i)
            {

                utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
                for(int j=0;j<initial_point.size();j++)
                {
                    current_initial_point[j]=initial_point[j][domain_index[j]];
                }        
                current_initial_point[2]=initial_t[g];

                if(newton(b, next_root, current_initial_point,maxIter,span_range,d_max_square,center,root_finding_grad_epsilon,hessian_det_epsilon))
                {
            
                    if(newRoot(next_root,root_in_original_domain,same_root_epsilon))
                    {        
                        root_in_original_domain.emplace_back(next_root);
                        // std::cout<<next_root_in_original_domain.transpose()<<std::endl;
                        // std::cout<<"record+++++"<<std::endl;
                        // if(total_root_num>=deg)
                        // {
                        //     break;
                        // }
                    } 
                }

            }

            if(!root_in_original_domain.empty())
            {
                root.insert(root.end(),root_in_original_domain.begin(),root_in_original_domain.end());
            }

        }

        return !root.empty();

    }



    template<typename T>
    bool root_finding(Block<real_t>* block, std::vector<std::vector<VectorXi>>& span_index, 
    std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
        int current_index,
        T root_finding_epsilon, T same_root_epsilon, int t_initial_number, T hessian_det_epsilon) //2^n+1 initial points) 
    {

        for(auto i=0;i<block->mfa->nvars();++i)
        {
            if(root_finding(block,span_index[i][current_index], root,
            root_finding_epsilon,same_root_epsilon,t_initial_number,hessian_det_epsilon))
            {
                return true;
            }
        }

        return false;
    }

    template<typename T>
    void test_root_finding(Block<T>* b,std::vector<std::vector<Eigen::VectorX<T>>>& points,T root_finding_epsilon)
    {
        std::cout<<root_finding_epsilon<<std::endl;
        for(auto i=0;i<points.size();++i)
        {
            int domain_dim=2;
            VectorX<T>f(domain_dim);
            for(auto j=0;j<points[i].size();++j)
            {
                VectorX<T> p = points[i][j];
               
                VectorXi deriv(3);
                VectorX<T> f_vector(1);
                for(int k=0;k<domain_dim;k++)
                {
                    deriv.setZero();
                    deriv[k]+=1;
                    mfa_extend::recover_mfa(b, p,f_vector, deriv);
                    f[k] = f_vector[0];
                }
                if(f.squaredNorm() > root_finding_epsilon*root_finding_epsilon)
                {
                    std::cout<<"not root "<<f.transpose()<<" "<<f.squaredNorm()<<std::endl;
                }
            }
        }

    }

}