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
#include "closed_form_function.h"
#include "degenerate_case.h"

namespace find_boundary_roots
{

    int remap(int idx, int remove_idx) {
        return idx - (idx > remove_idx);
    }
    
    //gradient exclude last dimension
    template<typename T>
    void compute_gradient(VectorX<T>& p, VectorX<T>& f, const int function_type=0, const Block<T>* b = nullptr)
    {
        VectorX<T> f_vector(1);
        int domain_dim = p.size()-1;
        VectorXi deriv(p.size());
        f.resize(domain_dim);
        for(int i=0;i<domain_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            switch (function_type)
            {
            case 0:
                mfa_extend::recover_mfa(b, p,f_vector, deriv);
                break;
            case 1:
                closed_form_function::closed_form_function(p, f_vector, function_type, deriv);
                break;
            default:
                break;
            }
            
            f[i] = f_vector[0];
        }
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
                    switch (function_type)
                    {
                    case 0:
                        mfa_extend::recover_mfa(b, p,dev_f_vector,deriv);
                        break;
                    case 1:
                        closed_form_function::quartic_potential(p, dev_f_vector, deriv);
                        break;
                    default:
                        break;
                    }
                    
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
                    switch (function_type)
                    {
                    case 0:
                        mfa_extend::recover_mfa(b, p,dev_f_vector,deriv);
                        break;
                    case 1:
                        closed_form_function::quartic_potential(p, dev_f_vector, deriv);
                        break;
                    default:
                        break;
                    }
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

    template<typename T>
    void compute_f_dev_f(VectorX<T>& p, VectorX<T>& f, MatrixX<T>& dev_f, int removed_dom, const int function_type=0, const Block<T>* b=nullptr)
    {
        compute_Hessian(p, dev_f, removed_dom, function_type, b);
        compute_gradient(p, f, function_type, b);
    }
    // newton method with single initial_point
    template<typename T>
    bool newton(VectorX<T>& result, VectorX<T>& p, int max_itr,
                    // T d_max_square, VectorX<T>& center,
                    T root_finding_epsilon, T hessian_det_epsilon, std::vector<int>& used_domain, int removed_dom, std::vector<T>& point_update_epsilon, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        int itr_num=0;

        VectorX<T> p_on_boundary(p.size()-1);
        for(int i=0;i<p_on_boundary.size();i++)
        {
            p_on_boundary[i]=p[used_domain[i]];
        }

        MatrixX<T> dev_f;
        VectorX<T> f;
        compute_f_dev_f(p,f,dev_f, removed_dom, function_type, b);

        if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon)
        {
            result = p;
            return true;
        }

        result=p;

        T temp_rec;

        VectorX<T> pre_point = p;
        while(itr_num<max_itr)
        {
            // T determinant = dev_f.determinant();

            pre_point=p;

            Eigen::ColPivHouseholderQR<MatrixX<T>> qr(dev_f);
            if(qr.rank() < dev_f.cols())
            {
                return false;
            }

            VectorX<T> tem = qr.solve(f);
            
            for(int i=0;i<tem.size();i++)
            {
                p[used_domain[i]]-= tem[i];
            }
            p_on_boundary -=tem;
            
            
            // if((p_on_boundary-center).squaredNorm()>d_max_square)
            // {
            //     return false;
            // }

             

            if(!utility::In_Domain(p,core_mins,core_maxs))
            {
                return false;
            }


            compute_f_dev_f(p,f,dev_f,removed_dom, function_type, b);


            if(itr_num>0){
                if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon
                && std::abs(p[p.size()-1]-pre_point[pre_point.size()-1])<point_update_epsilon.back()
                && (p.head(p.size()-1)-pre_point.head(pre_point.size()-1)).squaredNorm()<point_update_epsilon[0]*point_update_epsilon[0]
                ){           
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
    bool check_new_root(std::vector<T>& threshold, VectorX<T>& z, VectorX<T>& existed_point)
    {
        if((z.head(z.size()-1)-existed_point.head(z.size()-1)).squaredNorm()<threshold[0]*threshold[0])
        {
            if(std::abs(z[z.size()-1]-existed_point[existed_point.size()-1])<threshold.back())
            {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    bool newRoot(VectorX<T>& z, std::vector<VectorX<T>>& root_so_far, std::vector<T>& threshold)
    {
        for(int i=0;i<root_so_far.size();++i)
        {
            if(!check_new_root(threshold, z, root_so_far[i]))
            {
                return false;
            }
        }
        return true;
    }






//reset this using globle initial point
    template<typename T>
    void root_finding_on_one_boundary(std::vector<VectorX<T>>& root, T boundary_value, int boundary_dim_index,
        T root_finding_grad_epsilon, std::vector<T>& same_root_epsilon,
        T hessian_det_epsilon, std::vector<int>& used_dom,int maxIter, 
        std::vector<std::vector<T>>&initial_point, std::vector<std::array<int,2>>& initial_point_range,
        T point_itr_threshold, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr) // top plane is 1, bottom plane is -1 
    {

        VectorXi num_initial_point_every_domain(initial_point.size()-1);
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point_range[used_dom[i]][1]-initial_point_range[used_dom[i]][0];
        }
        int num_initial_point = num_initial_point_every_domain.prod();

        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size());
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);


        std::vector<VectorX<T>> root_in_original_domain;
        VectorX<T> next_root; 
        
        std::vector<T> point_update_epsilon(same_root_epsilon.size()-1);
        for(int i=0;i<point_update_epsilon.size();++i)
        {
            point_update_epsilon[i] =same_root_epsilon[used_dom[i]] *point_itr_threshold;
        }



        for(int i=0;i<num_initial_point;++i)
        {


            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<num_initial_point_every_domain.size();j++)
            {
                current_initial_point[used_dom[j]]=initial_point[used_dom[j]][initial_point_range[used_dom[j]][0]+domain_index[j]];
            }        
            current_initial_point[boundary_dim_index]=boundary_value;


            if(newton(next_root, current_initial_point,maxIter,root_finding_grad_epsilon,hessian_det_epsilon,used_dom,boundary_dim_index,point_update_epsilon,core_mins,core_maxs,function_type,b))
            {
        
                if(newRoot(next_root,root_in_original_domain,same_root_epsilon))
                {        
                    root_in_original_domain.emplace_back(next_root);
                } 
            }



        }

        if(!root_in_original_domain.empty())
        {
            root.insert(root.end(),root_in_original_domain.begin(),root_in_original_domain.end());
        }

    }


    // template<typename T>
    // void span_range_for_plane(std::vector<std::vector<T>>& span_range, const Block<T>* b, int skipeed_dim,VectorX<T>& center, VectorXi& span_index)
    // {
    //     VectorX<T> domain_range = b->core_maxs-b->core_mins;
    //     int j=0;
    //     for(int i=0;i<=span_range.size();++i)
    //     {
    //         if(i!=skipeed_dim)
    //         {
    //             span_range[j].clear();
    //             span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
    //             span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);
    //             center[j]=(span_range[j][0]+span_range[j][1])*0.5;
    //             j++;
    //         }
           
    //     }

    // }



    template<typename T>
    void root_finding_single_block(std::vector<VectorX<T>>& root, T root_finding_grad_epsilon, std::vector<T>& same_root_epsilon,
        T hessian_det_epsilon,int maxItr, T point_itr_threshold,const VectorX<T>& core_mins, const VectorX<T>& core_maxs,const VectorXi& point_num_in_block, const VectorXi& set_block_num,std::vector<vector<T>>& initial_points, const int function_type=0, const Block<T>* b=nullptr, const VectorXi& span_index= VectorXi())
        {
            root.clear();

            std::vector<T> fixed_value; //the plane of the fixed value
            std::vector<int> fixed_dim;
            std::vector<std::vector<int>> used_domain;
            std::vector<int> used_domain_temp;
            std::vector<std::vector<std::array<int,2>>> initial_point_range;

            for(int i=0;i<core_mins.size();++i)
            {
                if(span_index[i]==0)
                {
                    fixed_value.emplace_back(core_mins[i]);
                    fixed_dim.emplace_back(i);
                }
                if(i!=(core_mins.size()-1) && (span_index[i]==(set_block_num[i]-1)))
                {
                    fixed_value.emplace_back(core_maxs[i]);
                    fixed_dim.emplace_back(i);
                }
            }

            initial_point_range.resize(fixed_dim.size());

            for(int i=0;i<fixed_dim.size();++i)
            {
                used_domain_temp.clear();
                initial_point_range[i].resize(core_mins.size());
                for(int j=0;j<core_mins.size();++j)
                {
                    if(j!=fixed_dim[i])
                    {
                        used_domain_temp.emplace_back(j);
                    }
                    initial_point_range[i][j][0]= point_num_in_block[j]*span_index[j];
                    initial_point_range[i][j][1]= point_num_in_block[j]*(span_index[j]+1);
                }             
                used_domain.emplace_back(used_domain_temp);
            }


            for(int i=0;i<fixed_value.size();++i)
            {
                root_finding_on_one_boundary(root, fixed_value[i], fixed_dim[i], root_finding_grad_epsilon, same_root_epsilon, hessian_det_epsilon, used_domain[i],maxItr,initial_points,initial_point_range[i], point_itr_threshold, core_mins, core_maxs, function_type,b);
            }

        }




    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(std::vector<VectorXi>& span_index, std::vector<VectorX<T>>& root, T root_finding_grad_epsilon, std::vector<T>& same_root_epsilon,T hessian_det_epsilon,int maxItr, T point_itr_threshold, const VectorXi& set_block_num, VectorXi point_num_in_block, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {        
        // VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=100;

        int distance_stop_itr = 5;
        auto domain_range =core_maxs-core_mins;

        

     
        std::vector<vector<T>> initial_points;
        cp_tracking_degenerate_case::generate_initial_points(initial_points,core_mins,core_maxs,point_num_in_block,set_block_num,b);


        tbb::enumerable_thread_specific<std::vector<VectorX<T>>> local_root;
        tbb::affinity_partitioner ap;


        // tbb::parallel_for(tbb::blocked_range<size_t>(0,span_index.size()), //
        // [&](const tbb::blocked_range<size_t>& range)
        // {
        //     auto& root_thread = local_root.local();
            std::vector<VectorX<T>> block_root;

        //     for(int i=range.begin();i!=range.end();++i){

            for(int i=0;i<span_index.size();++i){
                
                root_finding_single_block(block_root,root_finding_grad_epsilon,same_root_epsilon,hessian_det_epsilon,maxItr,point_itr_threshold,core_mins,core_maxs,point_num_in_block,set_block_num,initial_points,function_type,b,span_index[i]);
                if(!block_root.empty())
                {
                    root.insert(root.end(), block_root.begin(), block_root.end());
                    // root_thread.insert(root_thread.end(), block_root.begin(), block_root.end());
                }
            }

        // },ap               
        // );

        // for (const auto& thread_vec : local_root) {
        //     root.insert(root.end(), thread_vec.begin(), thread_vec.end());
        // }

        return !root.empty();

    }



    template<typename T>
    void test_root_finding(Block<T>* b,std::vector<Eigen::VectorX<T>>& points,T root_finding_epsilon)
    {
        std::cout<<root_finding_epsilon<<std::endl;

        int domain_dim=2;
        VectorX<T>f(domain_dim);
        for(auto i=0;i<points.size();++i)
        {
            VectorX<T> p = points[i];
            
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

    // template<typename T>
    // void root_number(std::vector<std::vector<Eigen::VectorX<T>>>& points)
    // {
    //     int count=0;
    //     for(auto i=0;i<points.size();++i)
    //     {
    //         count+= points[i].size();
    //     }
    //     std::cout<<"total root number "<<count<<std::endl;
    // }

}