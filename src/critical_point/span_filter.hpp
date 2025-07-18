#pragma once

#include    <mfa/mfa.hpp>
#include    <mfa/block_base.hpp>


#include    <stdio.h>

#include    <Eigen/Dense>

#include    <random>
#include "utility_function.h"
#include <atomic>
#include <tbb/tbb.h>

using Index = MatrixXf::Index;


namespace span_filter
{

template <typename T>
bool check_cross_zero_from_sequence(T* value, size_t start_index, size_t windows_number)
{
    T max, min;
    size_t end = start_index + windows_number;

    max = value[start_index];
    min = value[start_index];
    for(size_t i = start_index + 1; i<end; ++i)
    {
        if (value[i] > max)
            max = value[i];
        if (value[i] < min)
            min = value[i];
    }

    if(max>=0 && min <=0)
    {
        return true;
    }

    return false;
}

//only for 1d domain
template <typename T>
void compute_valid_span(
    const MatrixX<T>&             control_points,
    Block<real_t>*              block,
    std::vector<int>&           valid_span
)
{
    int nvars           = block->mfa->nvars();   

    for (size_t i = 0; i < nvars; i++)                                      // science variables
    {
        auto& mfa_data          = block->mfa->var(i);
        auto& tmesh             = mfa_data.tmesh;
        auto& tensor_prods      = tmesh.tensor_prods;
        // auto& all_knots         = tmesh.all_knots;
        // auto& all_knot_levels   = tmesh.all_knot_levels;
        int ndom_dims           = block->dom_dim;  // number of geometry dims

        VectorX<T> value = control_points.col(control_points.cols()-1);

        for (auto t = 0; t < tensor_prods.size(); t++)                      // tensor products
        {
            auto& tc = mfa_data.tmesh.tensor_prods[t];

            for (auto k = 0; k < ndom_dims; k++)                            // domain dimensions
            {
                int degree = mfa_data.p(k)-1;
                int span_num =tc.nctrl_pts(k) -1- degree;

                for(int i=0;i<span_num;++i)
                {
                    if(check_cross_zero_from_sequence(value.data(),i,degree + 1))
                    {
                        valid_span.emplace_back(i+degree);
                    }

                }         

            }

        }
    }


}









template <typename T>
bool spanInRange(const mfa::MFA_Data<T>& mfa_data, VectorXi& span_index, std::vector<T>& interval)
{
    std::vector<std::vector<T>> span_range(span_index.size());
    for(int i=0;i<span_index.size();++i)
    {    
        span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]]);
        span_range[i].emplace_back(mfa_data.tmesh.all_knots[i][span_index[i]+1]);
    }   

    for(int i=0;i<span_index.size();++i)
    {
        if(span_range[i][1]<interval[2*i]){
            return false;
        }
        if(span_range[i][0]>interval[2*i+1]){
            return false;
        }
    }
    return true;
}




bool span_on_boundary(VectorXi& span_index, const VectorXi& ctrl_pts_num, const VectorXi& degree, bool exlucde_last_dim_max = false)
{
    for(int i=0;i<span_index.size();++i)
    {
        if(i==span_index.size()-1 && exlucde_last_dim_max)
        {
            if(span_index[i]==ctrl_pts_num[i]-1)
            {
                return false;
            }
        }
        if(span_index[i]==degree[i] || span_index[i]==ctrl_pts_num[i]-1)
        {
            return true;
        }
    }
    return false;
}


template <typename T> //exlucde_last_dim_max can determine if wnant to include boudary t=1
void compute_boundary_span(Block<T>*              block,std::vector<std::vector<VectorXi>>&           valid_span, bool exlucde_last_dim_max = false)
{
    int nvars           = block->mfa->nvars();   
    std::vector<std::vector<VectorXi>> boundary_span;
    boundary_span.resize(nvars);
    
    for (size_t i = 0; i < nvars; i++)                                      // science variables
    {
        auto& mfa_data          = block->mfa->var(i);
        const auto& tc = mfa_data.tmesh.tensor_prods[0];

        boundary_span[i].reserve(valid_span[i].size()/10);
        // int last_dim_max = tc.nctrl_pts[tc.nctrl_pts.size()-1]-mfa_data.p[mfa_data.p.size()-1]-1;
        for(int k = 0;k< valid_span[0].size();k++)
        {            
            if(span_on_boundary(valid_span[i][k],tc.nctrl_pts,mfa_data.p, exlucde_last_dim_max))
            {       
                boundary_span[i].emplace_back(valid_span[i][k]);       
            }            
        }
    }

    valid_span = boundary_span;

}


//for nd domain
template <typename T>
void compute_valid_span(
    std::vector<std::vector<MatrixX<T>>>&             control_points,
    Block<real_t>*              block,
    std::vector<std::vector<VectorXi>>&           valid_span, std::vector<T>& domain_limit, int limit_dimension_num = -1
)
{
    int nvars           = block->mfa->nvars();   
    valid_span.resize(nvars);

    if(limit_dimension_num == -1)
    {
        limit_dimension_num = block->dom_dim;
    }

    for (size_t i = 0; i < nvars; i++)                                      // science variables
    {
        auto& mfa_data          = block->mfa->var(i);
        auto& tmesh             = mfa_data.tmesh;
        auto& tensor_prods      = tmesh.tensor_prods;
        const auto& tc = tmesh.tensor_prods[0];
        int ndom_dims           = block->dom_dim;  // number of geometry dims
        // auto& all_knots         = tmesh.all_knots;
        // auto& all_knot_levels   = tmesh.all_knot_levels;

        std::vector<std::atomic<size_t>> count_in_original_function(tc.nctrl_pts.prod());//here,we project every block into the original funciton,
        //the common block means the number should be the same as the domain dimension
        // memset(count_in_original_function.data(),0,sizeof(size_t)*count_in_original_function.size());

        VectorXi number_in_every_domain_ori_func;
        utility::obtain_number_in_every_domain(tc.nctrl_pts,number_in_every_domain_ori_func);

        valid_span[i].reserve(tc.nctrl_pts.prod()/10);

        std::cout<<ndom_dims<<std::endl;


        int total_block=0;
    

        for(int k = 0;k< limit_dimension_num;k++)
        {            

            VectorXi p =mfa_data.p;
            p[k]-=1;

            VectorXi control_points_num = tc.nctrl_pts;
            control_points_num[k]-=1;

            VectorXi span_num = control_points_num-p;

            int spanned_block_num = span_num.prod();

            VectorXi number_in_every_domain; //span
            utility::obtain_number_in_every_domain(span_num,number_in_every_domain);

            VectorXi block_span_num = p+Eigen::VectorXi::Ones(p.size()); //every block
            VectorXi number_in_every_dmain_block;
            utility::obtain_number_in_every_domain(block_span_num,number_in_every_dmain_block);

            VectorXi number_in_every_domain_ori; //control_point 
            utility::obtain_number_in_every_domain(control_points_num,number_in_every_domain_ori);

            //loop over every block, 
            
            int block_span=block_span_num.prod();

            tbb::affinity_partitioner ap;
        
  

            tbb::parallel_for(tbb::blocked_range<size_t>(0,spanned_block_num),
            [&](const tbb::blocked_range<size_t>& range)
            {
                for (auto j = range.begin(); j!= range.end(); j++)
                {
                    VectorXi span_index(p.size());
                    //decode the span index in every dimension
                    utility::obtainDomainIndex(j,span_index,number_in_every_domain);


                    //get the index in control point for this block
                    std::vector<size_t>index_for_control_point(block_span);

                    utility::obtain_index_for_control_point(index_for_control_point,span_index,number_in_every_dmain_block,number_in_every_domain_ori);


                    // std::cout<<j<<" "<<index_for_control_point.size()<<std::endl;
                    // if(j==0)
                    // {
                    //     for(int l=0;l<index_for_control_point.size();++l)
                    //     {
                    //         std::cout<<index_for_control_point[l]<<" ";
                    //     }
                    //     std::cout<<std::endl;

                    // }

                    // if(span_index==test)
                    // {
                    //     std::cout<<"span_index "<<span_index.transpose()<<std::endl;
                    //     std::cout<<utility::check_valid_span(control_points[i][k].data(),index_for_control_point)<<std::endl;
                    // }
                    

                    if(utility::check_valid_span(control_points[i][k].data(),index_for_control_point))
                    {
                        // if(j==0)
                        // {
                        //     std::cout<<"find span true "<<std::endl;

                        // }

                        //project it to the original function span (fx->f)
                        span_index+=mfa_data.p;  //here the derivative is 
                        size_t ind = utility::obtain_index_from_domain_index(span_index,number_in_every_domain_ori_func);
                        count_in_original_function[ind].fetch_add(1);       
  
                    }
                }
                
            },ap
            );
            total_block = spanned_block_num;
            
        }

        VectorXi index(ndom_dims);

        for(int j=0;j<count_in_original_function.size();++j)
        {
            if(count_in_original_function[j].load()==limit_dimension_num)
            {
                utility::obtainDomainIndex(j,index,number_in_every_domain_ori_func);

                if(spanInRange(mfa_data,index,domain_limit))
                {
                    valid_span[i].emplace_back(index);
                }

            }
            // else if(count_in_original_function[j].load()>ndom_dims){
            //     std::cout<<"error "<<count_in_original_function[j].load()<<" "<<j<<std::endl;
            // }

        }


        std::cout<<"total span num "<<total_block<<" filtered span num "<< valid_span[i].size()<<std::endl;


    }
}

}
