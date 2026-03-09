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

#include "tracking_derivatives.h"
#include "degenerate_case.h"
#include "tracking_utility.h"

template<typename T>
class Find_boundary_roots
{
private:
    Block<T>* b;
    const int function_type;
    const VectorX<T> core_mins;
    const VectorX<T> core_maxs;
    T root_finding_epsilon;
    int max_itr;
    T point_itr_threshold;
    const VectorXi set_block_num;
    std::vector<T> same_root_epsilon;
    INRModel<T>* inr_model=nullptr;


    void compute_f_dev_f(VectorX<T>& p, VectorX<T>& f, MatrixX<T>& dev_f, int removed_dom)
    {
        if(inr_model!=nullptr)
        {
            inr_model->query_dim_reduced_grad_hessian(p, f, dev_f, removed_dom);
            return;
        }
        tracking_derivatives::compute_Hessian(p, dev_f, removed_dom, function_type, b);
        tracking_derivatives::compute_gradient(p, f, function_type, b);
    }

    
    // newton method with single initial_point
    bool newton(VectorX<T>& result, VectorX<T>& p,
                    // T d_max_square, VectorX<T>& center,
                    std::vector<int>& used_domain, int removed_dom, std::vector<T>& point_update_epsilon)
    {
        int itr_num=0;

        VectorX<T> p_on_boundary(p.size()-1);
        for(int i=0;i<p_on_boundary.size();i++)
        {
            p_on_boundary[i]=p[used_domain[i]];
        }

        MatrixX<T> dev_f;
        VectorX<T> f;
        compute_f_dev_f(p,f,dev_f, removed_dom);

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


            compute_f_dev_f(p,f,dev_f,removed_dom);


            if(itr_num>0){
                if(f.squaredNorm()<root_finding_epsilon*root_finding_epsilon
                // && std::abs(p[p.size()-1]-pre_point[pre_point.size()-1])<point_update_epsilon.back()
                // && (p.head(p.size()-1)-pre_point.head(pre_point.size()-1)).squaredNorm()<point_update_epsilon[0]*point_update_epsilon[0]
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


public:

    VectorXi point_num_in_block;

    Find_boundary_roots(T root_finding_epsi, const VectorX<T> domain_min, const VectorX<T> domain_max, VectorXi point_num_in_b, const VectorXi set_block_n,std::vector<T> same_root_epsi, int func_type=0, int max_it=50, T point_itr_thres=0.5, Block<T>* block=nullptr, INRModel<T>* inr_model_=nullptr): root_finding_epsilon(root_finding_epsi), core_mins(domain_min), core_maxs(domain_max), point_num_in_block(point_num_in_b), b(block), function_type(func_type),max_itr(max_it), point_itr_threshold(point_itr_thres), set_block_num(set_block_n), same_root_epsilon(same_root_epsi),inr_model(inr_model_) {}
    
    ~Find_boundary_roots() {}




//reset this using globle initial point
    void root_finding_on_one_boundary(std::vector<VectorX<T>>& root, T boundary_value, int boundary_dim_index,
        std::vector<int>& used_dom,
        std::vector<std::vector<T>>&initial_point, std::vector<std::array<int,2>>& initial_point_range, std::vector<T>& same_root_epsilon_) // top plane is 1, bottom plane is -1 
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
        
        std::vector<T> point_update_epsilon(same_root_epsilon_.size()-1);
        for(int i=0;i<point_update_epsilon.size();++i)
        {
            point_update_epsilon[i] =same_root_epsilon_[used_dom[i]] *point_itr_threshold;
        }



        for(int i=0;i<num_initial_point;++i)
        {


            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<num_initial_point_every_domain.size();j++)
            {
                current_initial_point[used_dom[j]]=initial_point[used_dom[j]][initial_point_range[used_dom[j]][0]+domain_index[j]];
            }        
            current_initial_point[boundary_dim_index]=boundary_value;


            if(newton(next_root, current_initial_point,used_dom,boundary_dim_index,point_update_epsilon))
            {
        
                if(tracking_utility::newRoot(next_root,root_in_original_domain,same_root_epsilon_))
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




    void root_finding_single_block(std::vector<VectorX<T>>& root,
        std::vector<vector<T>>& initial_points,const VectorXi& span_index)
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
                if(span_index[i]==(set_block_num[i]-1))
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
                root_finding_on_one_boundary(root, fixed_value[i], fixed_dim[i], used_domain[i],initial_points,initial_point_range[i],same_root_epsilon);
            }

        }




    // Function to find the roots of the polynomial using Newton's method
    bool root_finding(std::vector<VectorXi>& span_index, std::vector<VectorX<T>>& root)
    {        
        // VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int distance_stop_itr = 5;
        auto domain_range =core_maxs-core_mins;

        

     
        std::vector<vector<T>> initial_points;
        tracking_utility::generate_initial_points(initial_points,core_mins,core_maxs,point_num_in_block,set_block_num,b);


        tbb::enumerable_thread_specific<std::vector<VectorX<T>>> local_root;
        tbb::affinity_partitioner ap;

        tbb::parallel_for(tbb::blocked_range<size_t>(0, span_index.size()),
        [&](const tbb::blocked_range<size_t>& range)
        {
            auto& root_thread = local_root.local();
            std::vector<VectorX<T>> block_root;
            for(size_t i = range.begin(); i != range.end(); ++i){
                root_finding_single_block(block_root, initial_points, span_index[i]);
                if(!block_root.empty())
                {
                    root_thread.insert(root_thread.end(), block_root.begin(), block_root.end());
                }
            }
        }, ap);

        for (const auto& thread_vec : local_root) {
            root.insert(root.end(), thread_vec.begin(), thread_vec.end());
        }

        return !root.empty();

    }



};