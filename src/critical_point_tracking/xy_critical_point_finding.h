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

    int remap(int idx, int remove_idx) {
        return idx - (idx > remove_idx);
    }
    // compute derivative in a single span. here the third dimension is fixed. So only compute derivative in x and y directions
    template<typename T>
    void compute_f_dev_f(const Block<T>* b, VectorX<T>& p, VectorX<T>& f, MatrixX<T>& dev_f, int removed_dom)
    {
        int domain_dim = b->dom_dim-1;
        dev_f.resize(domain_dim,domain_dim);
        f.resize(domain_dim);

        VectorX<T> f_vector(1);
        VectorX<T> dev_f_vector(1);    

        // std::cout<<local_domain_range.transpose()<<std::endl;

        VectorXi deriv(b->dom_dim);

        if(removed_dom==domain_dim)
        {
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
        }
        else
        {
            for(int i=0;i<domain_dim;i++)
            {
                for(int j=i;j<b->dom_dim;j++)
                {
                    if(j==i && j==removed_dom)
                    {
                        continue;
                    }
                    deriv.setZero();
                    deriv[i]+=1;
                    deriv[j]+=1;
                    mfa_extend::recover_mfa(b, p,dev_f_vector,deriv);
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
                    T root_finding_epsilon, T hessian_det_epsilon, std::vector<int>& used_domain, int removed_dom)
    {
        int itr_num=0;

        VectorX<T> p_on_boundary(p.size()-1);
        for(int i=0;i<p_on_boundary.size();i++)
        {
            p_on_boundary[i]=p[used_domain[i]];
        }

        MatrixX<T> dev_f;
        VectorX<T> f;
        compute_f_dev_f(b,p,f,dev_f, removed_dom);

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

                for(int i=0;i<tem.size();i++)
                {
                    p[used_domain[i]]-= tem[i];
                }
                p_on_boundary -=tem;
            }
            
            

            if((p_on_boundary-center).squaredNorm()>d_max_square)
            {
                return false;
            }

             

            if(!utility::InBlock(span_range,p_on_boundary))
            {
                return false;
            }


            compute_f_dev_f(b,p,f,dev_f,removed_dom);


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

    template<typename T>
    void root_finding_on_one_boudnary(const Block<T>* b,std::vector<VectorX<T>>& root, std::vector<std::vector<T>>& span_range, T boundary_value, int boundary_dim_index,
        T root_finding_grad_epsilon, T same_root_epsilon,
        T hessian_det_epsilon, std::vector<int>& used_domain,int maxIter, T d_max_square,VectorX<T>& center,T step_size, T plane_step_size, T top_or_bottom) // top plane is 1, bottom plane is -1 
    {
   
 

        VectorXi degree(b->mfa->var(0).p.size()-1);

        std::vector<int> used_dom(degree.size());
        int j=0;
        for(int i=0;i<=degree.size();++i)
        {
            if(i!=boundary_dim_index)
            {
                used_dom[j]=i;
                degree[j]=b->mfa->var(0).p[i];
                j++;
            }
        }

        std::vector<std::vector<T>>initial_point;
        utility::compute_initial_points(initial_point,degree,span_range);

        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point[i].size();
        }
        int num_initial_point = num_initial_point_every_domain.prod();

        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size()+1);
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);


        std::vector<VectorX<T>> root_in_original_domain;
        VectorX<T> next_root; 
        


        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<initial_point.size();j++)
            {
                current_initial_point[used_dom[j]]=initial_point[j][domain_index[j]];
            }        
            current_initial_point[boundary_dim_index]=boundary_value;
            if(boundary_dim_index==degree.size())
            {
                current_initial_point[boundary_dim_index]+=step_size;
            }
            else
            {
                current_initial_point[boundary_dim_index]-=top_or_bottom*plane_step_size;
            }

            
            // if(boundary_dim_index==0)
            // {
            //     if(top_or_bottom==1.0)
            //     {
            //         std::cout<<current_initial_point.transpose()<<std::endl;

            //     }
            // }

            // std::cout<<"initial point "<<i<<" "<<  current_initial_point.transpose()<<std::endl;
            // VectorX<T> domain_range = b->core_maxs - b->core_mins;
            // VectorX<T> param = (current_initial_point-b->core_mins).cwiseQuotient(domain_range);
            // if(!utility::InDomain(param))
            // {
            //     std::cout<<param.transpose()<<std::endl;
            //     std::cout<<boundary_dim_index<<std::endl;
            //     std::cout<<"initial p "<<current_initial_point.transpose()<<std::endl;
            //     std::cout<<used_dom[0]<<" "<<used_dom[1]<<std::endl;
            // }

            if(newton(b, next_root, current_initial_point,maxIter,span_range,d_max_square,center,root_finding_grad_epsilon,hessian_det_epsilon,used_domain,boundary_dim_index))
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


    template<typename T>
    void span_range_for_plane(std::vector<std::vector<T>>& span_range, const Block<T>* b, int skipeed_dim,VectorX<T>& center, VectorXi& span_index)
    {
        VectorX<T> domain_range = b->core_maxs-b->core_mins;
        int j=0;
        for(int i=0;i<=span_range.size();++i)
        {
            if(i!=skipeed_dim)
            {
                span_range[j].clear();
                span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
                span_range[j].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);
                center[j]=(span_range[j][0]+span_range[j][1])*0.5;
                j++;
            }
           
        }

    }

    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
        T root_finding_grad_epsilon, T same_root_epsilon,
        T hessian_det_epsilon,int maxItr,T step_size, T plane_step_size) { 

        root.clear();
        
        // VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=100;

        int distance_stop_itr = 5;

        
        auto domain_range = b->core_maxs-b->core_mins;
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        //find all the boundary plane
        //
        std::vector<std::vector<std::vector<T>>> span_range;//
        std::vector<VectorX<T>> center;
        std::vector<T> fixed_value; //the plane of the fixed value
        std::vector<int> fixed_dim;
        auto& nctrl_pts = b->mfa->var(0).tmesh.tensor_prods[0].nctrl_pts;

        std::vector<std::vector<T>> span_range_for_one_plane(domain_range.size()-1);
        VectorX<T> temp_center(domain_range.size()-1);
        std::vector<std::vector<int>> used_domain;
        std::vector<int> used_domain_temp;
        std::vector<T> record_top_bottom_plane;
        for(int i=0;i<domain_range.size();++i)
        {
            if(span_index[i]==b->mfa->var(0).p[i])
            {
                fixed_value.emplace_back(b->core_mins[i]);

                span_range_for_plane(span_range_for_one_plane, b, i, temp_center, span_index);

                span_range.emplace_back(span_range_for_one_plane);
                center.emplace_back(temp_center);
                fixed_dim.emplace_back(i);

                record_top_bottom_plane.emplace_back(-1.0); // bottom plane is -1
                // if(i==2)
                // {
                //     if(span_range_for_one_plane[0][0]<0.675 && span_range_for_one_plane[0][1]>0.675 && span_range_for_one_plane[1][0]<0.0 && span_range_for_one_plane[1][1]>0.0)
                //     {
                //         std::cout<<"find the span "<<span_index.transpose()<<std::endl;
                //     }
                // }

                used_domain_temp.clear();
                for(int j=0;j<domain_range.size();++j)
                {
                    if(j!=i)
                    {
                        used_domain_temp.emplace_back(j);
                    }
                }
                used_domain.emplace_back(used_domain_temp);

            }
            if(i!=(domain_range.size()-1) && (span_index[i]==(nctrl_pts[i]-1)))
            {
                fixed_value.emplace_back(b->core_maxs[i]);
                span_range_for_plane(span_range_for_one_plane, b, i, temp_center, span_index);

                span_range.emplace_back(span_range_for_one_plane);
                center.emplace_back(temp_center);
                fixed_dim.emplace_back(i);
                record_top_bottom_plane.emplace_back(1.0); // top plane is 1

                used_domain_temp.clear();
                for(int j=0;j<domain_range.size();++j)
                {
                    if(j!=i)
                    {
                        used_domain_temp.emplace_back(j);
                    }
                }
                used_domain.emplace_back(used_domain_temp);
            }
        }

        std::vector<T> d_max_square;
        // compute distance to terminate iteration
        for(auto j=span_range.begin();j<span_range.end();++j)
        { 
            T d_max_square_temp=0;
            for(auto i=j->begin();i!=j->end();++i)
            {  
                d_max_square_temp+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
            }

            d_max_square_temp*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2
            d_max_square.emplace_back(d_max_square_temp);
        }
       
        for(int i=0;i<fixed_value.size();++i)
        {
            // std::cout<<"find root on boundary "<<fixed_value[i]<<std::endl;
            root_finding_on_one_boudnary(b, root, span_range[i], fixed_value[i], fixed_dim[i], root_finding_grad_epsilon, same_root_epsilon, hessian_det_epsilon, used_domain[i],maxItr,d_max_square[i],center[i], step_size,plane_step_size, record_top_bottom_plane[i]);
        }

        return !root.empty();

    }



    template<typename T>
    bool root_finding(Block<T>* block, std::vector<std::vector<VectorXi>>& span_index, 
    std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
        int current_index,
        T root_finding_epsilon, T same_root_epsilon, T hessian_det_epsilon,int maxItr, T step_size, T plane_step_size) //2^n+1 initial points) 
    {

        for(auto i=0;i<block->mfa->nvars();++i)
        {
            if(root_finding(block,span_index[i][current_index], root,
            root_finding_epsilon,same_root_epsilon,hessian_det_epsilon,maxItr, step_size, plane_step_size))
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