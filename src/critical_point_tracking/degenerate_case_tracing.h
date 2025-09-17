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
#include "CP_Trace.h"
#include "xy_critical_point_finding.h"
#include "particle_tracing.h"
#include "spatial_hashing_spatial_temporal.h"
#include "degenerate_case.h"

namespace degenerate_case_tracing
{
    template<typename T>
    void read_degenerate_point(std::string& filename, std::vector<VectorX<T>>& singular_points)
    {
        std::vector<Eigen::MatrixXd> root;

        std::ifstream file(filename.c_str());
        if (!file) {
            std::cerr << "File does not exist: "<< std::endl;
            return;
        }



        utility::loadMatrixVector(filename.c_str(),root);

        std::cout<<"read critical points 0"<<root[0].rows()<<std::endl;
        
        singular_points.resize(root[0].rows());
        for(auto i=0;i<root[0].rows();++i)
        {
            singular_points[i]=root[0].row(i).transpose();
        }
        

    }
  

    template<typename T>
    void set_one_dim_range(T degenerate_point, T step_size, T core_min, T core_max, T& min, T& max)
    {
        if(degenerate_point-step_size> core_min)
        {
            min = degenerate_point-step_size;
        }
        else
        {
            min = core_min;
        }

        if(degenerate_point+step_size < core_max)
        {
            max = degenerate_point+step_size;
        }
        else
        {
            max = core_max;
        }
    }

    template<typename T>
    void set_one_boudary_range(int dim, const VectorX<T>& degenerate_point,std::vector<std::vector<T>>& span_range_for_one_plane, std::vector<T>& step_size, const VectorX<T>& max, const VectorX<T>& min)
    {
        int j=0;
        for(int i=0;i<degenerate_point.size();++i)
        {
            if(i!=dim)
            {                
                span_range_for_one_plane[j].clear();
                if(degenerate_point[i]-step_size[i]> min[i])
                {
                    span_range_for_one_plane[j].emplace_back(degenerate_point[i]-step_size[i]);
                }
                else
                {
                    span_range_for_one_plane[j].emplace_back(min[i]);
                }
                if(degenerate_point[i]+step_size[i]<max[i])
                {
                    span_range_for_one_plane[j].emplace_back(degenerate_point[i]+step_size[i]);
                }
                else
                {
                    span_range_for_one_plane[j].emplace_back(max[i]);
                }
                j++;
            }
        }
    }


    template<typename T>
    void find_neighbor_critical_points(std::vector<T>& ori_step_size, std::vector<T>& step_size, const VectorX<T>& degenerate_point, std::vector<VectorX<T>>& start_points,T root_finding_grad_epsilon, T hessian_det_epsilon,int maxIter, T point_itr_threshold, VectorXi point_num_in_block, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {

        auto domain_range = core_maxs-core_mins;
        // std::vector<std::vector<T>> span_range_for_one_plane(domain_range.size()-1);
        std::vector<T> fixed_value; //the plane of the fixed value
        std::vector<int> fixed_dim;
        std::vector<std::vector<int>> used_domain(domain_range.size());
        // std::vector<T> d_max_square;
        int distance_stop_itr = 5;
        // std::vector<VectorX<T>> center;
        // VectorX<T> temp_center(domain_range.size()-1);
        VectorX<T> block_min(core_maxs.size());
        VectorX<T> block_max(core_maxs.size());
        for(int i=0;i<domain_range.size();++i)
        {
            for(int j=0;j<domain_range.size();++j)
            {
                if (i!=j)
                {
                    used_domain[i].emplace_back(j);
                }
            }
                //plane at degenerate_point[i]+ori_step_size
            set_one_dim_range(degenerate_point[i], ori_step_size[i], core_mins[i], core_maxs[i], block_min[i], block_max[i]);

            // for(int j=0;j<span_range_for_one_plane.size();++j)
            // {
            //     temp_center[j]=(span_range_for_one_plane[j][1]+span_range_for_one_plane[j][0])*0.5;
            // }
            // center.emplace_back(temp_center);
            // span_range.emplace_back(span_range_for_one_plane);
            fixed_dim.emplace_back(i);
            fixed_value.emplace_back(block_max[i]);
            fixed_value.emplace_back(block_min[i]);

        }


        // VectorXi point_num_in_block = b->mfa->var(0).p + VectorXi::Ones(b->mfa->var(0).p.size()); //n

        std::vector<std::vector<T>>initial_points;
        VectorXi set_block_num = VectorXi::Ones(point_num_in_block.size());
        cp_tracking_degenerate_case::generate_initial_points(initial_points,block_min,block_max,point_num_in_block,set_block_num,b);

        std::vector<std::array<int,2>> initial_point_range(domain_range.size());
        for(int i=0;i<domain_range.size();++i)
        {
            initial_point_range[i][0]=0;
            initial_point_range[i][1]=initial_points[i].size();
        }

        std::vector<std::vector<T>> span_range(domain_range.size());
        for(int i=0;i<span_range.size();++i)
        {
            span_range[i].emplace_back(block_min[i]);
            span_range[i].emplace_back(block_max[i]);
        }



        std::vector<VectorX<T>> start_points_on_one_boundary;
        for(int i=0;i<fixed_value.size();++i)
        {
            start_points_on_one_boundary.clear();

            find_boundary_roots::root_finding_on_one_boundary(start_points_on_one_boundary, fixed_value[i], fixed_dim[i/2], root_finding_grad_epsilon, step_size, hessian_det_epsilon, used_domain[i/2], maxIter, initial_points,  initial_point_range, point_itr_threshold,core_mins, core_maxs, function_type, b);



            for(auto& sp: start_points_on_one_boundary)
            {
                if(utility::InBlock(span_range,sp,fixed_dim[i/2]))
                {
                    start_points.emplace_back(sp);
                }
                
            }
        }

    }

    template<typename T>
    bool check_boundary_point_pass_degenerate(const VectorX<T>& degenerate_point, VectorX<T>& point, std::vector<T>& step_size,std::vector<T>& ori_step_size, T hessian_det_epsilon, T gradient_epsilon, int max_itr,const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        VectorX<T> step(point.size());
        for(int i=0;i<point.size();++i)
        {
            step[i] = ori_step_size[i];
        }

        VectorX<T> block_min = degenerate_point- step;
        VectorX<T> block_max = degenerate_point+ step;
        for(int i=0;i<point.size();++i)
        {
            if(block_min[i]<core_mins[i])
            {
                block_min[i] = core_mins[i];
            }
            if(block_max[i]>core_maxs[i])
            {
                block_max[i] = core_maxs[i];
            }
        }

        int distance_stop_itr = 2;
        T d_max_square = distance_stop_itr*distance_stop_itr*(block_max - block_min).squaredNorm();

        std::vector<VectorX<T>> trajectory;

        bool pass =  particle_tracing::trajectory_pass_degenerate_point(step_size.back(), step_size[0], point, point[point.size()-1]<degenerate_point[degenerate_point.size()-1], hessian_det_epsilon, gradient_epsilon, max_itr, d_max_square, block_min, block_max,degenerate_point,trajectory,core_mins, core_maxs,function_type,b);


        // std::cout<<"====="<<step_size.back()<<" "<<step_size[0]<<" "<<trajectory.size()<<std::endl;
 
        string name = "trajectory.csv";
        utility::saveToCSV(name, trajectory);
        
        //determine if the trajectory passes the degenerate point
        return pass;
    }

    template<typename T>
    void find_neighbor_start_points(const VectorX<T>& degenerate_point, std::vector<T>& ori_step_size, std::vector<T>& step_size, T root_finding_grad_epsilon, T hessian_det_epsilon,int maxIter, int correction_itr,std::vector<VectorX<T>>& start_points, T point_itr_threshold,VectorXi point_num_in_block, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {
        start_points.clear();

             
        std::vector<VectorX<T>> critical_points;
        find_neighbor_critical_points(ori_step_size, step_size, degenerate_point, critical_points, root_finding_grad_epsilon,hessian_det_epsilon, maxIter,point_itr_threshold, point_num_in_block, core_mins, core_maxs,function_type,b);

        start_points = critical_points;

        // for(auto& cp: critical_points)
        // {
        //     if(check_boundary_point_pass_degenerate(degenerate_point, cp, b, step_size, ori_step_size, hessian_det_epsilon, root_finding_grad_epsilon, correction_itr))
        //     {
        //         start_points.emplace_back(cp);
        //     }
        // }
    }



    template<typename T>
    void tracing_from_start_points(VectorX<T>& start_point, bool upper_tracing, std::vector<VectorX<T>>& trace, T time_step, T spatial_step, T hessian_det_epsilon, T gradient_epsilon, int correction_max_itr, T d_max_square, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {

        particle_tracing::tracing_one_direction(time_step, spatial_step, start_point, trace, upper_tracing, hessian_det_epsilon, gradient_epsilon, correction_max_itr, d_max_square,
        core_mins, core_maxs, core_mins, core_maxs, function_type,b);
        
        if(!upper_tracing)
        {
            std::reverse(trace.begin(), trace.end());
        }
    }

    
    template<typename T>
    void tracing_from_all_degenerate_points(std::vector<VectorX<T>>& degenerate_points, std::vector<CP_Trace<T>>& trace, T time_step, T spatial_step, T step_ratio, T hessian_det_epsilon, T gradient_epsilon, int maxIter, int correction_max_itr, T d_max_square, T point_itr_threshold,VectorXi point_num_in_block,  const VectorX<T>& core_mins, const VectorX<T>& core_maxs, const int function_type=0, const Block<T>* b=nullptr)
    {

        std::vector<VectorX<T>> raw_start_points;
        std::vector<int> upper_tracing; 


        std::vector<T> ori_step_size(core_mins.size(), spatial_step);
        ori_step_size.back() = time_step; // the last dimension is time

        std::vector<T> step_size(core_mins.size(),spatial_step* step_ratio);
        step_size.back() = time_step* step_ratio; // the last dimension is time
  
        // std::cout<<"ori step size "<<ori_step_size[0]<<" "<<ori_step_size.back()<<std::endl;
        // std::cout<<step_size[0]<<" "<<step_size.back()<<std::endl;

        for(int i=0;i<degenerate_points.size();++i)
        {

            // if((degenerate_points[i]-test_point).norm()>0.001 &&
            // (degenerate_points[i]-test_point2).norm()>0.001
            // )
            // {
            //     continue;
            // }
            // std::cout<<degenerate_points[i].transpose()<<std::endl;

            std::vector<VectorX<T>> temp_start_points;
            find_neighbor_start_points(degenerate_points[i], ori_step_size, step_size, gradient_epsilon, hessian_det_epsilon, maxIter, correction_max_itr,temp_start_points,point_itr_threshold, point_num_in_block, core_mins, core_maxs,function_type,b);

            // std::vector<int> temp_upper_tracing;
            for(auto& start_point: temp_start_points)
            {
                if(start_point[start_point.size()-1]>degenerate_points[i][degenerate_points[i].size()-1])
                {
                    // temp_upper_tracing.emplace_back(1);
                    raw_start_points.emplace_back(start_point);

                }
                // else
                // {
                //     temp_upper_tracing.emplace_back(0);
                // }
            }

            // raw_start_points.insert(raw_start_points.end(), temp_start_points.begin(), temp_start_points.end());
            // upper_tracing.insert(upper_tracing.end(), temp_upper_tracing.begin(), temp_upper_tracing.end());
        }

        std::vector<VectorX<T>> start_points;
        spatial_hashing_spatial_temporal::find_all_unique_root(raw_start_points, start_points,step_size[0],step_size.back());

        // VectorX<T> test_root(3);
        // test_root<< 1.8337,1.0408,2.91297;

        // for(auto sp:start_points)
        // {
            // if((sp-test_root).norm()<0.001)
            // {
            //     std::cout<<"find test root "<<sp.transpose()<<std::endl;
            // }
           
        // }


        int size = trace.size();
        
        trace.resize(size+start_points.size());
        for(int i=0;i<start_points.size();++i)
        {
            // if((start_points[i]-test_root).norm()>0.001)
            // {
            //     continue;
            // }
            
            // std::cout<<start_points[i].transpose()<<std::endl;
            // std::cout<<"upper tracing "<<upper_tracing[i]<<std::endl;

            tracing_from_start_points(start_points[i], true, trace[i+size].traces, time_step, spatial_step, hessian_det_epsilon, gradient_epsilon, correction_max_itr, d_max_square,core_mins, core_maxs,function_type,b);

        }

        // if (trace.size() == size)
        // {
        //     std::cout << "Error: trace size does not match start points size." << std::endl;
        // }
        

    }


}