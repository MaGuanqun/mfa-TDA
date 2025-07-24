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
    void find_neighbor_critical_points(std::vector<T>& ori_step_size, std::vector<T>& step_size, const Block<T>* b, const VectorX<T>& degenerate_point, std::vector<VectorX<T>>& start_points,T root_finding_grad_epsilon, T hessian_det_epsilon,int maxIter, T point_itr_threshold)
    {
        auto domain_range = b->core_maxs-b->core_mins;
        std::vector<std::vector<T>> span_range_for_one_plane(domain_range.size()-1);

        std::vector<std::vector<std::vector<T>>> span_range;
        std::vector<T> fixed_value; //the plane of the fixed value
        std::vector<int> fixed_dim;
        std::vector<std::vector<int>> used_domain(domain_range.size());
        std::vector<T> d_max_square;
        int distance_stop_itr = 5;

        std::vector<VectorX<T>> center;
        VectorX<T> temp_center(domain_range.size()-1);

        
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
            set_one_boudary_range(i, degenerate_point, span_range_for_one_plane, ori_step_size, b->core_maxs, b->core_mins);

            for(int j=0;j<span_range_for_one_plane.size();++j)
            {
                temp_center[j]=(span_range_for_one_plane[j][1]+span_range_for_one_plane[j][0])*0.5;
            }
            center.emplace_back(temp_center);

            span_range.emplace_back(span_range_for_one_plane);
            fixed_dim.emplace_back(i);
            if(b->core_maxs[i]-degenerate_point[i]>ori_step_size[i])
            {
                fixed_value.emplace_back(degenerate_point[i]+ori_step_size[i]);
            }
            else
            {
                fixed_value.emplace_back(b->core_maxs[i]);
            }
            if(degenerate_point[i]-b->core_mins[i]>ori_step_size[i])
            {
                fixed_value.emplace_back(degenerate_point[i]-ori_step_size[i]);
            }
            else
            {
                fixed_value.emplace_back(b->core_mins[i]);
            }

            T d_max_square_temp=0;
            for(auto j=span_range_for_one_plane.begin();j<span_range_for_one_plane.end();++j)
            {  
                d_max_square_temp+=((*j)[1]-(*j)[0])*((*j)[1]-(*j)[0]);
            }

            d_max_square_temp*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2
            d_max_square.emplace_back(d_max_square_temp);


        }

        for(int i=0;i<fixed_value.size();++i)
        {
            find_boundary_roots::root_finding_on_one_boundary(b, start_points, span_range[i/2], fixed_value[i], fixed_dim[i/2], root_finding_grad_epsilon, step_size, hessian_det_epsilon, used_domain[i/2], maxIter, d_max_square[i/2], center[i/2], point_itr_threshold);
        }

    }

    template<typename T>
    bool check_boundary_point_pass_degenerate(const VectorX<T>& degenerate_point, VectorX<T>& point, const Block<T>* b,std::vector<T>& step_size,std::vector<T>& ori_step_size, T hessian_det_epsilon, T gradient_epsilon, int max_itr)
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
            if(block_min[i]<b->core_mins[i])
            {
                block_min[i] = b->core_mins[i];
            }
            if(block_max[i]>b->core_maxs[i])
            {
                block_max[i] = b->core_maxs[i];
            }
        }

        int distance_stop_itr = 2;
        T d_max_square = distance_stop_itr*distance_stop_itr*(block_max - block_min).squaredNorm();

        std::vector<VectorX<T>> trajectory;

        bool pass =  particle_tracing::trajectory_pass_degenerate_point(step_size.back(), step_size[0], b, point, point[point.size()-1]<degenerate_point[degenerate_point.size()-1], hessian_det_epsilon, gradient_epsilon, max_itr, d_max_square, block_min, block_max,degenerate_point,trajectory);


        // std::cout<<"====="<<step_size.back()<<" "<<step_size[0]<<" "<<trajectory.size()<<std::endl;
        if(!pass){
        string name = "trajectory.csv";
        utility::saveToCSV(name, trajectory);
        }
        //determine if the trajectory passes the degenerate point
        return pass;
    }

    template<typename T>
    void find_neighbor_start_points(const Block<T>* b, const VectorX<T>& degenerate_point, T ori_time_step, T ori_spatial_step, T step_ratio, T root_finding_grad_epsilon, T hessian_det_epsilon,int maxIter, int correction_itr,std::vector<VectorX<T>>& start_points, T point_itr_threshold)
    {
        start_points.clear();

        std::vector<T> ori_step_size(degenerate_point.size(),ori_spatial_step);
        ori_step_size.back() = ori_time_step; // the last dimension is time

        std::vector<T> step_size(degenerate_point.size(),ori_spatial_step* step_ratio);
        step_size.back() = ori_time_step* step_ratio; // the last dimension is time
        
        std::vector<VectorX<T>> critical_points;
        find_neighbor_critical_points(ori_step_size, step_size, b, degenerate_point, critical_points, root_finding_grad_epsilon,hessian_det_epsilon, maxIter,point_itr_threshold);

        
        for(auto& cp: critical_points)
        {
            if(check_boundary_point_pass_degenerate(degenerate_point, cp, b, step_size, ori_step_size, hessian_det_epsilon, root_finding_grad_epsilon, correction_itr))
            {
                start_points.emplace_back(cp);
            }
        }
    }



    template<typename T>
    void tracing_from_start_points(const Block<T>* b, VectorX<T>& start_point, bool upper_tracing, std::vector<VectorX<T>>& trace, T time_step, T spatial_step, T hessian_det_epsilon, T gradient_epsilon, int correction_max_itr, T d_max_square)
    {

        particle_tracing::tracing_one_direction(time_step, spatial_step, b, start_point, trace, upper_tracing, hessian_det_epsilon, gradient_epsilon, correction_max_itr, d_max_square,
        b->core_mins, b->core_maxs);
        
        if(!upper_tracing)
        {
            std::reverse(trace.begin(), trace.end());
        }
    }

    
    template<typename T>
    void tracing_from_all_degenerate_points(const Block<T>* b, std::vector<VectorX<T>>& degenerate_points, std::vector<CP_Trace<T>>& trace, T time_step, T spatial_step, T step_ratio, T hessian_det_epsilon, T gradient_epsilon, int maxIter, int correction_max_itr, T d_max_square, T point_itr_threshold)
    {

        std::vector<VectorX<T>> start_points;
        std::vector<int> upper_tracing; 


        VectorX<T> test_point(3);
        test_point<<-0.0204463, 0.695221, 1.07636;

        for(int i=0;i<degenerate_points.size();++i)
        {

            if((degenerate_points[i]-test_point).norm()>0.001)
            {
                continue;
            }

            std::vector<VectorX<T>> temp_start_points;
            find_neighbor_start_points(b, degenerate_points[i], time_step, spatial_step, step_ratio, gradient_epsilon, hessian_det_epsilon, maxIter, correction_max_itr,temp_start_points,point_itr_threshold);

            std::vector<int> temp_upper_tracing;
            for(auto& start_point: temp_start_points)
            {
                if(start_point[start_point.size()-1]>degenerate_points[i][degenerate_points[i].size()-1])
                {
                    temp_upper_tracing.emplace_back(1);
                }
                else
                {
                    temp_upper_tracing.emplace_back(0);
                }
            }

            start_points.insert(start_points.end(), temp_start_points.begin(), temp_start_points.end());
            upper_tracing.insert(upper_tracing.end(), temp_upper_tracing.begin(), temp_upper_tracing.end());
        }

        if(start_points.empty())
        {
            std::cout<<"no start points found for degenerate points"<<std::endl;
        }

        int size = trace.size();
        
        trace.resize(size+start_points.size());
        for(int i=0;i<start_points.size();++i)
        {
            tracing_from_start_points(b, start_points[i], upper_tracing[i], trace[i+size].traces, time_step, spatial_step, hessian_det_epsilon, gradient_epsilon, correction_max_itr, d_max_square);
        }

        if (trace.size() == size)
        {
            std::cout << "Error: trace size does not match start points size." << std::endl;
        }
        

    }


}