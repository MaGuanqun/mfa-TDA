#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"

// #include "mfa_extend.h"
#include "CP_Trace.h"
#include "find_boundary_roots.h"
#include "particle_tracing.h"
#include "spatial_hashing_spatial_temporal.h"
#include "degenerate_case.h"
#include "tracking_utility.h"


template<typename T>
class Degenerate_case_tracing{

private:

    const Block<T>* b;
    const int function_type;
    const VectorX<T> core_mins;
    const VectorX<T> core_maxs;
    const VectorXi point_num_in_block;
    T time_step;
    T spatial_step;
    std::vector<T> ori_step_size;
    T root_finding_grad_epsilon;
    int correction_max_itr;

    Find_boundary_roots<T>* find_boundary_roots;

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

    

    void find_neighbor_critical_points(std::vector<T>& step_size, const VectorX<T>& degenerate_point, std::vector<VectorX<T>>& start_points)
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


        std::vector<std::vector<T>>initial_points;
        VectorXi set_block_num = VectorXi::Ones(point_num_in_block.size());
        tracking_utility::generate_initial_points(initial_points,block_min,block_max,point_num_in_block,set_block_num,b);

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

            find_boundary_roots->root_finding_on_one_boundary(start_points_on_one_boundary, fixed_value[i], fixed_dim[i/2], used_domain[i/2], initial_points,  initial_point_range, step_size);



            for(auto& sp: start_points_on_one_boundary)
            {
                if(utility::InBlock(span_range,sp,fixed_dim[i/2]))
                {
                    start_points.emplace_back(sp);
                }
                
            }
        }

    }



    
    // bool check_boundary_point_pass_degenerate(const VectorX<T>& degenerate_point, VectorX<T>& point, std::vector<T>& step_size, int max_itr)
    // {
    //     VectorX<T> step(point.size());
    //     for(int i=0;i<point.size();++i)
    //     {
    //         step[i] = ori_step_size[i];
    //     }

    //     VectorX<T> block_min = degenerate_point- step;
    //     VectorX<T> block_max = degenerate_point+ step;
    //     for(int i=0;i<point.size();++i)
    //     {
    //         if(block_min[i]<core_mins[i])
    //         {
    //             block_min[i] = core_mins[i];
    //         }
    //         if(block_max[i]>core_maxs[i])
    //         {
    //             block_max[i] = core_maxs[i];
    //         }
    //     }

    //     int distance_stop_itr = 2;
    //     T d_max_square = distance_stop_itr*distance_stop_itr*(block_max - block_min).squaredNorm();

    //     std::vector<VectorX<T>> trajectory;

    //     bool pass =  particle_tracing::trajectory_pass_degenerate_point(step_size.back(), step_size[0], point, point[point.size()-1]<degenerate_point[degenerate_point.size()-1],  root_finding_grad_epsilon, max_itr, d_max_square, block_min, block_max,degenerate_point,trajectory,core_mins, core_maxs,function_type,b);


    //     // std::cout<<"====="<<step_size.back()<<" "<<step_size[0]<<" "<<trajectory.size()<<std::endl;
 
    //     string name = "trajectory.csv";
    //     utility::saveToCSV(name, trajectory);
        
    //     //determine if the trajectory passes the degenerate point
    //     return pass;
    // }

    void find_neighbor_start_points(const VectorX<T>& degenerate_point, std::vector<T>& step_size, std::vector<VectorX<T>>& start_points)
    {
        start_points.clear();

             
        std::vector<VectorX<T>> critical_points;
        find_neighbor_critical_points(step_size, degenerate_point, critical_points);

        start_points = critical_points;

    }

    
    void tracing_from_start_points(VectorX<T>& start_point, bool upper_tracing, std::vector<VectorX<T>>& trace, T d_max_square)
    {

        particle_tracing::tracing_one_direction(time_step, spatial_step, start_point, trace, upper_tracing,  root_finding_grad_epsilon, correction_max_itr, d_max_square,
        core_mins, core_maxs, core_mins, core_maxs, function_type,b);
        
        if(!upper_tracing)
        {
            std::reverse(trace.begin(), trace.end());
        }
    }


public:

    Degenerate_case_tracing(const VectorX<T> domain_min, const VectorX<T> domain_max, const VectorXi point_num_in_b, Find_boundary_roots<T>* find_boundary_roots_, T time_step_, T spatial_step_size_, T root_finding_grad_epsilon_, int correction_max_it, int func_type=0, const Block<T>* block=nullptr): core_mins(domain_min), core_maxs(domain_max), b(block), function_type(func_type), point_num_in_block(point_num_in_b),find_boundary_roots(find_boundary_roots_), time_step(time_step_), spatial_step(spatial_step_size_), root_finding_grad_epsilon(root_finding_grad_epsilon_),correction_max_itr(correction_max_it)
    {
        ori_step_size.resize(domain_min.size(), spatial_step_size_);
        ori_step_size.back() = time_step_; // the last dimension is time
    }

    static void read_degenerate_point(std::string& filename, std::vector<VectorX<T>>& singular_points)
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
  

    
    
    void tracing_from_all_degenerate_points(std::vector<VectorX<T>>& degenerate_points, std::vector<CP_Trace<T>>& trace, T step_ratio,  T d_max_square)
    {

        std::vector<VectorX<T>> raw_start_points;
        std::vector<int> upper_tracing; 



        std::vector<T> step_size(core_mins.size(),spatial_step* step_ratio);
        step_size.back() = time_step* step_ratio; // the last dimension is time
  


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
            find_neighbor_start_points(degenerate_points[i], step_size, temp_start_points);

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

            tracing_from_start_points(start_points[i], true, trace[i+size].traces, d_max_square);

        }

        // if (trace.size() == size)
        // {
        //     std::cout << "Error: trace size does not match start points size." << std::endl;
        // }
        

    }


};