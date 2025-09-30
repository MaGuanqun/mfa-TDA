#pragma once
#include <vector>
#include <Eigen/Dense>
#include <mfa/mfa.hpp>
#include "utility_function.h"

namespace tracking_utility
{
    template<typename T>
    void convert_to_obj(const std::string& filename, std::vector<Eigen::VectorX<T>>& points)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        for(auto i=points.begin();i<points.end();++i)
        {

            outFile << std::setprecision(15) << "v " << (*i).data()[0] << " " << (*i).data()[1] << " " << (*i).data()[2] << "\n";
            
        }
        outFile.close();

    }

    
    template<typename T>
    void generate_initial_points(std::vector<vector<T>>& initial_points, const VectorX<T>& domain_min, const VectorX<T>& domain_max,
    const VectorXi& point_num_in_block,const VectorXi& set_block_num,const Block<T>* b=nullptr)
    {
        if(b==nullptr)
        {
            
            std::vector<std::vector<T>> domain_range(domain_min.size());
            for(int i=0;i<domain_min.size();++i)
            {
                domain_range[i].emplace_back(domain_min[i]);
                domain_range[i].emplace_back(domain_max[i]);
            }

            VectorXi initial_point_number = point_num_in_block.cwiseProduct(set_block_num);
           
            utility::compute_initial_points2(initial_points,initial_point_number,domain_range);
        }
        else
        {
            initial_points.resize(domain_min.size());
            for(int i=0;i<domain_min.size();++i)
            {
                initial_points[i].reserve(point_num_in_block[i]*set_block_num[i]);
                for(int j=0;j<set_block_num[i];++j)
                {
                    T span_min = b->mfa->var(0).tmesh.all_knots[i][j+b->mfa->var(0).p[i]]*(domain_max[i]-domain_min[i])+domain_min[i];
                    T span_max = b->mfa->var(0).tmesh.all_knots[i][j+b->mfa->var(0).p[i]+1]*(domain_max[i]-domain_min[i])+domain_min[i];

                    std::vector<T> current_span_range{span_min,span_max};
                    std::vector<T> initial_points_single_dim;
                    utility::compute_initial_points_single_dim(initial_points_single_dim,point_num_in_block[i],current_span_range);

                    initial_points[i].insert(initial_points[i].end(), initial_points_single_dim.begin(), initial_points_single_dim.end());
                }
            }
        }

    }

    template<typename T>
    bool newRoot(VectorX<T>& z, std::vector<VectorX<T>>& root_so_far, std::vector<T>& threshold)
    {
        for(int i=0;i<root_so_far.size();++i)
        {
            if(std::abs(z[z.size()-1]-root_so_far[i][z.size()-1])<threshold.back())
            {
                if((z.head(z.size()-1)-root_so_far[i].head(z.size()-1)).squaredNorm()<threshold[0]*threshold[0])
                {
                    return false;
                }
            }
        }
        return true;
    }


}