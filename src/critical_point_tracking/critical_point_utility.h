#pragma once
#include <vector>
#include <Eigen/Dense>
#include "CP_Trace.h"

#include "tracking_derivatives.h"
#include "INRModel.h"


namespace critical_point_utility
{
    template<typename T>
    int determine_type(MatrixX<T>& Hessian)
    {
        T D=Hessian.determinant();
        if(D>0)
        {
            if(Hessian(0,0)>0)
            {
                return 0; //minimum
            }
            else
            {
                return 3; //maximum
            }
        }
        else if(D<0)
        {
            return 1; //saddle point
        }
        else
        {
            return -1; //degenerate point
        }
    }

    template<typename T>
    void critical_point_type_single_point(VectorX<T>& p, std::vector<int>& critical_point_types, int function_type=0, Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        MatrixX<T> hessian(p.size()-1,p.size()-1);
        if(inr_model!=nullptr)
        {
            VectorX<T> gradient(2);
            inr_model->query_hessian_t(p,gradient, hessian);
        }
        else
        {
            tracking_derivatives::compute_Hessian(p, hessian, p.size()-1, function_type, b);
        }

        int cpt_type = determine_type(hessian);
        critical_point_types.push_back(cpt_type);
    }


    template<typename T>
    void compute_critical_point_type(std::vector<CP_Trace<T>>& traces,std::vector<VectorX<T>>& degenerate_points, std::vector<int>& critical_point_types,int function_type=0, Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        critical_point_types.clear();


        for(int i=0;i<degenerate_points.size();++i)
        {
            VectorX<T>& p = degenerate_points[i];
            critical_point_type_single_point(p, critical_point_types, function_type, b, inr_model);
        }
        

        for(int i=0;i<traces.size();++i)
        {
            if(!traces[i].duplicated){
                for(int j=0;j<traces[i].traces.size();++j)
                {
                    VectorX<T>& p = traces[i].traces[j];
                    critical_point_type_single_point(p, critical_point_types, function_type, b, inr_model);
                }
            }
        }
    }

    template<typename T>
    T compute_accuracy_single_point(VectorX<T>& point,const int function_type=0, const Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        VectorX<T> gradient;
        if(inr_model!=nullptr)
        {
            Eigen::MatrixX<T> Hessian;
            inr_model->query_dim_reduced_grad_hessian(point, gradient, Hessian, 2);
            // inr_model->query_grad_exclude_last(point, gradient);
        }
        else
        {
            tracking_derivatives::compute_gradient(point, gradient, function_type, b);
        }
        return gradient.norm();

    }

    template<typename T>
    void accuracy(std::vector<CP_Trace<T>>& traces,std::vector<VectorX<T>>& degenerate_points,int function_type=0, Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        T accuracy_value=0;
        int num_points=0;
        T max_accuracy = 0;
        T accuracy;
        for(int i=0;i<degenerate_points.size();++i)
        {
            accuracy= compute_accuracy_single_point(degenerate_points[i], function_type, b, inr_model);
            accuracy_value += accuracy;

            if(accuracy>max_accuracy)
            {
                max_accuracy = accuracy;
            }
        }   
        num_points += degenerate_points.size();

        std::cout<<"degenerat_accuracy "<<accuracy_value/ (T)num_points<<std::endl;
        std::cout<<"degenerate max accuracy "<<max_accuracy<<std::endl;

        //conmpute_gradient
        for(int i=0;i<traces.size();++i)
        {
            if(traces[i].duplicated)
            {
                continue;
            }

            for(int j=0;j<traces[i].traces.size();++j)
            {
                accuracy = compute_accuracy_single_point(traces[i].traces[j], function_type, b, inr_model);
                accuracy_value += accuracy;

                if(accuracy>max_accuracy)
                {
                    max_accuracy = accuracy;
                }
            }
            num_points += traces[i].traces.size();
           
        }

        if(num_points==0)
        {
            std::cout<<"No critical points found, cannot compute accuracy "<<std::endl;
            return;
        }
        accuracy_value = accuracy_value / (T)num_points;
        std::cout<<"Average critical point gradient norm: "<<accuracy_value<<std::endl;
        std::cout<<"Max critical point gradient norm: "<<max_accuracy<<std::endl;
        
    }
}