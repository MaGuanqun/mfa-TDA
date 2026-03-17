#pragma once
#include <vector>
#include <Eigen/Dense>
#include "CP_Trace.h"

#include "tracking_derivatives.h"
#include "INRModel.h"


namespace critical_point_utility
{

    template<typename T>
    int determine_type_2d(MatrixX<T>& Hessian)
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
    int determine_type_3d(MatrixX<T>& H)
    {
        // Use SelfAdjointEigenSolver for symmetric real matrices
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
        if (es.info() != Eigen::Success) {
            // fallback / indicate degenerate / failure
            return -1;
        }

        // const double eps = 1e-12; // adjust tolerance to your problem scale
        int pos = 0, neg = 0, zero = 0;
        auto vals = es.eigenvalues();

        for (int i = 0; i < vals.size(); ++i) {
            double v = vals[i];
            if (v > 0) ++pos;
            else if (v < 0) ++neg;
            else ++zero;
        }

        if (zero > 0) return -1;      // degenerate (zero eigenvalue)
        if (pos == (int)H.rows()) return 0;  // all positive -> minimum
        if (neg == (int)H.rows()) return 3;  // all negative -> maximum
        // Mixed signs -> saddle. For 3D specifically:
        if (neg == 1) return 1; // index-1 saddle (+,+,-)
        if (neg == 2) return 2; // index-2 saddle (+,-,-)

        return -1; // should not reach for non-degenerate real symmetric H
    }

    template<typename T>
    int determine_type(MatrixX<T>& Hessian)
    {
        if (Hessian.cols()==2)
        {
            return determine_type_2d(Hessian);
        }

        else if (Hessian.cols()==3)
        {
            return determine_type_3d(Hessian);
        }
        else
        {
            std::cout<<"error, only support 2D, 3D and 4D data"<<std::endl;
            return -1;
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
            inr_model->query_dim_reduced_grad_hessian(point, gradient, Hessian, point.size()-1);
            // inr_model->query_grad_exclude_last(point, gradient);
        }
        else
        {
            tracking_derivatives::compute_gradient(point, gradient, function_type, b);
        }
        return gradient.norm();
    }

    template<typename T>
    T compute_Hessian_single_point(VectorX<T>& point,int function_type=0, Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        MatrixX<T> Hessian(point.size()-1,point.size()-1);

        tracking_derivatives::compute_Hessian(point, Hessian, point.size()-1, function_type, b);
        
        std::cout<<Hessian.determinant() / Hessian.squaredNorm()<<std::endl;
        return Hessian.determinant();
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

    template<typename T>
    void test_accuracy(std::vector<VectorX<T>>& degenerate_points,int function_type=0, Block<T>* b=nullptr, INRModel<T>* inr_model=nullptr)
    {
        T accuracy_value=0;
        int num_points=0;
        T max_accuracy = 0;
        T accuracy;
        for(int i=0;i<degenerate_points.size();++i)
        {
            accuracy= compute_accuracy_single_point(degenerate_points[i], function_type, b, inr_model);
            compute_Hessian_single_point(degenerate_points[i], function_type, b, inr_model);
            std::cout<<"gradient norm "<<accuracy<<std::endl;
            accuracy_value += accuracy;

            if(accuracy>max_accuracy)
            {
                max_accuracy = accuracy;
            }
        }   
        num_points += degenerate_points.size();

        std::cout<<"degenerat_accuracy "<<accuracy_value/ (T)num_points<<std::endl;
        std::cout<<"degenerate max accuracy "<<max_accuracy<<std::endl;
        
    }
}