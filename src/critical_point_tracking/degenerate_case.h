//find the degenerate cases where the critical points are merging/splitting
//find the critical points of xy plane with singular hessian matrix
//J=[f_xxf_yy-f_xy^2, f_x, f_y]=0
#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <mfa/mfa.hpp>


#include "opts.h"

#include "block.hpp"
#include "tracking_utility.h"
#include "tracking_derivatives.h"
#include "INRModel.h"

template<typename T>
class Tracking_degenerate_case
{
private:
    const Block<T>* b;
    const int function_type;
    const VectorX<T> domain_min;
    const VectorX<T> domain_max;
    INRModel<T>* inr_model;
    
    T degenerate_finding_epsilon;
    T gradient_epsilon;
    int max_itr;

    std::vector<T> same_root_epsilon;

    void adjugate_of_symmetric(const MatrixX<T>& m, MatrixX<T>& adj)
    // m must be a symmetric matrix
    {
        int n = m.rows();
        adj.resize(n,n);
        MatrixX<T> sub_m(n-1,n-1);
        if(n==2)
        {   
            adj(0,0) = m(1,1);
            adj(0,1) = -m(0,1);
            adj(1,0) = -m(0,1);
            adj(1,1) = m(0,0);
            return;
        }
        
        for(int i=0;i<n;i++)
        {
            for(int j=i;j<n;j++)
            {
                sub_m.block(0,0,i,j) = m.block(0,0,i,j);
                if(i<n-1)
                {
                    sub_m.block(i,0,n-i-1,j) = m.block(i + 1, 0, n - i - 1, j);  
                }
                if(j<n-1)
                {
                    sub_m.block(0,j,i,n-j-1) = m.block(0, j + 1, i, n-j-1);
                }
                if(i<n-1 && j<n-1)
                {
                     sub_m.block(i,j,n-i-1,n-j-1) = m.block(i + 1, j + 1, n - i - 1, n - j - 1);
                }
               
                adj(j,i) = (i+j)%2==0 ? sub_m.determinant() : -sub_m.determinant();
                adj(i,j) = adj(j,i); // since m is symmetric, adj is also symmetric
            }
        }
    }

    
    void Hessian(VectorX<T>& p, MatrixX<T>& Hessian)
    {
        int n_var_dim = p.size() - 1;

        tracking_derivatives::compute_Hessian(p, Hessian, n_var_dim, function_type, b);
    }

    
    void partial_xt(VectorX<T>& p, VectorX<T>& result)
    {
        int n_var_dim = p.size() - 1;
        result.resize(n_var_dim);
        VectorXi deriv(p.size());
        VectorX<T> out(1);

        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[n_var_dim]+=1;

            query_function::query_function(p, out, function_type, b, deriv);

            result(i) = out[0];
        }
    }


    
    void partial_Hessian(VectorX<T>& p, MatrixX<T>& p_Hessian, int partial_deriv_index, VectorX<T>* third_deriv_INR=nullptr)
    {
        int n_var_dim = p.size() - 1;
        p_Hessian.resize(n_var_dim, n_var_dim);

        if(third_deriv_INR==nullptr)
        {
            VectorXi deriv(p.size());
            VectorX<T> result(1);

            for(int i=0;i<n_var_dim;i++)
            {
                for(int j=i;j<n_var_dim;j++)
                {
                    deriv.setZero();
                    deriv[i]+=1;
                    deriv[j]+=1;
                    deriv[partial_deriv_index] += 1; // partial derivative with respect to the partial_deriv_index

                    query_function::query_function(p, result, function_type, b, deriv);

                    p_Hessian(j,i) = result[0];
                    p_Hessian(i,j) = result[0];
                }
            }
        }
        else
        {
            for(int i=0;i<n_var_dim;i++)
            {
                for(int j=i;j<n_var_dim;j++)
                {
                    p_Hessian(j,i) = third_deriv_INR->data()[3*partial_deriv_index + i + j];//[_xx(), _xy(), _yy()]
                    p_Hessian(i,j) = p_Hessian(j,i);
                }
            }
        }
    }


    
    void diff_Hessian(VectorX<T>& p, std::vector<MatrixX<T>>& d_Hessian, VectorX<T>* third_deriv_INR=nullptr)
    {
        d_Hessian.resize(p.size());
        for(int i=0;i<p.size();i++)
        {
            partial_Hessian(p, d_Hessian[i], i,third_deriv_INR);       
        }
    }
    
    void diff_Hessian_determinant(VectorX<T>& p,const MatrixX<T>& hessian, VectorX<T>& diff, VectorX<T>* third_deriv_INR=nullptr)
    // third_deriv_INR: _xxx, _xyx, _yyx, _xxy, _xyy, _yyy, _xxt, _xyt, _yyt
    {
        MatrixX<T> adj;
        adjugate_of_symmetric(hessian, adj);
        std::vector<MatrixX<T>> diff_Hes;
        diff_Hessian(p, diff_Hes,third_deriv_INR);

        diff.resize(p.size());
        for(int i=0;i<p.size();i++)
        {
            diff[i]=adj.cwiseProduct(diff_Hes[i]).sum();
        }
    }




    
    // void compute_J(VectorX<T>& p, VectorX<T>& J)
    // {
    //     J.resize(p.size());
    //     VectorXi deriv(p.size());
    //     int n_var_dim=p.size()-1;
    //     VectorX<T> result(1); 
    //     for(int i=0;i<n_var_dim;i++)
    //     {
    //         deriv.setZero();
    //         deriv[i]+=1;
    //         query_function::query_function(p, result, function_type, b, deriv);
    //         J[i+1] = result[0];
    //     }

    //     MatrixX<T> Hessian_f;
    //     Hessian(p, Hessian_f);
    //     J[0] = Hessian_f.determinant();        
    // }

    void compute_J_dev_J_INR(VectorX<T>& p, VectorX<T>& J, MatrixX<T>& Jacobian_f)
    {
        J.resize(p.size());
        Jacobian_f.resize(p.size(),p.size());

        Eigen::MatrixX<T> Hessian_;
        VectorX<T> gradient;
        VectorX<T> third_spatial;
        VectorX<T> third_tmix;
               
        inr_model->query_up_to_third_derivative(p, gradient, Hessian_, third_spatial, third_tmix);
        J.tail(p.size()-1) = gradient.head(p.size()-1); // only xy gradient
        Eigen::MatrixX<T> Hessian_xy = Hessian_.block(0,0,p.size()-1,p.size()-1);
        J[0] = Hessian_xy.determinant();

            // third_deriv_INR: _xxx, _xyx, _yyx, _xxy, _xyy, _yyy, _xxt, _xyt, _yyt

        VectorX<T> third_deriv_INR(9);
        third_deriv_INR.head(3) = third_spatial.head(3); // _xxx, _xyx, _yyx
        third_deriv_INR.segment(3,3) = third_spatial.tail(3); // _xxy, _xyy, _yyy
        third_deriv_INR.tail(3) = third_tmix; // _xxt, _xyt, _yyt

        VectorX<T> diff;
        diff_Hessian_determinant(p, Hessian_xy, diff, &third_deriv_INR);

        Jacobian_f.row(0) = diff.transpose();
        Jacobian_f.block(1,0,p.size()-1,p.size()) = Hessian_.block(0,0,p.size()-1,p.size());


    //   Jacobian_f.row(0) = diff.transpose();
    //     Jacobian_f.block(1,0,n_var_dim,n_var_dim) = Hessian_f;

    //     VectorX<T> partial_xt_result;
    //     partial_xt(p, partial_xt_result);
    //     Jacobian_f.block(1,n_var_dim,n_var_dim,1) = partial_xt_result;

    }

    
    void compute_J_dev_J(VectorX<T>& p, VectorX<T>& J, MatrixX<T>& Jacobian_f)
    {

        if(inr_model!=nullptr)
        {
            compute_J_dev_J_INR(p,J,Jacobian_f);
            return;
        }

        J.resize(p.size());
        Jacobian_f.resize(p.size(),p.size());

        VectorXi deriv(p.size());

        int n_var_dim=p.size()-1;

        VectorX<T> result(1); 


        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;

            query_function::query_function(p, result, function_type, b, deriv);


            J[i+1] = result[0];
        }

        MatrixX<T> Hessian_f;
        Hessian(p, Hessian_f);
        J[0] = Hessian_f.determinant();     

        VectorX<T> diff;
        diff_Hessian_determinant(p, Hessian_f, diff);

        Jacobian_f.row(0) = diff.transpose();
        Jacobian_f.block(1,0,n_var_dim,n_var_dim) = Hessian_f;

        VectorX<T> partial_xt_result;
        partial_xt(p, partial_xt_result);
        Jacobian_f.block(1,n_var_dim,n_var_dim,1) = partial_xt_result;

    }

    
    bool newton(VectorX<T>& result, VectorX<T>& p)
    {
        int itr_num=0;
        MatrixX<T> dev_J;
        VectorX<T> J;
        compute_J_dev_J(p,J,dev_J);
        if(J.norm()<degenerate_finding_epsilon)
        {
            result = p;
            return true;
        }

        result=p;
        T temp_rec;

        VectorX<T> pre_point = p;
        while(itr_num<max_itr)
        {
            // T determinant = dev_J.determinant();
            Eigen::ColPivHouseholderQR<MatrixX<T>> qr(dev_J);

            pre_point=p;

            if(qr.rank() < dev_J.cols())
            {
                return false;
            }


            p -= qr.solve(J); 

            if(!utility::In_Domain(p,domain_min,domain_max))
            {
                return false;
            }

            compute_J_dev_J(p,J,dev_J);   

            if(itr_num>0){
                if(J.norm()< degenerate_finding_epsilon 
                && J.tail(J.size()-1).norm()<gradient_epsilon
                ){                    
                    result = p;
                    return true;
                }
            }

            itr_num++;

        }
        return false;

    }


    


   // Function to find the roots of the polynomial using Newton's method
    
    bool degenerate_finding_single_block(std::vector<VectorX<T>>& root,
        std::vector<std::vector<T>>&initial_point, std::vector<std::array<int,2>>& initial_point_range)
    {

        root.clear();
        VectorXi one = VectorXi::Ones(domain_min.size());       
        auto domain_range = domain_max-domain_min;

        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point_range[i][1]-initial_point_range[i][0];
        }

        int num_initial_point = num_initial_point_every_domain.prod();

        VectorX<T> next_root; 
        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size());
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);


        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<initial_point.size();j++)
            {
                current_initial_point[j]=initial_point[j][initial_point_range[j][0]+domain_index[j]];
            } 


            if(newton(next_root, current_initial_point))
            {
        
                if(tracking_utility::newRoot(next_root,root,same_root_epsilon))
                {        
                    root.emplace_back(next_root);
                } 
            }

        }        

        return !root.empty();

    }


public:

    Tracking_degenerate_case(const VectorX<T>& core_mins_, const VectorX<T>& core_maxs_, T degenerate_finding_epsilon_, T gradient_epsilon_,std::vector<T> same_root_epsilon_, int max_itr_=50,  const int function_type_=0, const Block<T>* b_=nullptr, INRModel<T>* inr_model_=nullptr)
    : domain_min(core_mins_), domain_max(core_maxs_), b(b_), function_type(function_type_), degenerate_finding_epsilon(degenerate_finding_epsilon_), gradient_epsilon(gradient_epsilon_), max_itr(max_itr_), same_root_epsilon(same_root_epsilon_), inr_model(inr_model_)
    {
        std::cout<<"degenerate case tracking initialized "<<std::endl;
        std::cout<<"domain min "<<domain_min.transpose()<<std::endl;
        std::cout<<"domain max "<<domain_max.transpose()<<std::endl;
    }
    ~Tracking_degenerate_case(){}
    
    void degenerate_finding(std::vector<VectorX<T>>& root,
    const VectorXi& point_num_in_block, const VectorXi& set_block_num, const std::vector<VectorXi>& selected_span_index = std::vector<VectorXi>())
    {
        std::vector<vector<T>> initial_points;
        tracking_utility::generate_initial_points(initial_points,domain_min,domain_max,point_num_in_block,set_block_num,b);
        

        // try to divide the initial points into different blocks. For MFA, every block is a span, for other functions, it is manually defined. So that it is easier to deduplicate

        //MFA. both spans and initial points are uniformly distributed, number of points = (degree+1)*span_num
        VectorXi number_in_every_dim;
        int num_block;
        if(selected_span_index!=std::vector<VectorXi>())
        {
            num_block = selected_span_index.size();
        }
        else
        {
            num_block = set_block_num.prod();
            utility::obtain_number_in_every_domain(set_block_num,number_in_every_dim);
        }


        tbb::enumerable_thread_specific<std::vector<VectorX<T>>> local_root;
        tbb::affinity_partitioner ap;


        tbb::parallel_for(tbb::blocked_range<size_t>(0,num_block), //
        [&](const tbb::blocked_range<size_t>& range)
        {
            auto& root_thread = local_root.local();
            std::vector<std::array<int,2>> initial_point_range(initial_points.size());
            VectorXi block_index;
            std::vector<VectorX<T>> block_root;
            for(int i=range.begin();i!=range.end();++i)
            {
                // if(i<64)
                // {
                //     continue;
                // }
                if(selected_span_index== std::vector<VectorXi>())
                {
                    utility::obtainDomainIndex(i,block_index,number_in_every_dim);
                }
                else
                {
                    block_index = selected_span_index[i];
                }
                for(int j=0;j<initial_points.size();j++)
                {
                    initial_point_range[j][0]= point_num_in_block[j]*block_index[j];
                    initial_point_range[j][1]= point_num_in_block[j]*(block_index[j]+1);
                }
                degenerate_finding_single_block(block_root,initial_points,initial_point_range);
                if(!block_root.empty())
                {
                    root_thread.insert(root_thread.end(), block_root.begin(), block_root.end());
                } 

                // if(i%10==0)
                // {
                    std::cout<<"degenerate case processing block "<<i<<std::endl;
                // }
            }

        },ap               
        );


        for (const auto& thread_vec : local_root) {
            root.insert(root.end(), thread_vec.begin(), thread_vec.end());
        }

    }
   

};