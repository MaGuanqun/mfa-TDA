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
#include "mfa_extend.h"
#include "utility_function.h"
#include "closed_form_function.h"

namespace cp_tracking_degenerate_case
{
    template<typename T>
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

    template<typename T>
    void Hessian(VectorX<T>& p, MatrixX<T>& Hessian,const Block<T>* b=nullptr, const int function_type=0)
    {
        int n_var_dim = p.size() - 1;
        Hessian.resize(n_var_dim, n_var_dim);
        VectorXi deriv(p.size());
        VectorX<real_t> result(1);

        for(int i=0;i<n_var_dim;i++)
        {
            for(int j=i;j<n_var_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                switch (function_type)
                {
                case 0:
                    mfa_extend::recover_mfa(b, p, result, deriv);
                    break;
                case 1:
                    closed_form_function::quartic_potential(p, result, deriv);
                    break;
                default:
                    std::cerr<<"invalid function type"<<std::endl;
                    exit(0);
                    break;
                }
                
                Hessian(j,i) = result[0];
                Hessian(i,j) = result[0];
            }
        }
    }

    template<typename T>
    void partial_xt(VectorX<T>& p, VectorX<T>& result, const Block<T>* b=nullptr, const int function_type=0)
    {
        int n_var_dim = p.size() - 1;
        result.resize(n_var_dim);
        VectorXi deriv(p.size());
        VectorX<real_t> out(1);

        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[n_var_dim]+=1;
            switch (function_type)
            {
            case 0:
                mfa_extend::recover_mfa(b, p, out, deriv);
                break;
            case 1:
                closed_form_function::quartic_potential(p, out, deriv);
                break;
            default:
                std::cerr<<"invalid function type"<<std::endl;
                exit(0);
                break;
            }
            result(i) = out[0];
        }
    }


    template<typename T>
    void partial_Hessian(VectorX<T>& p, MatrixX<T>& p_Hessian, int partial_deriv_index, const Block<T>* b=nullptr, const int function_type=0)
    {
        int n_var_dim = p.size() - 1;
        p_Hessian.resize(n_var_dim, n_var_dim);
        VectorXi deriv(p.size());
        VectorX<real_t> result(1);

        for(int i=0;i<n_var_dim;i++)
        {
            for(int j=i;j<n_var_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                deriv[partial_deriv_index] += 1; // partial derivative with respect to the partial_deriv_index
                switch (function_type)
                {
                case 0:
                    mfa_extend::recover_mfa(b, p, result, deriv);
                    break;
                case 1:
                    closed_form_function::quartic_potential(p, result, deriv);
                    break;
                default:
                    std::cerr<<"invalid function type"<<std::endl;
                    exit(0);
                    break;
                }
                p_Hessian(j,i) = result[0];
                p_Hessian(i,j) = result[0];
            }
        }
    }

    template<typename T>
    void diff_Hessian(VectorX<T>& p, std::vector<MatrixX<T>>& d_Hessian, const Block<T>* b=nullptr, const int function_type=0)
    {
        d_Hessian.resize(p.size());
        for(int i=0;i<p.size();i++)
        {
            partial_Hessian(p, d_Hessian[i], i, b, function_type);       
        }
    }
    template<typename T>
    void diff_Hessian_determinant(VectorX<T>& p,const MatrixX<T>& hessian, VectorX<T>& diff, const Block<T>* b=nullptr, const int function_type=0)
    {
        MatrixX<T> adj;
        adjugate_of_symmetric(hessian, adj);
        std::vector<MatrixX<T>> diff_Hes;
        diff_Hessian(p, diff_Hes, b, function_type);

        diff.resize(p.size());
        for(int i=0;i<p.size();i++)
        {
            diff[i]=adj.cwiseProduct(diff_Hes[i]).sum();
        }
    }

    template<typename T>
    void compute_J(VectorX<T>& p, VectorX<T>& J, const Block<T>* b=nullptr, const int function_type=0)
    {
        J.resize(p.size());
        VectorXi deriv(p.size());
        int n_var_dim=p.size()-1;
        VectorX<T> result(1); 
        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;

            switch (function_type)
            {   
            case 0:
                mfa_extend::recover_mfa(b, p,result, deriv);
                break;
            case 1:
                closed_form_function::quartic_potential(p, result, deriv);
                break;
            default:
                std::cerr<<"invalid function type"<<std::endl;
                exit(0);
                break;
            }

            J[i+1] = result[0];
        }

        MatrixX<T> Hessian_f;
        Hessian(b, p, Hessian_f);
        J[0] = Hessian_f.determinant();        
    }

    template<typename T>
    void compute_J_dev_J(VectorX<T>& p, VectorX<T>& J, MatrixX<T>& Jacobian_f,const Block<T>* b=nullptr, const int function_type=0)
    {
        J.resize(p.size());
        Jacobian_f.resize(p.size(),p.size());

        VectorXi deriv(p.size());

        int n_var_dim=p.size()-1;

        VectorX<T> result(1); 


        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            switch (function_type)
            {
            case 0:
                mfa_extend::recover_mfa(b, p,result, deriv);
                break;
            case 1:
                closed_form_function::quartic_potential(p, result, deriv);
                break;
            default:
                std::cerr<<"invalid function type"<<std::endl;
                exit(0);
                break;
            }

            J[i+1] = result[0];
        }

        MatrixX<T> Hessian_f;
        Hessian(p, Hessian_f, b, function_type);
        J[0] = Hessian_f.determinant();     

        VectorX<T> diff;
        diff_Hessian_determinant(p, Hessian_f, diff, b, function_type);

        Jacobian_f.row(0) = diff.transpose();
        Jacobian_f.block(1,0,n_var_dim,n_var_dim) = Hessian_f;

        VectorX<T> partial_xt_result;
        partial_xt(p, partial_xt_result,b, function_type);
        Jacobian_f.block(1,n_var_dim,n_var_dim,1) = partial_xt_result;

    }

    template<typename T>
    bool newton(VectorX<T>& result, VectorX<T>& p, int max_itr,
                    // T d_max_square, VectorX<T>& center,
                    T degenerate_finding_epsilon, T hessian_det_epsilon, T gradient_epsilon,
                    const VectorX<T>& domain_min, const VectorX<T>& domain_max,const Block<T>* b=nullptr, const int function_type=0)
    {
        int itr_num=0;
        MatrixX<T> dev_J;
        VectorX<T> J;
        compute_J_dev_J(p,J,dev_J, b, function_type);
        if(J.squaredNorm()<degenerate_finding_epsilon*degenerate_finding_epsilon)
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

            // if((p-center).squaredNorm()>d_max_square)
            // {
            //     return false;
            // }
            if(!utility::In_Domain(p,domain_min,domain_max))
            {
                return false;
            }

            compute_J_dev_J(p,J,dev_J, b, function_type);   

            if(itr_num>0){
                if(J.squaredNorm()< degenerate_finding_epsilon * degenerate_finding_epsilon 
                && J.tail(J.size()-1).squaredNorm()<gradient_epsilon*gradient_epsilon
                ){                    
                    result = p;
                    return true;
                }
            }

            itr_num++;

        }
        return false;

    }


    template<typename T>
    bool newRoot(VectorX<T>& z, std::vector<VectorX<T>>& root_so_far, std::vector<T>& threshold)
    {
        for(int i=0;i<root_so_far.size();++i)
        {
            if(z[z.size()-1]-root_so_far[i][z.size()-1]<threshold.back())
            {
                if((z.head(z.size()-1)-root_so_far[i].head(z.size()-1)).squaredNorm()<threshold[0]*threshold[0])
                {
                    return false;
                }
            }
        }
        return true;
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
            // if(b!=nullptr)
            // {
            //     VectorXi span_num = b->mfa->var(0).tmesh.tensor_prods[0].nctrl_pts-b->mfa->var(0).p;
            //     initial_point_number = span_num.cwiseProduct(b->mfa->var(0).p+VectorXi::Ones(span_num.size()));
            // }
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


 



   // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool degenerate_finding_single_block(std::vector<VectorX<T>>& root,
        T degenerate_finding_epsilon, std::vector<T>& same_root_epsilon,
        T hessian_det_epsilon, T point_itr_threshold, T gradient_epsilon,const VectorX<T>& domain_min, const VectorX<T>& domain_max,
        std::vector<std::vector<T>>initial_point, std::vector<std::array<int,2>>& initial_point_range,
        const Block<T>* b=nullptr, const int function_type=0)
    {

        root.clear();
        VectorXi one = VectorXi::Ones(domain_min.size());
        int maxIter=50;        
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

            if(!utility::In_Domain(current_initial_point,domain_min,domain_max))
            {
                std::cout<<current_initial_point.size()<<std::endl;
                std::cout<<domain_min.transpose()<<" "<<domain_max.transpose()<<std::endl;
                std::cerr<<"initial point is not in the domain"<<std::endl;
            }


            if(newton(next_root, current_initial_point,maxIter,degenerate_finding_epsilon,hessian_det_epsilon,gradient_epsilon,domain_min,domain_max,b,function_type))
            {
        
                if(newRoot(next_root,root,same_root_epsilon))
                {        
                    root.emplace_back(next_root);
                } 
            }

        }        

        return !root.empty();

    }


    
    template<typename T>
    void degenerate_finding(std::vector<VectorX<T>>& root,
    T root_finding_epsilon, std::vector<T>& same_root_epsilon, T hessian_det_epsilon, T point_itr_threshold, T gradient_epsilon, const VectorX<T>& domain_min, const VectorX<T>& domain_max,
    const VectorXi& point_num_in_block, const VectorXi& set_block_num, const Block<real_t>* block=nullptr,std::vector<VectorXi>& selected_span_index = std::vector<VectorXi>(),const int function_type=0)
    {
        std::vector<vector<T>> initial_points;
        generate_initial_points(initial_points,domain_min,domain_max,point_num_in_block,set_block_num,block);
        

        // try to divide the initial points into different blocks. For MFA, every block is a span, for other functions, it is manually defined. So that it is easier to deduplicate

        //MFA. both spans and initial points are uniformly distributed, number of points = (degree+1)*span_num
        VectorXi number_in_every_dim;
        int num_block;
        if(block!=nullptr)
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

        std::cout<<point_num_in_block.transpose()<<" initial points in each block"<<std::endl;
        std::cout<<initial_points[0].size()<<" "<<initial_points[1].size()<<" initial points in each dimension"<<std::endl;
        std::cout<<set_block_num.transpose()<<" blocks in each dimension"<<std::endl;



        tbb::parallel_for(tbb::blocked_range<size_t>(0,num_block), //
        [&](const tbb::blocked_range<size_t>& range)
        {
            auto& root_thread = local_root.local();
            std::vector<std::array<int,2>> initial_point_range(initial_points.size());
            VectorXi block_index;
            std::vector<VectorX<T>> block_root;
            for(int i=range.begin();i!=range.end();++i)
            {
                if(block==nullptr)
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
                degenerate_finding_single_block(block_root,root_finding_epsilon,same_root_epsilon,hessian_det_epsilon,point_itr_threshold,gradient_epsilon,domain_min,domain_max,initial_points,initial_point_range,block,function_type);
                if(!block_root.empty())
                {
                    root_thread.insert(root_thread.end(), block_root.begin(), block_root.end());
                } 
            }

        },ap               
        );


        for (const auto& thread_vec : local_root) {
            root.insert(root.end(), thread_vec.begin(), thread_vec.end());
        }

    }


    // template<typename T>
    // bool degenerate_finding(Block<real_t>* block, std::vector<std::vector<VectorXi>>& span_index, 
    // std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
    //     int current_index,
    //     T root_finding_epsilon, std::vector<T>& same_root_epsilon, T hessian_det_epsilon, T point_itr_threshold, T gradient_epsilon) //2^n+1 initial points) 
    // {

    //     for(auto i=0;i<block->mfa->nvars();++i)
    //     {
    //         if(degenerate_finding(block,span_index[i][current_index], root,
    //         root_finding_epsilon,same_root_epsilon,hessian_det_epsilon, point_itr_threshold,gradient_epsilon))
    //         {
    //             return true;
    //         }
    //     }

    //     return false;
    // }    

}