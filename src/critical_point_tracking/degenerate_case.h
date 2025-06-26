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
    void Hessian(const Block<T>* b, VectorX<T>& p, MatrixX<T>& Hessian)
    {
        int n_var_dim = b->dom_dim - 1;
        Hessian.resize(n_var_dim, n_var_dim);
        VectorXi deriv(b->dom_dim);
        VectorX<real_t> result(1);

        for(int i=0;i<n_var_dim;i++)
        {
            for(int j=i;j<n_var_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                mfa_extend::recover_mfa(b, p, result, deriv);
                Hessian(j,i) = result[0];
                Hessian(i,j) = result[0];
            }
        }
    }

    template<typename T>
    void partial_xt(const Block<T>* b, VectorX<T>& p, VectorX<T>& result)
    {
        int n_var_dim = b->dom_dim - 1;
        result.resize(n_var_dim);
        VectorXi deriv(b->dom_dim);
        VectorX<real_t> out(1);

        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            deriv[n_var_dim]+=1;
            mfa_extend::recover_mfa(b, p, out, deriv);
            result(i) = out[0];
        }
    }


    template<typename T>
    void partial_Hessian(const Block<T>* b, VectorX<T>& p, MatrixX<T>& p_Hessian, int partial_deriv_index)
    {
        int n_var_dim = b->dom_dim - 1;
        p_Hessian.resize(n_var_dim, n_var_dim);
        VectorXi deriv(b->dom_dim);
        VectorX<real_t> result(1);

        for(int i=0;i<n_var_dim;i++)
        {
            for(int j=i;j<n_var_dim;j++)
            {
                deriv.setZero();
                deriv[i]+=1;
                deriv[j]+=1;
                deriv[partial_deriv_index] += 1; // partial derivative with respect to the partial_deriv_index
                mfa_extend::recover_mfa(b, p, result, deriv);
                p_Hessian(j,i) = result[0];
                p_Hessian(i,j) = result[0];
            }
        }
    }

    template<typename T>
    void diff_Hessian(const Block<T>* b, VectorX<T>& p, std::vector<MatrixX<T>>& d_Hessian)
    {
        d_Hessian.resize(b->dom_dim);
        for(int i=0;i<b->dom_dim;i++)
        {
            partial_Hessian(b, p, d_Hessian[i], i);       
        }
    }
    template<typename T>
    void diff_Hessian_determinant(const Block<T>* b, VectorX<T>& p,const MatrixX<T>& hessian, VectorX<T>& diff)
    {
        MatrixX<T> adj;
        adjugate_of_symmetric(hessian, adj);
        std::vector<MatrixX<T>> diff_Hes;
        diff_Hessian(b, p, diff_Hes);

        diff.resize(b->dom_dim);
        for(int i=0;i<b->dom_dim;i++)
        {
            diff[i]=adj.cwiseProduct(diff_Hes[i]).sum();
        }
    }

    template<typename T>
    void compute_J(const Block<T>* b, VectorX<T>& p, VectorX<T>& J)
    {
        J.resize(b->dom_dim);
        VectorXi deriv(b->dom_dim);
        int n_var_dim=b->dom_dim-1;
        VectorX<T> result(1); 
        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            mfa_extend::recover_mfa(b, p,result, deriv);
            J[i+1] = result[0];
        }

        
        MatrixX<T> Hessian_f;
        Hessian(b, p, Hessian_f);
        J[0] = Hessian_f.determinant();        
    }

    template<typename T>
    void compute_J_dev_J(const Block<T>* b, VectorX<T>& p, VectorX<T>& J, MatrixX<T>& Jacobian_f)
    {
        J.resize(b->dom_dim);
        Jacobian_f.resize(b->dom_dim,b->dom_dim);

        VectorXi deriv(b->dom_dim);

        int n_var_dim=b->dom_dim-1;

        VectorX<T> result(1); 


        for(int i=0;i<n_var_dim;i++)
        {
            deriv.setZero();
            deriv[i]+=1;
            mfa_extend::recover_mfa(b, p,result, deriv);
            J[i+1] = result[0];
        }

        MatrixX<T> Hessian_f;
        Hessian(b, p, Hessian_f);
        J[0] = Hessian_f.determinant();     

        VectorX<T> diff;
        diff_Hessian_determinant(b, p, Hessian_f, diff);

        Jacobian_f.row(0) = diff.transpose();
        Jacobian_f.block(1,0,n_var_dim,n_var_dim) = Hessian_f;

        VectorX<T> partial_xt_result;
        partial_xt(b, p, partial_xt_result);
        Jacobian_f.block(1,n_var_dim,n_var_dim,1) = partial_xt_result;
        // VectorX<T> Hessian_f(b->dom_dim*(b->dom_dim+1)/2-1);
        // //2d: f_xx,xy,xt,yy,yt
        // //[0 1 2]
        // //[1 3 4]
        // //[2 4 5]
        // //3d: f_xx,xy,xz,xt,yy,yz,yt,zz,zt
        // //[0 1 2 3]
        // //[1 4 5 6]
        // //[2 5 7 8]

        // int k=0;
        // for(int i=0;i<n_var_dim;i++)
        // {
        //     for(int j=i;j<b->dom_dim;j++)
        //     {
        //         deriv.setZero();
        //         deriv[i]+=1;
        //         deriv[j]+=1;
        //         mfa_extend::recover_mfa(b, p, result, deriv);
        //         Hessian_f[k] = result[0];
        //         k++;
        //     }
        // }

        // VectorX<T> f(n_var_dim);
        // for(int i=0;i<n_var_dim;i++)
        // {
        //     deriv.setZero();
        //     deriv[i]+=1;
        //     mfa_extend::recover_mfa(b, p,result, deriv);
        //     f[i] = result[0];
        // }

        // if(n_var_dim==2)
        // {
        //     J[0]=Hessian_f[0]*Hessian_f[3]-Hessian_f[1]*Hessian_f[1];
        // }
        // else if(n_var_dim==3)
        // {
        //     J[0]= Hessian_f[0]*Hessian_f[4]*Hessian_f[7] + 2.0 * Hessian_f[1]*Hessian_f[5]*Hessian_f[2] - Hessian_f[2]*Hessian_f[4]*Hessian_f[2] - Hessian_f[1]*Hessian_f[1]*Hessian_f[7] - Hessian_f[0]*Hessian_f[5]*Hessian_f[5];
        // }
        // else
        // {
        //     std::cerr<<"Error: only support 2D and 3D case"<<std::endl;
        // }

        // J.segment(1,n_var_dim)=f;

        // VectorX<T> f_3(4);
        // //F_xxx,xxy,xyy,yyy
        // deriv.setZero();
        // deriv[0]=3;
        // for(int i=0;i<4;i++)
        // {
        //     mfa_extend::recover_mfa(b, p,result, deriv);
        //     f_3[i] = result[0];
        //     deriv[0]--;
        //     deriv[1]++;
        // }

        // deriv.setZero();
        // deriv[2]=1;
        // VectorX<T> f_t_3(3);
        // //F_txx,txy,tyy
        // for(int i=0;i<n_var_dim;i++)
        // {
        //     for(int j=i;j<n_var_dim;j++)
        //     {
        //         deriv[i]+=1;
        //         deriv[j]+=1;
        //         mfa_extend::recover_mfa(b, p,result, deriv);
        //         f_t_3[n_var_dim*i+j] = result[0];
        //         deriv[i]--;
        //         deriv[j]--;
        //     }
        // }

        // //f_xx,xy,xt,yy,yt
        // //F_xxx,xxy,xyy,yyy
        // //F_txx,txy,tyy
        // Jacobian_f(0,0) = Hessian_f[3]*f_3[0]+Hessian_f[0]*f_3[2]-2.0*Hessian_f[1]*f_3[1];
        // Jacobian_f(0,1) = Hessian_f[0]*f_3[3]+Hessian_f[3]*f_3[1]-2.0*Hessian_f[1]*f_3[2];
        // Jacobian_f(0,2) = Hessian_f[3]*f_t_3[0]+Hessian_f[0]*f_t_3[2]-2.0*Hessian_f[1]*f_t_3[1];

        // Jacobian_f.row(1) = Hessian_f.segment(0,3).transpose();
        // Jacobian_f(2,0) = Hessian_f[1];
        // Jacobian_f(2,1) = Hessian_f[3];
        // Jacobian_f(2,2) = Hessian_f[4];

    }

    template<typename T>
    bool newton(const Block<T>* b,VectorX<T>& result, VectorX<T>& p, int max_itr, std::vector<std::vector<T>>& span_range,
                    T d_max_square, VectorX<T>& center,
                    T degenerate_finding_epsilon, T hessian_det_epsilon)
    {
        int itr_num=0;
        MatrixX<T> dev_J;
        VectorX<T> J;
        compute_J_dev_J(b,p,J,dev_J);
        if(J.squaredNorm()<degenerate_finding_epsilon*degenerate_finding_epsilon)
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
            T determinant = dev_J.determinant();
            if(std::abs(determinant) < hessian_det_epsilon)
            {
                return false;                
            }
            else
            {  
                p -= dev_J.colPivHouseholderQr().solve(J); 
            }

            if((p-center).squaredNorm()>d_max_square)
            {
                return false;
            }
            if(!utility::In_Domain(p,b->core_mins,b->core_maxs))
            {
                return false;
            }

            compute_J_dev_J(b,p,J,dev_J);   

            if(itr_num>0){
                if(J.squaredNorm()< degenerate_finding_epsilon * degenerate_finding_epsilon){//|| (result-p).squaredNorm()<ROOT_FINDING_EPSILON*ROOT_FINDING_EPSILON
              
                    if(!utility::InBlock(span_range,p))
                    {
                        return false;
                    }
                    
                    result = p;
                    return true;
                }
            }

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


   // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool degenerate_finding(const Block<T>* b, VectorXi& span_index, std::vector<VectorX<T>>& root,
        T degenerate_finding_epsilon, T same_root_epsilon,
        T hessian_det_epsilon) { 

        root.clear();
        
        VectorXi one = VectorXi::Ones(b->mfa->var(0).p.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=100;

        int distance_stop_itr = 5;

        
        auto domain_range = b->core_maxs-b->core_mins;
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        std::vector<std::vector<T>> span_range(b->dom_dim);


        VectorX<T> center(3);
        for(int i=0;i<span_range.size();++i)
        {    
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+b->core_mins[i]);
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+b->core_mins[i]);

            center[i]=(span_range[i][0]+span_range[i][1])*0.5;
        }   



        // compute distance to terminate iteration
        T d_max_square=0;
        for(auto i=span_range.begin();i!=span_range.end();++i)
        {  
            d_max_square+=((*i)[1]-(*i)[0])*((*i)[1]-(*i)[0]);
        }

        d_max_square*=distance_stop_itr * distance_stop_itr; //d^2=(2*diagonal of span)^2

        std::vector<std::vector<T>>initial_point;


        utility::compute_initial_points(initial_point,b->mfa->var(0).p,span_range);

        VectorXi num_initial_point_every_domain(initial_point.size());
        for(int i=0;i<num_initial_point_every_domain.size();i++)
        {
            num_initial_point_every_domain[i]=initial_point[i].size();
        }

        int num_initial_point = num_initial_point_every_domain.prod();

        VectorX<T> next_root; 


        VectorXi domain_index;
        VectorXi number_in_every_domain;
        VectorX<T> current_initial_point(initial_point.size());
        utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);

        std::vector<VectorX<T>> root_in_original_domain;
        // std::cout<<"num_initial_point "<<num_initial_point<<std::endl;

        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<initial_point.size();j++)
            {
                current_initial_point[j]=initial_point[j][domain_index[j]];
            } 

            // std::cout<<"initial point "<<i<<" "<<  current_initial_point.transpose()<<std::endl;    

            if(newton(b, next_root, current_initial_point,maxIter,span_range,d_max_square,center,degenerate_finding_epsilon,hessian_det_epsilon))
            {
        
                if(newRoot(next_root,root_in_original_domain,same_root_epsilon))
                {        
                    root_in_original_domain.emplace_back(next_root);
                } 
            }

        }

        if(!root_in_original_domain.empty())
        {
            root.insert(root.end(),root_in_original_domain.begin(),root_in_original_domain.end());
        }

        

        return !root.empty();

    }

    template<typename T>
    bool degenerate_finding(Block<real_t>* block, std::vector<std::vector<VectorXi>>& span_index, 
    std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
        int current_index,
        T root_finding_epsilon, T same_root_epsilon, T hessian_det_epsilon) //2^n+1 initial points) 
    {

        for(auto i=0;i<block->mfa->nvars();++i)
        {
            if(degenerate_finding(block,span_index[i][current_index], root,
            root_finding_epsilon,same_root_epsilon,hessian_det_epsilon))
            {
                return true;
            }
        }

        return false;
    }    

}