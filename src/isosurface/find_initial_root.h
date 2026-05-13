// find roots in the initial domain
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <map>

#include <mfa/mfa.hpp>



#include "opts.h"

#include "block.hpp"
#include "utility_function.h"
#include "query_function.h"
#include "INRModel.h"

#include"find_unique_root.h"

template<typename T>
class find_initial_root
{   
    
    private:

    const Block<T>* b;
    const int function_type;
    const VectorX<T> domain_min;
    const VectorX<T> domain_max;
    INRModel<T>* inr_model;
    
    T root_finding_epsilon;
    int max_itr;

    std::vector<T> same_root_epsilon;

    public:
    
    Tracking_degenerate_case(const VectorX<T>& core_mins_, const VectorX<T>& core_maxs_, T degenerate_finding_epsilon_, T gradient_epsilon_,std::vector<T> same_root_epsilon_, int max_itr_=50,  const int function_type_=0, const Block<T>* b_=nullptr, INRModel<T>* inr_model_=nullptr)
    : domain_min(core_mins_), domain_max(core_maxs_), b(b_), function_type(function_type_), degenerate_finding_epsilon(degenerate_finding_epsilon_), gradient_epsilon(gradient_epsilon_), max_itr(max_itr_), same_root_epsilon(same_root_epsilon_), inr_model(inr_model_)
    {
        std::cout<<"degenerate case tracking initialized "<<std::endl;
        std::cout<<"domain min "<<domain_min.transpose()<<std::endl;
        std::cout<<"domain max "<<domain_max.transpose()<<std::endl;
    }
    ~Tracking_degenerate_case(){}

    template<typename T>
    void compute_gradient(VectorX<T>& p, VectorX<T>& grad)
    {
        grad.resize(p.size());
        VectorX<T> f_vector(1);
        VectorXi deriv = VectorXi::Zero(p.size());
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            query_function::query_function(p, f_vector, function_type, b, deriv, inr_model);
            // mfa_extend::recover_mfa_selected(mfa,mfa_data,span_index, p,deriv,weights,f_vector,domain_min,domain_range);
            grad[i] = f_vector[0];///local_domain_range[i];
            deriv[i]=0;
        }
    }

    template<typename T>
    //point_func_value: the value of the point we want to extract
    bool gradient_descent(VectorX<T>& result, VectorX<T>& p, int max_itr,
                //std::vector<std::vector<T>>& span_range,
                //T d_max_square, //VectorX<T>& center,
                T root_finding_epsilon, VectorX<T>& function_value, T point_func_value)
    {

        query_function::query_function(p, function_value, function_type, b, VectorXi(), inr_model);
        function_value[0]-=point_func_value;

        if(std::abs(function_value[0])<root_finding_epsilon)
        {
            result = p;
            return true;
        }

        VectorX<T> grad;
        int i=0;
        while(i<max_itr)
        {
            compute_gradient(p, grad);
            T step_size=function_value[0]/grad.squaredNorm();
            // printf("step_size %f\n",step_size);
            p-= step_size*grad;
            // if((p-center).squaredNorm()>d_max_square)
            // {
            //     return false;
            // }
            if(!utility::In_Domain(p,domain_min,domain_max))
            {
                return false;
            }
            query_function::query_function(p, function_value, function_type, b, VectorXi(), inr_model);

            function_value[0]-=point_func_value;
            

            if(std::abs(function_value[0])<root_finding_epsilon)
            {
                result = p;
                function_value[0]+=point_func_value;
                return true;
            }
            i++;

        }

        return false;

    }

    template<typename T>
    bool root_finding(std::vector<VectorX<T>>& root,
        T root_finding_grad_epsilon, std::vector<T>& function_value, T point_func_value, T same_root_epsilon) { 

        if(b==nullptr)
        {
            
        }
        else
        {
            //compute initial points all at once.
            root_finding_for_entire_domain(root, num_of_initial_point, root_finding_epsilon, function_value, point_func_value, same_root_epsilon);

        }


    }


       // Function to find the roots of the polynomial using Newton's method
       template<typename T>
       bool root_finding_for_entire_domain(std::vector<VectorX<T>>& root, VectorXi& num_of_initial_point,
           T root_finding_epsilon, std::vector<T>& function_value, T point_func_value, T same_root_epsilon) { 
   
           function_value.clear();
           root.clear();
           
           VectorXi one = VectorXi::Ones(domain_min.size());
           // int deg = (mfa_data->p-one).prod();
   
           int maxIter=30;
   
   
           std::vector<VectorX<T>> root_in_original_domain;
           
           // std::cout<<"max_iteration--"<<maxIter<<std::endl;
   
           std::vector<std::vector<T>> span_range(span_index.size());
           
           std::vector<std::vector<T>> domain_range(num_of_initial_point.size());

           for(int i=0;i<num_of_initial_point.size();i++)
           {
               domain_range[i].emplace_back(domain_min[i]);
               domain_range[i].emplace_back(domain_max[i]);
           }
        
           std::vector<std::vector<T>>initial_point;
   
           utility::compute_initial_points2(initial_point,num_of_initial_point,domain_range);

   
           VectorXi num_initial_point_every_domain(initial_point.size());
           for(int i=0;i<num_initial_point_every_domain.size();i++)
           {
               num_initial_point_every_domain[i]=initial_point[i].size();
           }
   
           int num_initial_point = num_initial_point_every_domain.prod();
   
           VectorX<T> next_root; 
   
   

           VectorXi number_in_every_domain;
           VectorX<T> current_initial_point(initial_point.size());
           utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);
   



        //    for(int i=0;i<num_initial_point;++i)
           tbb::enumerable_thread_specific<std::vector<VectorX<T>>> local_root;
           tbb::affinity_partitioner ap;
           tbb::parallel_for(tbb::blocked_range<size_t>(0,num_initial_point), //
           [&](const tbb::blocked_range<size_t>& range)
           {

            auto& root_thread = local_root.local();
            VectorX<T> func_value;
            VectorXi domain_index;
            std::vector<VectorX<T>> temp_root;
            for(size_t i=range.begin();i!=range.end();++i)
            {
               utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
               for(int j=0;j<current_initial_point.size();j++)
               {
                   current_initial_point[j]=initial_point[j][domain_index[j]];
               }        

   
               if(gradient_descent(next_root, current_initial_point,maxIter,root_finding_epsilon,func_value,point_func_value))
               {

                    temp_root.emplace_back(next_root);
                   // std::cout<<"is a new root "<<std::endl;                
               }    
            }         
                      
            
            spatial_hashing::find_all_unique_root(temp_root, root_thread,same_root_epsilon);
            temp_root.clear();
            temp_root.shrink_to_fit();

           },ap               
        );

        std::vector<VectorX<T>> combined_root;
        for (const auto& thread_vec : local_root) {
            combined_root.insert(combined_root.end(), thread_vec.begin(), thread_vec.end());
        }
        spatial_hashing::find_all_unique_root(combined_root, root,same_root_epsilon);
        combined_root.clear();
        combined_root.shrink_to_fit();
   
        return !root.empty();
   
       }



    // Function to find the roots of the polynomial using Newton's method
    template<typename T>
    bool root_finding(VectorXi& span_index, std::vector<VectorX<T>>& root,
        T root_finding_grad_epsilon, std::vector<T>& function_value, T point_func_value, T same_root_epsilon) { 

        function_value.clear();
        root.clear();
        
        VectorXi one = VectorXi::Ones(domain_min.size());
        // int deg = (mfa_data->p-one).prod();

        int maxIter=30;


        std::vector<VectorX<T>> root_in_original_domain;
        
        // std::cout<<"max_iteration--"<<maxIter<<std::endl;

        std::vector<std::vector<T>> span_range(span_index.size());
        
        auto domain_range = domain_max-domain_min;

        VectorX<T> center(span_index.size());
        for(int i=0;i<span_index.size();++i)
        {    
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]]*domain_range[i]+domain_min[i]);
            span_range[i].emplace_back(b->mfa->var(0).tmesh.all_knots[i][span_index[i]+1]*domain_range[i]+domain_min[i]);
        }   



        std::vector<std::vector<T>>initial_point;

        // degree + 1 initial points
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

        VectorX<T> func_value;

        for(int i=0;i<num_initial_point;++i)
        {

            utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
            for(int j=0;j<current_initial_point.size();j++)
            {
                current_initial_point[j]=initial_point[j][domain_index[j]];
            }        
            // current_initial_point=ini_p;


            // std::cout<<"intial point "<< current_initial_point.transpose()<<std::endl;
            // std::cout<< "initial_point "<<i<<" "<<  current_initial_point.transpose()<<std::endl;

            if(gradient_descent(b, next_root, current_initial_point,maxIter,span_range,root_finding_grad_epsilon,func_value,point_func_value))
            {
                bool duplicate = false;
                for(auto i=root.begin();i!=root.end();++i)
                {
                    if(((*i)-next_root).squaredNorm()<same_root_epsilon*same_root_epsilon)
                    {
                        duplicate = true;
                        break;
                    }
                }
                if(!duplicate)
                {
                    root.emplace_back(next_root);
                }
                // std::cout<<"is a new root "<<std::endl;                
            }                        
        }


        return !root.empty();

    }



    template<typename T>
    bool root_finding(std::vector<std::vector<VectorXi>>& span_index, 
    std::vector<VectorX<T>>& root,//std::vector<int>& multi_of_root,
        //MatrixX<T>&             ctrl_pts,   //control points of first derivative
        int current_index,
        T root_finding_epsilon, std::vector<T>& function_value, T point_func_value,T same_root_epsilon) //2^n+1 initial points) 
    {

        for(auto i=0;i<block->mfa->nvars();++i)
        {
            VectorXi span_index_local = span_index[i][current_index]+block->mfa->var(0).p;

            if(root_finding(block,span_index_local, root,
            root_finding_epsilon,function_value,point_func_value,same_root_epsilon))
            {
                return true;
            }
        }

        return false;
    }
}
