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

#include "../critical_point/find_unique_root.h"

template<typename T>
class Find_initial_root
{   
    
    private:
    const Block<T>* b;
    const int function_type;
    const VectorX<T> domain_min;
    const VectorX<T> domain_max;
    INRModel<T>* inr_model;
    
    T root_finding_epsilon;
    int max_itr;

    T same_root_epsilon;
    VectorXi num_of_initial_point;

    public:
    
    Find_initial_root(const VectorX<T>& core_mins_, const VectorX<T>& core_maxs_,T same_root_epsilon_, T root_finding_epsilon_, int max_itr_=50,  const int function_type_=0, const Block<T>* b_=nullptr, INRModel<T>* inr_model_=nullptr, const VectorXi& num_of_initial_point_=VectorXi())
    : domain_min(core_mins_), domain_max(core_maxs_), b(b_), function_type(function_type_), root_finding_epsilon(root_finding_epsilon_), max_itr(max_itr_), same_root_epsilon(same_root_epsilon_), inr_model(inr_model_), num_of_initial_point(num_of_initial_point_) {}
    ~Find_initial_root(){}


    void compute_gradient(VectorX<T>& p, VectorX<T>& grad)
    {
        grad.resize(p.size());
        VectorX<T> f_vector(1);
        VectorXi deriv = VectorXi::Zero(p.size());
        for(int i=0;i<p.size();i++)
        {
            deriv[i]=1;
            query_function::query_function(p, f_vector, function_type, b, deriv, inr_model);
            grad[i] = f_vector[0];///local_domain_range[i];
            deriv[i]=0;
        }
    }

    //point_func_value: the value of the point we want to extract
    bool gradient_descent(VectorX<T>& result, VectorX<T>& p, int max_itr,
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

    void test_function_value(std::vector<VectorX<T>>& root)
    {
        T min_value = std::numeric_limits<T>::max();
        T max_value = std::numeric_limits<T>::min();
        T mean_value = 0.0;
        for(int i=0;i<root.size();i++)
        {
            VectorX<T> func_value;
            query_function::query_function(root[i], func_value, function_type, b, VectorXi(), inr_model);
            if(func_value[0]<min_value)
            {
                min_value = func_value[0];
            }
            if(func_value[0]>max_value)
            {
                max_value = func_value[0];
            }
            mean_value += func_value[0];
        }
        mean_value /= root.size();
        std::cout<<"min_value "<<min_value<<std::endl;
        std::cout<<"max_value "<<max_value<<std::endl;
        std::cout<<"mean_value "<<mean_value<<std::endl;
    }


    bool root_finding(std::vector<VectorX<T>>& root, T point_func_value, std::vector<VectorXi>* span_index=nullptr) { 

        if(b!=nullptr)
        {
            std::vector<VectorX<T>> root_all;
            for(int i=0;i<span_index->size();i++)
            {
                VectorXi span_index_local = (*span_index)[i]+b->mfa->var(0).p;

                std::vector<VectorX<T>> root_local; 

                root_finding_mfa(span_index_local, root_local,point_func_value);

                root_all.insert(root_all.end(), root_local.begin(), root_local.end());


            }

            spatial_hashing::find_all_unique_root(root_all, root,same_root_epsilon);
            
        }
        else
        {
            std::cout<<"explicit root_finding_for_entire_domain"<<std::endl;
            //compute initial points all at once.
            return root_finding_for_entire_domain(root, num_of_initial_point, point_func_value);

            

        }

        return false;


    }


       // Function to find the roots of the polynomial using Newton's method
       bool root_finding_for_entire_domain(std::vector<VectorX<T>>& root, VectorXi& num_of_initial_point, T point_func_value) 
       { 
   
           root.clear();
           
           VectorXi one = VectorXi::Ones(domain_min.size());
           // int deg = (mfa_data->p-one).prod();
   
   
           std::vector<VectorX<T>> root_in_original_domain;
           
           // std::cout<<"max_iteration--"<<maxIter<<std::endl;
   
           std::vector<std::vector<T>> domain_range(num_of_initial_point.size());

           std::cout<<"num_of_initial_point "<<num_of_initial_point.transpose()<<std::endl;

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

           std::cout<<"initial_point size "<<num_initial_point_every_domain.transpose()<<std::endl;

            std::cout<<"same_root_epsilon "<<same_root_epsilon<<std::endl;


           int num_initial_point = num_initial_point_every_domain.prod();
   
   
   

           VectorXi number_in_every_domain;
           utility::obtain_number_in_every_domain(num_initial_point_every_domain,number_in_every_domain);
   



        //    for(int i=0;i<num_initial_point;++i)
           struct RootSlot {
               bool found = false;
               VectorX<T> position;
           };
           std::vector<RootSlot> root_per_initial_point(static_cast<size_t>(num_initial_point));

           tbb::parallel_for(tbb::blocked_range<size_t>(0,num_initial_point), //
           [&](const tbb::blocked_range<size_t>& range)
           {
                VectorX<T> current_initial_point(initial_point.size());
                VectorX<T> next_root;
                VectorX<T> func_value;
                VectorXi domain_index;
                for(size_t i=range.begin();i!=range.end();++i)
                {
                    utility::obtainDomainIndex(i,domain_index,number_in_every_domain);
                    for(int j=0;j<current_initial_point.size();j++)
                    {
                        current_initial_point[j]=initial_point[j][domain_index[j]];
                    }
                    if(gradient_descent(next_root, current_initial_point,max_itr,root_finding_epsilon,func_value,point_func_value))
                    {
                        root_per_initial_point[i].found = true;
                        root_per_initial_point[i].position = next_root;
                    }
                }
           });

           std::vector<VectorX<T>> combined_root;
           combined_root.reserve(static_cast<size_t>(num_initial_point));
           for(size_t i=0;i<root_per_initial_point.size();++i)
           {
               if(root_per_initial_point[i].found)
               {
                   combined_root.emplace_back(std::move(root_per_initial_point[i].position));
               }
           }
           root_per_initial_point.clear();
           root_per_initial_point.shrink_to_fit();

      
        spatial_hashing::find_all_unique_root(combined_root, root,same_root_epsilon);
        combined_root.clear();
        combined_root.shrink_to_fit();
   
        return !root.empty();
   
       }



    // Function to find the roots of the polynomial using Newton's method
    bool root_finding_mfa(VectorXi& span_index, std::vector<VectorX<T>>& root,T point_func_value) { 

        root.clear();
        
        VectorXi one = VectorXi::Ones(domain_min.size());
        // int deg = (mfa_data->p-one).prod()


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
                
            if(gradient_descent(next_root, current_initial_point,max_itr,root_finding_epsilon,func_value,point_func_value))
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

};
