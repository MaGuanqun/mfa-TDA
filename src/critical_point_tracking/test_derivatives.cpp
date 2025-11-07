// #include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

// #include <diy/master.hpp>
// #include <diy/reduce-operations.hpp>
// #include <diy/decomposition.hpp>
// #include <diy/assigner.hpp>
// #include <diy/io/block.hpp>

#include <chrono>


#include <tbb/tbb.h>
#include <atomic>

// #include "opts.h"

// #include "block.hpp"


#include <fstream>
#include <iostream>
#include <Eigen/Dense>

#include "degenerate_case.h"

// #include "tracking_utility.h"

#include "spatial_hashing_spatial_temporal.h"
#include "closed_form_function.h"
#include "../INRModel.h"

using namespace std;


namespace {
    tbb::global_control globalControl(tbb::global_control::max_allowed_parallelism, 1);
}



struct ThirdOrder {
    // T[i][j][k] = ∂^3 f / (∂x_i ∂x_j ∂x_k)  (fully symmetric)
    std::vector<std::vector<std::vector<float>>> T;
};

struct Derivatives {
    Eigen::VectorXf  grad;   // size = n
    Eigen::MatrixXf  H;      // n x n
    ThirdOrder       third;  // n x n x n
};

// After:
template<typename T>
inline float eval_func(INRModel<T>& model, const Eigen::Matrix<T,-1,1>& p) {
    Eigen::Matrix<T,-1,1> out;
    model.query(p, out);            // query is non-const
    return static_cast<float>(out(0));
}

// Numerical derivatives; vectors are T-typed to match INRModel<T>
template<typename T>
Derivatives numerical_derivatives(
    INRModel<T>& model,
    const Eigen::Matrix<T,-1,1>& p,
    const Eigen::Matrix<T,-1,1>& h)
{
    const int n = static_cast<int>(p.size());
    Derivatives d;
    d.grad = Eigen::VectorXf::Zero(n);
    d.H    = Eigen::MatrixXf::Zero(n, n);
    d.third.T.resize(n, std::vector<std::vector<float>>(n, std::vector<float>(n, 0.0)));

    const float f0 = eval_func(model, p);

    // First derivatives
    for (int i = 0; i < n; ++i) {
        Eigen::Matrix<T,-1,1> pp = p, pm = p;
        pp(i) += h(i);
        pm(i) -= h(i);
        const float fp = eval_func(model, pp);
        const float fm = eval_func(model, pm);
        d.grad(i) = (fp - fm) / (2.0 * static_cast<float>(h(i)));
        std::cout<<"grad "<<i<<" : "<<d.grad(i)<<std::endl;
    }

    // Second derivatives (pure + mixed)
    for (int i = 0; i < n; ++i) {
        Eigen::Matrix<T,-1,1> pp = p, pm = p;
        pp(i) += h(i);
        pm(i) -= h(i);
        const float fp = eval_func(model, pp);
        const float fm = eval_func(model, pm);
        d.H(i,i) = (fp - 2.0*f0 + fm) / std::pow(static_cast<float>(h(i)), 2);

        for (int j = i+1; j < n; ++j) {
            Eigen::Matrix<T,-1,1> ppp = p, ppm = p, pmp = p, pmm = p;
            ppp(i) += h(i); ppp(j) += h(j);
            ppm(i) += h(i); ppm(j) -= h(j);
            pmp(i) -= h(i); pmp(j) += h(j);
            pmm(i) -= h(i); pmm(j) -= h(j);
            const float val =
                ( eval_func(model, ppp) - eval_func(model, ppm)
                - eval_func(model, pmp) + eval_func(model, pmm))
                / (4.0 * static_cast<float>(h(i)) * static_cast<float>(h(j)));
            d.H(i,j) = d.H(j,i) = val;
        }
    }

    // Third derivatives (pure)
    for (int i = 0; i < n; ++i) {
        Eigen::Matrix<T,-1,1> p2p = p, p1p = p, p1m = p, p2m = p;
        p2p(i) += 2*h(i);
        p1p(i) +=   h(i);
        p1m(i) -=   h(i);
        p2m(i) -= 2*h(i);
        const float f2p = eval_func(model, p2p);
        const float f1p = eval_func(model, p1p);
        const float f1m = eval_func(model, p1m);
        const float f2m = eval_func(model, p2m);
        d.third.T[i][i][i] = (f2m - 2*f1m + 2*f1p - f2p)
                           / std::pow(static_cast<float>(h(i)), 3)
                           / 2.0;
    }

    // Third derivatives (two-same-one-different)
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        if (i == j) continue;

        auto second_i = [&](const Eigen::Matrix<T,-1,1>& q)->float {
            Eigen::Matrix<T,-1,1> qp = q, qm = q;
            qp(i) += h(i);
            qm(i) -= h(i);
            const float fp = eval_func(model, qp);
            const float fm = eval_func(model, qm);
            const float f0q= eval_func(model, q);
            return (fp - 2.0*f0q + fm) / std::pow(static_cast<float>(h(i)), 2);
        };

        Eigen::Matrix<T,-1,1> q_plus = p, q_minus = p;
        q_plus(j)  += h(j);
        q_minus(j) -= h(j);
        const float val = (second_i(q_plus) - second_i(q_minus))
                         / (2.0 * static_cast<float>(h(j)));

        d.third.T[i][i][j] = val;
        d.third.T[i][j][i] = val;
        d.third.T[j][i][i] = val;
    }

    // Fully mixed d^3/(dx dy dz) for 3D
    if (n == 3) {
        float sum = 0.0;
        for (int sx : {-1, 1})
            for (int sy : {-1, 1})
                for (int sz : {-1, 1}) {
                    Eigen::Matrix<T,-1,1> q = p;
                    q(0) += static_cast<T>(sx) * h(0);
                    q(1) += static_cast<T>(sy) * h(1);
                    q(2) += static_cast<T>(sz) * h(2);
                    sum += (sx*sy*sz) * eval_func(model, q);
                }
        const float denom = 8.0 * static_cast<float>(h(0))*static_cast<float>(h(1))*static_cast<float>(h(2));
        const float v = sum / denom;
        d.third.T[0][1][2] = d.third.T[0][2][1] = d.third.T[1][0][2] =
        d.third.T[1][2][0] = d.third.T[2][0][1] = d.third.T[2][1][0] = v;
    }

    // d.grad.reverseInPlace();

    // d.H.row(0).swap(d.H.row(2));
    // d.H.col(0).swap(d.H.col(2));  // do both to preserve symmetry/semantics

    // model.convert_gradient_to_domain(d.grad);
    // model.convert_hessian_to_domain(d.H);

    return d;
}



int main(int argc, char** argv)
{

    // setenv("CUDA_VISIBLE_DEVICES", "", /*overwrite=*/1);


    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    // default command line arguments
    //int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    //int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help
    // get command line arguments
    opts::Options ops;


    float shrink_factor = 0.5; // shrink factor for RKF45 as the minimum shrink factor
    // string input_sample_point_number = "100-100";

    string degenerate_point_file = "degenerate_point.dat";
    
    float J_threshold = 1e-5;

    int correction_max_itr = 30;
    float time_step = 1e-3;
    float spatial_step_size = 1.0;

    real_t hessian_threshold = 1e-20;

    string input_shrink_ratio = "0-1-0-1-0-1";

    real_t point_itr_threshold = 0.5;

    real_t grad_epsilon = 1e-8;
    string input_function_name="quartic_potential";
    
    int max_itr=50;
    string input_model="";

        string trace_file="";
    string output_csv_path = "gradients.csv";

    ops >> opts::Option('f', "input_function_name",  input_function_name,  " diy input function name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('m', "input_model", input_model, " input INR model");
    ops >> opts::Option('t', "trace_file", trace_file, " file name of critical point trajectory");
    ops >> opts::Option('o', "output_csv", output_csv_path, " output CSV file for gradients");
  
    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }


    std::istringstream iss(input_shrink_ratio);
    std::vector<float> shrink_ratio;
    float number;
    std::string token;
    while (std::getline(iss, token, '-')) {
        std::istringstream tokenStream(token);
        if (tokenStream >> number) {
            shrink_ratio.push_back(number);
        }
    }  

    printf("%s\n", input_model.c_str());
    std::ifstream model_file(input_model);
    if (!model_file) {
        std::cerr << "Error: Input model file '" << input_model << "' does not exist." << std::endl;
        return 1;
    }

    INRModel<float> inr_model(input_function_name,input_model);
    int function_type=-1; 

    auto start_time = std::chrono::high_resolution_clock::now();
        

        Eigen::VectorXf local_domain_range=inr_model.domain_max-inr_model.domain_min;
        VectorXf core_maxs = inr_model.domain_max;
        VectorXf core_mins = inr_model.domain_min;

        VectorXi span_num = inr_model.block_num;


        std::vector<VectorX<float>> root; //the inner vector store the root in a span


        VectorXi point_num_in_block = inr_model.point_num_in_block; //number of initial points in a block


        Eigen::VectorXf p(3);
        // p<<-0.292091, 0.30837, 1.11537;
        // p<<0.288695, 0.00458166,0.016144;
        // p<<30, 59,20;
        p << 0.5,0.5,0.5;

        Eigen::VectorXf h(3);
        h << 1e-2, 1e-2, 1e-2;
        h[0]*=inr_model.domain_range[2]/2;
        h[1]*=inr_model.domain_range[1]/2;
        h[2]*=inr_model.domain_range[0]/2;
        // // Step sizes: choose small but not too small; scale by |p_i|+1 to reduce cancellation.
        // float eps = std::cbrt(std::numeric_limits<float>::epsilon()); // ~1.5e-5
        // float hx = eps * std::max(1.0, std::abs(p.x()));
        // float hy = eps * std::max(1.0, std::abs(p.y()));
        // float hz = eps * std::max(1.0, std::abs(p.z()));

        std::cout<< "query point in [-1,1]: "<< p.transpose() <<std::endl;
        Derivatives d = numerical_derivatives(inr_model, p, h);

        // d.H /= 4.0;
        // d



        VectorXf grad(3);
        MatrixXf hessian(3,3);
        VectorXf third_derivative_spatial(4);//[fxxx, fxxy, fxyy, fyyy]
        VectorXf third_derivative_time(3);//[fxxt, fxyt, fyyt]

        inr_model.query_up_to_third_derivative(p,grad,hessian,third_derivative_spatial,third_derivative_time);

        std::cout << "grad:\n" << d.grad.transpose()<<" with "<< grad.transpose() << "\n\n";
        std::cout << "Hessian:\n" << d.H << "\n\n";
        std::cout << "Hessian from INRModel:\n" << hessian << "\n\n";
        std::cout<< d.H - hessian <<std::endl;

        VectorXf third_derivative_spatial_test(4);
        third_derivative_spatial_test << d.third.T[0][0][0], d.third.T[0][0][1], d.third.T[0][1][1], d.third.T[1][1][1];
        third_derivative_spatial_test /= 8.0;
        VectorXf third_derivative_time_test(3);
        third_derivative_time_test << d.third.T[0][0][2], d.third.T[0][1][2], d.third.T[1][1][2];
        third_derivative_time_test /= 8.0;

        std::cout<< "third derivative spatial:\n" << third_derivative_spatial.transpose()<< "\n";
        std::cout<< third_derivative_spatial_test.transpose()<< "\n\n";
        std::cout<< "third derivative time:\n" << third_derivative_time.transpose()<< "\n";
        std::cout<< third_derivative_time_test.transpose()<< "\n\n";


        // VectorXf p_test(3);
        // p_test<<0.5,0.5,0.5;
        // VectorXf result;
        // inr_model.query(p_test,result);
        // std::cout<<"test query "<<result.transpose()<<std::endl;
        // VectorXi deriv(3);
        // deriv<<1,1,0;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test second derivative "<<result.transpose()<<std::endl;
        // deriv<<1,0,0;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test first derivative " <<result.transpose()<<std::endl;
        // deriv<<1,1,1;
        // inr_model.query(p_test,result,deriv);
        // std::cout<<"test third derivative " <<result.transpose()<<std::endl;

        // inr_model.derivative(p_test);


}