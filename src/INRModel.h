#pragma once
#include <torch/torch.h>       // <- brings torch::cuda, autograd helpers, etc.
#include <torch/script.h>
#include <torch/autograd.h>
#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <type_traits>

template<typename T>
class INRModel {
public:
    const string function_name;
    VectorX<T> domain_min;
    VectorX<T> domain_max;
    VectorX<T> function_range;
    VectorXi block_num;
    VectorXi point_num_in_block;

    INRModel(const string& func_name, const std::string& model_path = "inr_base.pt"): function_name(func_name),
    device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU),
    loaded(false) {
        try {
            module = torch::jit::load(model_path);
            module.to(device);
            module.eval();
            loaded = true;
            std::cout << "Loaded INR base model (float32) on "
                      << (device.is_cuda() ? "CUDA" : "CPU") << std::endl;
        } catch (const c10::Error& e) {
            std::cerr << "Error loading INR model: " << e.what() << std::endl;
            loaded = false;
        }
        this->domain_min = domain_min_(func_name);
        this->domain_max = domain_max_(func_name);
        this->function_range = function_range_(func_name);
        this->block_num = block_num_(func_name);
        this->point_num_in_block = point_num_in_block_(func_name);

    }

    bool isLoaded() const { return loaded; }

static VectorX<T> domain_min_(const string& func_name)
    {
        VectorX<T> result(3);
        if(func_name=="quartic_potential_2")
        {
            result << -2.0,-2.0,0.0;
        }
        return result;
    } 

static VectorX<T> domain_max_(const string& func_name)
    {
        VectorX<T> result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 2.0,2.0,4.0;
        }
        return result;
    } 

static VectorX<T> function_range_(const string& func_name)
    {
        VectorX<T> result(2);
        if(func_name=="quartic_potential_2")
        {
            result <<  -2.0, 2.0;
        }
        return result;
    }

static VectorXi block_num_(const string& func_name) //number of blocks that splits the domain in each dimension
    {
        VectorXi result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 10,10,10;
        }

        return result;
    }


static VectorXi point_num_in_block_(const string& func_name) //number of initial points in each block in each dimension
    {
        VectorXi result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 5,5,5;
        }
        return result;
    }


    void query(const VectorX<T>& p_, VectorX<T>& out, const VectorXi& derivs = VectorXi()) {
        out.resize(1);
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // Convert input (float or double) -> float32 tensor
        Eigen::VectorXf p = p_.template cast<float>();

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, p_.size()},
            torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
        ).clone();
        input = input.to(device, /*non_blocking=*/false, /*copy=*/true);


        if(derivs.size()==0 || derivs.sum()==0)
        {
           // Forward pass
            torch::Tensor output = module.forward({input}).toTensor().reshape({}); // [1,1]
            out(0) = static_cast<T>(output.detach().to(torch::kCPU).item<float>());
            return;
        }

        input.set_requires_grad(true);
        // Forward pass
        torch::Tensor output = module.forward({input}).toTensor().reshape({}); // [1,1]



        if(derivs.sum()==1)
        {
            int coord = (derivs(0)==1) ? 0 : (derivs(1)==1 ? 1 : 2);
            auto g_list = torch::autograd::grad(
            /*outputs=*/{output},
            /*inputs=*/{input},
            /*grad_outputs=*/{},
            /*retain_graph=*/false,
            /*create_graph=*/false,
            /*allow_unused=*/true
            );

            torch::Tensor g_full = g_list[0];  // may be undefined if unused
            float val = 0.0f;
            if (g_full.defined()) {
                // g_full is [1,3]
                val = g_full.index({0, coord}).detach().to(torch::kCPU).item<float>();
            } else {
                // input truly not used in graph => derivative is 0
                val = 0.0f;
            }

            out(0) = static_cast<T>(val);
            return;
        }


            // Helper to map (nx,ny,nt) -> which coord (0,1,2) is still “active” at each step
        auto pick_coord = [](int nx, int ny, int nt) {
            if (nx > 0) return 0;
            if (ny > 0) return 1;
            return 2;
        };


        if(derivs.sum()==2)
        {
            int nx = derivs(0), ny = derivs(1), nt = derivs(2);
            int c1 = pick_coord(nx, ny, nt);

            // full gradient of f wrt input
            auto g_full = torch::autograd::grad({output}, {input}, {}, /*retain_graph=*/true, /*create_graph=*/true)[0]; // [1,3]
            auto g1 = g_full.index({0, c1});  // scalar: ∂f/∂x or ∂f/∂y or ∂f/∂t
            // consume that count
            if      (c1==0) --nx;
            else if (c1==1) --ny;
            else            --nt;

            // Second derivative: grad of that scalar wrt full input, then index the remaining coord
            auto g2_full = torch::autograd::grad({g1}, {input}, {}, /*retain_graph=*/false, /*create_graph=*/false)[0]; // [1,3]
            int c2 = pick_coord(nx, ny, nt);
            auto g2 = g2_full.index({0, c2});
            out(0) = g2.detach().to(torch::kCPU).item<float>();
            return;
        }

        if(derivs.sum()==3)
        {
            int nx=derivs(0), ny=derivs(1), nt=derivs(2);
            // Helper lambda to take grad of a scalar tensor 's' w.r.t. a single coordinate
            int c1 = pick_coord(nx, ny, nt);
            auto g1_full = torch::autograd::grad({output}, {input}, {}, /*retain_graph=*/true, /*create_graph=*/true)[0]; // [1,3]
            auto d1 = g1_full.index({0, c1}); // scalar
            if      (c1==0) --nx;
            else if (c1==1) --ny;
            else            --nt;

            // 2nd derivative
            int c2 = pick_coord(nx, ny, nt);
            auto g2_full = torch::autograd::grad({d1}, {input}, {}, /*retain_graph=*/true, /*create_graph=*/true)[0]; // [1,3]
            auto d2 = g2_full.index({0, c2}); // scalar
            if      (c2==0) --nx;
            else if (c2==1) --ny;
            else            --nt;

            // 3rd derivative
            int c3 = pick_coord(nx, ny, nt);
            auto g3_full = torch::autograd::grad({d2}, {input}, {}, /*retain_graph=*/false, /*create_graph=*/false)[0]; // [1,3]
            auto d3 = g3_full.index({0, c3}); // scalar

            out(0) = d3.detach().to(torch::kCPU).item<float>();
            return;
        }

        
    }


    void derivative(const VectorX<T>& p_) {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // Convert input (float or double) -> float32 tensor
        Eigen::VectorXf p = p_.template cast<float>();

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, p_.size()},
            torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
        ).clone();
        input = input.to(device, /*non_blocking=*/false, /*copy=*/true);
        input.set_requires_grad(true);


        // Forward pass
        torch::Tensor output = module.forward({input}).toTensor(); // [1,1]
        torch::Tensor cur = output.reshape({});

        // Compute first-order derivatives
        std::vector<torch::Tensor> first_order = torch::autograd::grad(
            {cur}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);

        torch::Tensor g1 = first_order.at(0);

        torch::Tensor dx = g1.index({0, 0}); // ∂f/∂x
        torch::Tensor dy = g1.index({0, 1}); // ∂f/∂y
        torch::Tensor dz = g1.index({0, 2}); // ∂f/∂z

        // Compute second-order derivatives
        std::vector<torch::Tensor> second_order_x = torch::autograd::grad(
            {dx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor g2=second_order_x.at(0);
        torch::Tensor dxx = g2.index({0, 0}); // ∂²f/∂x²
        torch::Tensor dxy = g2.index({0, 1}); // ∂²f/∂x∂y
        torch::Tensor dxt = g2.index({0, 2}); // ∂²f/∂x∂t

        std::vector<torch::Tensor> second_order_y = torch::autograd::grad(
            {dy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor gy2=second_order_y.at(0);
        torch::Tensor dyy = gy2.index({0, 1}); // ∂²f/∂y²
        torch::Tensor dyt = gy2.index({0, 2}); // ∂²f/∂y∂t

        // Compute third-order derivatives
        std::vector<torch::Tensor> third_order_xx = torch::autograd::grad(
            {dxx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor g3=third_order_xx.at(0);

        torch::Tensor dxxx = g3.index({0, 0}); // ∂³f/∂x³
        torch::Tensor dxxy = g3.index({0, 1}); // ∂³f/∂x²∂y
        torch::Tensor dxxt = g3.index({0, 2}); // ∂³f/∂x²∂t

        std::vector<torch::Tensor> third_order_yy = torch::autograd::grad(
            {dyy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor gyy3=third_order_yy.at(0);
        torch::Tensor dyyx = gyy3.index({0, 0}); // ∂³f/∂y²∂x
        torch::Tensor dyyy = gyy3.index({0, 1}); // ∂³f/∂y³
        torch::Tensor dyyt = gyy3.index({0, 2}); // ∂³f/∂y²∂t

        std::vector<torch::Tensor> third_order_xy = torch::autograd::grad(
            {dxy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/false);
        torch::Tensor gxy3=third_order_xy.at(0);
        torch::Tensor dxyt = gxy3.index({0, 2}); // ∂³f/∂x∂y∂t

        // Return the required derivatives as a custom structure or vector
        VectorX<T> grad(p_.size());

        grad(0) = static_cast<T>(dx.item<float>());
        grad(1) = static_cast<T>(dy.item<float>());
        grad(2) = static_cast<T>(dz.item<float>());

        std::cout<<"grad "<<grad.transpose()<<std::endl;
    }

private:
    torch::jit::script::Module module;
    torch::Device device{torch::kCPU};
    bool loaded = false;
};