#pragma once
#include <torch/torch.h>       // <- brings torch::cuda, autograd helpers, etc.
#include <torch/script.h>
#include <torch/autograd.h>
#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <type_traits>

class INRModel {
public:
    const string function_name;

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
    }

    bool isLoaded() const { return loaded; }

static VectorXd domain_min(const string& func_name)
    {
        VectorXd result(3);
        if(func_name=="quartic_potential_2")
        {
            result << -2.0,-2.0,0.0;
        }
        return result;
    } 

static VectorXd domain_max(const string& func_name)
    {
        VectorXd result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 2.0,2.0,4.0;
        }
        return result;
    } 

static VectorXd function_range(const string& func_name)
    {
        VectorXd result(2);
        if(func_name=="quartic_potential_2")
        {
            result <<  -2.0, 2.0;
        }
        return result;
    }

static VectorXi block_num(const string& func_name) //number of blocks that splits the domain in each dimension
    {
        VectorXi result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 10,10,10;
        }

        return result;
    }


static VectorXi point_num_in_block(const string& func_name) //number of initial points in each block in each dimension
    {
        VectorXi result(3);
        if(func_name=="quartic_potential_2")
        {
            result << 5,5,5;
        }
        return result;
    }

    template <typename T>
    void derivative(const VectorX<T>& p_) {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // Convert input (float or double) -> float32 tensor
        Eigen::VectorXf p = p_.template cast<float>();

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, 3},
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
            {dxy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
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