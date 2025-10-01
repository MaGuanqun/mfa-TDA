#pragma once
#include <torch/script.h>
#include <Eigen/Dense>
#include <iostream>
#include <tuple>
#include <type_traits>

class INRModel {
public:
    INRModel(const std::string& model_path = "inr_base.pt") {
        try {
            device = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
            module = torch::jit::load(model_path, device);
            module.eval();
            loaded = true;
            std::cout << "Loaded INR base model (float32) on "
                      << (device.is_cuda() ? "CUDA" : "CPU") << std::endl;
        }
        catch (const c10::Error& e) {
            std::cerr << "Error loading INR model: " << e.what() << std::endl;
            loaded = false;
        }
    }

    bool isLoaded() const { return loaded; }

    template <typename T>
    void query_derivative(const VectorX<T>& p_) {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // Convert input (float or double) -> float32 tensor
        Eigen::VectorXf p = p_.template cast<float>();

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, 3},
            torch::TensorOptions().dtype(torch::kFloat32).device(device)
        ).clone();
        input.set_requires_grad(true);

        // Forward pass
        torch::Tensor output = module.forward({input}).toTensor(); // [1,1]
        torch::Tensor cur = output[0][0];

        // Compute first-order derivatives
        std::vector<torch::Tensor> first_order = torch::autograd::grad(
            {cur}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);

        torch::Tensor dx = first_order_vec[0]; // ∂f/∂x
        torch::Tensor dy = first_order_vec[1]; // ∂f/∂y
        torch::Tensor dz = first_order_vec[2]; // ∂f/∂z

        // Compute second-order derivatives
        std::vector<torch::Tensor> second_order_x = torch::autograd::grad(
            {dx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor dxx = second_order_x[0][0][0]; // ∂²f/∂x²
        torch::Tensor dxy = second_order_x[0][0][1]; // ∂²f/∂x∂y
        torch::Tensor dxt = second_order_x[0][0][2]; // ∂²f/∂x∂t

        std::vector<torch::Tensor> second_order_y = torch::autograd::grad(
            {dy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor dyy = second_order_y[0][0][1]; // ∂²f/∂y²
        torch::Tensor dyt = second_order_y[0][0][2]; // ∂²f/∂y∂t

        // Compute third-order derivatives
        std::vector<torch::Tensor> third_order_xx = torch::autograd::grad(
            {dxx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor dxxx = third_order_xx[0][0][0]; // ∂³f/∂x³
        torch::Tensor dxxy = third_order_xx[0][0][1]; // ∂³f/∂x²∂y
        torch::Tensor dxxt = third_order_xx[0][0][2]; // ∂³f/∂x²∂t

        std::vector<torch::Tensor> third_order_yy = torch::autograd::grad(
            {dyy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor dyyx = third_order_yy[0][0][0]; // ∂³f/∂y²∂x
        torch::Tensor dyyy = third_order_yy[0][0][1]; // ∂³f/∂y³
        torch::Tensor dyyt = third_order_yy[0][0][2]; // ∂³f/∂y²∂t

        std::vector<torch::Tensor> third_order_xy = torch::autograd::grad(
            {dxy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor dxyt = third_order_xy[0][0][2]; // ∂³f/∂x∂y∂t

        // Return the required derivatives as a custom structure or vector
        VectorX<T> grad(p_.size());

        grad(0) = dx.item<float>();
        grad(1) = dy.item<float>();
        grad(2) = dz.item<float>();

        std::cout<<"grad "<<grad.transpose()<<std::endl;
    }

private:
    torch::jit::script::Module module;
    torch::Device device;
    bool loaded = false;
};