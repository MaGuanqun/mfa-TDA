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
    VectorX<T> domain_range;
    VectorX<T> function_range;
    VectorXi block_num;
    VectorXi point_num_in_block;
    T delta_h = 1e-3;

    INRModel(const string& func_name, const std::string& model_path = "inr_base.pt"): function_name(func_name),
    device(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU),
    // device(torch::kCPU),
    loaded(false) {
        try {
            module = torch::jit::load(model_path, device);
            // module.to(device);
            module.eval();
            loaded = true;
            std::cout << "Loaded INR base model on "
                      << (device.is_cuda() ? "CUDA" : "CPU") << std::endl;
        } catch (const c10::Error& e) {
            std::cerr << "Error loading INR model: " << e.what() << std::endl;
            loaded = false;
        }
        this->domain_min = domain_min_(func_name);
        this->domain_max = domain_max_(func_name);
        this->domain_range = this->domain_max - this->domain_min;
        this->function_range = function_range_(func_name);
        this->block_num = block_num_(func_name);
        this->point_num_in_block = point_num_in_block_(func_name);

        std::cout<<"INR model domain min: "<<this->domain_min.transpose()<<std::endl;
        std::cout<<"INR model domain max: "<<this->domain_max.transpose()<<std::endl;
        std::cout<<"INR model domain_range range: "<<this->domain_range.transpose()<<std::endl;

        init_buffers();
        torch::jit::setGraphExecutorOptimize(false);
        torch::jit::getProfilingMode() = false;

        hx = delta_h * this->domain_range(0);
        hy = delta_h * this->domain_range(1);
        ht = delta_h * this->domain_range(2);
        inv2hx = static_cast<T>(0.5) / hx;
        inv2hy = static_cast<T>(0.5) / hy;
        inv2ht = static_cast<T>(0.5) / ht;

        init_stencil_offsets();
    }

    bool isLoaded() const { return loaded; }

static VectorX<T> domain_min_(const string& func_name)
    {
        VectorX<T> result(3);
        if(func_name=="quartic_potential_2")
        {
            result << -2.0,-2.0,0.0;
        }
        else if (func_name=="vortex_street")
        {
            result << 0,0,0;
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
        else if (func_name=="vortex_street")
        {
            result << 99,79,49;
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
        else if (func_name=="vortex_street")
        {
            result << 10,8,5;
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
        else if (func_name=="vortex_street")
        {
            result << 5,5,5;
        }
        return result;
    }

    // data is in [-1,1]
    VectorX<T> convert_point_to_domain_reverse_order(const VectorX<T>& p_) 
    {
        VectorX<T> param = (p_-domain_min).cwiseQuotient(domain_range) * 2.0 - VectorX<T>::Constant(p_.size(),1.0);
        param.reverseInPlace();
        return param;
    }

    void convert_gradient_to_domain(VectorX<T>& grad)
    {
        grad = grad.cwiseQuotient(domain_range)*2.0;
    }

    void convert_hessian_to_domain(Eigen::MatrixX<T>& Hessian)
    {
        for(int i=0;i<Hessian.rows();++i)
        {
            for(int j=i;j<Hessian.cols();++j)
            {
                Hessian(i,j)=4.0*Hessian(i,j)/(domain_range(i)*domain_range(j));
                Hessian(j,i)=Hessian(i,j);
            }
        }
    }


    // VectorX<T>& third_spatial_out,            // size 4: [fxxx, fxxy, fxyy, fyyy]
    // VectorX<T>& third_tmix_out)               // size 3: [fxxt, fxyt, fyyt]
    void convert_third_order_to_domain(VectorX<T>& third_spatial, VectorX<T>& third_tmix)
    {
        // Spatial third
        third_spatial(0)=third_spatial(0)/(domain_range(0)*domain_range(0)*domain_range(0)); //fxxx
        third_spatial(1)=third_spatial(1)/(domain_range(0)*domain_range(0)*domain_range(1)); //fxxy
        third_spatial(2)=third_spatial(2)/(domain_range(0)*domain_range(1)*domain_range(1)); //fxyy
        third_spatial(3)=third_spatial(3)/(domain_range(1)*domain_range(1)*domain_range(1)); //fyyy

        third_spatial = third_spatial * 8.0;

        // Time-mixed third
        third_tmix(0)=third_tmix(0)/(domain_range(0)*domain_range(0)*domain_range(2)); //fxxt
        third_tmix(1)=third_tmix(1)/(domain_range(0)*domain_range(1)*domain_range(2)); //fxyt
        third_tmix(2)=third_tmix(2)/(domain_range(1)*domain_range(1)*domain_range(2)); //fyyt

        third_tmix = third_tmix * 8.0;
    }

    // hessian and \partial_xt, \partial_yt
    void query_hessian_t(const VectorX<T>& point,
                VectorX<T>& grad_t,                     // size 3: [fx, fy]
                Eigen::MatrixX<T>& Hessian)
    {
            VectorX<T> p = convert_point_to_domain_reverse_order(point);
            input_.detach_();                 // drop previous graph

            if (input_.is_cpu()) {
                // CPU path: write directly via pointer (no accessor overhead)
                float* buf = input_.data_ptr<float>();   // contiguous [1,3]
                buf[0] = static_cast<float>(p(0));
                buf[1] = static_cast<float>(p(1));
                buf[2] = static_cast<float>(p(2));
            } else {
                // CUDA path: never use accessor or data_ptr to write from host.
                float* h = input_host_.data_ptr<float>();
                h[0] = static_cast<float>(p(0));
                h[1] = static_cast<float>(p(1));
                h[2] = static_cast<float>(p(2));
                input_.copy_(input_host_, /*non_blocking=*/false);
            }


            input_.requires_grad_(true);      // we need autograd
            torch::Tensor f = module.forward({input_}).toTensor().reshape({}); // scalar
            torch::Tensor g_full = torch::autograd::grad(
                /*outputs=*/{f},
                /*inputs=*/{input_},
                /*grad_outputs=*/{},
                /*retain_graph=*/true,
                /*create_graph=*/true
            )[0];  // shape [1,3]


        // Extract grad components
        auto fy = g_full.index({0, 1});
        auto fx = g_full.index({0, 2});

        // // Grad (fx, fy)
        // grad_out.resize(2);
        // grad_out(0) = static_cast<T>(fx.detach().to(torch::kCPU).item<float>());
        // grad_out(1) = static_cast<T>(fy.detach().to(torch::kCPU).item<float>());
        // grad_out = grad_out.cwiseQuotient(domain_range.head(2))*2.0;

        // ---- Second order: rows of Hessian are grads of fx and fy wrt FULL input ----
        // We exclude any t-columns for the Hessian output (xy-block only).
        torch::Tensor Hx_full = torch::autograd::grad(
            /*outputs=*/{fx},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/false
        )[0];  // [1,3] = [fxx, fxy, fxt]

        torch::Tensor Hy_full = torch::autograd::grad(
            /*outputs=*/{fy},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/false
        )[0];  // [1,3] = [fyx, fyy, fyt]

        Hessian.resize(2,2);

        Hessian(0,0)= static_cast<T>(Hx_full.index({0,2}).detach().to(torch::kCPU).item<float>()); //fxx
        Hessian(0,1)= static_cast<T>(Hx_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fxy
        Hessian(1,0)= Hessian(0,1); //fyx
        Hessian(1,1)= static_cast<T>(Hy_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fyy
        convert_hessian_to_domain(Hessian);


        grad_t.resize(2);
        grad_t(0) = static_cast<T>(Hx_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fxt
        grad_t(0) = 4.0 * grad_t(0)/(domain_range(0)*domain_range(2));
        grad_t(1) = static_cast<T>(Hy_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fyt
        grad_t(1) = 4.0 * grad_t(1)/(domain_range(1)*domain_range(2));

    }

    void query_dim_reduced_grad_hessian(const VectorX<T>& point,
                VectorX<T>& grad_out,                     // size 3: [fx, fy]
                Eigen::MatrixX<T>& Hessian, int removed_dim)
    {
        VectorX<T> p = convert_point_to_domain_reverse_order(point);
        input_.detach_();                 // drop previous graph

        if (input_.is_cpu()) {
            // CPU path: write directly via pointer (no accessor overhead)
            float* buf = input_.data_ptr<float>();   // contiguous [1,3]
            buf[0] = static_cast<float>(p(0));
            buf[1] = static_cast<float>(p(1));
            buf[2] = static_cast<float>(p(2));
        } else {
            // CUDA path: never use accessor or data_ptr to write from host.
            float* h = input_host_.data_ptr<float>();
            h[0] = static_cast<float>(p(0));
            h[1] = static_cast<float>(p(1));
            h[2] = static_cast<float>(p(2));
            input_.copy_(input_host_, /*non_blocking=*/false);
        }


        input_.requires_grad_(true);      // we need autograd
        torch::Tensor f = module.forward({input_}).toTensor().reshape({}); // scalar
        torch::Tensor g_full = torch::autograd::grad(
            /*outputs=*/{f},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/true
        )[0];  // shape [1,3]

        // Extract grad components
        // auto ft = g_full.index({0, 0});
        auto fy = g_full.index({0, 1});
        auto fx = g_full.index({0, 2});

        // Grad (fx, fy)
        grad_out.resize(2);
        grad_out(0) = static_cast<T>(fx.detach().to(torch::kCPU).item<float>());
        grad_out(1) = static_cast<T>(fy.detach().to(torch::kCPU).item<float>());


        grad_out = grad_out.cwiseQuotient(domain_range.head(2))*2.0;

        // ---- Second order: rows of Hessian are grads of fx and fy wrt FULL input ----
        // We exclude any t-columns for the Hessian output (xy-block only).
        torch::Tensor Hx_full = torch::autograd::grad(
            /*outputs=*/{fx},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/false
        )[0];  // [1,3] = [fxx, fxy, fxt]

        torch::Tensor Hy_full = torch::autograd::grad(
            /*outputs=*/{fy},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/false
        )[0];  // [1,3] = [fyx, fyy, fyt]

        Hessian.resize(2,2);
        if(removed_dim==2) //remove t
        {

            Hessian(0,0)= static_cast<T>(Hx_full.index({0,2}).detach().to(torch::kCPU).item<float>()); //fxx
            Hessian(0,1)= static_cast<T>(Hx_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fxy
            Hessian(1,0)= Hessian(0,1); //fyx
            Hessian(1,1)= static_cast<T>(Hy_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fyy
            convert_hessian_to_domain(Hessian);

        }
        else if(removed_dim==1) //remove y
        {
            Hessian(0,0)= static_cast<T>(Hx_full.index({0,2}).detach().to(torch::kCPU).item<float>()); //fxx
            Hessian(0,0)= 4.0*Hessian(0,0)/(domain_range(0)*domain_range(0));

            Hessian(0,1)= static_cast<T>(Hx_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fxt
            Hessian(0,1)= 4.0*Hessian(0,1)/(domain_range(0)*domain_range(2));

            Hessian(1,0)= static_cast<T>(Hy_full.index({0,2}).detach().to(torch::kCPU).item<float>()); //fyx
            Hessian(1,0)= 4.0*Hessian(1,0)/(domain_range(1)*domain_range(0));

            Hessian(1,1)= static_cast<T>(Hy_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fyt
            Hessian(1,1)= 4.0*Hessian(1,1)/(domain_range(1)*domain_range(2));
        }
        else if(removed_dim==0) //remove x
        {
            Hessian(0,0)= static_cast<T>(Hx_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fxy
            Hessian(0,0)= 4.0*Hessian(0,0)/(domain_range(0)*domain_range(1));

            Hessian(0,1)= static_cast<T>(Hx_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fxt
            Hessian(0,1)= 4.0*Hessian(0,1)/(domain_range(0)*domain_range(2));

            Hessian(1,0)= static_cast<T>(Hy_full.index({0,1}).detach().to(torch::kCPU).item<float>()); //fyy
            Hessian(1,0)= 4.0*Hessian(1,0)/(domain_range(1)*domain_range(1));

            Hessian(1,1)= static_cast<T>(Hy_full.index({0,0}).detach().to(torch::kCPU).item<float>()); //fyt
            Hessian(1,1)= 4.0*Hessian(1,1)/(domain_range(1)*domain_range(2));
        }

    }

   // Batch first-order gradients in *domain* coordinates (x,y,t)
    void eval_grad_batch_in_domain(
        const std::vector<VectorX<T>>& points_domain,   // N points in (x,y,t)
        std::vector<VectorX<T>>&       grads_domain_out // N gradients [fx, fy, ft]
    ) {
        if (points_domain.empty()) {
            grads_domain_out.clear();
            return;
        }

        const int N = static_cast<int>(points_domain.size());
        constexpr int D = 3;

        grads_domain_out.clear();
        grads_domain_out.resize(N);

        torch::Tensor input_batch;

        // ---- Build [N,3] input on CPU, then move to GPU if needed ----
        if (device.is_cuda()) {
            // 1) build on CPU
            auto host_opts = torch::TensorOptions()
                                .dtype(torch::kFloat32)
                                .device(torch::kCPU);
            torch::Tensor input_host = torch::empty({N, D}, host_opts);

            {
                auto acc = input_host.accessor<float, 2>();
                for (int i = 0; i < N; ++i) {
                    VectorX<T> p_model =
                        convert_point_to_domain_reverse_order(points_domain[i]);
                    acc[i][0] = static_cast<float>(p_model(0));
                    acc[i][1] = static_cast<float>(p_model(1));
                    acc[i][2] = static_cast<float>(p_model(2));
                }
            }

            // 2) move to GPU
            input_batch = input_host.to(device, /*non_blocking=*/false, /*copy=*/true);
        } else {
            // Pure CPU path: we can write directly into the tensor
            auto opts = torch::TensorOptions()
                            .dtype(torch::kFloat32)
                            .device(device);
            input_batch = torch::empty({N, D}, opts);

            {
                auto acc = input_batch.accessor<float, 2>();
                for (int i = 0; i < N; ++i) {
                    VectorX<T> p_model =
                        convert_point_to_domain_reverse_order(points_domain[i]);
                    acc[i][0] = static_cast<float>(p_model(0));
                    acc[i][1] = static_cast<float>(p_model(1));
                    acc[i][2] = static_cast<float>(p_model(2));
                }
            }
        }

        input_batch.set_requires_grad(true);

        // Forward: f for all points; shape [N] or [N,1]
        torch::Tensor out = module.forward({input_batch}).toTensor();
        out = out.view({N});  // ensure shape [N]

        // grad(sum_i f_i) wrt each row is just grad(f_i) for that row.
        std::vector<torch::Tensor> grads = torch::autograd::grad(
            /*outputs=*/{out.sum()},
            /*inputs=*/{input_batch},
            /*grad_outputs=*/{},
            /*retain_graph=*/false,
            /*create_graph=*/false
        );

        // Bring gradients back to CPU for Eigen
        torch::Tensor g_all = grads[0].to(torch::kCPU).contiguous();
        auto g_acc = g_all.accessor<float, 2>();

        for (int i = 0; i < N; ++i) {
            VectorX<T> g_model(D);

            // model order is [t, y, x] -> domain [x, y, t]
            g_model(0) = static_cast<T>(g_acc[i][2]); // fx
            g_model(1) = static_cast<T>(g_acc[i][1]); // fy
            g_model(2) = static_cast<T>(g_acc[i][0]); // ft

            // rescale from [-1,1] model coords to physical domain
            convert_gradient_to_domain(g_model);      // divides by domain_range * 2

            grads_domain_out[i] = g_model;
        }

        input_batch.set_requires_grad(false);
    }
    // Second-order: Hessian in *domain* coordinates via batched FD on gradient
    void compute_hessian_numeric_in_domain(
        const VectorX<T>& point_domain,
        Eigen::MatrixX<T>& H_out          // 3x3, H(i,j) = ∂²f/∂u_i∂u_j, u=(x,y,t)
    ) {
        constexpr int D = 3;
        H_out.resize(D, D);
        H_out.setZero();

        // We need 2D points: for each dimension j, point +/- h_j e_j
        std::vector<VectorX<T>> points_domain;
        points_domain.reserve(2 * D);

        std::array<T, D> h_step;

        for (int j = 0; j < D; ++j) {
            T range_j = domain_range(j);
            T h = delta_h * range_j;
            h_step[j] = h;

            VectorX<T> p_plus  = point_domain;
            VectorX<T> p_minus = point_domain;
            p_plus(j)  += h;
            p_minus(j) -= h;
            points_domain.push_back(p_plus);   // index 2*j
            points_domain.push_back(p_minus);  // index 2*j+1
        }

        // Compute all gradients at once
        std::vector<VectorX<T>> grads_domain;
        eval_grad_batch_in_domain(points_domain, grads_domain); // size = 2D

        // Fill Hessian columns
        for (int j = 0; j < D; ++j) {
            const VectorX<T>& grad_plus  = grads_domain[2 * j + 0];
            const VectorX<T>& grad_minus = grads_domain[2 * j + 1];

            H_out.col(j) = (grad_plus - grad_minus) / (static_cast<T>(2) * h_step[j]);
        }

        // Symmetrize for numerical robustness
        H_out = static_cast<T>(0.5) * (H_out + H_out.transpose());

        // Optional: if you want f_tt = 0 to match previous behavior:
        // H_out(2, 2) = static_cast<T>(0);
    }

    void query_up_to_third_derivative(
        const VectorX<T>& point_domain,   // (x,y,t) in physical domain
        VectorX<T>& grad_out,             // [fx, fy, ft]
        Eigen::MatrixX<T>& Hessian_out,   // 3x3
        VectorX<T>& third_spatial_out,    // [fxxx, fxxy, fxyy, fyyy]
        VectorX<T>& third_tmix_out       // [fxxt, fxyt, fyyt]
    ) {
        if (!loaded)
            throw std::runtime_error("INR model not loaded");

        constexpr int D = 3; // (x,y,t)
      // ---------------------------------------------------------------------
        // 3) Build physical points for all precomputed integer offsets
        // ---------------------------------------------------------------------
        const int N = static_cast<int>(stencil_offsets_.size());
        std::vector<VectorX<T>> points_domain;
        points_domain.reserve(N);

        for (const auto& off : stencil_offsets_) {
            VectorX<T> p = point_domain;
            p(0) += static_cast<T>(off[0]) * hx; // x + ix * hx
            p(1) += static_cast<T>(off[1]) * hy; // y + iy * hy
            p(2) += static_cast<T>(off[2]) * ht; // t + it * ht
            points_domain.push_back(p);
        }

        // ---------------------------------------------------------------------
        // 4) Single batched autograd: f and ∇f for ALL stencil points
        // ---------------------------------------------------------------------
        std::vector<VectorX<T>> grads; // grads[i] = [fx, fy, ft] in domain coords
        eval_grad_batch_in_domain(points_domain, grads);

        // Helper to fetch grad at a given integer offset
        auto grad_at = [&](int ix, int iy, int it) -> const VectorX<T>& {
            Offset key{ix, iy, it};
            auto it_f = stencil_offset_to_index_.find(key);
            if (it_f == stencil_offset_to_index_.end()) {
                throw std::runtime_error("Missing gradient for requested offset");
            }
            return grads[it_f->second];
        };


        auto indi = [&](int ix, int iy, int it) -> const int& {
            Offset key{ix, iy, it};
            auto it_f = stencil_offset_to_index_.find(key);
            if (it_f == stencil_offset_to_index_.end()) {
                throw std::runtime_error("Missing gradient for requested offset");
            }
            return it_f->second;
        };
        // Base index (0,0,0)
        Offset base_off{0,0,0};
        int base_idx = stencil_offset_to_index_[base_off];

        // Base gradient
        grad_out = grads[base_idx];  // [fx, fy, ft]
        // ---------------------------------------------------------------------
        // 5) Hessian components Hxx, Hxy, Hyy at the 7 centers:
        //
        // centers (integer offsets):
        //   c0: (0,0,0)
        //   c1: (1,0,0)
        //   c2: (-1,0,0)
        //   c3: (0,1,0)
        //   c4: (0,-1,0)
        //   c5: (0,0,1)
        //   c6: (0,0,-1)
        //
        // Using:
        //   Hxx(c) = [ fx(c+e_x) - fx(c-e_x) ] / (2hx)
        //   Hyy(c) = [ fy(c+e_y) - fy(c-e_y) ] / (2hy)
        //   Hxy(c) = 0.5 * ( [fx(c+e_y) - fx(c-e_y)]/(2hy)
        //                  + [fy(c+e_x) - fy(c-e_x)]/(2hx) )
        // ---------------------------------------------------------------------
        struct Center { int ix, iy, it; };
        const int NUM_CENTERS = 7;
        Center centers[NUM_CENTERS] = {
            { 0,  0,  0},  // c0
            { 1,  0,  0},  // c1
            {-1,  0,  0},  // c2
            { 0,  1,  0},  // c3
            { 0, -1,  0},  // c4
            { 0,  0,  1},  // c5
            { 0,  0, -1}   // c6
        };

        T Hxx[NUM_CENTERS];
        T Hxy[NUM_CENTERS];
        T Hyy[NUM_CENTERS];

        for (int ci = 0; ci < NUM_CENTERS; ++ci) {
            int cx = centers[ci].ix;
            int cy = centers[ci].iy;
            int ct = centers[ci].it;

            // neighbors along x at this center
            const VectorX<T>& g_px = grad_at(cx + 1, cy, ct);
            const VectorX<T>& g_mx = grad_at(cx - 1, cy, ct);

            // neighbors along y at this center
            const VectorX<T>& g_py = grad_at(cx, cy + 1, ct);
            const VectorX<T>& g_my = grad_at(cx, cy - 1, ct);

            // Hxx = ∂/∂x (fx)
            T Hxx_from_x = (g_px(0) - g_mx(0)) * inv2hx;

            // Hxy via ∂/∂y (fx)
            T Hxy_from_y = (g_py(0) - g_my(0)) * inv2hy;

            // Optionally symmetrize using ∂/∂x (fy)
            T Hyx_from_x = (g_px(1) - g_mx(1)) * inv2hx;
            T Hxy_sym    = static_cast<T>(0.5) * (Hxy_from_y + Hyx_from_x);

            // Hyy = ∂/∂y (fy)
            T Hyy_from_y = (g_py(1) - g_my(1)) * inv2hy;

            Hxx[ci] = Hxx_from_x;
            Hxy[ci] = Hxy_sym;
            Hyy[ci] = Hyy_from_y;

        }

        // Time-mixed Hessian entries at base (0,0,0):
        // Hxt = ∂/∂t (fx) ≈ [ fx(0,0,1) - fx(0,0,-1) ] / (2ht)
        // Hyt = ∂/∂t (fy) ≈ [ fy(0,0,1) - fy(0,0,-1) ] / (2ht)
        const VectorX<T>& g_t_plus  = grad_at(0, 0, +1);
        const VectorX<T>& g_t_minus = grad_at(0, 0, -1);

        T Hxt = (g_t_plus(0) - g_t_minus(0)) * inv2ht; // f_xt
        T Hyt = (g_t_plus(1) - g_t_minus(1)) * inv2ht; // f_yt

        // ---------------------------------------------------------------------
        // 6) Assemble Hessian at base point (center c0 = (0,0,0)):
        //
        //     H =
        //     [ f_xx  f_xy  f_xt ]
        //     [ f_xy  f_yy  f_yt ]
        //     [ f_xt  f_yt   0   ]   // no f_tt by design
        // ---------------------------------------------------------------------
        Hessian_out.setZero(D, D);

        Hessian_out(0,0) = Hxx[0];
        Hessian_out(0,1) = Hxy[0];
        Hessian_out(1,0) = Hxy[0];
        Hessian_out(1,1) = Hyy[0];

        Hessian_out(0,2) = Hxt;
        Hessian_out(2,0) = Hxt;
        Hessian_out(1,2) = Hyt;
        Hessian_out(2,1) = Hyt;
        Hessian_out(2,2) = static_cast<T>(0); // explicitly no f_tt

        // ---------------------------------------------------------------------
        // 7) Third derivatives from Hxx/Hxy/Hyy at centers:
        //
        // Spatial:
        //   fxxx = d/dx (fxx) = (Hxx(c1) - Hxx(c2)) / (2hx)
        //   fxxy = d/dy (fxx) = (Hxx(c3) - Hxx(c4)) / (2hy)
        //   fxyy = d/dy (fxy) = (Hxy(c3) - Hxy(c4)) / (2hy)
        //   fyyy = d/dy (fyy) = (Hyy(c3) - Hyy(c4)) / (2hy)
        //
        // Time-mixed:
        //   fxxt = d/dt (fxx) = (Hxx(c5) - Hxx(c6)) / (2ht)
        //   fxyt = d/dt (fxy) = (Hxy(c5) - Hxy(c6)) / (2ht)
        //   fyyt = d/dt (fyy) = (Hyy(c5) - Hyy(c6)) / (2ht)
        // ---------------------------------------------------------------------
        third_spatial_out.resize(4);
        third_tmix_out.resize(3);

        constexpr int c0 = 0;
        constexpr int c1 = 1;
        constexpr int c2 = 2;
        constexpr int c3 = 3;
        constexpr int c4 = 4;
        constexpr int c5 = 5;
        constexpr int c6 = 6;

        // Spatial 3rd derivatives
        third_spatial_out(0) = (Hxx[c1] - Hxx[c2]) * inv2hx; // fxxx
        third_spatial_out(1) = (Hxx[c3] - Hxx[c4]) * inv2hy; // fxxy
        third_spatial_out(2) = (Hxy[c3] - Hxy[c4]) * inv2hy; // fxyy
        third_spatial_out(3) = (Hyy[c3] - Hyy[c4]) * inv2hy; // fyyy

        // Time-mixed 3rd derivatives
        third_tmix_out(0) = (Hxx[c5] - Hxx[c6]) * inv2ht; // fxxt
        third_tmix_out(1) = (Hxy[c5] - Hxy[c6]) * inv2ht; // fxyt
        third_tmix_out(2) = (Hyy[c5] - Hyy[c6]) * inv2ht; // fyyt
    }
    


//     void query(const VectorX<T>& point,
//                         VectorX<T>&       out)
//     {
//         if (!loaded)
//             throw std::runtime_error("INR model not loaded");

//         out.resize(1);
//         input_.detach_(); // drop previous graph

//         // Convert to model input domain and reverse order
//         VectorX<T> p = convert_point_to_domain_reverse_order(point);

//         std::cout<<"query point in model domain: "<<p.transpose()<<std::endl;

//         if (input_.is_cpu()) {
//             // CPU path: write directly via pointer (no accessor overhead)
//             float* buf = input_.data_ptr<float>();   // contiguous [1,3]
//             buf[0] = static_cast<float>(p(0));
//             buf[1] = static_cast<float>(p(1));
//             buf[2] = static_cast<float>(p(2));
//         } else {
//             // CUDA path: never use accessor or data_ptr to write from host.
//             float* h = input_host_.data_ptr<float>();
//             h[0] = static_cast<float>(p(0));
//             h[1] = static_cast<float>(p(1));
//             h[2] = static_cast<float>(p(2));
//             input_.copy_(input_host_, /*non_blocking=*/false);
//         }

//         input_.requires_grad_(false);      // we need autograd
//                 // auto in_acc = input_.accessor<T,2>();
//                 // in_acc[0][0] = p[0];
//                 // in_acc[0][1] = p[1];
//                 // in_acc[0][2] = p[2];

//         torch::Tensor f = module.forward({input_}).toTensor().reshape({}); // scalar
//    // Helper to map (nx, ny, nt) → coordinate index in reversed order

//         out(0) = static_cast<T>(f.detach().to(torch::kCPU).item<float>());
       

//     }
    void query(const VectorX<T>& point, VectorX<T>& out, const VectorXi& derivs = VectorXi()) {
        out.resize(1);
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // Convert input (float or float) -> float32 tensor
        input_.detach_();                 // drop previous graph


        VectorX<T> p = convert_point_to_domain_reverse_order(point);
        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, p.size()},
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
            int coord = (derivs(0)==1) ? 2 : (derivs(1)==1 ? 1 : 0);
            int coord_in_domain = 2 - coord; // because of reverse order
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
                val = val / domain_range(coord_in_domain)*2.0;
            } else {
                // input truly not used in graph => derivative is 0
                val = 0.0f;
            }

            out(0) = static_cast<T>(val);
            return;
        }


            // Helper to map (nx,ny,nt) -> which coord (0,1,2) is still “active” at each step
        auto pick_coord = [](int nx, int ny, int nt) {
            if (nx > 0) return 2;
            if (ny > 0) return 1;
            return 0;
        };


        if(derivs.sum()==2)
        {
            int nx = derivs(0), ny = derivs(1), nt = derivs(2);
            int c1 = pick_coord(nx, ny, nt);

            // full gradient of f wrt input
            auto g_full = torch::autograd::grad({output}, {input}, {}, /*retain_graph=*/true, /*create_graph=*/true)[0]; // [1,3]
            auto g1 = g_full.index({0, c1});  // scalar: ∂f/∂x or ∂f/∂y or ∂f/∂t
            // consume that count
            if      (c1==2) --nx;
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
            if      (c1==2) --nx;
            else if (c1==1) --ny;
            else            --nt;

            // 2nd derivative
            int c2 = pick_coord(nx, ny, nt);
            auto g2_full = torch::autograd::grad({d1}, {input}, {}, /*retain_graph=*/true, /*create_graph=*/true)[0]; // [1,3]
            auto d2 = g2_full.index({0, c2}); // scalar
            if      (c2==2) --nx;
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


// Computes value, grad, Hessian (xy-only), and selected 3rd-order partials.
    // Excludes any Hessian entries with t and any 3rd-order terms with t^2 or t^3.
    void query_up_to_third_derivative2(const VectorX<T>& point,
                VectorX<T>& grad_out,                     // size 3: [fx, fy]
                Eigen::MatrixX<T>& Hessian,            // [ [fxx, fxy, fxt], [fyx, fyy, fyt], [fxt, fyt]], exclude ftt
                VectorX<T>& third_spatial_out,            // size 4: [fxxx, fxxy, fxyy, fyyy]
                VectorX<T>& third_tmix_out)               // size 3: [fxxt, fxyt, fyyt]
    {

        
        VectorX<T> p = convert_point_to_domain_reverse_order(point);

        if (!loaded) throw std::runtime_error("INR model not loaded");
        // TORCH_CHECK(p_.size() == 3, "Input must be (x,y,t)");

        // 0) Build input (CPU -> clone -> to(device))
        // Eigen::VectorXf p = p_.template cast<float>();

        input_.detach_();                 // drop previous graph

        // torch::NoGradGuard nog;


        if (input_.is_cpu()) {
            // CPU path: write directly via pointer (no accessor overhead)
            float* buf = input_.data_ptr<float>();   // contiguous [1,3]
            buf[0] = static_cast<float>(p(0));
            buf[1] = static_cast<float>(p(1));
            buf[2] = static_cast<float>(p(2));
        } else {
            // CUDA path: never use accessor or data_ptr to write from host.
            float* h = input_host_.data_ptr<float>();
            h[0] = static_cast<float>(p(0));
            h[1] = static_cast<float>(p(1));
            h[2] = static_cast<float>(p(2));
            input_.copy_(input_host_, /*non_blocking=*/false);
        }


        input_.requires_grad_(true);      // we need autograd
        // auto in_acc = input_.accessor<T,2>();
        // in_acc[0][0] = p[0];
        // in_acc[0][1] = p[1];
        // in_acc[0][2] = p[2];

        torch::Tensor f = module.forward({input_}).toTensor().reshape({}); // scalar


        // torch::Tensor input = torch::from_blob(
        //     (void*)p.data(), {1, 3},
        //     torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU)
        // ).clone();//.to(device, /*non_blocking=*/false, /*copy=*/true);

        // // 1) Forward (scalar)
        // input.set_requires_grad(true);
        // torch::Tensor f = module.forward({input}).toTensor().reshape({});  // scalar

        // ---- First order gradient wrt FULL input (needed once) ----
        // Keep graph because we'll need higher-order derivatives.
        torch::Tensor g_full = torch::autograd::grad(
            /*outputs=*/{f},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/true
        )[0];  // shape [1,3]

        // Extract grad components
        auto ft = g_full.index({0, 0});
        auto fy = g_full.index({0, 1});
        auto fx = g_full.index({0, 2});

        // ---- Second order: rows of Hessian are grads of fx and fy wrt FULL input ----
        // We exclude any t-columns for the Hessian output (xy-block only).
        torch::Tensor Hx_full = torch::autograd::grad(
            /*outputs=*/{fx},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/true
        )[0];  // [1,3] = [fxx, fxy, fxt]

        torch::Tensor Hy_full = torch::autograd::grad(
            /*outputs=*/{fy},
            /*inputs=*/{input_},
            /*grad_outputs=*/{},
            /*retain_graph=*/true,
            /*create_graph=*/true
        )[0];  // [1,3] = [fyx, fyy, fyt]

        // ---- Third order (no t^2 or t^3): take grads of second-order spatial terms ----
        // Spatial third: from fxx and fyy and fxy
        torch::Tensor fxx = Hx_full.index({0, 2}); // fxx
        torch::Tensor fxy = Hx_full.index({0, 1}); // fxy
        // torch::Tensor fyx = Hy_full.index({0, 2}); // fyx (== fxy ideally)
        torch::Tensor fyy = Hy_full.index({0, 1}); // fyy

        // Grad of fxx wrt FULL input -> [fxxx, fxxy, fxxt]
        torch::Tensor G_fxx = torch::autograd::grad(
            {fxx}, {input_}, {}, /*retain_graph=*/true, /*create_graph=*/false
        )[0]; // [1,3]

        // Grad of fxy wrt FULL input -> [fxyx, fxyy, fxyt]
        torch::Tensor G_fxy = torch::autograd::grad(
            {fxy}, {input_}, {}, /*retain_graph=*/true, /*create_graph=*/false
        )[0]; // [1,3]

        // Grad of fyy wrt FULL input -> [fyyx, fyyy, fyyt]
        torch::Tensor G_fyy = torch::autograd::grad(
            {fyy}, {input_}, {}, /*retain_graph=*/true, /*create_graph=*/false
        )[0]; // [1,3]

        // Grad (fx, fy)
        grad_out.resize(3);
        grad_out(0) = static_cast<T>(fx.detach().to(torch::kCPU).item<float>());
        grad_out(1) = static_cast<T>(fy.detach().to(torch::kCPU).item<float>());
        grad_out(2) = static_cast<T>(ft.detach().to(torch::kCPU).item<float>());

        convert_gradient_to_domain(grad_out);


        // 
        Hessian.resize(3,3);
        Hessian(0,0) = static_cast<T>(Hx_full.index({0,2}).detach().to(torch::kCPU).item<float>()); // fxx
        Hessian(0,1) = static_cast<T>(Hx_full.index({0,1}).detach().to(torch::kCPU).item<float>()); // fxy
        Hessian(1,0) = Hessian(0,1); 
        Hessian(1,1) = static_cast<T>(Hy_full.index({0,1}).detach().to(torch::kCPU).item<float>()); // fyy
        Hessian(0,2) = static_cast<T>(Hx_full.index({0,0}).detach().to(torch::kCPU).item<float>()); // fxt
        Hessian(2,0) = Hessian(0,2);
        Hessian(1,2) = static_cast<T>(Hy_full.index({0,0}).detach().to(torch::kCPU).item<float>()); // fyt
        Hessian(2,1) = Hessian(1,2);
        Hessian(2,2) = 0.0; // ftt excluded

        convert_hessian_to_domain(Hessian);

        // Third (spatial only): [fxxx, fxxy, fxyy, fyyy]
        third_spatial_out.resize(4);
        third_spatial_out(0) = static_cast<T>(G_fxx.index({0,2}).detach().to(torch::kCPU).item<float>()); // fxxx
        third_spatial_out(1) = static_cast<T>(G_fxx.index({0,1}).detach().to(torch::kCPU).item<float>()); // fxxy
        third_spatial_out(2) = static_cast<T>(G_fxy.index({0,1}).detach().to(torch::kCPU).item<float>()); // fxyy
        third_spatial_out(3) = static_cast<T>(G_fyy.index({0,1}).detach().to(torch::kCPU).item<float>()); // fyyy

        // Third (with exactly one t): [fxxt, fxyt, fyyt]
        third_tmix_out.resize(3);
        third_tmix_out(0) = static_cast<T>(G_fxx.index({0,0}).detach().to(torch::kCPU).item<float>()); // fxxt
        third_tmix_out(1) = static_cast<T>(G_fxy.index({0,0}).detach().to(torch::kCPU).item<float>()); // fxyt
        third_tmix_out(2) = static_cast<T>(G_fyy.index({0,0}).detach().to(torch::kCPU).item<float>()); // fyyt

        convert_third_order_to_domain(third_spatial_out, third_tmix_out);
    }


    void derivative(const VectorX<T>& point) {
        if (!loaded) throw std::runtime_error("INR model not loaded");


        VectorX<T> p = convert_point_to_domain_reverse_order(point);
        // Convert input (float or float) -> float32 tensor

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, p.size()},
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

        torch::Tensor dx = g1.index({0, 2}); // ∂f/∂x
        torch::Tensor dy = g1.index({0, 1}); // ∂f/∂y
        torch::Tensor dz = g1.index({0, 0}); // ∂f/∂z

        // Compute second-order derivatives
        std::vector<torch::Tensor> second_order_x = torch::autograd::grad(
            {dx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor g2=second_order_x.at(0);
        torch::Tensor dxx = g2.index({0, 2}); // ∂²f/∂x²
        torch::Tensor dxy = g2.index({0, 1}); // ∂²f/∂x∂y
        torch::Tensor dxt = g2.index({0, 0}); // ∂²f/∂x∂t

        std::vector<torch::Tensor> second_order_y = torch::autograd::grad(
            {dy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor gy2=second_order_y.at(0);
        torch::Tensor dyy = gy2.index({0, 1}); // ∂²f/∂y²
        torch::Tensor dyt = gy2.index({0, 0}); // ∂²f/∂y∂t

        // Compute third-order derivatives
        std::vector<torch::Tensor> third_order_xx = torch::autograd::grad(
            {dxx}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor g3=third_order_xx.at(0);

        torch::Tensor dxxx = g3.index({0, 2}); // ∂³f/∂x³
        torch::Tensor dxxy = g3.index({0, 1}); // ∂³f/∂x²∂y
        torch::Tensor dxxt = g3.index({0, 0}); // ∂³f/∂x²∂t

        std::vector<torch::Tensor> third_order_yy = torch::autograd::grad(
            {dyy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/true);
        torch::Tensor gyy3=third_order_yy.at(0);
        torch::Tensor dyyx = gyy3.index({0, 2}); // ∂³f/∂y²∂x
        torch::Tensor dyyy = gyy3.index({0, 1}); // ∂³f/∂y³
        torch::Tensor dyyt = gyy3.index({0, 0}); // ∂³f/∂y²∂t

        std::vector<torch::Tensor> third_order_xy = torch::autograd::grad(
            {dxy}, {input}, /*grad_outputs=*/{}, /*retain_graph=*/true, /*create_graph=*/false);
        torch::Tensor gxy3=third_order_xy.at(0);
        torch::Tensor dxyt = gxy3.index({0, 0}); // ∂³f/∂x∂y∂t

        // Return the required derivatives as a custom structure or vector
        VectorX<T> grad(p.size());

        grad(0) = static_cast<T>(dx.item<float>());
        grad(1) = static_cast<T>(dy.item<float>());
        grad(2) = static_cast<T>(dz.item<float>());
    }

private:
    torch::jit::script::Module module;
    torch::Device device{torch::kCPU};
    bool loaded = false;

    torch::Tensor input_;        // [1,3], cached
    torch::Tensor input_host_;   // [1,3], pinned CPU (only used when device is CUDA)
    torch::Tensor grad_buf_;     // [1,3], cached

    T hx, hy, ht;               // grid spacing in domain coords
    T inv2hx, inv2hy, inv2ht;     // 1/(2*hx), etc.

    using Offset = std::array<int, 3>;  // (ix, iy, it)

    // Precomputed integer offsets for stencil
    std::vector<Offset> stencil_offsets_;               // all offsets we use
    std::map<Offset, int> stencil_offset_to_index_;     // offset -> local index

    // Integer offsets of the 7 "centers"
    std::array<Offset, 7> center_offsets_;              // c0..c6

    // Call this once after loading module
    void init_buffers() {
        auto opts = torch::TensorOptions().dtype(torch::kFloat32).device(device);
        input_    = torch::zeros({1,3}, opts);
        grad_buf_ = torch::empty({1,3}, opts);

        if (device.is_cuda()) {
            input_host_ = torch::empty({1,3},
                torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCPU).pinned_memory(true));
        } else {
            input_host_ = torch::Tensor();
        }

        at::set_num_threads(std::max(1u, std::thread::hardware_concurrency()));
        at::set_num_interop_threads(2);
    }

    void init_stencil_offsets() {
        if (!stencil_offsets_.empty())
            return; // already initialized

        using std::abs;

        // Build offsets:
        //   it = 0: |ix| + |iy| <= 2
        //   it = +1: |ix| + |iy| <= 1
        //   it = -1: |ix| + |iy| <= 1
        auto add_layer = [&](int it, int max_r) {
            for (int ix = -max_r; ix <= max_r; ++ix) {
                for (int iy = -max_r; iy <= max_r; ++iy) {
                    if (abs(ix) + abs(iy) <= max_r) {
                        stencil_offsets_.push_back(Offset{ix, iy, it});
                    }
                }
            }
        };

        stencil_offsets_.clear();
        stencil_offset_to_index_.clear();

        // t = 0 layer: |ix| + |iy| <= 2
        add_layer(/*it=*/0, /*max_r=*/2);
        // t = +1 layer: |ix| + |iy| <= 1
        add_layer(/*it=*/+1, /*max_r=*/1);
        // t = -1 layer: |ix| + |iy| <= 1
        add_layer(/*it=*/-1, /*max_r=*/1);

        // Dedup just in case (should be unique already with this logic)
        stencil_offset_to_index_.clear();
        std::vector<Offset> unique_offsets;
        unique_offsets.reserve(stencil_offsets_.size());

        for (const auto& off : stencil_offsets_) {
            if (stencil_offset_to_index_.find(off) == stencil_offset_to_index_.end()) {
                int idx = static_cast<int>(unique_offsets.size());
                unique_offsets.push_back(off);
                stencil_offset_to_index_[off] = idx;
            }
        }
        stencil_offsets_.swap(unique_offsets);

        // Set up center offsets c0..c6 in integer space
        // c0: (0,0,0)
        // c1: (1,0,0)
        // c2: (-1,0,0)
        // c3: (0,1,0)
        // c4: (0,-1,0)
        // c5: (0,0,1)
        // c6: (0,0,-1)
        center_offsets_[0] = Offset{ 0,  0,  0};
        center_offsets_[1] = Offset{ 1,  0,  0};
        center_offsets_[2] = Offset{-1,  0,  0};
        center_offsets_[3] = Offset{ 0,  1,  0};
        center_offsets_[4] = Offset{ 0, -1,  0};
        center_offsets_[5] = Offset{ 0,  0,  1};
        center_offsets_[6] = Offset{ 0,  0, -1};
    }

};