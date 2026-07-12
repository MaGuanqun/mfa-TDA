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

    INRModel(const string& func_name, const std::string& model_path = "inr_base.pt",int initial_point_num_in_a_block_=-1, torch::Device input_device = torch::kCPU): function_name(func_name),
    device(input_device), //torch::cuda::is_available() ? torch::kCUDA : torch::kCPU
    // device(torch::kCPU),
    loaded(false) 
    
    {
        try {
            module = torch::jit::load(model_path, device);
            // module.to(device);
            module.eval();
            loaded = true;
            // std::cout << "Loaded INR base model on "
                    //   << (device.is_cuda() ? "CUDA" : "CPU") << std::endl;
        } catch (const c10::Error& e) {
            std::cerr << "Error loading INR model: " << e.what() << std::endl;
            loaded = false;
        }
        this->domain_min = domain_min_(func_name);
        this->domain_max = domain_max_(func_name);
        this->domain_range = this->domain_max - this->domain_min;
        this->function_range = function_range_(func_name);
        this->block_num = block_num_(func_name);
        if(initial_point_num_in_a_block_==-1)
        {
            this->point_num_in_block = point_num_in_block_(func_name);
        }
        else
        {
            this->point_num_in_block = VectorXi::Constant(3, initial_point_num_in_a_block_);
            // std::cout<<"initial point num set in a block "<<this->point_num_in_block.transpose()<<std::endl;
        }

        // std::cout<<"INR model domain min: "<<this->domain_min.transpose()<<std::endl;
        // std::cout<<"INR model domain max: "<<this->domain_max.transpose()<<std::endl;
        //  std::cout<<"INR model domain_range range: "<<this->domain_range.transpose()<<std::endl;

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
        else if (func_name=="vortex_street_3d")
        {
            result << -0.5,-0.5,13.5;
        }
        else if (func_name=="boussinesq_3d")
        {
            result << -0.5,-0.5,0.0;
        }
        else if (func_name=="fluid")
        {
            result << 0,0,0;
        }
        else if (func_name=="cylinder")
        {
            result << 1.5,0.5,0.0;
        }
        else if (func_name=="cylinder2")
        {
            result << 3.2,0.5,0.0;
        }
        else if (func_name=="cylinder3")
        {
            result << 3.5,0.5,0.0;
        }
        else
        {
            result.resize(4);
            if (func_name=="quartic_potential_3d")
            {
                result << -2.0,-2.0,-2.0,0.0;
            }
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
        else if (func_name=="vortex_street_3d")
        {
            result << 7.5,0.5,15.0;
        }
        else if (func_name=="boussinesq_3d")
        {
            result << 0.5,2.5,1.5;
        }
        else if (func_name=="fluid")
        {
            result << 1.0,1.0,1.0;
        }
        else if(func_name=="cylinder")
        {
            result << 5.5,1.5,1.0;
        }
        else if(func_name=="cylinder2")
        {
            result << 5.5,1.5,1.0;
        }
        else if(func_name=="cylinder3")
        {
            result << 5.5,1.5,1.0;
        }
        else
        {
            result.resize(4);
            if (func_name=="quartic_potential_3d")
            {
                result << 2.0,2.0,2.0,4.0;
            }
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
        else if (func_name=="vortex_street_3d")
        {
            result <<  0, 1.8358269;
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
        else if (func_name=="vortex_street_3d")
        {
            result << 80,10,15;
        }
        else if (func_name=="boussinesq_3d")
        {
            result << 10,30,15;
        }
        else if (func_name=="fluid")
        {
            result << 10,10,10;
        }
        else if (func_name=="cylinder")
        {
            result << 40,10,10;
        }   
        else if (func_name=="cylinder2")
        {
            result << 23,10,10;
        }   
        else if (func_name=="cylinder3")
        {
            result << 20,10,10;
        }
        else
        {
            result.resize(4);
            if (func_name=="quartic_potential_3d")
            {
                result << 10,10,10,10;
            }
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
        else if (func_name=="vortex_street_3d")
        {
            result << 2,2,2;
        }
        else if (func_name=="boussinesq_3d")
        {
            result << 2,2,2;
        }
        else if (func_name=="fluid")
        {
            result << 4,4,4;
        }
        else if (func_name=="cylinder")
        {
            result << 4,4,4;
        }
        else if (func_name=="cylinder2")
        {
            result << 4,4,4;
        }
        else if (func_name=="cylinder3")
        {
            result << 4,4,4;
        }
        else
        {
            result.resize(4);
            if (func_name=="quartic_potential_3d")
            {
                result << 4,4,4,4;
            }
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
                VectorX<T>& grad_t,                     //[fxt, fyt]
                Eigen::MatrixX<T>& Hessian)
    {
            // ------------------------------------------------------------------
        // 1) Build stencil points in domain coordinates
        //     center: (x, y, t)
        //     x± : (x ± hx, y, t)
        //     y± : (x, y ± hy, t)
        //     t± : (x, y, t ± ht)
        // ------------------------------------------------------------------
        std::vector<VectorX<T>> points_domain;
        points_domain.reserve(5);

        // center
        points_domain.push_back(point);          // idx 0

        // x+ , x-
        VectorX<T> p_xp = point;
        VectorX<T> p_xm = point;
        p_xp(0) += hx;
        p_xm(0) -= hx;
        points_domain.push_back(p_xp);                  // idx 1
        points_domain.push_back(p_xm);                  // idx 2

        // y+ , y-
        VectorX<T> p_yp = point;
        VectorX<T> p_ym = point;
        p_yp(1) += hy;
        p_ym(1) -= hy;
        points_domain.push_back(p_yp);                  // idx 3
        points_domain.push_back(p_ym);                  // idx 4


        // ------------------------------------------------------------------
        // 2) Evaluate gradients at all stencil points (domain coords)
        //     grads[i] = [fx, fy, ft] in physical domain
        // ------------------------------------------------------------------
        std::vector<VectorX<T>> grads;
        eval_grad_batch_in_domain(points_domain, grads);

        const VectorX<T>& g_c  = grads[0]; // center
        const VectorX<T>& g_xp = grads[1];
        const VectorX<T>& g_xm = grads[2];
        const VectorX<T>& g_yp = grads[3];
        const VectorX<T>& g_ym = grads[4];

        // ------------------------------------------------------------------
        // 3) Numeric Hessian in xy-plane from gradient field [fx, fy]
        //
        // Hxx = d/dx (fx) ≈ [ fx(x+hx) - fx(x-hx) ] / (2hx)
        // Hyy = d/dy (fy) ≈ [ fy(y+hy) - fy(y-hy) ] / (2hy)
        //
        // Hxy = Hyx ≈ 0.5 * ( d/dy fx + d/dx fy )
        //             = 0.5 * ( [ fx(y+hy) - fx(y-hy) ] / (2hy)
        //                      + [ fy(x+hx) - fy(x-hx) ] / (2hx) )
        // ------------------------------------------------------------------
        Hessian.resize(2,2);
        // Hxx
        T Hxx = (g_xp(0) - g_xm(0)) * inv2hx;  // fx difference along x

        // Hyy
        T Hyy = (g_yp(1) - g_ym(1)) * inv2hy;  // fy difference along y

        // Hxy from d/dy fx
        T Hxy_from_y = (g_yp(0) - g_ym(0)) * inv2hy;

        // Hyx from d/dx fy
        T Hyx_from_x = (g_xp(1) - g_xm(1)) * inv2hx;

        T Hxy = static_cast<T>(0.5) * (Hxy_from_y + Hyx_from_x);

        Hessian(0,0) = Hxx;
        Hessian(0,1) = Hxy;
        Hessian(1,0) = Hxy;
        Hessian(1,1) = Hyy;

        // ------------------------------------------------------------------
        // 4) Time derivatives of spatial gradients: f_xt, f_yt
        //
       // f_xt = d/dx (f_t) ≈ [ f_t(x+hx) - f_t(x-hx) ] / (2hx)
        // f_yt = d/dy (f_t) ≈ [ f_t(y+hy) - f_t(y-hy) ] / (2hy)
        // ------------------------------------------------------------------
        grad_t.resize(2);
        grad_t(0) = (g_xp(2) - g_xm(2)) * inv2hx;  // f_xt
        grad_t(1) = (g_yp(2) - g_ym(2)) * inv2hy;  // f_yt
    }


    void query_grad_exclude_last(const VectorX<T>& point,
                VectorX<T>& grad_out)
    {
        std::vector<VectorX<T>> pts;
        pts.push_back(point);
        std::vector<VectorX<T>> grads;
        eval_grad_batch_in_domain(pts, grads);
        grad_out = grads[0].head(point.size()-1);
    }

    void query_dim_reduced_grad_hessian(const VectorX<T>& point,
                VectorX<T>& grad_out,                     // size 3: [fx, fy]
                Eigen::MatrixX<T>& Hessian, int removed_dim)
    {
        std::array<int, 2> kept_dims;
        int k = 0;
        for (int d = 0; d < point.size(); ++d)
            if (d != removed_dim) kept_dims[k++] = d;
        
        T hs[3]={hx,hy,ht};
        T inv2h[3]={inv2hx,inv2hy,inv2ht};
        
        std::vector<VectorX<T>> pts;
        pts.reserve(5);
        pts.push_back(point); // center

        VectorX<T> p_plus = point, p_minus = point;
        // along first kept dim
        p_plus(kept_dims[0]) += hs[kept_dims[0]];
        p_minus(kept_dims[0]) -= hs[kept_dims[0]];
        pts.push_back(p_plus);
        pts.push_back(p_minus);

        // along second kept dim
        p_plus = point;
        p_minus = point;
        p_plus(kept_dims[1]) += hs[kept_dims[1]];
        p_minus(kept_dims[1]) -= hs[kept_dims[1]];
        pts.push_back(p_plus);
        pts.push_back(p_minus);


        // -----------------------------------------------
        // Evaluate gradients at all stencil points
        // -----------------------------------------------
        std::vector<VectorX<T>> grads;
        eval_grad_batch_in_domain(pts, grads);
        // grads[i] = [fx, fy, ft] in domain coords

        const VectorX<T>& g_center = grads[0];
        grad_out.resize(2);
        grad_out(0) = g_center(0);
        grad_out(1) = g_center(1);


  
         // -----------------------------------------------
        // 4) Numeric Hessian in the reduced plane
        // -----------------------------------------------
        Hessian.resize(2, 2);
        Hessian.setZero();

        // std::cout<<"kept_dims: "<<kept_dims[0]<<" , "<<kept_dims[1]<<std::endl;
        // std::cout<<pts[0].transpose()<<std::endl;
        // std::cout<<pts[1].transpose()<<std::endl;
        // std::cout<<pts[2].transpose()<<std::endl;
        // std::cout<<pts[3].transpose()<<std::endl;
        // std::cout<<pts[4].transpose()<<std::endl;

        // central finite differences on [fx, fy]
        // Keep only x,y partials (kept_dims order may vary)
        const VectorX<T>& g_p1 = grads[1]; // +d1
        const VectorX<T>& g_m1 = grads[2]; // -d1
        const VectorX<T>& g_p2 = grads[3]; // +d2
        const VectorX<T>& g_m2 = grads[4]; // -d2


        // std::cout<<"g_center: "<<g_center.transpose()<<std::endl;
        // std::cout<<"g_p1: "<<g_p1.transpose()<<std::endl;
        // std::cout<<"g_m1: "<<g_m1.transpose()<<std::endl;
        // std::cout<<"g_p2: "<<g_p2.transpose()<<std::endl;
        // std::cout<<"g_m2: "<<g_m2.transpose()<<std::endl;


        //fx(+kept_dims[0])-fx(-kept_dims[0])
        // d(fx,fy)/d first kept dim
        VectorX<T> dfd1 = (g_p1.head(2) - g_m1.head(2)) * inv2h[kept_dims[0]];
        // d(fx,fy)/d second kept dim
        VectorX<T> dfd2 = (g_p2.head(2) - g_m2.head(2)) * inv2h[kept_dims[1]];

        // Diagonal terms: ∂fx/∂x, ∂fy/∂y
        Hessian(0, 0) = dfd1(0);
        Hessian(1, 1) = dfd2(1);

        // Off-diagonal ∂fx/∂y ≈ ∂fy/∂x = 0.5*(dfd2(0)+dfd1(1))
        Hessian(0, 1) = dfd2(0);
        Hessian(1, 0) = dfd1(1);
        

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
                                .dtype(torch::kFloat64)
                                .device(torch::kCPU);
            torch::Tensor input_host = torch::empty({N, D}, host_opts);

            {
                auto acc = input_host.accessor<double, 2>();
                for (int i = 0; i < N; ++i) {
                    VectorX<T> p_model =
                        convert_point_to_domain_reverse_order(points_domain[i]);
                    acc[i][0] = static_cast<double>(p_model(0));
                    acc[i][1] = static_cast<double>(p_model(1));
                    acc[i][2] = static_cast<double>(p_model(2));
                }
            }

            // 2) move to GPU
            input_batch = input_host.to(device, /*non_blocking=*/false, /*copy=*/true);
        } else {
            // Pure CPU path: we can write directly into the tensor
            auto opts = torch::TensorOptions()
                            .dtype(torch::kFloat64)
                            .device(device);
            input_batch = torch::empty({N, D}, opts);

            {
                auto acc = input_batch.accessor<double, 2>();
                for (int i = 0; i < N; ++i) {
                    VectorX<T> p_model =
                        convert_point_to_domain_reverse_order(points_domain[i]);
                    acc[i][0] = static_cast<double>(p_model(0));
                    acc[i][1] = static_cast<double>(p_model(1));
                    acc[i][2] = static_cast<double>(p_model(2));
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
        auto g_acc = g_all.accessor<double, 2>();

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
   
    void query_up_to_third_derivative(
        const VectorX<T>& point_domain,   // (x,y,t) in physical domain
        VectorX<T>& grad_out,             // [fx, fy, ft]
        Eigen::MatrixX<T>& Hessian_out,   // 3x3
        VectorX<T>& third_spatial_out,    // [fxxx, fxxy, fxyy, fyyy]
        VectorX<T>& third_tmix_out,       // [fxxt, fxyt, fyyt]
        T* fxtt_out = nullptr,            // optional: f_xtt (one-spatial + two-time)
        T* fytt_out = nullptr             // optional: f_ytt
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

        // Optional one-spatial + two-time third derivatives, from a central
        // second time difference of the spatial gradient (reuses the t-stencil
        // gradients already evaluated above; no extra model evaluation):
        //   f_xtt = ( f_x(t+ht) - 2 f_x(t) + f_x(t-ht) ) / ht^2   (and f_ytt)
        if (fxtt_out || fytt_out)
        {
            const T inv_ht2 = static_cast<T>(1) / (ht * ht);
            if (fxtt_out)
                *fxtt_out = (g_t_plus(0) - static_cast<T>(2) * grad_out(0) + g_t_minus(0)) * inv_ht2;
            if (fytt_out)
                *fytt_out = (g_t_plus(1) - static_cast<T>(2) * grad_out(1) + g_t_minus(1)) * inv_ht2;
        }

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
//             double* buf = input_.data_ptr<double>();   // contiguous [1,3]
//             buf[0] = static_cast<double>(p(0));
//             buf[1] = static_cast<double>(p(1));
//             buf[2] = static_cast<double>(p(2));
//         } else {
//             // CUDA path: never use accessor or data_ptr to write from host.
//             double* h = input_host_.data_ptr<double>();
//             h[0] = static_cast<double>(p(0));
//             h[1] = static_cast<double>(p(1));
//             h[2] = static_cast<double>(p(2));
//             input_.copy_(input_host_, /*non_blocking=*/false);
//         }

//         input_.requires_grad_(false);      // we need autograd
//                 // auto in_acc = input_.accessor<T,2>();
//                 // in_acc[0][0] = p[0];
//                 // in_acc[0][1] = p[1];
//                 // in_acc[0][2] = p[2];

//         torch::Tensor f = module.forward({input_}).toTensor().reshape({}); // scalar
//    // Helper to map (nx, ny, nt) → coordinate index in reversed order

//         out(0) = static_cast<T>(f.detach().to(torch::kCPU).item<double>());
       

//     }
    void query(const VectorX<T>& point, VectorX<T>& out, const VectorXi& derivs = VectorXi()) {
        out.resize(1);
        if (!loaded) throw std::runtime_error("INR model not loaded");

        int order = 0;
        if (derivs.size() > 0)
                order = derivs.sum();

        if (order > 3) {
            throw std::runtime_error("query(): derivative order > 3 is not supported");
        }

    // For convenience, extract (nx, ny, nt); assume missing entries are 0.
        int nx = (derivs.size() > 0) ? derivs(0) : 0;
        int ny = (derivs.size() > 1) ? derivs(1) : 0;
        int nt = (derivs.size() > 2) ? derivs(2) : 0;

        if (order == 0) {
            input_.detach_();

            VectorX<T> p = convert_point_to_domain_reverse_order(point);
            torch::Tensor input = torch::from_blob(
                (void*)p.data(),
                {1, p.size()},
                torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
            ).clone();
            input = input.to(device, /*non_blocking=*/false, /*copy=*/true);

            torch::Tensor output = module.forward({input}).toTensor().reshape({});
            out(0) = static_cast<T>(output.detach().to(torch::kCPU).item<double>());
            return;
        }


        if (order == 1) {
            input_.detach_();

            VectorX<T> p = convert_point_to_domain_reverse_order(point);
            torch::Tensor input = torch::from_blob(
                (void*)p.data(),
                {1, p.size()},
                torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
            ).clone();
            input = input.to(device, /*non_blocking=*/false, /*copy=*/true);

            input.set_requires_grad(true);
            torch::Tensor output = module.forward({input}).toTensor().reshape({});

            // Which coordinate? (x,y,t) in domain corresponds to (2,1,0) in model
            int coord_model;
            if (nx == 1 && ny == 0 && nt == 0) {
                coord_model = 2; // x
            } else if (nx == 0 && ny == 1 && nt == 0) {
                coord_model = 1; // y
            } else if (nx == 0 && ny == 0 && nt == 1) {
                coord_model = 0; // t
            } else {
                throw std::runtime_error("query(): unsupported first-order multiindex");
            }

            int coord_domain = 2 - coord_model; // reverse order mapping

            auto g_list = torch::autograd::grad(
                /*outputs=*/{output},
                /*inputs=*/{input},
                /*grad_outputs=*/{},
                /*retain_graph=*/false,
                /*create_graph=*/false,
                /*allow_unused=*/true
            );

            torch::Tensor g_full = g_list[0];
            double val = 0.0;
            if (g_full.defined()) {
                val = g_full.index({0, coord_model}).detach().to(torch::kCPU).item<double>();
                // scale from [-1,1] to domain
                val = val / domain_range(coord_domain) * 2.0f;
            } else {
                val = 0.0f;
            }

            out(0) = static_cast<T>(val);
            return;
        }

        // =========================================================
    // 2nd order: numeric from gradients on *minimal* stencil
    //            using eval_grad_batch_in_domain()
    // =========================================================
    if (order == 2) {
        // We support: f_xx, f_yy, f_xy, f_xt, f_yt
        T val = static_cast<T>(0);
        bool ok = false;

        // ---- f_xx: ∂/∂x (fx) ----
        if (nx == 2 && ny == 0 && nt == 0) {
            std::vector<VectorX<T>> pts(2);
            pts[0]=point;
            pts[1]=point;
            pts[0](0) += hx;   // x + hx
            pts[1](0) -= hx;   // x - hx

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads); // [fx, fy, ft]

            const VectorX<T>& g_xp = grads[0];
            const VectorX<T>& g_xm = grads[1];

            val = (g_xp(0) - g_xm(0)) * inv2hx; // ∂/∂x (fx)
            ok  = true;
        }
        // ---- f_yy: ∂/∂y (fy) ----
        else if (nx == 0 && ny == 2 && nt == 0) {
            std::vector<VectorX<T>> pts(2);
            pts[0]=point;
            pts[1]=point;
            pts[0](1) += hy;   // y + hy
            pts[1](1) -= hy;   // y - hy
            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            const VectorX<T>& g_yp = grads[0];
            const VectorX<T>& g_ym = grads[1];

            val = (g_yp(1) - g_ym(1)) * inv2hy; // ∂/∂y (fy)
            ok  = true;
        }
        // ---- f_xy: mixed, symmetric central difference ----
        else if (nx == 1 && ny == 1 && nt == 0) {
            // Use 4 points:
            //   (x, y±hy) for d/dy fx
            //   (x±hx, y) for d/dx fy
            std::vector<VectorX<T>> pts;
            pts.reserve(4);

            VectorX<T> p_yp = point; p_yp(1) += hy;
            VectorX<T> p_ym = point; p_ym(1) -= hy;
            VectorX<T> p_xp = point; p_xp(0) += hx;
            VectorX<T> p_xm = point; p_xm(0) -= hx;

            pts.push_back(p_yp); // 0
            pts.push_back(p_ym); // 1
            pts.push_back(p_xp); // 2
            pts.push_back(p_xm); // 3

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            const VectorX<T>& g_yp = grads[0];
            const VectorX<T>& g_ym = grads[1];
            const VectorX<T>& g_xp = grads[2];
            const VectorX<T>& g_xm = grads[3];

            T Hxy_from_y = (g_yp(0) - g_ym(0)) * inv2hy; // d/dy (fx)
            T Hyx_from_x = (g_xp(1) - g_xm(1)) * inv2hx; // d/dx (fy)

            val = static_cast<T>(0.5) * (Hxy_from_y + Hyx_from_x);
            ok  = true;
        }
        // ---- f_xt: ∂/∂x (f_t) ----
        else if (nx == 1 && ny == 0 && nt == 1) {
            std::vector<VectorX<T>> pts(2);
            pts[0]=point;
            pts[1]=point;
            pts[0](0) += hx;
            pts[1](0) -= hx;

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            const VectorX<T>& g_xp = grads[0];
            const VectorX<T>& g_xm = grads[1];

            val = (g_xp(2) - g_xm(2)) * inv2hx; // d/dx (f_t)
            ok  = true;
        }
        // ---- f_yt: ∂/∂y (f_t) ----
        else if (nx == 0 && ny == 1 && nt == 1) {
            std::vector<VectorX<T>> pts(2);
            pts[0]=point;
            pts[1]=point;
            pts[0](1) += hy;
            pts[1](1) -= hy;

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            const VectorX<T>& g_yp = grads[0];
            const VectorX<T>& g_ym = grads[1];

            val = (g_yp(2) - g_ym(2)) * inv2hy; // d/dy (f_t)
            ok  = true;
        }

        if (!ok)
            throw std::runtime_error("query(): unsupported second-order multiindex for numeric approximation");

        out(0) = val;
        return;
    }

    // =========================================================
    // 3rd order: numeric via eval_grad_batch_in_domain()
    //            (per-derivative, but still using a stencil of grads)
    // =========================================================
    if (order == 3) {
        // We support: f_xxx, f_xxy, f_xyy, f_yyy, f_xxt, f_xyt, f_yyt
        T val = static_cast<T>(0);
        bool ok = false;

        using Offset = std::array<int,3>;
        auto make_points_from_offsets =
            [&](const std::vector<Offset>& offs,
                std::vector<VectorX<T>>&   pts_out) {
                pts_out.clear();
                pts_out.reserve(offs.size());
                for (const auto& o : offs) {
                    VectorX<T> p = point;
                    p(0) += static_cast<T>(o[0]) * hx;
                    p(1) += static_cast<T>(o[1]) * hy;
                    p(2) += static_cast<T>(o[2]) * ht;
                    pts_out.push_back(p);
                }
            };

        // ---------------- f_xxx ----------------
        if (nx == 3 && ny == 0 && nt == 0) {
            // Use 1D 3rd-derivative stencil on fx along x:
            // f_xxx ≈ ( -0.5 fx(x-2h) + fx(x-h) - fx(x+h) + 0.5 fx(x+2h) ) / h^3
            std::vector<Offset> offs = {
                Offset{-2, 0, 0},
                Offset{-1, 0, 0},
                Offset{ 1, 0, 0},
                Offset{ 2, 0, 0}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fx_m2 = grads[0](0);
            T fx_m1 = grads[1](0);
            T fx_p1 = grads[2](0);
            T fx_p2 = grads[3](0);

            T h3 = hx * hx * hx;
            val = (static_cast<T>(-0.5) * fx_m2
                 + fx_m1
                 - fx_p1
                 + static_cast<T>(0.5) * fx_p2) / h3;
            ok = true;
        }
        // ---------------- f_xxy ----------------
        else if (nx == 2 && ny == 1 && nt == 0) {
            // f_xxy = ∂/∂y (f_xx)
            // f_xx(y) ≈ [ fx(x+hx,y) - fx(x-hx,y) ] / (2hx)
            // so:
            // f_xxy ≈ ( f_xx(y+hy) - f_xx(y-hy) ) / (2hy)
            //
            // Need fx at:
            //  (x±hx, y+hy), (x±hx, y-hy)
            std::vector<Offset> offs = {
                Offset{1, 1, 0},
                Offset{-1, 1, 0},
                Offset{1, -1, 0},
                Offset{-1, -1, 0}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            // Hxx(y+hy)
            T fx_xp_yp = grads[0](0); // (1,1)
            T fx_xm_yp = grads[1](0); // (-1,1)
            T Hxx_yp   = (fx_xp_yp - fx_xm_yp) * inv2hx;

            // Hxx(y-hy)
            T fx_xp_ym = grads[2](0); // (1,-1)
            T fx_xm_ym = grads[3](0); // (-1,-1)
            T Hxx_ym   = (fx_xp_ym - fx_xm_ym) * inv2hx;

            val = (Hxx_yp - Hxx_ym) * inv2hy; // d/dy (f_xx)
            ok  = true;
        }
        // ---------------- f_xyy ----------------
        else if (nx == 1 && ny == 2 && nt == 0) {
            // f_xyy = ∂/∂y (f_xy)
            // Use f_xy ≈ d/dy (fx) at y±hy and then central diff in y.
            //
            // f_xy(y) ≈ [ fx(x, y+hy) - fx(x, y-hy) ] / (2hy)
            //
            // Need fx at:
            //   (x, y+2hy), (x, y), (x, y-2hy)
            std::vector<Offset> offs = {
                Offset{0, 2, 0},
                Offset{0,  0, 0},
                Offset{0, -2, 0}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fx_p2 = grads[0](0); // (0,2)
            T fx_0  = grads[1](0); // (0,0)
            T fx_m2 = grads[2](0); // (0,-2)

            // central second derivative in y of fx:
            // f_xyy ≈ (fx(y+2h) - 2 fx(y) + fx(y-2h)) / (4 h^2)
            // (this is a 2h step variant; error O(h^2))
            T four_h2 = static_cast<T>(4) * hy * hy;
            val = (fx_p2 - static_cast<T>(2) * fx_0 + fx_m2) / four_h2;
            ok  = true;
        }
        // ---------------- f_yyy ----------------
        else if (nx == 0 && ny == 3 && nt == 0) {
            // 1D 3rd-derivative stencil on fy along y:
            // f_yyy ≈ ( -0.5 fy(y-2h) + fy(y-h) - fy(y+h) + 0.5 fy(y+2h) ) / h^3
            std::vector<Offset> offs = {
                Offset{0, -2, 0},
                Offset{0, -1, 0},
                Offset{0, 1, 0},
                Offset{0, 2, 0}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fy_m2 = grads[0](1);
            T fy_m1 = grads[1](1);
            T fy_p1 = grads[2](1);
            T fy_p2 = grads[3](1);

            T h3 = hy * hy * hy;
            val = (static_cast<T>(-0.5) * fy_m2
                 + fy_m1
                 - fy_p1
                 + static_cast<T>(0.5) * fy_p2) / h3;
            ok = true;
        }
        // ---------------- f_xxt ----------------
        else if (nx == 2 && ny == 0 && nt == 1) {
            // f_xxt = ∂/∂t (f_xx)
            // f_xx(t) ≈ [ fx(x+hx) - fx(x-hx) ] / (2hx) at t±ht
            //
            // Need fx at:
            //  (x±hx, t+ht), (x±hx, t-ht)
            std::vector<Offset> offs = {
                Offset{1, 0, 1},
                Offset{-1, 0, 1},
                Offset{1, 0, -1},
                Offset{-1, 0, -1}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fx_xp_tp = grads[0](0);
            T fx_xm_tp = grads[1](0);
            T fx_xp_tm = grads[2](0);
            T fx_xm_tm = grads[3](0);

            T Hxx_tp = (fx_xp_tp - fx_xm_tp) * inv2hx;
            T Hxx_tm = (fx_xp_tm - fx_xm_tm) * inv2hx;

            val = (Hxx_tp - Hxx_tm) * inv2ht; // d/dt (f_xx)
            ok  = true;
        }
        // ---------------- f_xyt ----------------
        else if (nx == 1 && ny == 1 && nt == 1) {
            // f_xyt = ∂/∂t (f_xy)
            // Use f_xy ≈ d/dy (fx) at t±ht and central diff in t.
            //
            // Need fx at:
            //   (x, y±hy, t+ht), (x, y±hy, t-ht)
            std::vector<Offset> offs = {
                Offset{0, 1, 1},
                Offset{0, -1, 1},
                Offset{0, 1, -1},
                Offset{0, -1, -1}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fx_yp_tp = grads[0](0);
            T fx_ym_tp = grads[1](0);
            T fx_yp_tm = grads[2](0);
            T fx_ym_tm = grads[3](0);

            T fxy_tp = (fx_yp_tp - fx_ym_tp) * inv2hy;
            T fxy_tm = (fx_yp_tm - fx_ym_tm) * inv2hy;

            val = (fxy_tp - fxy_tm) * inv2ht; // d/dt (f_xy)
            ok  = true;
        }
        // ---------------- f_yyt ----------------
        else if (nx == 0 && ny == 2 && nt == 1) {
            // f_yyt = ∂/∂t (f_yy)
            // f_yy(t) ≈ [ fy(y+hy) - fy(y-hy) ] / (2hy) at t±ht
            //
            // Need fy at:
            //   (x, y±hy, t+ht), (x, y±hy, t-ht)
            std::vector<Offset> offs = {
                Offset{0, 1, 1},
                Offset{0, -1, 1},
                Offset{0, 1, -1},
                Offset{0, -1, -1}
            };
            std::vector<VectorX<T>> pts;
            make_points_from_offsets(offs, pts);

            std::vector<VectorX<T>> grads;
            eval_grad_batch_in_domain(pts, grads);

            T fy_yp_tp = grads[0](1);
            T fy_ym_tp = grads[1](1);
            T fy_yp_tm = grads[2](1);
            T fy_ym_tm = grads[3](1);

            T Hyy_tp = (fy_yp_tp - fy_ym_tp) * inv2hy;
            T Hyy_tm = (fy_yp_tm - fy_ym_tm) * inv2hy;

            val = (Hyy_tp - Hyy_tm) * inv2ht; // d/dt (f_yy)
            ok  = true;
        }

        if (!ok)
            throw std::runtime_error("query(): unsupported third-order multiindex for numeric approximation");

        out(0) = val;
        return;
    }

    // Should never reach here
    throw std::runtime_error("query(): unreachable state");
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
        // Eigen::VectorXf p = p_.template cast<double>();

        input_.detach_();                 // drop previous graph

        // torch::NoGradGuard nog;


        if (input_.is_cpu()) {
            // CPU path: write directly via pointer (no accessor overhead)
            double* buf = input_.data_ptr<double>();   // contiguous [1,3]
            buf[0] = static_cast<double>(p(0));
            buf[1] = static_cast<double>(p(1));
            buf[2] = static_cast<double>(p(2));
        } else {
            // CUDA path: never use accessor or data_ptr to write from host.
            double* h = input_host_.data_ptr<double>();
            h[0] = static_cast<double>(p(0));
            h[1] = static_cast<double>(p(1));
            h[2] = static_cast<double>(p(2));
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
        //     torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
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
        grad_out(0) = static_cast<T>(fx.detach().to(torch::kCPU).item<double>());
        grad_out(1) = static_cast<T>(fy.detach().to(torch::kCPU).item<double>());
        grad_out(2) = static_cast<T>(ft.detach().to(torch::kCPU).item<double>());

        convert_gradient_to_domain(grad_out);


        auto Hx_cpu = Hx_full.detach().to(torch::kCPU);
        auto Hy_cpu = Hy_full.detach().to(torch::kCPU);

        Hessian.resize(3,3);
        Hessian(0,0) = Hx_cpu.index({0,2}).item<T>(); // fxx
        Hessian(0,1) = Hx_cpu.index({0,1}).item<T>(); // fxy
        Hessian(1,0) = Hessian(0,1);

        Hessian(1,1) = Hy_cpu.index({0,1}).item<T>(); // fyy

        Hessian(0,2) = Hx_cpu.index({0,0}).item<T>(); // fxt
        Hessian(2,0) = Hessian(0,2);

        Hessian(1,2) = Hy_cpu.index({0,0}).item<T>(); // fyt
        Hessian(2,1) = Hessian(1,2);

        Hessian(2,2) = 0.0; // ftt excluded

        convert_hessian_to_domain(Hessian);

        // Third (spatial only): [fxxx, fxxy, fxyy, fyyy]
        third_spatial_out.resize(4);

        auto fxx_cpu = G_fxx.detach().to(torch::kCPU);
        auto fxy_cpu = G_fxy.detach().to(torch::kCPU);
        auto fyy_cpu = G_fyy.detach().to(torch::kCPU);

        third_spatial_out(0) = fxx_cpu.index({0,2}).item<T>(); // fxxx
        third_spatial_out(1) = fxx_cpu.index({0,1}).item<T>(); // fxxy
        third_spatial_out(2) = fxy_cpu.index({0,1}).item<T>(); // fxyy
        third_spatial_out(3) = fyy_cpu.index({0,1}).item<T>(); // fyyy  // Third (with exactly one t): [fxxt, fxyt, fyyt]
        
        third_tmix_out.resize(3);

        third_tmix_out(0) = fxx_cpu.index({0,0}).item<T>(); // fxxt
        third_tmix_out(1) = fxy_cpu.index({0,0}).item<T>(); // fxyt
        third_tmix_out(2) = fyy_cpu.index({0,0}).item<T>(); // fyyt

        convert_third_order_to_domain(third_spatial_out, third_tmix_out);
    }


    void derivative(const VectorX<T>& point) {
        if (!loaded) throw std::runtime_error("INR model not loaded");


        VectorX<T> p = convert_point_to_domain_reverse_order(point);
        // Convert input (float or float) -> float32 tensor

        torch::Tensor input = torch::from_blob(
            (void*)p.data(),
            {1, p.size()},
            torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
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

        grad(0) = static_cast<T>(dx.item<double>());
        grad(1) = static_cast<T>(dy.item<double>());
        grad(2) = static_cast<T>(dz.item<double>());
    }

    // Accurate Hessian via nested autograd (second-order), returned in physical
    // domain coordinates. Follows the autograd scheme of derivative() above and
    // rescales with convert_hessian_to_domain. The full 3x3 Hessian is filled,
    // including the time-time term (f_tt), so the leading 2x2 block is the exact
    // spatial Hessian:
    //   [ f_xx  f_xy  f_xt ]
    //   [ f_xy  f_yy  f_yt ]
    //   [ f_xt  f_yt  f_tt ]
    void query_hessian_autograd(const VectorX<T>& point, Eigen::MatrixX<T>& Hessian)
    {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // model input is in [-1,1] with reversed axis order (t, y, x)
        VectorX<T> p = convert_point_to_domain_reverse_order(point);

        torch::Tensor input = torch::from_blob(
            (void*)p.data(), {1, p.size()},
            torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
        ).clone();
        input = input.to(device, /*non_blocking=*/false, /*copy=*/true);
        input.set_requires_grad(true);

        torch::Tensor output = module.forward({input}).toTensor().reshape({});

        // first-order gradient wrt full input (col 0 = t, 1 = y, 2 = x)
        torch::Tensor g1 = torch::autograd::grad(
            {output}, {input}, /*grad_outputs=*/{},
            /*retain_graph=*/true, /*create_graph=*/true)[0];
        torch::Tensor dx = g1.index({0, 2});
        torch::Tensor dy = g1.index({0, 1});
        torch::Tensor dt = g1.index({0, 0});

        // second-order: differentiate each first-order component again
        torch::Tensor g2x = torch::autograd::grad({dx}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fxt, fxy, fxx]
        torch::Tensor g2y = torch::autograd::grad({dy}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fyt, fyy, fyx]
        torch::Tensor g2t = torch::autograd::grad({dt}, {input}, {}, /*retain*/false, /*create*/false)[0]; // [ftt, fty, ftx]

        auto val = [](const torch::Tensor& t) {
            return static_cast<T>(t.detach().to(torch::kCPU).item<double>());
        };

        const T fxx = val(g2x.index({0, 2}));
        const T fxy = val(g2x.index({0, 1}));
        const T fxt = val(g2x.index({0, 0}));
        const T fyy = val(g2y.index({0, 1}));
        const T fyt = val(g2y.index({0, 0}));
        const T ftt = val(g2t.index({0, 0}));

        Hessian.resize(3, 3);
        Hessian(0, 0) = fxx; Hessian(0, 1) = fxy; Hessian(0, 2) = fxt;
        Hessian(1, 0) = fxy; Hessian(1, 1) = fyy; Hessian(1, 2) = fyt;
        Hessian(2, 0) = fxt; Hessian(2, 1) = fyt; Hessian(2, 2) = ftt;

        // [-1,1] model coords -> physical domain: H_ij *= 4/(range_i*range_j)
        convert_hessian_to_domain(Hessian);
    }

    // All derivatives needed for the spatial-acceleration (d^2x/dt^2) of a
    // critical point, computed via nested autograd and returned in physical
    // domain coordinates. Domain axis order is (x, y, t); the network input is
    // [-1,1] with reversed axis order (t, y, x).
    //
    // Outputs:
    //   H     : 2x2 spatial Hessian   [[f_xx, f_xy], [f_xy, f_yy]]
    //   gt    : d/dt of spatial grad  [f_xt, f_yt]
    //   third : size 9, in the order
    //           [f_xxx, f_xxy, f_xyy, f_yyy,   // pure spatial thirds
    //            f_xxt, f_xyt, f_yyt,          // two-spatial + one-time
    //            f_xtt, f_ytt]                 // one-spatial + two-time
    void query_accel_derivs(const VectorX<T>& point,
                            Eigen::MatrixX<T>& H,
                            VectorX<T>& gt,
                            VectorX<T>& third)
    {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // model input is in [-1,1] with reversed axis order (t, y, x)
        VectorX<T> p = convert_point_to_domain_reverse_order(point);

        torch::Tensor input = torch::from_blob(
            (void*)p.data(), {1, p.size()},
            torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
        ).clone();
        input = input.to(device, /*non_blocking=*/false, /*copy=*/true);
        input.set_requires_grad(true);

        torch::Tensor output = module.forward({input}).toTensor().reshape({});

        // L1: first-order gradient (col 0 = t, 1 = y, 2 = x)
        torch::Tensor g1 = torch::autograd::grad(
            {output}, {input}, {}, /*retain*/true, /*create*/true)[0];
        torch::Tensor fx = g1.index({0, 2});
        torch::Tensor fy = g1.index({0, 1});

        // L2: rows of the Hessian (keep graph for third order)
        torch::Tensor Hx = torch::autograd::grad({fx}, {input}, {}, /*retain*/true, /*create*/true)[0]; // [fxt, fxy, fxx]
        torch::Tensor Hy = torch::autograd::grad({fy}, {input}, {}, /*retain*/true, /*create*/true)[0]; // [fyt, fyy, fyx]

        torch::Tensor fxx = Hx.index({0, 2});
        torch::Tensor fxy = Hx.index({0, 1});
        torch::Tensor fxt = Hx.index({0, 0});
        torch::Tensor fyy = Hy.index({0, 1});
        torch::Tensor fyt = Hy.index({0, 0});

        // L3: differentiate each needed second-order term once more
        torch::Tensor Gxx = torch::autograd::grad({fxx}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fxxt, fxxy, fxxx]
        torch::Tensor Gxy = torch::autograd::grad({fxy}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fxyt, fxyy, fxyx]
        torch::Tensor Gyy = torch::autograd::grad({fyy}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fyyt, fyyy, fyyx]
        torch::Tensor Gxt = torch::autograd::grad({fxt}, {input}, {}, /*retain*/true,  /*create*/false)[0]; // [fxtt, fxty, fxtx]
        torch::Tensor Gyt = torch::autograd::grad({fyt}, {input}, {}, /*retain*/false, /*create*/false)[0]; // [fytt, fyty, fytx]

        auto val = [](const torch::Tensor& t) {
            return static_cast<T>(t.detach().to(torch::kCPU).item<double>());
        };

        // model-coordinate values
        const T m_fxx = val(fxx), m_fxy = val(fxy), m_fyy = val(fyy);
        const T m_fxt = val(fxt), m_fyt = val(fyt);
        const T m_fxxx = val(Gxx.index({0, 2})), m_fxxy = val(Gxx.index({0, 1})), m_fxxt = val(Gxx.index({0, 0}));
        const T m_fxyy = val(Gxy.index({0, 1})), m_fxyt = val(Gxy.index({0, 0}));
        const T m_fyyy = val(Gyy.index({0, 1})), m_fyyt = val(Gyy.index({0, 0}));
        const T m_fxtt = val(Gxt.index({0, 0}));
        const T m_fytt = val(Gyt.index({0, 0}));

        // [-1,1] model coords -> physical domain: each derivative is scaled by
        // prod_axis (2/range_axis)^(order in that axis).
        const T rx = domain_range(0), ry = domain_range(1), rt = domain_range(2);
        const T sx = T(2) / rx, sy = T(2) / ry, st = T(2) / rt;

        H.resize(2, 2);
        H(0, 0) = m_fxx * sx * sx;
        H(0, 1) = m_fxy * sx * sy;
        H(1, 0) = H(0, 1);
        H(1, 1) = m_fyy * sy * sy;

        gt.resize(2);
        gt(0) = m_fxt * sx * st;
        gt(1) = m_fyt * sy * st;

        third.resize(9);
        third(0) = m_fxxx * sx * sx * sx;   // fxxx
        third(1) = m_fxxy * sx * sx * sy;   // fxxy
        third(2) = m_fxyy * sx * sy * sy;   // fxyy
        third(3) = m_fyyy * sy * sy * sy;   // fyyy
        third(4) = m_fxxt * sx * sx * st;   // fxxt
        third(5) = m_fxyt * sx * sy * st;   // fxyt
        third(6) = m_fyyt * sy * sy * st;   // fyyt
        third(7) = m_fxtt * sx * st * st;   // fxtt
        third(8) = m_fytt * sy * st * st;   // fytt
    }

    // Finite-difference counterpart of query_accel_derivs. All derivatives are
    // obtained from central finite differences of the (autograd) gradient field
    // (step delta_h*range) instead of exact nested autograd, and returned in
    // physical-domain coordinates. Same outputs and ordering as
    // query_accel_derivs:
    //   H     : 2x2 spatial Hessian   [[f_xx, f_xy], [f_xy, f_yy]]
    //   gt    : d/dt of spatial grad  [f_xt, f_yt]
    //   third : size 9
    //           [f_xxx, f_xxy, f_xyy, f_yyy,   f_xxt, f_xyt, f_yyt,   f_xtt, f_ytt]
    // Unlike the exact autograd version, the finite differences low-pass smooth
    // the INR's high-frequency curvature, typically yielding smaller, less noisy
    // velocities/accelerations.
    void query_accel_derivs_fd(const VectorX<T>& point,
                               Eigen::MatrixX<T>& H,
                               VectorX<T>& gt,
                               VectorX<T>& third)
    {
        if (!loaded) throw std::runtime_error("INR model not loaded");

        // 1) Full 3x3 FD Hessian (f_tt = 0), spatial and time-mixed third
        //    derivatives, and the two-time thirds (f_xtt, f_ytt) -- all from the
        //    single central-difference stencil of the gradient field.
        VectorX<T>        grad, third_spatial, third_tmix;
        Eigen::MatrixX<T> Hfull;
        T fxtt = T(0), fytt = T(0);
        query_up_to_third_derivative(point, grad, Hfull, third_spatial, third_tmix,
                                     &fxtt, &fytt);

        // 2) Assemble outputs in the query_accel_derivs layout
        H.resize(2, 2);
        H(0, 0) = Hfull(0, 0);
        H(0, 1) = Hfull(0, 1);
        H(1, 0) = Hfull(1, 0);
        H(1, 1) = Hfull(1, 1);

        gt.resize(2);
        gt(0) = Hfull(0, 2);   // f_xt
        gt(1) = Hfull(1, 2);   // f_yt

        third.resize(9);
        third(0) = third_spatial(0);   // fxxx
        third(1) = third_spatial(1);   // fxxy
        third(2) = third_spatial(2);   // fxyy
        third(3) = third_spatial(3);   // fyyy
        third(4) = third_tmix(0);      // fxxt
        third(5) = third_tmix(1);      // fxyt
        third(6) = third_tmix(2);      // fyyt
        third(7) = fxtt;               // fxtt
        third(8) = fytt;               // fytt
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
        auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(device);
        input_    = torch::zeros({1,3}, opts);
        grad_buf_ = torch::empty({1,3}, opts);

        if (device.is_cuda()) {
            input_host_ = torch::empty({1,3},
                torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU).pinned_memory(true));
        } else {
            input_host_ = torch::Tensor();
        }

        at::set_num_threads(1);
        at::set_num_interop_threads(1);
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