#include "coordnet_field.h"
#include <torch/script.h>
#include <torch/csrc/autograd/autograd.h>

namespace fiber_surface {
namespace {
torch::Tensor point_tensor(const Vec& point) {
    return torch::tensor({point.x(), point.y(), point.z()}, torch::kFloat64).view({1, 3});
}

class CoordNetPair {
    std::array<torch::jit::Module, 2> modules_;
    torch::Tensor origin_, scale_;

public:
    CoordNetPair(const std::string& first, const std::string& second, const Box& box)
        : modules_{torch::jit::load(first, torch::kCPU), torch::jit::load(second, torch::kCPU)},
          origin_(point_tensor(box.min)), scale_(point_tensor(2.0 * (box.max-box.min).cwiseInverse())) {
        for (auto& module : modules_) module.eval();
    }

    torch::Tensor normalized(const torch::Tensor& xyz) const {
        return ((xyz - origin_) * scale_ - 1).flip({1});
    }

    torch::Tensor scalar(int component, const torch::Tensor& input) {
        auto output = modules_[component].forward({input});
        if (!output.isTensor() || output.toTensor().numel() != 1)
            throw std::invalid_argument("Each CoordNet must return one scalar per input point");
        return output.toTensor().reshape({});
    }
};

// Affine and constant fields have exact zero derivatives at higher orders.
torch::Tensor gradient(const torch::Tensor& scalar, const torch::Tensor& point,
                       bool create_graph = false, bool retain_graph = false) {
    if (scalar.requires_grad()) {
        auto grad = torch::autograd::grad({scalar}, {point}, {}, retain_graph || create_graph, create_graph, true).front();
        if (grad.defined()) return grad;
    }
    return torch::zeros_like(point);
}

Vec vector_from_tensor(const torch::Tensor& tensor) {
    auto data = tensor.contiguous();
    return Eigen::Map<const Vec>(data.data_ptr<double>());
}
} // namespace

Field coordnet_field(const std::string& first, const std::string& second, const Box& training_box) {
    if (!training_box.min.allFinite() || !training_box.max.allFinite() ||
        (training_box.max.array() <= training_box.min.array()).any())
        throw std::invalid_argument("CoordNet requires a finite training domain with positive extents");
    auto models = std::make_shared<CoordNetPair>(first, second, training_box);
    Field field;
    field.box = training_box;
    field.value = [models](const Vec& x) {
        torch::NoGradGuard guard;
        auto input = models->normalized(point_tensor(x));
        return Range(models->scalar(0, input).item<double>(), models->scalar(1, input).item<double>());
    };
    // Differentiate through normalization and axis reversal with respect to
    // physical (x,y,t). Each call owns its input and autograd graph.
    field.jacobian = [models](const Vec& x) {
        torch::AutoGradMode guard(true);
        auto point = point_tensor(x).requires_grad_(true);
        Jacobian result;
        for (int i = 0; i < 2; ++i) {
            auto value = models->scalar(i, models->normalized(point));
            result.row(i) = vector_from_tensor(gradient(value, point)).transpose();
        }
        return result;
    };
    field.hessians = [models](const Vec& x) {
        torch::AutoGradMode guard(true);
        auto point = point_tensor(x).requires_grad_(true);
        std::array<Hessian, 2> result;
        for (int i = 0; i < 2; ++i) {
            auto value = models->scalar(i, models->normalized(point));
            auto grad = gradient(value, point, true).view({3});
            for (int d = 0; d < 3; ++d) {
                auto derivative = gradient(grad[d], point, false, d < 2);
                result[i].row(d) = vector_from_tensor(derivative).transpose();
            }
        }
        return result;
    };
    return field;
}
} // namespace fiber_surface
