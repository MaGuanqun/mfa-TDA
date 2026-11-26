#include "coordnet_field.h"
#include <torch/script.h>
#include <torch/csrc/autograd/autograd.h>
#include <mutex>

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

// Adjacent triangles, Newton corrections and prospective front edges often
// visit identical points. Cache numerical results only, scoped to one model
// pair and training box. Concurrent misses may compute the same sample.
class SampleCache {
    struct Entry {
        Vec point;
        std::optional<FirstOrderSample> sample;
    };
    std::array<Entry,65536> entries_;
    std::array<std::mutex,256> locks_;
    size_t slot(const Vec& point) const {
        size_t h=0;
        for (int d=0; d<3; ++d)
            h ^= std::hash<double>{}(point[d])+0x9e3779b9+(h<<6)+(h>>2);
        return h % entries_.size();
    }
public:
    std::optional<FirstOrderSample> find(const Vec& point) {
        const size_t i=slot(point);
        std::lock_guard<std::mutex> lock(locks_[i % locks_.size()]);
        const auto& entry=entries_[i];
        if (entry.sample && entry.point==point) return entry.sample;
        return {};
    }
    void save(const Vec& point,const FirstOrderSample& sample) {
        const size_t i=slot(point);
        std::lock_guard<std::mutex> lock(locks_[i % locks_.size()]);
        entries_[i]={point,sample};
    }
};
} // namespace

Field coordnet_field(const std::string& first, const std::string& second, const Box& training_box) {
    if (!training_box.min.allFinite() || !training_box.max.allFinite() ||
        (training_box.max.array() <= training_box.min.array()).any())
        throw std::invalid_argument("CoordNet requires a finite training domain with positive extents");
    auto models = std::make_shared<CoordNetPair>(first, second, training_box);
    auto cache = std::make_shared<SampleCache>();
    Field field;
    field.box = training_box;
    field.value = [models,cache](const Vec& x) {
        if (auto cached=cache->find(x)) return cached->value;
        torch::NoGradGuard guard;
        auto input = models->normalized(point_tensor(x));
        return Range(models->scalar(0, input).item<double>(), models->scalar(1, input).item<double>());
    };
    // Differentiate through normalization and axis reversal with respect to
    // physical (x,y,t). Each call owns its input and autograd graph.
    field.value_jacobian = [models,cache](const Vec& x) {
        if (auto cached=cache->find(x)) return *cached;
        torch::AutoGradMode guard(true);
        auto point = point_tensor(x).requires_grad_(true);
        FirstOrderSample result;
        for (int i = 0; i < 2; ++i) {
            auto value = models->scalar(i, models->normalized(point));
            result.value[i] = value.item<double>();
            result.jacobian.row(i) = vector_from_tensor(gradient(value, point)).transpose();
        }
        cache->save(x,result);
        return result;
    };
    field.jacobian = [sample=field.value_jacobian](const Vec& x) { return sample(x).jacobian; };
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
