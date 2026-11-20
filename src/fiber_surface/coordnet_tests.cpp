#include "pipeline.h"
#include <torch/script.h>
#include <iostream>

using namespace fiber_surface;

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("Expected a fixture directory");
        const fs::path dir = argv[1];
        fs::create_directories(dir);
        RunConfig config;
        config.model = "coordnet";
        config.first = (dir / "quadratic.pt").string();
        config.second = (dir / "affine.pt").string();
        torch::jit::Module quadratic("Quadratic"), affine("Affine");
        quadratic.define(R"(
def forward(self, p: Tensor) -> Tensor:
    t, y, x = p[:,0], p[:,1], p[:,2]
    return x*x + 2*y*y + 3*t*t + 4*x*y + 5*x*t + 6*y*t
)");
        affine.define(R"(
def forward(self, p: Tensor) -> Tensor:
    return 7*p[:,2] - 3*p[:,1] + 2*p[:,0]
)");
        quadratic.save(config.first);
        affine.save(config.second);
        config.first_digest = file_digest(config.first);
        config.second_digest = file_digest(config.second);
        config.training_box = {Vec(-2, 1, 10), Vec(6, 3, 14)};
        config.box = {Vec(1, 1.2, 11), Vec(3, 2, 12)};
        write_config(dir, config);
        auto loaded = load_field(read_config(dir), true);
        const Field& field = loaded.field;
        if (field.box.min != config.box.min || field.box.max != config.box.max)
            throw std::runtime_error("Saved extraction box changed");

        // In normalized (x,y,t) coordinates the first field is q^T H q / 2.
        Hessian h;
        h << 2,4,5, 4,4,6, 5,6,6;
        const Vec scale(0.25, 1, 0.5), linear(7, -3, 2);
        const Hessian physical_h = scale.asDiagonal() * h * scale.asDiagonal();
        for (const Vec& x : discovery_grid(field.box, 3)) {
            const Vec q = (x-config.training_box.min).cwiseProduct(scale)-Vec::Ones();
            Range expected(0.5*q.dot(h*q), linear.dot(q));
            Jacobian j;
            j.row(0) = (h*q).cwiseProduct(scale).transpose();
            j.row(1) = linear.cwiseProduct(scale).transpose();
            const auto second = field.hessians(x);
            if (!field.value(x).isApprox(expected, 1e-12) || !field.jacobian(x).isApprox(j, 1e-12) ||
                !second[0].isApprox(physical_h, 1e-12) || !second[1].isZero(1e-12))
                throw std::runtime_error("Incorrect CoordNet values or physical derivatives in cropped domain");
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
