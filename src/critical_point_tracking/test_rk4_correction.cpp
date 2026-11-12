#include <mfa/mfa.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include "RK4.h"

int main()
{
    const Eigen::VectorXd mins = closed_form_function::domain_min(2);
    const Eigen::VectorXd maxs = closed_form_function::domain_max(2);
    Eigen::VectorXd seed(3);
    seed << 1.0, std::sqrt(-std::cos(2.0)), 2.0;
    auto residual = [](const Eigen::VectorXd& p) {
        return std::hypot(p[0] * (p[0] * p[0] + 1.0 - p[2]),
                          p[1] * (p[1] * p[1] + std::cos(p[2])));
    };
    auto require = [](bool condition, const char* message) {
        if (!condition)
            std::cerr << message << '\n';
        return condition;
    };

    // Exercise temporal and spatial predictors in both tracking directions.
    for (bool forward : {false, true})
    {
        for (double spatial_step : {0.02, 0.2})
        {
            Eigen::VectorXd predicted(3), corrected(3);
            if (!require(RK4::RK4_correction(seed, predicted, 0.1, spatial_step,
                    1e-12, forward, 0, 1.0, mins, maxs, 2), "Prediction failed"))
                return 1;
            if (!require(residual(predicted) > 1e-12,
                         "Uncorrected prediction unexpectedly satisfies Newton tolerance"))
                return 1;
            Eigen::VectorXd unchanged = predicted;
            if (!require(RK4::correction_newton(unchanged, 0, 1.0, 1e-12, mins, maxs, 2)
                         && (unchanged - predicted).norm() == 0.0,
                         "Disabled correction changed the prediction"))
                return 1;
            if (!require(RK4::RK4_correction(seed, corrected, 0.1, spatial_step,
                    1e-12, forward, 30, 1.0, mins, maxs, 2)
                         && residual(corrected) < 1e-12,
                         "Enabled correction failed to converge"))
                return 1;
        }
    }

    Eigen::VectorXd boundary(3), stale = seed;
    boundary << std::sqrt(3.0), std::sqrt(-std::cos(4.0)), 4.0;
    if (!require(!RK4::RK4_choose_direction(boundary, stale, 0.1, 0.1, 1e-12,
                    true, 0, 1.0, false, mins, maxs, 2),
                 "Failed RK4 predictors accepted stale output"))
        return 1;
    Eigen::VectorXd invalid = seed;
    invalid[0] = 3.0;
    if (!require(!RK4::correction_newton(invalid, 0, 1.0, 1e-12, mins, maxs, 2),
                 "Accepted an out-of-domain prediction"))
        return 1;
    invalid[0] = std::numeric_limits<double>::quiet_NaN();
    if (!require(!RK4::correction_newton(invalid, 0, 1.0, 1e-12, mins, maxs, 2),
                 "Accepted a nonfinite prediction"))
        return 1;
    std::cout << "RK4 correction ablation checks passed\n";
    return 0;
}
