#include <mfa/mfa.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include "RK4.h"
#include "particle_tracing.h"
#include "trace_deduplication.h"

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
    // Reproduce boundary-seeded rotating traces in both directions, including
    // the stationary center and curves using the spatial predictor. Every trace
    // must reach BOTH time faces before deduplication chooses one direction.
    std::vector<CP_Trace<double>> traces;
    for (int u : {-1, 0, 1})
    {
        for (int v : {-1, 0, 1})
        {
            for (bool forward : {false, true})
            {
                const double start_time = forward ? 0.0 : 4.0;
                Eigen::VectorXd start(3);
                start << u*std::cos(start_time)-v*std::sin(start_time),
                         u*std::sin(start_time)+v*std::cos(start_time), start_time;
                CP_Trace<double> trace;
                particle_tracing::tracing_one_direction(0.1, 0.1, start, trace.traces,
                    forward, 1e-12, 30, 1.0, mins, maxs, mins, maxs, 3);
                if (!require(trace.traces.size() > 2 && (trace.traces.front()-start).norm() == 0.0,
                             "Boundary-seeded trace lost its initial point")
                    || !require(trace.traces.back()[2] == (forward ? 4.0 : 0.0),
                                "Boundary-seeded trace stopped before the opposite time face"))
                    return 1;
                double previous_time = start_time;
                for (size_t i = 0; i < trace.traces.size(); ++i)
                {
                    const auto& p = trace.traces[i];
                    const double uu = p[0]*std::cos(p[2])+p[1]*std::sin(p[2]);
                    const double vv = -p[0]*std::sin(p[2])+p[1]*std::cos(p[2]);
                    if (!require(p.allFinite() && (p.array() >= mins.array()).all()
                                 && (p.array() <= maxs.array()).all()
                                 && std::hypot(uu-u, vv-v) < 1e-9,
                                 "Completed boundary trace left its analytic component"))
                        return 1;
                    if (i && !require(forward ? p[2] > previous_time : p[2] < previous_time,
                                      "Boundary completion did not advance monotonically"))
                        return 1;
                    previous_time = p[2];
                }
                if (!forward)
                    std::reverse(trace.traces.begin(), trace.traces.end());
                traces.push_back(std::move(trace));
            }
        }
    }
    std::vector<Eigen::VectorXd> events;
    deduplication::deduplicate_traces(traces, events, 0.1, 0.1, mins, maxs);
    int kept = 0;
    for (const auto& trace : traces)
    {
        if (trace.duplicated)
            continue;
        ++kept;
        if (!require(trace.traces.front()[2] == 0.0 && trace.traces.back()[2] == 4.0,
                     "Deduplication retained an incomplete component"))
            return 1;
    }
    if (!require(kept == 9, "Boundary completion changed the nine-component count"))
        return 1;

    // A final interval shorter than the stagnation threshold still needs to be
    // emitted. Test both directions, including Newton-disabled predictions.
    for (bool forward : {false, true})
    {
        for (int corrections : {0, 30})
        {
            const double time = forward ? 4.0-1e-4 : 1e-4;
            Eigen::VectorXd start(3);
            start << std::cos(time), std::sin(time), time;
            std::vector<Eigen::VectorXd> result;
            particle_tracing::tracing_one_direction(0.1, 0.1, start, result,
                forward, 1e-12, corrections, 1.0, mins, maxs, mins, maxs, 3);
            if (!require(result.size() == 2 && result.back()[2] == (forward ? 4.0 : 0.0),
                         "Short terminal step was discarded as stagnation"))
                return 1;
            Eigen::VectorXd endpoint = result.back(), unused(3);
            if (!require(!RK4::RK4_correction(endpoint, unused, 0.1, 0.1, 1e-12,
                         forward, corrections, 1.0, mins, maxs, 3),
                         "An outward step from the time boundary should terminate"))
                return 1;
        }
    }
    std::cout << "RK4 correction and boundary completion checks passed\n";
    return 0;
}
