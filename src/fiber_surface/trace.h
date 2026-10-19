#pragma once
#include "corrector.h"

namespace fiber_surface {
class FiberTracer {
    const Corrector& c_;
    Result& result_;
    struct Trace { std::vector<Vec> points; bool closed=false, complete=false; };
    void event(const Vec& x, const std::string& why);
    bool contact(const Vec& p, const Vec& trial, const Range& q, Vec& x) const;
    Trace direction(const Vec& seed, const Range& q, double sign);
    bool covered(const Vec& seed, int endpoint) const;
public:
    FiberTracer(const Corrector& corrector, Result& result) : c_(corrector), result_(result) {}
    void run(const std::vector<Vec>& grid);
    void run_seeds(const std::array<std::vector<Vec>,2>& endpoints);
};
} // namespace fiber_surface
