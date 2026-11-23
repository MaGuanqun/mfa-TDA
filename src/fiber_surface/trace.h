#pragma once
#include "corrector.h"

namespace fiber_surface {
class FiberTracer {
    const Corrector& c_;
    Result& result_;
    struct Trace { std::vector<Vec> points; bool closed=false, complete=false; };
    struct Sample { Vec point, tangent=Vec::Zero(); bool on_box=false; };
    void event(const Vec& x, const std::string& why);
    bool contact(const Vec& p, const Vec& trial, const Range& q, double h, Vec& x) const;
    const char* try_step(const Sample& start, const Range& q, double sign, double h, Sample& next) const;
    const char* advance(const Sample& start, const Range& q, double sign, double h, Sample& next) const;
    Trace direction(const Vec& seed, const Range& q, double sign);
    bool covered(const Vec& seed, int endpoint) const;
public:
    FiberTracer(const Corrector& corrector, Result& result) : c_(corrector), result_(result) {}
    void run(const std::vector<Vec>& grid);
    void run_seeds(const std::array<std::vector<Vec>,2>& endpoints);
};
} // namespace fiber_surface
