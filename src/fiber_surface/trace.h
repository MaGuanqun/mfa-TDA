#pragma once
#include "corrector.h"

namespace fiber_surface {
// The starting point lies on line[edge..edge+1]. Both polylines must supply the requested length.
bool agrees_over_arclength(const std::vector<Vec>& line, size_t edge, Vec start,
                           const std::vector<Vec>& reference, double length, double tolerance);
class FiberTracer {
    const Corrector& c_;
    Result& result_;
    struct Trace { std::vector<Vec> points; bool closed=false, complete=false; };
    struct Sample { Vec point, tangent=Vec::Zero(); bool on_box=false; };
    struct Arc { int sign=1; std::vector<Vec> points; };
    struct Witnesses { std::vector<Arc> tips; Arc seed; std::vector<Vec> contacts; };
    enum class Recovery { unresolved, complete, ambiguous };
    struct Fragment {
        size_t seed_index;
        Trace trace;
        std::vector<Event> events;
        std::optional<Witnesses> evidence;
        Recovery recovery=Recovery::unresolved;
    };
    bool contact(const Vec& p, const Vec& trial, const Range& q, double h, Vec& x) const;
    const char* try_step(const Sample& start, const Range& q, double sign, double h, Sample& next) const;
    const char* advance(const Sample& start, const Range& q, double sign, double h, Sample& next) const;
    Trace direction(const Vec& seed, const Range& q, double sign, std::vector<Event>& events) const;
    Fragment discover(const Vec& seed, const Range& q, size_t seed_index) const;
    Trace recovery_walk(const Vec& seed, const Range& q, int sign, double spacing, double length) const;
    bool shared_arc(const std::vector<Vec>& line, const std::vector<Vec>& reference,
                    const Range& q, int sign, double spacing, size_t start_edges) const;
    Witnesses witnesses(const Trace& fragment, const Vec& seed, const Range& q,
                        bool stalled_backward, bool stalled_forward) const;
    bool contains_fragment(const Trace& candidate, const Range& q, const Witnesses& witnesses) const;
    Recovery recover(Trace& fragment, const Range& q, const std::vector<Vec>& seeds,
                     size_t seed_index, const Witnesses& witnesses) const;
    bool same_component(const Trace& a, const Trace& b, const Range& q) const;
    bool covered(const Vec& seed, int endpoint) const;
public:
    FiberTracer(const Corrector& corrector, Result& result) : c_(corrector), result_(result) {}
    void run(const std::vector<Vec>& grid);
    void run_seeds(const std::array<std::vector<Vec>,2>& endpoints);
};
} // namespace fiber_surface
