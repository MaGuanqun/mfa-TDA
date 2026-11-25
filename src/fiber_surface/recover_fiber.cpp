#include "trace.h"

namespace fiber_surface {
// Compare ordered curves at equal arclength, including every chord boundary and
// the exact end of the required overlap. Neither input may run out early.
bool agrees_over_arclength(const std::vector<Vec>& line, size_t edge, Vec a,
                           const std::vector<Vec>& reference, double length, double tolerance) {
    if (edge>=line.size() || reference.size()<2 || length<=0) return false;
    Vec b=reference.front();
    size_t j=0;
    const double roundoff=1e-12*length;
    if ((a-b).norm()>tolerance) return false;
    while (length>roundoff && edge+1<line.size() && j+1<reference.size()) {
        const Vec u=line[edge+1]-a, v=reference[j+1]-b;
        const double la=u.norm(), lb=v.norm();
        if (la<=roundoff) { a=line[++edge]; continue; }
        if (lb<=roundoff) { b=reference[++j]; continue; }
        if (u.dot(v)<std::sqrt(.5)*la*lb) return false;
        const double ds=std::min({la,lb,length});
        a+=u*(ds/la); b+=v*(ds/lb);
        if ((a-b).norm()>tolerance) return false;
        length-=ds;
        if (ds==la) a=line[++edge];
        if (ds==lb) b=reference[++j];
    }
    return length<=roundoff;
}

namespace {
// Extend a closed curve through its seam once to allow an overlap that wraps.
std::vector<Vec> ordered_curve(std::vector<Vec> points, bool closed, int sign) {
    if (sign<0) std::reverse(points.begin(),points.end());
    if (closed) {
        const auto lap=points;
        points.insert(points.end(),lap.begin(),lap.end());
    }
    return points;
}

bool follows_line(const std::vector<Vec>& points, const std::vector<Vec>& line, bool closed, double tolerance) {
    const size_t edges=line.size()-(closed ? 0 : 1);
    size_t near=0;
    for (const Vec& p:points) {
        bool found=false;
        for (size_t j=0; j<edges; ++j) {
            const size_t i=(near+j)%edges;
            const Vec closest=marching_triangles::closest_on_segment(p,line[i],line[(i+1)%line.size()]);
            if ((closest-p).norm()<=tolerance) { near=i; found=true; break; }
        }
        if (!found) return false;
    }
    return true;
}
}

bool FiberTracer::shared_arc(const std::vector<Vec>& line, const std::vector<Vec>& reference,
                             const Range& q, int sign, double spacing, size_t start_edges) const {
    if (line.size()<2 || reference.size()<2) return false;
    const Vec& anchor=reference.front();
    Vec normal;
    if (!c_.tangent(anchor,normal)) return false;
    normal*=sign;
    const double eps=c_.options.fiber_plane_tolerance();
    const double tube=spacing/4;
    for (size_t i=0; i<std::min(start_edges,line.size()-1); ++i) {
        const double a=normal.dot(line[i]-anchor), b=normal.dot(line[i+1]-anchor);
        if (a>eps || b< -eps || b<=a) continue;
        const Vec crossing=line[i]+std::clamp(-a/(b-a),0.,1.)*(line[i+1]-line[i]);
        if ((crossing-anchor).norm()>tube) continue;
        // Distinct nearby fibers may be almost parallel for the whole overlap.
        // Correct their crossings in the same transverse plane before comparing.
        Vec root=crossing;
        if (!c_.fiber(root,q,normal,anchor,tube) || (root-anchor).norm()>4*eps) continue;
        if (agrees_over_arclength(line,i,crossing,reference,c_.options.step,tube)) return true;
    }
    return false;
}

FiberTracer::Trace FiberTracer::recovery_walk(const Vec& seed, const Range& q, int sign,
                                               double spacing, double length) const {
    Trace trace;
    trace.points.push_back(seed);
    Sample current{seed,Vec::Zero()};
    if (!c_.tangent(seed,current.tangent)) return trace;
    current.tangent*=sign;
    const Vec seed_tangent=current.tangent;
    const double h=c_.options.step;
    double traveled=0;
    for (int k=0; k<c_.options.max_trace_steps; ++k) {
        Sample next;
        if (advance(current,q,sign,spacing,next)) return trace;
        const double ds=(next.point-current.point).norm();
        if (next.on_box) {
            if (ds>h*1e-7) trace.points.push_back(next.point);
            trace.complete=true;
            return trace;
        }
        if (ds<spacing*.1 || (next.point-current.point).dot(current.tangent)<=0) return trace;
        if (!std::isfinite(length) && k>=2 &&
            (next.point-seed).norm()<=1.5*spacing &&
            seed_tangent.dot(current.point-seed)<=c_.options.fiber_plane_tolerance() &&
            seed_tangent.dot(next.point-seed)>=-c_.options.fiber_plane_tolerance()) {
            // A sub-resolution loop must not accumulate the overlap by taking extra laps.
            if (traveled<h) return trace;
            auto overlap=recovery_walk(next.point,q,sign,spacing,h+spacing);
            overlap.points.insert(overlap.points.begin(),current.point);
            if (shared_arc(overlap.points,trace.points,q,sign,spacing,1)) {
                // The return arc duplicates the beginning. Keep just one lap.
                if ((current.point-seed).norm()<=4*c_.options.fiber_plane_tolerance()) trace.points.pop_back();
                trace.closed=true; trace.complete=true;
                return trace;
            }
        }
        trace.points.push_back(next.point);
        traveled+=ds;
        if (traveled>=length) return trace;
        current=next;
    }
    return trace;
}

FiberTracer::Witnesses FiberTracer::witnesses(const Trace& fragment, const Vec& seed, const Range& q,
                                              bool stalled_backward, bool stalled_forward) const {
    const double h=c_.options.step, spacing=std::max(c_.options.min_step,h/8);
    Witnesses w;
    auto arc=[&](const Vec& p,int sign) { return Arc{sign,recovery_walk(p,q,sign,spacing,h+spacing).points}; };
    // Reconstruct the overlap into the fragment, opposite the stalled direction.
    if (stalled_backward) w.tips.push_back(arc(fragment.points.front(),1));
    if (stalled_forward) w.tips.push_back(arc(fragment.points.back(),-1));
    for (const Vec& end:{fragment.points.front(),fragment.points.back()})
        if (c_.field.box.faces(end,h*1e-8)) w.contacts.push_back(end);
    w.seed=arc(seed,1);
    double length=0;
    for (size_t i=1;i<w.seed.points.size();++i) length+=(w.seed.points[i]-w.seed.points[i-1]).norm();
    if (length<h) w.seed=arc(seed,-1);
    return w;
}

bool FiberTracer::contains_fragment(const Trace& candidate, const Range& q, const Witnesses& w) const {
    const double spacing=std::max(c_.options.min_step,c_.options.step/8);
    for (const Vec& end:w.contacts) {
        if (candidate.closed || std::min((end-candidate.points.front()).norm(),(end-candidate.points.back()).norm())
            >4*c_.options.fiber_plane_tolerance()) return false;
    }
    const auto positive=ordered_curve(candidate.points,candidate.closed,1);
    const auto negative=ordered_curve(candidate.points,candidate.closed,-1);
    auto matches=[&](const Arc& arc) {
        const auto& path=arc.sign>0 ? positive : negative;
        return shared_arc(path,arc.points,q,arc.sign,spacing,candidate.points.size());
    };
    return matches(w.seed) && std::all_of(w.tips.begin(),w.tips.end(),matches);
}

FiberTracer::Recovery FiberTracer::recover(Trace& fragment, const Range& q, const std::vector<Vec>& seeds,
                          size_t seed_index, const Witnesses& w) const {
    const double h=c_.options.step, spacing=std::max(c_.options.min_step,h/8);
    std::optional<Trace> replacement;
    for (const auto& tip:w.tips) {
        for (size_t i=seed_index+1; i<seeds.size(); ++i) {
            const double distance=(seeds[i]-tip.points.front()).norm();
            if (distance>h || distance<=4*c_.options.fiber_plane_tolerance()) continue;
            const auto approach=recovery_walk(seeds[i],q,tip.sign,spacing,3*h);
            if (!shared_arc(approach.points,tip.points,q,tip.sign,spacing,approach.points.size())) continue;

            if (replacement) {
                const auto line=ordered_curve(replacement->points,replacement->closed,tip.sign);
                if (shared_arc(line,approach.points,q,tip.sign,spacing,replacement->points.size())) continue;
            }
            // Reconstruct from the independent seed, trimming the duplicated arc.
            // Publish it atomically only after accounting for the old seed and
            // every stalled end, including the opposite end of a partial trace.
            auto candidate=recovery_walk(seeds[i],q,1,spacing,std::numeric_limits<double>::infinity());
            if (!candidate.complete) continue;
            if (!candidate.closed) {
                auto backward=recovery_walk(seeds[i],q,-1,spacing,std::numeric_limits<double>::infinity());
                if (!backward.complete || backward.closed) continue;
                std::reverse(backward.points.begin(),backward.points.end());
                backward.points.pop_back();
                backward.points.insert(backward.points.end(),candidate.points.begin(),candidate.points.end());
                candidate.points=std::move(backward.points);
            }
            if (!contains_fragment(candidate,q,w)) continue;
            if (!replacement) replacement=std::move(candidate);
            else if (!same_component(*replacement,candidate,q)) return Recovery::ambiguous;
        }
    }
    if (!replacement) return Recovery::unresolved;
    fragment=std::move(*replacement);
    return Recovery::complete;
}

bool FiberTracer::same_component(const Trace& a, const Trace& b, const Range& q) const {
    if (a.closed!=b.closed || a.points.size()<2 || b.points.size()<2) return false;
    const double spacing=std::max(c_.options.min_step,c_.options.step/8);
    if (!a.closed && ((a.points.front()-b.points.front()).norm()>4*c_.options.fiber_plane_tolerance() ||
                      (a.points.back()-b.points.back()).norm()>4*c_.options.fiber_plane_tolerance())) return false;
    const auto line=ordered_curve(b.points,b.closed,1);
    return shared_arc(line,a.points,q,1,spacing,b.points.size()) &&
           follows_line(a.points,b.points,b.closed,spacing/4) && follows_line(b.points,a.points,a.closed,spacing/4);
}
} // namespace fiber_surface
