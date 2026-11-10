#include "sampling_grid.h"
#include <iostream>

using namespace local_reference;
void check(bool ok, const char* message) {if(!ok) throw std::runtime_error(message);}
int main() {
    try {
        for(const auto& name : {"boussinesq_3d", "vortex_street_3d"}) {
            const auto d = dataset(name);
            Triple lo{-.25,-.25,d.min[2]+.125}, hi{.25,.25,d.min[2]+.375};
            Grid g(d,lo,hi,{4,4,8},{0,0,0});
            check(g.dimensions == std::array<std::size_t,3>{21,21,21},"physical ratio/count mismatch");
            check(g.point(0) == lo && g.point(g.size-1) == hi,"box endpoints not preserved");
            check(std::abs(g.point(1)[0]-(lo[0]+.025)) < 1e-12,"x must vary fastest");
            check(g.point(21)[1] == lo[1]+g.spacing[1],"y stride mismatch");
            check(g.point(441)[2] == lo[2]+g.spacing[2],"time stride mismatch");
            Grid padded(d,lo,hi,{4,4,8},{100,100,100});
            check(padded.min == d.min && padded.max == d.max,"padding must clip to dataset bounds");
            bool rejected = false;
            try {Grid bad(d,hi,lo,{1,1,1},{0,0,0});} catch(const std::invalid_argument&) {rejected=true;}
            check(rejected,"reversed box accepted");
            rejected=false;
            try {Grid bad(d,lo,hi,{0,1,1},{0,0,0});} catch(const std::invalid_argument&) {rejected=true;}
            check(rejected,"zero ratio accepted");
            auto outside=lo;outside[0]=d.min[0]-.01;rejected=false;
            try {Grid bad(d,outside,hi,{1,1,1},{0,0,0});} catch(const std::invalid_argument&) {rejected=true;}
            check(rejected,"outside box accepted");
            auto tiny_hi=lo;
            for(auto& x:tiny_hi) x+=.001;
            Grid tiny(d,lo,tiny_hi,{1,1,1},{0,0,0});
            check(tiny.dimensions == std::array<std::size_t,3>{2,2,3},"minimum grid for legacy TTK");
        }
        check(uniform_span_count({0,0,0,.25,.5,.75,1,1,1}) == 4,
              "Repeated end knots must not count as spans");
        check(uniform_span_count({0,0,.5,.5,1,1}) == 2,
              "Repeated interior knots must not count as spans");
        bool rejected = false;
        try {uniform_span_count({0,.25,1});} catch(const std::invalid_argument&) {rejected=true;}
        check(rejected,"Nonuniform knots accepted for fixed points per span");
        const Dataset d{{-.5,-.5,13.5},{7.5,.5,15},{4,2,3}};
        Grid a(d,{.51,-.1,14.3},{.56,.1,14.6},{160,160,160},{0,0,0},true);
        Grid b(d,{.36,-.1,13.6},{.76,.42,14.},{160,160,160},{0,0,0},true);
        for(int axis=0;axis<3;++axis) {
            check(a.spacing[axis] == b.spacing[axis],"Spacing changes across boxes");
            check(std::abs(a.spacing[axis]-(d.max[axis]-d.min[axis])/d.base_intervals[axis]/160)<1e-14,
                  "160 subdivisions per span must have span_width/160 spacing");
        }
        check(a.min[0]<=.51 && a.max[0]>=.56 && .51-a.min[0]<a.spacing[0] && a.max[0]-.56<a.spacing[0],
              "Partial spans must snap outward by less than one step");
        Grid full(d,d.min,d.max,{160,160,160},{0,0,0},true);
        check(full.dimensions == std::array<std::size_t,3>{641,321,481},"Full span counts are wrong");
        check(full.point(full.size-1)==d.max,"Snapped grid must preserve full-domain upper bounds");
        Grid edge(d,{7.499,.499,14.999},d.max,{160,160,160},{0,0,0},true);
        check(edge.dimensions == std::array<std::size_t,3>{2,2,3},"Tiny end box needs three valid time steps");
        std::cout << "Local sampling grid checks passed\n";
    } catch(const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
