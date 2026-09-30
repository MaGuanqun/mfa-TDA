#pragma once
#include "surface_field.h"
#include <array>
#include <map>
#include <vector>

namespace marching_triangles {
// Sampling finds starting points, not a proof that every component was found.
// Reuse the same safeguarded projector as triangle growth.
inline std::vector<Point<double>> find_initial_roots(const SurfaceField<double>& field,
    const MeshingOptions<double>& options,int samples_per_axis,double separation) {
    options.validate();
    if(samples_per_axis<2||samples_per_axis>500||!std::isfinite(separation)||separation<=0)
        throw std::invalid_argument("Invalid root sampling resolution or separation");
    if((options.domain_max-options.domain_min).maxCoeff()/separation > std::numeric_limits<int>::max()-2.0)
        throw std::invalid_argument("Root separation is too small for the sampling domain");
    SurfaceProjector<double> projector(field,options);
    std::vector<Point<double>> roots;
    std::map<std::array<int,3>,std::vector<size_t>> cells;
    for(int i=0;i<samples_per_axis;++i)for(int j=0;j<samples_per_axis;++j)for(int k=0;k<samples_per_axis;++k) {
        Point<double> p=options.domain_min+(options.domain_max-options.domain_min).cwiseProduct(Point<double>(i,j,k)/double(samples_per_axis-1));
        if(!projector.project(p))continue;
        std::array<int,3> cell;
        for(int d=0;d<3;++d)cell[d]=static_cast<int>(std::floor((p[d]-options.domain_min[d])/separation));
        bool duplicate=false;
        for(int x=-1;x<=1&&!duplicate;++x)for(int y=-1;y<=1&&!duplicate;++y)for(int z=-1;z<=1&&!duplicate;++z) {
            auto it=cells.find({cell[0]+x,cell[1]+y,cell[2]+z});
            if(it!=cells.end())for(size_t index:it->second)if((p-roots[index]).norm()<separation){duplicate=true;break;}
        }
        if(!duplicate){cells[cell].push_back(roots.size());roots.push_back(p);}
    }
    return roots;
}
} // namespace marching_triangles
