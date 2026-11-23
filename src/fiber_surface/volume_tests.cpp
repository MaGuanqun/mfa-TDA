#include "volume_examples.h"
#include <filesystem>

using namespace fiber_surface;
static void require(bool ok,const std::string& why) {if(!ok)throw std::runtime_error(why);}
int main(int argc,char** argv) {
    try {
        const std::string name=argc>1?argv[1]:"prism";
        auto e=volume_example(name=="validation"||name=="box"?"prism":name); Options o;
        if(name=="curved-prism" || name=="sphere")o.residual_tolerance=.005;
        if(name=="hollow-cylinder") {o.step=.25;o.min_step=o.step/8;o.spatial_tolerance=.075*o.step;}
        if(name=="box")e.triangle=RangeTriangle(Range(-10,-10),Range(20,-10),Range(-10,20));
        auto r=extract_volume(e.field,e.triangle,o);
        write_volume((std::filesystem::current_path()/("volume_test_"+name)/"volume").string(),r,e.field,e.triangle,o);
        require(r.complete(),"Incomplete "+name+"; inspect volume_test_"+name+"/volume_report.json");
        const double expected=name=="box"?9:e.exact_volume;
        const bool hollow=name=="hollow-cylinder" || name=="hollow-sphere" || name=="sphere";
        require(std::abs(r.report.volume-expected)<((hollow || name=="disconnected")?.01*expected:(name=="curved-prism"?.002:1e-9)),"Incorrect reference volume");
        const bool cavity=name=="hollow-sphere" || name=="sphere";
        require(r.boundary_report.holes.size()==(cavity?1:0) &&
                r.boundary_report.signed_volumes.size()==(cavity || name=="disconnected"?2:1),"Unexpected shells/cavities");
        if(name=="validation") {
            bool rejected=false;
            try {RangeTriangle invalid(Range(0,0),Range(0,1),Range(1,0));}catch(const std::invalid_argument&){rejected=true;}
            require(rejected,"Clockwise triangle accepted");
            auto open=r.boundary;open.mesh.triangles.pop_back();open.patches.pop_back();
            require(!validate_volume_boundary(open,e.field,e.triangle,o).valid,"Open boundary accepted");
            auto missing=r.mesh;missing.cells.pop_back();
            require(!assess_volume(missing,r.boundary,r.boundary_report,e.field,e.triangle,o).valid,"Missing cell accepted");
            auto duplicate=r.mesh;duplicate.cells.push_back(duplicate.cells.front());
            require(!assess_volume(duplicate,r.boundary,r.boundary_report,e.field,e.triangle,o).valid,"Duplicate cell accepted");
            auto mislabeled=r.mesh;mislabeled.parent_facets[0]=(mislabeled.parent_facets[0]+17)%r.boundary.mesh.triangles.size();
            require(!assess_volume(mislabeled,r.boundary,r.boundary_report,e.field,e.triangle,o).valid,"Incorrect boundary marker accepted");
        }
        std::cout<<name<<": "<<r.mesh.cells.size()<<" tetrahedra, volume "<<r.report.volume<<" (expected "<<expected<<")\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
