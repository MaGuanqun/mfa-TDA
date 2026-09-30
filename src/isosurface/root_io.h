#pragma once
#include "surface_field.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace marching_triangles {
// Compatible with utility::writeMatrixVector's one column-major double matrix.
// Also accept text xyz rows (and OBJ v rows) for small reproducible examples.
inline std::vector<Point<double>> read_roots(const std::string& filename) {
    std::ifstream in(filename,std::ios::binary);
    if(!in) throw std::runtime_error("Cannot read root file: "+filename);
    int count=0;in.read(reinterpret_cast<char*>(&count),sizeof(count));
    in.clear();in.seekg(0);
    std::vector<Point<double>> roots;
    if(count==1) {
        size_t rows=0,cols=0;
        in.read(reinterpret_cast<char*>(&count),sizeof(count));
        in.read(reinterpret_cast<char*>(&rows),sizeof(rows));
        in.read(reinterpret_cast<char*>(&cols),sizeof(cols));
        if(!in||cols!=3||rows>10000000) throw std::runtime_error("Invalid 3D root matrix header: "+filename);
        const auto offset=in.tellg();in.seekg(0,std::ios::end);
        if(in.tellg()-offset!=static_cast<std::streamoff>(rows*cols*sizeof(double))) throw std::runtime_error("Truncated or unexpected root matrix data: "+filename);
        in.seekg(offset);roots.resize(rows);
        for(size_t c=0;c<3;++c) for(auto& p:roots) in.read(reinterpret_cast<char*>(&p[c]),sizeof(double));
    } else {
        std::string line;
        while(std::getline(in,line)) {
            const auto comment=line.find('#');if(comment!=std::string::npos)line.resize(comment);
            std::istringstream row(line);row>>std::ws;if(row.eof())continue;
            if(row.peek()=='v')row.get();
            Point<double> p;std::string extra;
            if(!(row>>p[0]>>p[1]>>p[2])||(row>>extra)) throw std::runtime_error("Expected three coordinates per root: "+filename);
            roots.push_back(p);
        }
    }
    for(const auto& p:roots) if(!p.allFinite())throw std::runtime_error("Non-finite root coordinate: "+filename);
    if(roots.empty())throw std::runtime_error("Root file is empty: "+filename);
    return roots;
}
inline void write_roots(const std::string& filename,const std::vector<Point<double>>& roots) {
    std::ofstream out(filename,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot write root file: "+filename);
    const int count=1;const size_t rows=roots.size(),cols=3;
    out.write(reinterpret_cast<const char*>(&count),sizeof(count));
    out.write(reinterpret_cast<const char*>(&rows),sizeof(rows));
    out.write(reinterpret_cast<const char*>(&cols),sizeof(cols));
    for(int c=0;c<3;++c)for(const auto& p:roots)out.write(reinterpret_cast<const char*>(&p[c]),sizeof(double));
    out.close();if(!out)throw std::runtime_error("Failed writing root file: "+filename);
}
} // namespace marching_triangles
