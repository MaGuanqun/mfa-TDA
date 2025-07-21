#pragma once
#include <vector>
#include <Eigen/Dense>

namespace tracking_utility
{
    template<typename T>
    void convert_to_obj(const std::string& filename, std::vector<Eigen::VectorX<T>>& points)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        for(auto i=points.begin();i<points.end();++i)
        {

            outFile << std::setprecision(15) << "v " << (*i).data()[0] << " " << (*i).data()[1] << " " << (*i).data()[2] << "\n";
            
        }
        outFile.close();

    }

}