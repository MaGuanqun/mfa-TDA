#pragma once
#include <vector>
#include <Eigen/Dense>

namespace tracking_utility
{
    template<typename T>
    void convert_to_obj(const std::string& filename, std::vector<std::vector<Eigen::VectorX<T>>>& points)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        for(auto i=points.begin();i<points.end();++i)
        {
            for(auto j=0;j<i->size();++j)
            {
                outFile << std::setprecision(15) << "v " << (*i)[j].data()[0] << " " << (*i)[j].data()[1] << " " << (*i)[j].data()[2] << "\n";
            }
        }
        outFile.close();

    }

}