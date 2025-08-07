#pragma once

template<typename T>
struct CP_Trace
{
    // size_t span_index;
    std::vector<VectorX<T>> traces;

    
    std::array<int,2> connect_info{-1,-1}; //store the degenerate_point_index
    size_t prefix_start; //store the start index of every trace in the prefix sum of the number of points in every trace

};

namespace CP_Trace_fuc
{
    template<typename T>
    void convert_to_obj(const std::string& filename, std::vector<CP_Trace<T>>& traces)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }



        for (auto i=traces.begin();i<traces.end();++i)
        {
            for(auto j=0;j<i->traces.size();++j)
            {
                    outFile << std::setprecision(15) << "v " << i->traces[j].data()[0] << " " << i->traces[j].data()[1] << " " << i->traces[j].data()[2] << "\n";
            }
            
        }


        int obj_index=1;
        for (auto i=traces.begin();i<traces.end();++i)
        {
            // std::cout<<"write start "<<obj_index<< std::endl;
            for(auto j=1;j<i->traces.size();++j)
            {
                // std::cout<<"write trace "<<obj_index<< std::endl;
                outFile <<  "l " << obj_index << " " << obj_index+1 << "\n";
                obj_index++;
                
            }
            if(!i->traces.empty())
            {
                obj_index++;
            }


        }
        outFile.close();


    }

}