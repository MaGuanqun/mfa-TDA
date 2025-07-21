#pragma once

template<typename T>
struct CP_Trace
{
    // size_t span_index;
    std::vector<VectorX<T>> traces;

    
    std::array<int,8> connect_info; //for every trace, store the connect information of two end points. for one end point, store [x,y,z, m]. x: 0,1 represents critical points or other trace. y: span index. z: critical point index in the span or trace index in the span. m: only for trace. 0,1 represents left or right end of the trace. So every trace has 8 int.

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
            for(auto j=1;j<i->traces.size();++j)
            {
                outFile <<  "l " << obj_index << " " << obj_index+1 << "\n";
                obj_index++;
                
            }
            obj_index++;

        }
        outFile.close();


    }

}