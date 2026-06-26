#pragma once

template<typename T>
struct CP_Trace
{
    // size_t span_index;
    std::vector<VectorX<T>> traces;

    
    std::array<int,2> connect_info{-1,-1}; //store the degenerate_point_index
    size_t prefix_start; //store the start index of every trace in the prefix sum of the number of points in every trace
    bool duplicated{false}; //if the trace is duplicated, it will not be written to the output file
    bool is_loop{false}; //if the trace closed back onto its seed (closed feature line / loop)

};

namespace CP_Trace_fuc
{

    void save_edge_type(const std::string& filename, std::vector<int>& edge_type)
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

       for (const auto& value : edge_type) {
            outFile << value << "\n";  // Each value on a new line
        }
        outFile.close();
        std::cout << "Edge values saved to " << filename << std::endl;
    }

        

    // for 4d, add t to y. So that y dim can go to 3y in the end.  
    template<typename T>
    void convert_to_obj(const std::string& filename, std::vector<CP_Trace<T>>& traces, std::vector<VectorX<T>>& degenerate_points,  VectorX<T>& domain_min, VectorX<T>& domain_range, std::vector<int>* critical_point_types=nullptr, std::string edge_type_filename="")
    {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }

        if(domain_min.size()==3)
        {
            for(auto i=0;i<degenerate_points.size();++i)
            {
                outFile << std::setprecision(15) << "v " << degenerate_points[i].data()[0] << " " << degenerate_points[i].data()[1] << " " << degenerate_points[i].data()[2] << "\n";
            }

            for (auto i=traces.begin();i<traces.end();++i)
            {
                if(i->duplicated)
                {
                    continue;
                }
                for(auto j=0;j<i->traces.size();++j)
                {
                    outFile << std::setprecision(15) << "v " << i->traces[j].data()[0] << " " << i->traces[j].data()[1] << " " << i->traces[j].data()[2] << "\n";
                }
                
            }
        }
        else if(domain_min.size()==4)
        {
            // T ratio = 2.0*domain_range(1) / domain_range(domain_range.size()-1);
            for(auto i=0;i<degenerate_points.size();++i)
            {
                outFile << std::setprecision(15) << "v " << degenerate_points[i].data()[0] << " " << degenerate_points[i].data()[1] << " " << degenerate_points[i].data()[2] << " " << degenerate_points[i].data()[3] << "\n";
                // outFile << std::setprecision(15) << "v " << degenerate_points[i].data()[0] << " " << degenerate_points[i].data()[1] << " "  << degenerate_points[i].data()[3] << "\n";
            }


            for (auto i=traces.begin();i<traces.end();++i)
            {
                if(i->duplicated)
                {
                    continue;
                }
                for(auto j=0;j<i->traces.size();++j)
                {
                    // outFile << std::setprecision(15) << "v " << i->traces[j].data()[0] << " " << i->traces[j].data()[1] << " "  << i->traces[j].data()[3] << "\n";
                    outFile << std::setprecision(15) << "v " << i->traces[j].data()[0] << " " << i->traces[j].data()[1] << " " << i->traces[j].data()[2]<<" "<<i->traces[j].data()[3] << "\n";
                }
                
            }

        }

        std::vector<int> edge_type;

        int obj_index=1+degenerate_points.size();
        for (auto i=traces.begin();i<traces.end();++i)
        {
            if(i->duplicated)
            {
                continue;
            }

            if(i->connect_info[0]>-1)
            {
                outFile <<  "l " << i->connect_info[0]+1 << " " << obj_index << "\n";
                if(critical_point_types!=nullptr)
                {
                    edge_type.emplace_back((*critical_point_types)[obj_index-1]);
                }
            }
            // std::cout<<"write start "<<obj_index<< std::endl;
            for(auto j=1;j<i->traces.size();++j)
            {
                // std::cout<<"write trace "<<obj_index<< std::endl;
                outFile <<  "l " << obj_index << " " << obj_index+1 << "\n";

                if(critical_point_types!=nullptr)
                {
                    edge_type.emplace_back((*critical_point_types)[obj_index-1]);
                }

                obj_index++;
            }
            if(i->connect_info[1]>-1)
            {
                outFile <<  "l " << obj_index << " " << i->connect_info[1]+1 << "\n";

                if(critical_point_types!=nullptr)
                {
                    edge_type.emplace_back((*critical_point_types)[obj_index-1]);
                }
            }

            obj_index++;
            
        }
        outFile.close();

        if(critical_point_types!=nullptr)
        {
            save_edge_type(edge_type_filename, edge_type);
        }

    }

}