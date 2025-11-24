#pragma once
#include <vector>
#include <Eigen/Dense>
#include "CP_Trace.h"

namespace error_deduplication 
{

template<typename T>
struct TraceTimeRange
{
    T start_time;
    T end_time;
    int trace_id;
};

template<typename T>
struct Hash
{
    size_t operator()(const Eigen::VectorX<T>& v) const
    {
        std::hash<T> hasher;
        size_t seed = 0;
        int index=(v.size()>>1)+2;
        for (int i =3;i<index;++i) {
            seed ^= hasher(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


template<typename T>
struct Equal {
    bool operator()(const Eigen::VectorX<T>& a, const Eigen::VectorX<T>& b) const {
        int size=(a.size()>>1)-1;
        T temporal_distance = std::abs(a[a.size()-1] - b[b.size()-1]);

        T squaredDistance = (a.segment(size+3,size-1)-b.segment(size+3,size-1)).squaredNorm();

        return (squaredDistance <= a[0] * a[0])&&(temporal_distance <= a[1]);
    }
};

//true: close enough
template<typename T>
bool check_compute_distance_from_P_to_A_B(const VectorX<T>& P, const VectorX<T>& A, const VectorX<T>& B, const T spatial_step_size)
{
    if((P-A).head(P.size()-1).norm() < spatial_step_size)
    {
        return true; //point P is close enough to point A
    }
    if((P-B).head(P.size()-1).norm() < spatial_step_size)
    {
        return true; //point P is close enough to point B
    }

    return false; //point P is not close enough to point A or B

    // VectorX<T> AP= (P- A).head(P.size()-1);
    // VectorX<T> AB = (B - A).head(P.size()-1);

    // VectorX<T> closest_point;
    // if(AB.squaredNorm()>1e-20)
    // {
    //     T t = AP.dot(AB) / AB.squaredNorm();
    //     t = std::clamp(t,0.0,1.0);
    //     closest_point = A + t * AB;
    // }
    // else
    // {
    //     closet_point = A;
    // }

    // T squared_distance = ((P-closest_point).head(P.size()-1)).squaredNorm();
    // if(squared_distance < spatial_step_size*spatial_step_size)
    // {
    //     return true;
    // }
    // return false;
    
}

template<typename T>
bool find_matching_traces(const VectorX<T>& degenerate_point, std::vector<TraceTimeRange<T>>& time_ranges, std::vector<int>& matching_trace_ids, std::vector<int>& closest_trace_point_index, std::vector<CP_Trace<T>>& traces, T spatial_step_size)
{
    T time = degenerate_point[degenerate_point.size()-1];
    auto it=std::upper_bound(time_ranges.begin(), time_ranges.end(), time, [](T t, const TraceTimeRange<T>& a) {
        return t < a.start_time;
    });
    
    for (auto jt = it; jt != time_ranges.begin();) {
        --jt;
        if (jt->end_time > time) {
            //find first index that is larger than t
            auto it2 = std::upper_bound(
            traces[jt->trace_id].traces.begin(), traces[jt->trace_id].traces.end(), time,
            [](T t, const VectorX<T>& p) {
                return t < p[p.size()-1];
            });
           
            if(it2==traces[jt->trace_id].traces.begin() || it2==traces[jt->trace_id].traces.end())
            {
                continue; //continue to the next trace
            }
        
            const auto& pre=*(it2-1);
            if(check_compute_distance_from_P_to_A_B(degenerate_point, *it2,pre,spatial_step_size))
            {
                auto idx = std::distance(traces[jt->trace_id].traces.begin(), it2);//the first point that is larger than t
                closest_trace_point_index.emplace_back(idx);
                matching_trace_ids.emplace_back(jt->trace_id);
            }
            
        }   
    }

    if(matching_trace_ids.empty())
    {
        return false; 
    }
    return true;
}

template<typename T>
void splitting(std::vector<CP_Trace<T>>& traces,std::vector<VectorX<T>>& degenerate_points,T spatial_step_size)
{
    std::cout<<"start splitting traces"<<std::endl;
    // for(auto& degerate_point: degenerate_points)
    // {
    //     if(degerate_point.size()<3)
    //     {
    //         std::cout<<degerate_point<<std::endl;
    //     }
    // }

    std::cout<<"start sorting "<<degenerate_points.size()<<" degenerate points"<<std::endl;
    std::sort(degenerate_points.begin(), degenerate_points.end(), [](const VectorX<T>& a, const VectorX<T>& b) {
        return a[a.size()-1] < b[b.size()-1]; //sort by time
    });

    //check if traces pass through the degenerate points
    std::vector<TraceTimeRange<T>> time_ranges(traces.size());
    for(int i=0;i<traces.size();++i)
    {
        if(traces[i].traces.empty())
        {
            time_ranges[i].start_time = std::numeric_limits<T>::max();
            time_ranges[i].end_time = std::numeric_limits<T>::lowest();
            time_ranges[i].trace_id = i;
            continue;
        }
        time_ranges[i].start_time = traces[i].traces.front()[traces[i].traces.front().size()-1];
        time_ranges[i].end_time = traces[i].traces.back()[traces[i].traces.back().size()-1];
        time_ranges[i].trace_id = i;
    }

    std::sort(time_ranges.begin(), time_ranges.end(), [](const auto& a, const auto& b) {
        return a.start_time < b.start_time;
    });

    std::vector<std::vector<int>> matching_trace_ids(degenerate_points.size());
    std::vector<std::vector<int>> closest_trace_point_index(degenerate_points.size());//the first point index that is larger than t

    for(auto i=0;i< degenerate_points.size();i++)
    {
        find_matching_traces(degenerate_points[i], time_ranges, matching_trace_ids[i], closest_trace_point_index[i], traces, spatial_step_size);
    }


    //record info of degenerate points into traces

    std::vector<std::vector<int>> matching_degenerate_point_id(traces.size());
    std::vector<std::vector<int>> matching_degenerate_trace_point_id(traces.size()); //record which point in the trace that the degenerate point is added before


    for(auto i=0;i<degenerate_points.size();++i)
    {
        if(matching_trace_ids[i].empty())
        {
            continue;
        }
        //add degenerate point to the trace
        for(auto j=0;j<matching_trace_ids[i].size();++j)
        {
            int trace_id = matching_trace_ids[i][j];
            int point_index = closest_trace_point_index[i][j];
            matching_degenerate_point_id[trace_id].emplace_back(i);
            matching_degenerate_trace_point_id[trace_id].emplace_back(point_index);   
        }
    }
    //split traces
    int trace_size = traces.size();
    for(auto i=0;i<trace_size;++i)
    {
        if(matching_degenerate_point_id[i].empty())
        {
            continue;
        }

        traces[i].connect_info[1] = matching_degenerate_point_id[i][0];
        for(auto j=0;j<(matching_degenerate_point_id[i].size()-1);++j)
        {
            int size = traces.size();
            traces.resize(size+1);
            traces.back().traces.insert(traces.back().traces.end(), traces[i].traces.begin()+matching_degenerate_trace_point_id[i][j], traces[i].traces.begin()+matching_degenerate_trace_point_id[i][j+1]);
            // traces.back().connect_info[0] = matching_degenerate_point_id[i][j];
            // traces.back().connect_info[1] = matching_degenerate_point_id[i][j+1];
        }
        //add the last part of the trace
        int size = traces.size();
        traces.resize(size+1);
        traces.back().traces.insert(traces.back().traces.end(), traces[i].traces.begin()+matching_degenerate_trace_point_id[i].back(), traces[i].traces.end());
        // traces.back().connect_info[0] = matching_degenerate_point_id[i].back();
        traces[i].traces.erase(traces[i].traces.begin()+matching_degenerate_trace_point_id[i][0], traces[i].traces.end());
    }

}


template<typename T>
bool close_to_point_A(VectorX<T>& p, T x, T y, T z)
{
    T dx = p[0]-x;
    T dy = p[1]-y;
    T dz = p[2]-z;

    T dist_square = sqrt(dx*dx + dy*dy);

    if(dist_square < 1e-5 && std::abs(dz) < 1e-5)
    {
        return true;
    }
    return false;
}

template<typename T>
//for two traces with same start, test if they are the same. Test mid point and end point
bool two_traces_are_equal(std::vector<CP_Trace<T>>& traces, int trace_id_0, int trace_id_1, T spatial_step_size, T time_step)
{

    int dim = traces[trace_id_0].traces[0].size();
    if(std::abs(traces[trace_id_0].traces.back()[dim-1]-traces[trace_id_1].traces.back()[dim-1]) > time_step)
    {
        // std::cout<<"error "<<std::abs(traces[trace_id_0].traces.back()[dim-1]-traces[trace_id_1].traces.back()[dim-1]) <<std::endl;

        return false;
    }
    if((traces[trace_id_0].traces.back().head(dim-1)-traces[trace_id_1].traces.back().head(dim-1)).squaredNorm() > spatial_step_size*spatial_step_size)
    {
        // std::cout<<"error 2 "<<(traces[trace_id_0].traces.back().head(dim-1)-traces[trace_id_1].traces.back().head(dim-1)).squaredNorm()<<std::endl;
        // std::cout<<trace_id_0<<" "<<trace_id_1<<std::endl;
        // std::cout<< traces[trace_id_0].traces[0].transpose()<<std::endl;
        // std::cout<< traces[trace_id_1].traces[1].transpose()<<std::endl;
        return false;
    }
    auto& mid = traces[trace_id_0].traces[traces[trace_id_0].traces.size()>>1];

    // std::cout<<"mid point test "<<mid.transpose()<<std::endl;
    

    auto it2 = std::upper_bound(
            traces[trace_id_1].traces.begin(), traces[trace_id_1].traces.end(), mid[dim-1],
            [](T t, const VectorX<T>& p) {
                return t < p[p.size()-1];
            });
    
    if(it2==traces[trace_id_1].traces.end())
    {
        if((mid.head(dim-1)-traces[trace_id_1].traces.back().head(dim-1)).squaredNorm() > spatial_step_size*spatial_step_size)
        {
            return false;
        }

        return true;
    }

    if((mid.head(dim-1)-it2->head(dim-1)).squaredNorm() < spatial_step_size*spatial_step_size)
    {
        return true;
    }

    if(it2 != traces[trace_id_1].traces.begin())
    {
        if((mid.head(dim-1)-(it2-1)->head(dim-1)).squaredNorm() < spatial_step_size*spatial_step_size)
        {
            return true;
        }
    }

    // std::cout<< it2-traces[trace_id_1].traces.begin()<< " mid point test failed "<<std::endl;
    
    return false;

}

//the point should be 2n+3 <spatial threshold, temporal threshold, trace_index, current hash index, original position>
template<typename T>
bool registerCell(Eigen::VectorX<T>& Position, T spatial_epsilon, T temporal_epsilon, std::vector<std::vector<size_t>>& temp_index,std::unordered_set<Eigen::VectorX<T>, Hash<T>, Equal<T>>& points_step_1, const VectorX<T>& domain_min,std::vector<CP_Trace<T>>& traces)
{
    size_t k1;

    T k0, t0;
    T cell_size=10.0;
    T left_range=1.0 / cell_size; T right_range=1.0 - left_range;

    int j=0;
    for(int i=((Position.size()>>1)+2);i<Position.size()-1;++i)
    {
        temp_index[j].clear();
        k0=(Position[i]-domain_min[j]) / (cell_size*spatial_epsilon);
        t0=k0-std::floor(k0);
        k1=std::floor(k0);
        temp_index[j].emplace_back(k1); 
        if(t0<left_range)
        {            
            temp_index[j].emplace_back(k1-1);
        }
        else if(t0>right_range)
        {
            temp_index[j].emplace_back(k1+1);
        }
        j++;
    }

    temp_index[j].clear();
    k0=(Position[Position.size()-1]-domain_min[j]) / (cell_size*temporal_epsilon);
    t0=k0-std::floor(k0);
    k1=std::floor(k0);
    temp_index[j].emplace_back(k1);
    if(t0<left_range)
    {            
        temp_index[j].emplace_back(k1-1);
    }
    else if(t0>right_range)
    {
        temp_index[j].emplace_back(k1+1);
    }

   if(Position.size()==7) //dimension = 2
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                Position[3]=i;
                Position[4]=j;

                auto range = points_step_1.equal_range(Position);

                for (auto it = range.first; it != range.second; ++it) {
                    const auto& elem = *it;
                    // elem is equal to Position according to VectorEqual<T>
                    if(two_traces_are_equal(traces, elem[2], Position[2], spatial_epsilon, temporal_epsilon))
                    {
                        return false; // deplicated root
                    }
                }
            }
        }
    }
    else if(Position.size()==9) //dimension = 3
    {
        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    Position[3]=i;
                    Position[4]=j;
                    Position[5]=k;

                    auto range = points_step_1.equal_range(Position);

                    for (auto it = range.first; it != range.second; ++it) {
                        const auto& elem = *it;
                        // elem is equal to Position according to VectorEqual<T>
                        if(two_traces_are_equal(traces, elem[2], Position[2], spatial_epsilon, temporal_epsilon))
                        {
                            //  std::cout<<"test deplicated root "<<elem[2] << " "<<Position[2]<< std::endl;
                            return false; // deplicated root
                        }
                        if(Position[2]==31)
                        {
                            std::cout<<"check trace 31 "<<elem[2]<<std::endl;
                        }
                    }

                    // auto status =  points_step_1.insert(Position);
                    // if(!status.second)
                    // {
                    //     duplicated=true;

                    // }
                }
            }
        }
    }
    else if(Position.size()==11)
    {

        for(auto i:temp_index[0])
        {
            for(auto j:temp_index[1])
            {
                for(auto k:temp_index[2])
                {
                    for(auto l:temp_index[3])
                    {
                        Position[3]=i;
                        Position[4]=j;
                        Position[5]=k;
                        Position[6]=l;

                        auto range = points_step_1.equal_range(Position);

                        for (auto it = range.first; it != range.second; ++it) {
                            const auto& elem = *it;
                            if(two_traces_are_equal(traces, elem[2], Position[2], spatial_epsilon, temporal_epsilon))
                            {
                                return false; // deplicated root
                            }
                        }
                       
                    }
                }
            }
        }
    
    }
    else
    {
        int total_case=1;
        for(auto i:temp_index)
        {
            total_case*=i.size();
        }
        std::vector<int> domain_store_size(temp_index.size());
        domain_store_size.back()=1; //domain_store_size [...,d1*d2*d3,d1*d2,d1,1] dn is 1 or 2
        for(int i=temp_index.size()-2;i>=0;--i)
        {
            domain_store_size[i]=domain_store_size[i+1]*temp_index[i+1].size();
        }

        for(int i=0;i<total_case;++i)
        {
            int temp=i;
            for(int j=0;j<temp_index.size();++j)
            {
                Position[j+2]=temp_index[j][temp/domain_store_size[j]];
                temp=temp%domain_store_size[j];
            }

            auto range = points_step_1.equal_range(Position);

            for (auto it = range.first; it != range.second; ++it) {
                const auto& elem = *it;
                if(two_traces_are_equal(traces, elem[2], Position[2], spatial_epsilon, temporal_epsilon))
                {
                    return false; // deplicated root
                }
            }
        }
    }



    for(int i=0;i<temp_index.size();++i)
    {
        Position[i+3]=temp_index[i][0]; 
    }

    if(Position[2]==31)
    {
        std::cout<<"insert trace 31 "<<std::endl;
        std::cout<<Position.transpose()<<std::endl;
    }

    points_step_1.insert(Position); //insert the point into the hash table
    
    return true;
    
}

template<typename T>
void deduplicate_traces(std::vector<CP_Trace<T>>& traces,std::vector<VectorX<T>>& degenerate_points, T spatial_step_size, T time_step,const VectorX<T>& domain_min)
{
    
    //sort degenerate points by time in splitting()

    if(!degenerate_points.empty())
    {
        splitting(traces, degenerate_points, spatial_step_size);

        // std::cout<<"end splitting "<<std::endl;
    }
    
    std::vector<std::vector<size_t>> temp_index;
    
    temp_index.resize(domain_min.size());

    std::unordered_set<Eigen::VectorX<T>, Hash<T>, Equal<T>> points_step_1; 

    for(int i=0; i<traces.size();++i)
    {
        if(traces[i].traces.empty())
        {
            traces[i].duplicated=true; //skip empty traces
            continue;  
        }

        // if(i!=11 && i!=31 && i!=33 && i!=46)
        // {
        //     traces[i].duplicated=true;
        //     continue;
        // }

        VectorX<T> Pos(2*domain_min.size()+3); //<spatial threshold, temporal threshold, trace_index, current hash index, original position>
        Pos[0]= spatial_step_size; Pos[1]= time_step; Pos[2]=i;
        Pos.tail(domain_min.size()) = traces[i].traces[0];

        if(close_to_point_A(traces[i].traces[0],0.0,0.0,1.02151) || close_to_point_A(traces[i].traces[0],0.0,0.0,1.05) || close_to_point_A(traces[i].traces[0],0.0,0.0,1.0208))
        {
            std::cout<<"find close to point A "<<i<<std::endl;
        }
        if(close_to_point_A(traces[i].traces.back(),0.0,0.0,1.02151) || close_to_point_A(traces[i].traces.back(),0.0,0.0,1.05) || close_to_point_A(traces[i].traces.back(),0.0,0.0,1.0208))
        {
            std::cout<<"find close to point A at the end "<<i<<std::endl;
        }
        //function true means not duplicated
        if(!registerCell(Pos, spatial_step_size, time_step, temp_index, points_step_1, domain_min, traces))
        {
                traces[i].duplicated=true;

            
        }


        if(i == 11 || i == 33 || i ==31 || i==46)
        {
            std::cout<<"trace "<<i<<" is kept "<<std::endl;
            std::cout<< traces[i].duplicated<<std::endl;
        }
    }

}





}