#pragma once
#include <vector>
#include <Eigen/Dense>
#include "CP_Trace.h"

namespace deduplication 
{

template<typename T>
struct TraceTimeRange
{
    T start_time;
    T end_time;
    int trace_id;
};

//true: close enough
template<typename T>
bool check_compute_distance_from_point_to_segment(const VectorX<T>& P, const VectorX<T>& A, const VectorX<T>& B, const T spatial_step_size)
{
    VectorX<T> AP= (P- A).head(P.size()-1);
    VectorX<T> AB = (B - A).head(P.size()-1);

    VectorX<T> closest_point;
    if(AB.squaredNorm()>1e-20)
    {
        T t = AP.dot(AB) / AB.squaredNorm();
        t = std::clamp(t,0.0,1.0);
        closest_point = A + t * AB;
    }
    else
    {
        closet_point = A;
    }

    T squared_distance = ((P-closest_point).head(P.size()-1)).squaredNorm();
    if(squared_distance < spatial_step_size*spatial_step_size)
    {
        return true;
    }
    return false;
    
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
            // compute distance from the point to the corresponding segment
            const auto& pre=*(it2-1);
            if(check_compute_distance_from_point_to_segment(degenerate_point, *it2,pre,spatial_step_size))
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
    //check if traces pass through the degenerate points
    std::vector<TraceTimeRange<T>> time_ranges(traces.size());
    for(int i=0;i<traces.size();++i)
    {
        time_ranges[i].start_time = traces[i].traces.front()[traces[i].traces.front().size()-1];
        time_ranges[i].end_time = traces[i].traces.back()[traces[i].traces.back().size()-1];
        time_ranges[i].trace_id = i;
    }

    std::sort(time_ranges.begin(), time_ranges.end(), [](const auto& a, const auto& b) {
        return a.start_time < b.start_time;
    });

    std::vector<std::vector<int>> matching_trace_ids(degenerate_points.size());
    std::vector<std::vector<int>> closest_trace_point_index(degenerate_points.size());

    for(auto i=0;i< degenerate_points.size();i++)
    {
        find_matching_traces(degenerate_points[i], time_ranges, matching_trace_ids[i], closest_trace_point_index[i], traces, spatial_step_size);
    }
    //record info of degenerate points into traces
    std::vector<std::vector<int>> matching_degenerate_point_id(traces.size());
    std::vector<std::vector<int>> matching_degenerate_trace_point_id(traces.size()); //record which point in the trace that hte degenerate point is added to


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
    for(auto i=0;i<traces.size();++i)
    {
        if(matching_degenerate_point_id[i].empty())
        {
            continue;
        }

        traces[i].connect_info[1] = matching_degenerate_point_id[i][0];
        for(auto j=0;j<matching_degenerate_point_id[i].size()-1;++j)
        {
            int size = traces.size();
            traces.resize(size+1);
            traces.back().traces.insert(traces.back().traces.end(), traces[i].traces.begin()+matching_degenerate_trace_point_id[i][j], traces[i].traces.begin()+matching_degenerate_trace_point_id[i][j+1]);
            traces.back().connect_info[0] = matching_degenerate_point_id[i][j];
            traces.back().connect_info[1] = matching_degenerate_point_id[i][j+1];
        }
        //add the last part of the trace
        int size = traces.size();
        traces.resize(size+1);
        traces.back().traces.insert(traces.back().traces.end(), traces[i].traces.begin()+matching_degenerate_trace_point_id[i].back(), traces[i].traces.end());
        traces.back().connect_info[0] = matching_degenerate_point_id[i].back();
        traces[i].traces.erase(traces[i].traces.begin()+matching_degenerate_trace_point_id[i][0], traces[i].traces.end());
        
    }

}


template<typename T>
void deduplicate_traces(std::vector<CP_Trace<T>>& traces, T spatial_step_size, T time_step)
{
    //only check if both of the two ends are close enough




}