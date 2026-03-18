#pragma once

#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <stdexcept>
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

    std::cout<<"start time range sorting "<<std::endl;
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

    // std::cout<<"end time range sorting "<<std::endl;

    std::vector<std::vector<int>> matching_trace_ids(degenerate_points.size());
    std::vector<std::vector<int>> closest_trace_point_index(degenerate_points.size());//the first point index that is larger than t

    for(auto i=0;i< degenerate_points.size();i++)
    {
        find_matching_traces(degenerate_points[i], time_ranges, matching_trace_ids[i], closest_trace_point_index[i], traces, spatial_step_size);
    }

    // std::cout<<"end finding matching traces "<<std::endl;
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
// ============================================================================
// 1. Cell key & hash for spatial hashing
// ============================================================================

struct CellKey
{
    // One integer index per dimension (spatial dims + time dim)
    std::vector<long long> idx;

    bool operator==(const CellKey& other) const noexcept
    {
        return idx == other.idx;
    }
};

struct CellKeyHash
{
    size_t operator()(CellKey const& key) const noexcept
    {
        std::hash<long long> hasher;
        size_t seed = 0;
        for (auto v : key.idx)
        {
            seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Map: cell -> list of trace IDs whose *first point* lies in that cell
using CellMap = std::unordered_map<CellKey, std::vector<int>, CellKeyHash>;

// ============================================================================
// 2. You already have this in your old code: we just forward-declare it here.
//    (Keep your existing definition; remove this declaration if you include
//     a header that already declares it.)
// ============================================================================

template<typename T>
bool two_traces_are_equal(std::vector<CP_Trace<T>>& traces,
                          int trace_id_0,
                          int trace_id_1,
                          T spatial_step_size,
                          T time_step)
{

    int dim = traces[trace_id_0].traces[0].size();

    if(std::abs(traces[trace_id_0].traces[0][dim-1]-traces[trace_id_1].traces[0][dim-1]) > time_step)
    {

        return false;
    }



    if((traces[trace_id_0].traces[0].head(dim-1)-traces[trace_id_1].traces[0].head(dim-1)).squaredNorm() > spatial_step_size*spatial_step_size)
    {

        // if((traces[trace_id_0].traces[0].head(dim-1)-traces[trace_id_1].traces[0].head(dim-1)).squaredNorm() < 9*spatial_step_size*spatial_step_size)
        // {
            // T distance = (traces[trace_id_0].traces[0].head(dim-1)-traces[trace_id_1].traces[0].head(dim-1)).norm();
            // if(distance/spatial_step_size < 1.5)
        //     std::cout<<distance<<" "<<spatial_step_size<< " "<<distance/spatial_step_size<<std::endl;
        // }
        return false;
    }



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

// ============================================================================
// 3. Helper: compute candidate cell indices for one coordinate
//    (central cell + possible neighbor cells if near boundary)
// ============================================================================

template <typename T>
inline void build_candidate_indices_1d(
    T coord,              // coordinate value
    T domain_min,         // domain min in this dimension
    T epsilon,            // spatial_step_size (for x/y/z) or time_step (for t)
    T cell_factor,        // e.g. 10.0  -> cell_size = cell_factor * epsilon
    std::vector<long long>& out_indices)
{
    out_indices.clear();

    T cell_size = cell_factor * epsilon;
    if (cell_size <= T(0))
    {
        // Avoid division by zero; treat as single "cell 0"
        out_indices.push_back(0);
        return;
    }

    T k = (coord - domain_min) / cell_size;

    auto k_floor = static_cast<long long>(std::floor(k));
    T frac = k - static_cast<T>(k_floor);

    const T left_range  = T(1.0) / cell_factor;
    const T right_range = T(1.0) - left_range;

    // Central cell
    out_indices.push_back(k_floor);

    // Neighbor cells if near cell boundary
    if (frac < left_range)
    {
        out_indices.push_back(k_floor - 1);
    }
    else if (frac > right_range)
    {
        out_indices.push_back(k_floor + 1);
    }
}

// ============================================================================
// 4. Core: register one trace by its first point
//    - Searches nearby cells for a first point close enough
//    - For each candidate trace, calls two_traces_are_equal()
//    - If any equal: return false (duplicate)
//    - Else: insert this trace into the hash and return true
// ============================================================================

template <typename T>
bool register_trace_by_first_point(
    int trace_id,
    std::vector<CP_Trace<T>>& traces,            // non-const for your existing two_traces_are_equal
    CellMap& cell_map,                           // stateful structure reused across calls
    const Eigen::VectorX<T>& domain_min,         // size = dim (spatial + time)
    T spatial_step_size,
    T time_step,
    T cell_factor = T(10.0))                     // controls cell size: cell_size = cell_factor * epsilon
{
    // Skip empty traces
    if (traces[trace_id].traces.empty())
        return true; // nothing to compare; you can alternatively mark duplicated=true

    const auto& first_point = traces[trace_id].traces.front();
    const int dim = static_cast<int>(first_point.size()); // number of coordinates: space_dim + 1 (time)

    if (domain_min.size() != dim)
    {
        throw std::runtime_error("domain_min.size() != first_point.size() in register_trace_by_first_point");
    }

    // ------------------------------------------------------------------------
    // 4.1 Build candidate indices per dimension (including time)
    // ------------------------------------------------------------------------
    std::vector<std::vector<long long>> candidates(dim);
    for (int d = 0; d < dim; ++d)
    {
        const bool is_time_dim = (d == dim - 1);
        T eps = is_time_dim ? time_step : spatial_step_size;

        build_candidate_indices_1d(
            first_point[d],
            domain_min[d],
            eps,
            cell_factor,
            candidates[d]);
    }

    // ------------------------------------------------------------------------
    // 4.2 Iterate over the cartesian product of candidate indices
    //     (= neighbor cells) and look for equal traces
    // ------------------------------------------------------------------------
    std::vector<int> sizes(dim);
    int total_combinations = 1;
    for (int d = 0; d < dim; ++d)
    {
        sizes[d] = static_cast<int>(candidates[d].size());
        total_combinations *= sizes[d];
    }

    for (int combo = 0; combo < total_combinations; ++combo)
    {
        int tmp = combo;

        CellKey key;
        key.idx.resize(dim);

        // Decode combo index into indices per dimension
        for (int d = 0; d < dim; ++d)
        {
            int idx_in_dim = tmp % sizes[d];
            tmp /= sizes[d];

            key.idx[d] = candidates[d][idx_in_dim];
        }

        auto it = cell_map.find(key);
        if (it == cell_map.end())
            continue;

        // For all traces in this cell, test equality by full trace comparison
        for (int other_trace_id : it->second)
        {
            if (two_traces_are_equal(traces,
                                     other_trace_id,
                                     trace_id,
                                     spatial_step_size,
                                     time_step))
            {
                // Found an equal trace -> treat this trace as duplicate
                return false;
            }
        }
    }

    // ------------------------------------------------------------------------
    // 4.3 No equal trace found, insert this trace into its *central* cell
    // ------------------------------------------------------------------------
    CellKey central_key;
    central_key.idx.resize(dim);

    for (int d = 0; d < dim; ++d)
    {
        const bool is_time_dim = (d == dim - 1);
        T eps = is_time_dim ? time_step : spatial_step_size;
        T cell_size = cell_factor * eps;
        if (cell_size <= T(0))
        {
            central_key.idx[d] = 0;
        }
        else
        {
            T k = (first_point[d] - domain_min[d]) / cell_size;
            auto k_floor = static_cast<long long>(std::floor(k));
            central_key.idx[d] = k_floor;
        }
    }

    cell_map[central_key].push_back(trace_id);
    return true;
}

// ============================================================================
// 5. Top-level dedup function: loop over all traces
//    - For each trace, try to register it by first point
//    - If register_trace_by_first_point() returns false, mark duplicated=true
// ==
// 
// 
// ==========================================================================


template<typename T>
bool is_point_on_boundary(const VectorX<T>& point, const VectorX<T>& core_mins, const VectorX<T>& core_maxs)
{
    for(int i=0;i<point.size();++i)
    {
        if(abs(point[i]-core_mins[i]) < 1e-8 || abs(point[i]-core_maxs[i]) < 1e-8)
        {
            return true;
        }
    }
    return false;
}

template<typename T>
void check_point_on_boundary(std::vector<CP_Trace<T>>& traces, const VectorX<T>& core_mins, const VectorX<T>& core_maxs, std::vector<int>& has_boundary)
{
    has_boundary.resize(traces.size(),0);
    for(int i=0;i<traces.size();++i)
    {
        if(traces[i].traces.empty())
        {
            // has_boundary[i]=0;
            continue;
        }
        if (is_point_on_boundary(traces[i].traces[0], core_mins, core_maxs))
        {
            has_boundary[i]=1;
            continue;
        }
        if (is_point_on_boundary(traces[i].traces.back(), core_mins, core_maxs))
        {
            has_boundary[i]=1;
            continue;
        }
    }
}


//put all traces that contain boundary to the front, and then put all traces that do not contain boundary to the back
template<typename T>
void sort_all_traces(std::vector<CP_Trace<T>>& traces, std::vector<int>& has_boundary, std::vector<size_t>& ordered_trace_id)
{
    ordered_trace_id.reserve(traces.size());
    for(auto i=0;i<traces.size();++i)
    {
        if(has_boundary[i]==1)
        {
            ordered_trace_id.emplace_back(i);
        }
    }
    for(auto i=0;i<traces.size();++i)
    {
        if(has_boundary[i]==0)
        {
            ordered_trace_id.emplace_back(i);
        }
    }
}


template <typename T>
void deduplicate_traces(
    std::vector<CP_Trace<T>>& traces,
    std::vector<VectorX<T>>& degenerate_points, 
    T spatial_step_size,
    T time_step,
    const Eigen::VectorX<T>& domain_min,
    const Eigen::VectorX<T>& domain_max)
{

    if(!degenerate_points.empty())
    {
        splitting(traces, degenerate_points, spatial_step_size);

        std::cout<<"end splitting "<<std::endl;
    }

    std::vector<int> has_boundary;
    check_point_on_boundary(traces, domain_min, domain_max, has_boundary);

    // std::cout<<"end checking point on boundary "<<std::endl;
    std::vector<size_t> ordered_trace_id;
    sort_all_traces(traces, has_boundary, ordered_trace_id);


    // std::cout<<"end sorting traces "<<std::endl;
    CellMap cell_map;

    for(auto k=0; k<ordered_trace_id.size();++k)
    // for (int i = 0; i < static_cast<int>(traces.size()); ++i)
    {
        size_t i = ordered_trace_id[k];
        if (traces[i].traces.empty())
        {
            // You can decide whether empty traces are considered duplicated.
            traces[i].duplicated = true;
            continue;
        }

        bool kept = register_trace_by_first_point(
            i,
            traces,
            cell_map,
            domain_min,
            spatial_step_size,
            time_step);

        if (!kept)
        {
            traces[i].duplicated = true;
        }
        else
        {
            traces[i].duplicated = false;
        }
    }
}

} // namespace deduplication
