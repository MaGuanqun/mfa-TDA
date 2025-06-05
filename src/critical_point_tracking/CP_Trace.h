#pragma once

template<typename T>
struct CP_Trace
{
    // size_t span_index;
    std::vector<std::vector<VectorX<T>>> traces;

    
    std::vector<std::array<int,8>> connect_info; //for every trace, store the connect information of two end points. for one end point, store [x,y,z, m]. x: 0,1 represents critical points or other trace. y: span index. z: critical point index in the span or trace index in the span. m: only for trace. 0,1 represents left or right end of the trace. So every trace has 8 int.

    std::vector<size_t> prefix_start; //store the start index of every trace in the prefix sum of the number of points in every trace

};