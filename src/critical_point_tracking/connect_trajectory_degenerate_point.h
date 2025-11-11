#pragma once
#include <vector>
#include <Eigen/Dense>
#include "CP_Trace.h"

namespace connect_trajectory_degenerate_point
{

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

template<typename T>
void registerPoints(Eigen::VectorX<T>& Position, T spatial_epsilon, T temporal_epsilon, std::unordered_set<Eigen::VectorX<T>, Hash<T>, Equal<T>>& points_step_1, const VectorX<T>& domain_min)
{

    size_t k1;

    T k0, t0;
    T cell_size=10.0;
    T left_range=1.0 / cell_size; T right_range=1.0 - left_range;

    std::vector<std::vector<size_t>> temp_index(domain_min.size());

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


    if(Position.size()==9) //dimension = 3
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

                    points_step_1.insert(Position);
                }
            }
        }
    }
    else 
    {
        std::cout<<"error, only support 3D data"<<std::endl;
    }


}

template<typename T>
int check_distance(Eigen::VectorX<T>& Position, T spatial_epsilon, T temporal_epsilon, std::unordered_set<Eigen::VectorX<T>, Hash<T>, Equal<T>>& points_step_1, const VectorX<T>& domain_min, int start_point=1) // end_point=-1
{
    size_t k1;
    T k0, t0;
    T cell_size=10.0;
    T left_range=1.0 / cell_size; T right_range=1.0 - left_range;

    std::vector<size_t> temp_index(domain_min.size());
    int j=0;
    for(int i=((Position.size()>>1)+2);i<Position.size()-1;++i)
    {
        k0=(Position[i]-domain_min[j]) / (cell_size*spatial_epsilon);
        k1=std::floor(k0);
        temp_index[j]=k1; 
        j++;
    }

    k0=(Position[Position.size()-1]-domain_min[j]) / (cell_size*temporal_epsilon);
    k1=std::floor(k0);
    temp_index[j]=k1;

    for(int i=0;i<domain_min.size();++i)
    {
        Position[i+3]=temp_index[i]; 
    }

    auto range = points_step_1.equal_range(Position);
    T distance =std::numeric_limits<T>::max();
    int degenerate_point_index=-1;
    for (auto it = range.first; it != range.second; ++it) {
        const auto& elem = *it;
        if(start_point*elem[elem.size()-1]<start_point*Position[Position.size()-1])
        {
            T temp_distance = (elem.segment(3+domain_min.size(),domain_min.size()-1)-Position.segment(3+domain_min.size(),domain_min.size()-1)).norm();
            if(temp_distance<distance)
            {
                distance=temp_distance;
                degenerate_point_index=elem[2];
            }
        }
    }

    return degenerate_point_index;

}

template<typename T>
void connect_trajectory(std::vector<CP_Trace<T>>& traces,std::vector<VectorX<T>>& degenerate_points, T spatial_step_size, T time_step,const VectorX<T>& domain_min)
{
    // register all degenerate points into the hash table
    std::unordered_set<Eigen::VectorX<T>, Hash<T>, Equal<T>> points_step_1; 
    VectorX<T> Pos(2*domain_min.size()+3); //<spatial threshold, temporal threshold, point_index, current hash index, original position>
    for(int i=0;i<degenerate_points.size();++i)
    {   
        Pos[0]= spatial_step_size*1.732; Pos[1]= time_step*1.732; Pos[2]=i;
        Pos.tail(domain_min.size()) = degenerate_points[i];
        registerPoints(Pos, spatial_step_size, time_step, points_step_1, domain_min);
    }

    //connect traces through degenerate points
    for(int i=0;i<traces.size();++i)
    {
        VectorX<T> Position(2*domain_min.size()+3); //<spatial threshold, 
        Position[0]= spatial_step_size*1.732; Position[1]= time_step*1.732; Position[2]=i;

        if(traces[i].connect_info[0]==-1)
        {
            Position.tail(domain_min.size()) = traces[i].traces[0];
            //find_nearest degenerate point
            int index = check_distance(Position, spatial_step_size, time_step, points_step_1, domain_min, 1);
            if(index!=-1)
            {
                traces[i].connect_info[0]=index;
            }
        }
        if(traces[i].connect_info[1]==-1)
        {
            Position.tail(domain_min.size()) = traces[i].traces.back();
            //find_nearest degenerate point
            int index = check_distance(Position, spatial_step_size, time_step, points_step_1, domain_min, -1);
            if(index!=-1)
            {
                traces[i].connect_info[1]=index;
            }
        }
    }

}


}