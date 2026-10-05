#pragma once

#include "../block.hpp"
#include <fstream>
#include <memory>
#include <stdexcept>

namespace marching_triangles {
// Use the same current MFA/DIY serialization as the other project tools.
class MfaModel {
public:
    MfaModel(const std::string& filename,const diy::mpi::communicator& world,int block_index=-1)
        : master_(world,1,-1,&Block<double>::create,&Block<double>::destroy,nullptr,&Block<double>::save,&Block<double>::load) {
        if(world.size()!=1)throw std::invalid_argument("Run MFA isosurfacing with one MPI rank (invoke the executable directly)");
        std::ifstream input(filename,std::ios::binary);
        if(!input)throw std::runtime_error("Cannot read MFA model: "+filename);
        input.seekg(0,std::ios::end);
        if(input.tellg()<static_cast<std::streamoff>(sizeof(size_t)))throw std::runtime_error("Truncated MFA file");
        const auto size=input.tellg();size_t footer_size=0;input.seekg(-static_cast<std::streamoff>(sizeof(size_t)),std::ios::end);
        input.read(reinterpret_cast<char*>(&footer_size),sizeof(size_t));
        if(footer_size<sizeof(size_t) || footer_size>static_cast<size_t>(size)-sizeof(size_t))throw std::runtime_error("Invalid DIY MFA file footer");
        diy::ContiguousAssigner assigner(1,-1);
        diy::io::read_blocks(filename,world,assigner,master_,&Block<double>::load);
        if(master_.size()!=1 && block_index<0)throw std::invalid_argument("Multiple MFA blocks: select a block with --block; cross-block stitching is not implemented");
        if(block_index<-1)throw std::invalid_argument("MFA block index must be nonnegative");
        if(block_index<0)block_index=0;
        if(block_index>=master_.size())throw std::invalid_argument("MFA block index out of range");
        block_=master_.block<Block<double>>(block_index);
        std::cout<<"Read "<<master_.size()<<" MFA block(s); selected block "<<block_index<<std::endl;
    }
    const Block<double>& block() const {return *block_;}
private:
    diy::Master master_;
    Block<double>* block_=nullptr;
};
} // namespace marching_triangles
