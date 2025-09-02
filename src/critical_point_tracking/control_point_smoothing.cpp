//-----------
//Compute derivative control points of an MFA for an entrie block
#include <mfa/mfa.hpp>

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <diy/master.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/io/block.hpp>

#include "opts.h"

#include "block.hpp"


using namespace std;



#include <Eigen/Dense>
#include <iostream>
#include <cmath>

// Generate 1D Gaussian kernel
template<typename T>
Eigen::VectorXd makeGaussianKernel(T sigma) {
    int radius = std::ceil(3 * sigma);
    int size = 2 * radius + 1;
    Eigen::VectorX<T> kernel(size);
    T sum = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        T value = std::exp(-(i*i)/(2*sigma*sigma));
        kernel(i + radius) = value;
        sum += value;
    }
    kernel /= sum; // normalize
    return kernel;
}


// Convert flat index -> N-dimensional index
inline std::vector<int> unravelIndex(int flatIndex, const VectorXi& shape) {
    int dims = shape.size();
    std::vector<int> idx(dims);
    for (int d = dims - 1; d >= 0; --d) {
        idx[d] = flatIndex % shape[d];
        flatIndex /= shape[d];
    }
    return idx;
}

// Convert N-dimensional index -> flat index
inline int ravelIndex(const std::vector<int>& idx, const VectorXi& shape) {
    int flat = 0;
    int stride = 1;
    for (int d = shape.size() - 1; d >= 0; --d) {
        flat += idx[d] * stride;
        stride *= shape[d];
    }
    return flat;
}


// Convolve along a given dimension
template<typename T>
void convolve1D_along_dim(
    const MatrixX<T>& data,
    MatrixX<T>& out,
    const Eigen::VectorX<T>& kernel,
    const VectorXi& shape,
    int dim)
{
    int radius = (kernel.size() - 1) / 2;
    int totalSize = data.rows();
    int dims = shape.size();

    out.setZero();

    for (int flat = 0; flat < totalSize; ++flat) {
        auto idx = unravelIndex(flat, shape);
        T accum = 0;
        int base = idx[dim];
        for (int k = -radius; k <= radius; ++k) {
            int pos = base + k;
            // clamp
            if (pos < 0) pos = 0;
            if (pos >= shape[dim]) pos = shape[dim] - 1;

            auto nidx = idx;
            nidx[dim] = pos;
            int nFlat = ravelIndex(nidx, shape);
            accum += kernel(k + radius) * data.data()[nFlat];
        }
        out.data()[flat] = accum;
    }
}

template<typename T>
void gaussianBlurND(
    const Eigen::MatrixX<T>& data,
    const VectorXi& shape,
    T sigma, Eigen::MatrixX<T>& out)
{

    std::cout<<"the size of ctrl point mat "<<data.rows()<<" "<<data.cols()<<std::endl;

    out.resize(data.rows(),data.cols());
    auto kernel = makeGaussianKernel(sigma);

    Eigen::MatrixX<T> tmp = data;

    int dims = shape.size()-1; //for space-time data, only smooth spatial dimensions
    for (int dim = 0; dim < dims; ++dim) {
        convolve1D_along_dim(tmp, out, kernel, shape, dim);
        std::cout<<"finsh smoothing dim "<<dim<<std::endl;
        tmp = out; // feed result into next dimension
    }
}



int main(int argc, char** argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    string infile = "approx.mfa";               // diy input file
    string outfile = "smoothed.mfa";

    // default command line arguments
    int  deriv     = 1;                         // which derivative to take (1st, 2nd, ...)
    int  partial   = -1;                        // limit derivatives to one partial in this dimension
    bool help;                                  // show help

    double sigma = 1.0; // smoothing parameter in row units
    // get command line arguments
    opts::Options ops;
    ops >> opts::Option('f', "infile",  infile,  " diy input file name");
    ops >> opts::Option('h', "help",    help,    " show help");
    ops >> opts::Option('o', "outfile", outfile,    " smoothed control points");
    ops >> opts::Option('s', "sigma",   sigma,   " Gaussian sigma in row units");


    if (!ops.parse(argc, argv) || help)
    {
        if (world.rank() == 0)
            std::cout << ops;
        return 1;
    }

    // echo args
    fprintf(stderr, "\n--------- Input arguments ----------\n");
    cerr <<
        "deriv = "    << deriv << endl;
    #ifdef MFA_TBB
        cerr << "threading: TBB" << endl;
    #endif
    #ifdef MFA_KOKKOS
        cerr << "threading: Kokkos" << endl;
    #endif
    #ifdef MFA_SYCL
        cerr << "threading: SYCL" << endl;
    #endif
    #ifdef MFA_SERIAL
        cerr << "threading: serial" << endl;
    #endif
        fprintf(stderr, "-------------------------------------\n\n");

    
    // initialize DIY
    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks to be moved out of core
    diy::Master      master(world,
            -1,
            -1,
            &Block<real_t>::create,
            &Block<real_t>::destroy,
            &storage,
            &Block<real_t>::save,
            &Block<real_t>::load);
    diy::ContiguousAssigner   assigner(world.size(), -1);   // number of blocks set by read_blocks()

    // read MFA model
    diy::io::read_blocks(infile.c_str(), world, assigner, master, &Block<real_t>::load);
    int nblocks = master.size();
    std::cout << nblocks << " blocks read from file "<< infile << "\n";

    //compute control point

    std::vector<std::vector<std::vector<std::vector<real_t>>>> geo_control_point;
    std::vector<std::vector<MatrixX<real_t>>> sci_deriv_control_points;



    Eigen::VectorXd kernel = makeGaussianKernel(sigma);

    Eigen::MatrixXd smoothed;



    master.foreach([&](Block<real_t>* b, const diy::Master::ProxyWithLink& cp)
            {

                const MatrixXd& controlPoints=b->mfa->var(0).tmesh.tensor_prods[0].ctrl_pts;
                const VectorXi& shape = b->mfa->var(0).tmesh.tensor_prods[0].nctrl_pts;

                std::cout<<"num of ctrl point "<< shape.transpose()<<std::endl;

                gaussianBlurND(controlPoints, shape, sigma, smoothed);
                
                std::cout<<"finish smoothing "<<std::endl;

                const auto& tp_const = b->mfa->var(0).tmesh.tensor_prods[0];
                const_cast<std::decay_t<decltype(tp_const)>&>(tp_const).ctrl_pts=smoothed;          

                });


    diy::io::write_blocks(outfile, world, master);
    

}

