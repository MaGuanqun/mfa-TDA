#!/bin/bash

# activate the environment
export SPACKENV=mfa-TDA-env
spack env deactivate > /dev/null 2>&1
spack env activate $SPACKENV
echo "activated spack environment $SPACKENV"


echo "setting flags for building mfa-TDA"
export MFA_PATH=$(spack location -i mfa)
export TBB_PATH=$(spack location -i tbb)
export EIGEN_PATH=$(spack location -i eigen)
export TORCH_PATH=$(spack location -i libtorch)


export Torch_DIR="$TORCH_PATH/share/cmake/Torch"
export TBB_DIR="$TBB_PATH/lib/cmake/TBB"
export TBB_ROOT="$TBB_PATH"

export CMAKE_PREFIX_PATH="$TORCH_PATH${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"

# Optional: extend your library path if you plan to run compiled executables
# if [[ "$OSTYPE" == "darwin"* ]]; then
#     export DYLD_LIBRARY_PATH=$HDF5_PATH/lib:$TBB_PATH/lib:$DYLD_LIBRARY_PATH
# else
#     export LD_LIBRARY_PATH=$HDF5_PATH/lib:$TBB_PATH/lib:$LD_LIBRARY_PATH
# fi

export CXXFLAGS="-I${MFA_PATH}/include -I${EIGEN_PATH}/include"
export LDFLAGS="-L${TORCH_PATH}/lib -L${TBB_PATH}/lib"

if [[ "$OSTYPE" == "darwin"* ]]; then
    export DYLD_LIBRARY_PATH=$HDF5_PATH/$TORCH_PATH/lib:$TBB_PATH/lib:$DYLD_LIBRARY_PATH
else
    export CUDA_PATH=/usr/local/cuda-12.6
    export LD_LIBRARY_PATH="$TORCH_PATH/lib:$TBB_PATH/lib:$CUDA_PATH/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    # export LD_LIBRARY_PATH=$HDF5_PATH/$TORCH_PATH/lib:$CUDA_PATH/lib:$TBB_PATH/lib:$LD_LIBRARY_PATH
fi

# Optional: echo them for debugging
echo "MFA_PATH        = $MFA_PATH"
echo "TBB_PATH        = $TBB_PATH"
echo "EIGEN_PATH      = $EIGEN_PATH"
echo "TORCH_PATH      = $TORCH_PATH"
echo "Torch_DIR = $Torch_DIR"
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "CUDA_PATH       = $CUDA_PATH"
fi
# give openMP 1 core for now to prevent using all cores for threading
# could set a more reasonable number to distribute cores between mpi + openMP
# export OMP_NUM_THREADS=1
