#!/bin/bash

# Exports must affect the calling shell so CMake can find Spack packages.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf 'Load this environment with: source ./load-env.sh\n' >&2
    exit 1
fi

# activate the environment
export SPACKENV=mfa-TDA-env
spack env deactivate > /dev/null 2>&1 || true
spack env activate "$SPACKENV" || return 1
echo "activated spack environment $SPACKENV"


echo "setting flags for building mfa-TDA"
if ! MFA_PATH=$(spack location -i mfa) ||
   ! TBB_PATH=$(spack location -i tbb) ||
   ! EIGEN_PATH=$(spack location -i eigen) ||
   ! TORCH_PATH=$(spack location -i libtorch); then
    printf 'Could not locate the build dependencies in %s. Run: source ./create-env.sh\n' "$SPACKENV" >&2
    return 1
fi
export MFA_PATH TBB_PATH EIGEN_PATH TORCH_PATH

if ! TETGEN_PATH=$(spack location -i 'mfa_tda.tetgen@1.6.0') ||
   ! EARCUT_PATH=$(spack location -i 'mfa_tda.earcut-hpp@2.2.4'); then
    printf 'Volume dependencies are missing. Run: source ./create-env.sh\n' >&2
    return 1
fi
export TETGEN_PATH EARCUT_PATH


export Torch_DIR="$TORCH_PATH/share/cmake/Torch"
export TBB_DIR="$TBB_PATH/lib/cmake/TBB"
export TBB_ROOT="$TBB_PATH"

export CMAKE_PREFIX_PATH="$TETGEN_PATH:$EARCUT_PATH:$TORCH_PATH${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"

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
echo "TETGEN_PATH     = $TETGEN_PATH"
echo "EARCUT_PATH     = $EARCUT_PATH"
echo "Torch_DIR = $Torch_DIR"
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "CUDA_PATH       = $CUDA_PATH"
fi
# give openMP 1 core for now to prevent using all cores for threading
# could set a more reasonable number to distribute cores between mpi + openMP
# export OMP_NUM_THREADS=1
