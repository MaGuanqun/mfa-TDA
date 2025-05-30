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

# Optional: extend your library path if you plan to run compiled executables
if [[ "$OSTYPE" == "darwin"* ]]; then
    export DYLD_LIBRARY_PATH=$HDF5_PATH/lib:$TBB_PATH/lib:$DYLD_LIBRARY_PATH
else
    export LD_LIBRARY_PATH=$HDF5_PATH/lib:$TBB_PATH/lib:$LD_LIBRARY_PATH
fi

# Optional: echo them for debugging
echo "MFA_PATH        = $MFA_PATH"
echo "TBB_PATH        = $TBB_PATH"
echo "EIGEN_PATH      = $EIGEN_PATH"
# give openMP 1 core for now to prevent using all cores for threading
# could set a more reasonable number to distribute cores between mpi + openMP
# export OMP_NUM_THREADS=1
