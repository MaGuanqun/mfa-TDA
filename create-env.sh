#!/bin/bash

export SPACKENV=mfa-TDA-env
export YAML=$PWD/env.yaml

# create spack environment
echo "creating spack environment $SPACKENV"
spack env deactivate > /dev/null 2>&1
spack env remove -y $SPACKENV > /dev/null 2>&1
spack env create $SPACKENV $YAML

# activate environment
echo "activating spack environment"
spack env activate $SPACKENV

spack add mpich@4
spack add tbb
spack add eigen
spack add mfa~examples~tests thread=tbb


# Optional: add other dependencies you may need
# spack add zlib
# spack add fmt
# spack add diy

# install
echo "installing dependencies in environment"
spack install mfa~examples~tests thread=tbb   # install separately so that MFA_PATH is set for later packages
export MFA_PATH=`spack location -i mfa`


# Install everything else
echo "Installing remaining dependencies"
spack install


spack env deactivate

