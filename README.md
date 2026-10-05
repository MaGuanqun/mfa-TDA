# mfa-TDA
Workflow example using mfa-TDA

# Instructions for Building and Running

Installation is done through Spack. If you don't have Spack installed or if Spack is new to you, go [here](https://spack.readthedocs.io/en/latest/) first.

-----

## Adding the following Spack repositories to your local Spack installation

MFA
```
git clone https://github.com/tpeterka/mfa.git
spack repo add mfa
```
Check the beginning of `mfa/packages/mfa/package.py`, make sure the file begins with 
```
from spack.package import *
```

-----

## Setting up Spack environment

### First time: create and load the Spack environment

```
git clone https://github.com/MaGuanqun/mfa-TDA.git
cd /path/to/mfa-TDA
source ./create-env.sh     # requires being in the same directory to work properly
source ./load-env.sh
```

### Subsequent times: load the Spack environment

```
source /path/to/mfa-TDA/load-env.sh
```

Use `source` in the same terminal where you run CMake. Executing
`./load-env.sh` cannot update the parent shell, and a separate IDE CMake process
may not inherit these variables. Configure once from the activated terminal to
cache the dependency paths in `build`, and use that build directory in the IDE.

----

## Building mfa-TDA

Run from the repository root after sourcing `load-env.sh`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
cmake --build build -j 3
```

The main build includes `critical_point_tracking`, `feature_flow_fields`, and
`isosurface`. Executables are placed in `build/src/<folder>/`.

To run the isosurface tests and example:

```sh
ctest --test-dir build -R '^isosurface_' --output-on-failure
./isosurface_explicit.sh
```

For extraction directly from an MFA model, run `./isosurface_mfa.sh`. Its default
input is `build/src/boussinesq_3d_iso/boussinesq_3d.mfa` at isovalue 0.5. The shared
mesher handles both MFA and explicit scalar fields.

See [the isosurface README](src/isosurface/README.md) for parameters, MFA model
compatibility, validation results, and the optional standalone build.
