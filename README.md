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

Also register the LibTorch package using [LibTorch-setting.md](LibTorch-setting.md),
with a LibTorch build suitable for the new machine. The setup below uses Spack
1.x and automatically registers this checkout's `spack-repo` in the environment.
It contains checksum-pinned recipes for TetGen 1.6.0 and earcut.hpp 2.2.4; no
separate checkout or manual header/library installation is needed for them.

-----

## Setting up Spack environment

### First time: create and load the Spack environment

```
git clone https://github.com/MaGuanqun/mfa-TDA.git
cd /path/to/mfa-TDA
source ./create-env.sh
source ./load-env.sh
```

`env.yaml` lists all root dependencies. `create-env.sh` creates `mfa-TDA-env` if
needed, adds missing specs to an existing environment, and installs dependencies.
It preserves the existing environment and installed packages, so the same command
also upgrades an older checkout's environment to include the volume dependencies:

```bash
source ./create-env.sh
source ./load-env.sh
```

On another computer, clone this repository, register MFA and LibTorch as above,
and run these commands again. Spack downloads and verifies the pinned archives
and builds TetGen with that machine's compiler. Keep the checkout available: the
environment records its local `spack-repo` path. The manifest pins these two
volume libraries; the other dependencies are resolved by the local Spack setup.

### Subsequent times: load the Spack environment

```
source /path/to/mfa-TDA/load-env.sh
```

Use `source` in the same terminal where you run CMake. Executing
`./load-env.sh` cannot update the parent shell, and a separate IDE CMake process
may not inherit these variables. Configure once from the activated terminal to
cache the dependency paths in `build`, and use that build directory in the IDE.
`load-env.sh` exports `TETGEN_PATH` and `EARCUT_PATH` and adds both installation
prefixes to `CMAKE_PREFIX_PATH`. Volume targets then use the installed Spack
headers and TetGen library without downloading sources during CMake configure.

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

Critical-curve extraction is available through `critical_curves_explicit`,
`critical_values_mfa`, and `critical_curves_mfa`. See the
[critical-curve workflow](src/isosurface/README.md#critical-curves) for known-value
analytic examples and MFA candidate-value discovery followed by tracing.

See [the isosurface README](src/isosurface/README.md) for parameters, MFA model
compatibility, validation results, and the optional standalone build.
