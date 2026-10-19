# Build instructions for Codex

These instructions apply to the entire repository.

## Activate the environment

Run build commands from the repository root in Bash. Before configuring,
compiling, testing, or running C++ executables, activate the existing Spack
environment:

```bash
source ./load-env.sh
```

The script activates `mfa-TDA-env` and exports the dependency paths for MFA,
Eigen, TBB, and LibTorch, along with compiler, linker, and runtime library paths.
Use these settings for the build.

Source the script in the **same shell** that runs CMake or the executable.
Running `bash ./load-env.sh` or `./load-env.sh` does not activate the calling
shell. Separate Codex shell tool calls do not preserve exported variables, so
source the script again in each call that needs the environment.

Spack must already be initialized with shell support. If `spack` is unavailable
or activation reports missing shell support, source the installed Spack's
`share/spack/setup-env.sh` before sourcing `load-env.sh`.

## Configure and compile

Run this sequence in one shell invocation from the repository root:

```bash
source ./load-env.sh &&
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb &&
cmake --build build -j 3
```

For an incremental build in a fresh shell invocation:

```bash
source ./load-env.sh && cmake --build build -j 3
```

Build outputs belong in `build/`. Check the active `add_subdirectory(...)`
entries in `src/CMakeLists.txt` when determining which targets are available.
For code changes, build the affected targets and run their available tests
with the environment loaded.

## CMake compatibility for new files

New source files and targets must be compatible with the existing root
`CMakeLists.txt` and its subdirectory build configuration. Follow the project's
C++17 standard and existing dependency, include, and linking conventions.
Register new sources or targets in the appropriate `CMakeLists.txt` and ensure
their directories are reachable through the root's `add_subdirectory(...)`
chain. New code must build through the root configure and build commands above
without requiring a separate standalone build. Verify this integration after
adding new compiled code.

## Environment troubleshooting

If CMake cannot find a dependency, check activation and the paths printed by
`load-env.sh` before changing CMake files or hardcoding installation paths.

`create-env.sh` is for environment provisioning. It removes and recreates
`mfa-TDA-env`; do not run it as part of a normal build. Use the existing
environment through `source ./load-env.sh`.
