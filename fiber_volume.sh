#!/usr/bin/env bash
set -eo pipefail

# Edit these settings, then run ./fiber_volume.sh with no arguments.
default_model=coordnet        # coordnet, mfa, or analytical
default_stage=all             # all, boundary, or fill
coordnet_dataset=vortex_street_3d_2_variables
model_dir="CoordNet/Exp/$coordnet_dataset"
first_coordnet="$model_dir/u-super-spatial-temporal-64-5-float64.pt"
second_coordnet="$model_dir/v-super-spatial-temporal-64-5-float64.pt"
domain_min=-0.5,-0.5,13.5     # Full training box in physical (x,y,time).
domain_max=7.5,0.5,15
box_min=0.38,0.133,13.988     # Local extraction box in the same coordinates.
box_max=0.42,0.167,14.012
q0=0.402,-0.344               # Counterclockwise triangle in normalized (u,v).
q1=0.429,-0.320
q2=0.402,-0.296
coordnet_mesh=(--step 0.006 --fiber-step 0.004 --spatial-tolerance 0.00045
               --residual-tolerance 0.0002 --seed-grid 3 --max-triangles 4000)

usage() {
    cat <<'EOF'
Usage: ./fiber_volume.sh [coordnet|analytical|mfa [all|boundary|fill]]
  coordnet     Vortex street u and v TorchScript models, local (x,y,time) box
  analytical   Planar triangular prism, exact volume 0.5
  mfa          Hurricane Isabel pressure and temperature, local spatial box
  all          Run boundary then fill in separate processes (default)
  boundary     Extract and save the closed triangle mesh
  fill         Read the saved PLY, validate it, and tetrahedralize

FIBER_RESULTS_ROOT overrides build/fiber-volume-results.
FIBER_SKIP_BUILD=1 reuses the existing build.
FIRST_MFA and SECOND_MFA override the two Hurricane model paths.
EOF
}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if (( $# > 2 )); then usage >&2; exit 1; fi
stage=${2:-$default_stage}
case "$stage" in all|boundary|fill) ;; *) usage >&2; exit 1;; esac
fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$fiber_repo"
case "${1:-$default_model}" in
    coordnet)
        target=fiber_volume_coordnet
        dataset=$coordnet_dataset
        options=(--first "$first_coordnet" --second "$second_coordnet"
                 --domain-min "$domain_min" --domain-max "$domain_max"
                 --min "$box_min" --max "$box_max" --q0 "$q0" --q1 "$q1" --q2 "$q2"
                 "${coordnet_mesh[@]}")
        # Seed discovery already runs in parallel across independent model queries.
        export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
        ;;
    analytical)
        target=fiber_volume_explicit
        dataset=analytical_prism
        options=(--example prism)
        ;;
    mfa)
        target=fiber_volume_mfa
        dataset=hurricane_isabel_3d
        first=${FIRST_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/Pf30_3d.mfa}
        second=${SECOND_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/TCf30_3d.mfa}
        options=(--first "$first" --second "$second" --variable1 0 --variable2 0
                 --q0 -1100,-40 --q1 10000,-40 --q2 -1100,0
                 --min 240,240,45 --max 260,260,55
                 --step 1 --fiber-step 0.6 --spatial-tolerance 0.08
                 --residual-tolerance 0.005 --seed-grid 7 --max-triangles 8000)
        ;;
    *) usage >&2; exit 1;;
esac

if [[ $(type -t spack || true) != function ]]; then
    spack_setup=${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}
    source "$spack_setup"
fi
source ./load-env.sh || exit 1
if [[ ${FIBER_SKIP_BUILD:-0} != 1 ]]; then
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
    cmake --build build --target "$target" -j 3
fi
prefix="${FIBER_RESULTS_ROOT:-build/fiber-volume-results}/$dataset/volume"
command=("build/src/fiber_surface/$target" "${options[@]}" --output "$prefix")
if [[ $stage == all || $stage == boundary ]]; then
    "${command[@]}" --stage boundary
fi
if [[ $stage == all || $stage == fill ]]; then
    "${command[@]}" --stage fill --boundary-input "${prefix}_boundary.ply"
fi
