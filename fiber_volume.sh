#!/usr/bin/env bash
set -eo pipefail

default_model=coordnet        # coordnet, mfa, or analytical
default_stage=all             # all, boundary, or fill

usage() {
    cat <<'EOF'
Usage: ./fiber_volume.sh [MODEL [STAGE [QUERY]]]
  coordnet     Full vortex street space-time domain, query up (default) or down
  analytical   Planar triangular prism, exact volume 0.5
  mfa          Full Hurricane Isabel grid, query 1200 (default) or 1300
  all          Run boundary then fill in separate processes (default)
  boundary     Extract and save the closed triangle mesh
  fill         Read the saved PLY, validate it, and tetrahedralize

FIBER_RESULTS_ROOT overrides build/fiber-volume-results.
FIBER_SKIP_BUILD=1 reuses the existing build.
FIRST_MFA and SECOND_MFA override the two Hurricane model paths.
EOF
}

if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if (( $# > 3 )); then usage >&2; exit 1; fi
stage=${2:-$default_stage}
case "$stage" in all|boundary|fill) ;; *) usage >&2; exit 1;; esac
fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$fiber_repo"
source "$fiber_repo/src/fiber_surface/model_experiments.sh"
case "${1:-$default_model}" in
    coordnet)
        target=fiber_volume_coordnet
        vortex_velocities
        query=${3:-up}
        case "$query" in
            up)   triangle=(--q0 "$edge_a" --q1 "$edge_b" --q2 "$wake_up");;
            down) triangle=(--q0 "$edge_b" --q1 "$edge_a" --q2 "$wake_down");;
            *) usage >&2; exit 1;;
        esac
        dataset=vortex_street_3d_2_variables_full/$query
        options=("${model_args[@]}" "${triangle[@]}" "${mesh_args[@]}")
        # Seed discovery already runs in parallel across independent model queries.
        export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
        ;;
    analytical)
        if (( $# == 3 )); then usage >&2; exit 1; fi
        target=fiber_volume_explicit
        dataset=analytical_prism
        options=(--example prism)
        ;;
    mfa)
        target=fiber_volume_mfa
        isabel_pressure_temperature
        query=${3:-1200}
        case "$query" in 1200|1300) ;; *) usage >&2; exit 1;; esac
        dataset=hurricane_isabel_3d_full/$query
        options=("${model_args[@]}" --q0 "$edge_a" --q1 "$edge_b" --q2 "-$query,-22" "${mesh_args[@]}")
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
