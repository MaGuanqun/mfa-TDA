#!/usr/bin/env bash
set -eo pipefail

usage() {
    cat <<'EOF'
Usage: ./fiber_volume.sh analytical|mfa [all|boundary|fill]
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
if (( $# < 1 || $# > 2 )); then usage >&2; exit 1; fi
stage=${2:-all}
case "$stage" in all|boundary|fill) ;; *) usage >&2; exit 1;; esac
fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$fiber_repo"
case "$1" in
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
