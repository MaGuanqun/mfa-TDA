#!/usr/bin/env bash
set -eo pipefail

usage() {
    cat <<'EOF'
Usage: run_vortex_street.sh [up|down]
Extract and audit the full-domain CoordNet triangle surface (default: up).
VORTEX_STEP sets the step (default 0.005) and scales the spatial and surface range tolerances.
VORTEX_SEED_GRID sets the shortest-axis seed count (default 40).
VORTEX_MAX_TRIANGLES sets the per-side budget (default 500000).
FIBER_RESULTS_ROOT sets the output parent (default build/fiber-volume-results).
EOF
}
if (( $# > 1 )); then usage >&2; exit 1; fi
query=${1:-up}
case "$query" in
    up|down) ;;
    --help|-h) usage; exit 0 ;;
    *) usage >&2; exit 1 ;;
esac

experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
fiber_repo=$(cd -- "$experiment_dir/../../.." && pwd)
cd "$fiber_repo"
source "$fiber_repo/src/fiber_surface/model_experiments.sh"
source "$experiment_dir/run_triangle_stages.sh"
vortex_velocities
for model in "${model_args[1]}" "${model_args[3]}"; do
    [[ -r "$model" ]] || { printf 'Missing CoordNet model: %s\n' "$model" >&2; exit 1; }
done
case "$query" in
    up)   polyline="$edge_a;$edge_b;$wake_up;$edge_a" ;;
    down) polyline="$edge_b;$edge_a;$wake_down;$edge_b" ;;
esac
step=${VORTEX_STEP:-0.005}
if ! [[ $step =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
    printf 'VORTEX_STEP must be a positive number.\n' >&2
    exit 1
fi
# Isabel's surface range tolerance is 0.005 at h=0.5.
residual_tolerance=$(LC_ALL=C awk -v h="$step" 'BEGIN {
    if (h+0 <= 0) exit 1
    printf "%.17g", 0.01 * h
}') || {
    printf 'VORTEX_STEP must be positive.\n' >&2; exit 1;
}
seed_grid=${VORTEX_SEED_GRID:-40}
max_triangles=${VORTEX_MAX_TRIANGLES:-500000}

if [[ $(type -t spack) != function ]]; then
    source "${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}"
fi
source ./load-env.sh || exit 1
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
cmake --build build --target fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface \
    fiber_stitch_surface fiber_volume_from_patches fiber_grid_recovery -j 3

results_root=${FIBER_RESULTS_ROOT:-"$fiber_repo/build/fiber-volume-results"}
mkdir -p -- "$results_root"
run_dir=$(mktemp -d -- "$results_root/vortex_street_${query}.XXXXXXXX")
printf 'Results: %s\n' "$run_dir"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run_triangle_stages "$run_dir" surface "$step" --model coordnet "${model_args[@]}" \
    --polyline "$polyline" --residual-tolerance "$residual_tolerance" \
    --fiber-field-tolerance 1e-9 --surface-solve-tolerance 1e-10 --derivative-tolerance 1e-10 \
    --seed-grid "$seed_grid" --max-seed-starts 1632160 \
    --max-normal-angle 85 --max-triangles "$max_triangles"
printf 'Assembled boundary and surface audit: %s\n' "$run_dir"
