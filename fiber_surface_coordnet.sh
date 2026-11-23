#!/usr/bin/env bash
set -eo pipefail

dataset=vortex_street_3d_2_variables_full
stage=all                     # all, seeds, fibers, march, stitch, or evaluate

fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$fiber_repo/src/fiber_surface/run_pipeline.sh"
source "$fiber_repo/src/fiber_surface/model_experiments.sh"
vortex_velocities
fiber_check_stage "$stage"
fiber_setup
# TBB parallelizes seed discovery. Each individual network query is one point.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
fiber_run_pipeline "$fiber_results_root/$dataset" "$stage" \
    --model coordnet "${model_args[@]}" --a "$edge_a" --b "$edge_b" "${mesh_args[@]}"
