#!/usr/bin/env bash
set -eo pipefail

# Edit this block, then run ./fiber_surface_coordnet.sh (no arguments).
dataset=vortex_street_3d_2_variables
model_dir="CoordNet/Exp/$dataset"
first="$model_dir/u-super-spatial-temporal-64-5-float64.pt"
second="$model_dir/v-super-spatial-temporal-64-5-float64.pt"
domain_min=-0.5,-0.5,13.5     # Full training box, in physical (x,y,time).
domain_max=7.5,0.5,15
box_min=0.38,0.133,13.988     # Local extraction box in the same coordinates.
box_max=0.42,0.167,14.012
range_args=(--a 0.402,-0.344 --b 0.402,-0.296) # Normalized (u,v) outputs.
# For a polyline, replace range_args with e.g. (--polyline 'u0,v0;u1,v1;u2,v2').
mesh_args=(--step 0.006 --fiber-step 0.004 --spatial-tolerance 0.00045
           --residual-tolerance 0.0002 --seed-grid 3 --max-triangles 4000)
stage=all                     # all, seeds, fibers, march, stitch, or evaluate

fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$fiber_repo/src/fiber_surface/run_pipeline.sh"
fiber_check_stage "$stage"
fiber_setup
# TBB parallelizes seed discovery. Each individual network query is one point.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
fiber_run_pipeline "$fiber_results_root/$dataset" "$stage" \
    --model coordnet --first "$first" --second "$second" \
    --domain-min "$domain_min" --domain-max "$domain_max" \
    --min "$box_min" --max "$box_max" "${range_args[@]}" "${mesh_args[@]}"
