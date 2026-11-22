#!/usr/bin/env bash
set -eo pipefail

# Usage: ./fiber_surface_mfa.sh [all|seeds|fibers|march|stitch|evaluate]
# These are the validated local Hurricane Isabel pressure/temperature parameters.
# Changing parameters takes effect when the seeds stage is rerun.
# --seed-grid sets the count on the shortest axis; other counts scale with length.
parameters=(--a -1100,-27 --b -1100,-22
            --min 240,240,45 --max 260,260,55
            --step 1 --fiber-step 0.6 --spatial-tolerance 0.08
            --residual-tolerance 0.005 --seed-grid 7 --max-triangles 8000)

fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$fiber_repo/src/fiber_surface/run_pipeline.sh"
stage="${1:-all}"
if [[ "$stage" == --help || "$stage" == -h ]]; then
    printf 'Usage: %s [all|seeds|fibers|march|stitch|evaluate]\n' "$0"
    printf 'FIBER_SKIP_BUILD=1 reuses a build; FIBER_RESULTS_ROOT overrides the output root.\n'
    printf 'FIRST_MFA, SECOND_MFA, and FIBER_DATASET override input files and the dataset folder.\n'
    exit 0
fi
if [[ $# -gt 1 ]]; then printf 'Too many arguments; use --help.\n' >&2; exit 1; fi
fiber_check_stage "$stage"
dataset="${FIBER_DATASET:-hurricane_isabel_3d}"
if [[ ! "$dataset" =~ ^[a-zA-Z0-9_-]+$ ]]; then printf 'FIBER_DATASET must be a simple folder name.\n' >&2; exit 1; fi
first="${FIRST_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/Pf30_3d.mfa}"
second="${SECOND_MFA:-$fiber_repo/build/src/hurricane_isabel_3d/TCf30_3d.mfa}"
fiber_setup
# No additional pressure_temperature_local directory: all stages live here.
fiber_run_pipeline "$fiber_results_root/$dataset" "$stage" \
    --model mfa --first "$first" --second "$second" --variable1 0 --variable2 0 "${parameters[@]}"
