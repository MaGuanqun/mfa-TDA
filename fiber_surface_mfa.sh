#!/usr/bin/env bash
set -eo pipefail

# Usage: ./fiber_surface_mfa.sh [all|seeds|fibers|march|stitch|evaluate]
fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$fiber_repo/src/fiber_surface/run_pipeline.sh"
source "$fiber_repo/src/fiber_surface/model_experiments.sh"
stage="${1:-all}"
if [[ "$stage" == --help || "$stage" == -h ]]; then
    printf 'Usage: %s [all|seeds|fibers|march|stitch|evaluate]\n' "$0"
    printf 'FIBER_SKIP_BUILD=1 reuses a build; FIBER_RESULTS_ROOT overrides the output root.\n'
    printf 'FIRST_MFA, SECOND_MFA, and FIBER_DATASET override input files and the dataset folder.\n'
    exit 0
fi
if [[ $# -gt 1 ]]; then printf 'Too many arguments; use --help.\n' >&2; exit 1; fi
fiber_check_stage "$stage"
dataset="${FIBER_DATASET:-hurricane_isabel_3d_full}"
if [[ ! "$dataset" =~ ^[a-zA-Z0-9_-]+$ ]]; then printf 'FIBER_DATASET must be a simple folder name.\n' >&2; exit 1; fi
isabel_pressure_temperature
fiber_setup
fiber_run_pipeline "$fiber_results_root/$dataset" "$stage" \
    --model mfa "${model_args[@]}" --a "$edge_a" --b "$edge_b" "${mesh_args[@]}"
