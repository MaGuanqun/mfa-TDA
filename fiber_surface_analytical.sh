#!/usr/bin/env bash
set -eo pipefail

# Usage: ./fiber_surface_analytical.sh [sphere|closed-sphere|strip|disconnected|all] [stage]
# Stage: all (default), seeds, fibers, march, stitch, evaluate.
# Edit these seed-stage parameters, then rerun all stages to apply the change.
# --seed-grid sets the count on the shortest axis; other counts scale with length.
parameters=(--step 0.2 --fiber-step 0.15 --spatial-tolerance 0.015
            --residual-tolerance 1e-9 --seed-grid 7 --max-triangles 30000)

fiber_repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$fiber_repo/src/fiber_surface/run_pipeline.sh"
example="${1:-sphere}"
stage="${2:-all}"
if [[ "$example" == --help || "$example" == -h ]]; then
    printf 'Usage: %s [sphere|closed-sphere|strip|disconnected|all] [all|seeds|fibers|march|stitch|evaluate]\n' "$0"
    printf 'FIBER_SKIP_BUILD=1 reuses an existing build; FIBER_RESULTS_ROOT overrides the output root.\n'
    exit 0
fi
if [[ $# -gt 2 ]]; then printf 'Too many arguments; use --help.\n' >&2; exit 1; fi
case "$example" in sphere|closed-sphere|strip|disconnected|all) ;; *) printf 'Unknown example: %s\n' "$example" >&2; exit 1;; esac
fiber_check_stage "$stage"
fiber_setup
examples=("$example")
if [[ "$example" == all ]]; then examples=(sphere closed-sphere strip disconnected); fi
overall=0
for example in "${examples[@]}"; do
    status=0
    fiber_run_pipeline "$fiber_results_root/analytical_${example//-/_}" "$stage" \
        --model analytical --example "$example" "${parameters[@]}" || status=$?
    if [[ "$status" != 0 && "$status" != 2 ]]; then exit "$status"; fi
    if [[ "$status" == 2 ]]; then overall=2; fi
done
exit "$overall"
