#!/usr/bin/env bash
set -eo pipefail

usage() {
    printf 'Usage: %s [outer|narrow|wedge|core|storm|loop] (default: outer)\nISABEL_OUTER_STEP sets h (default 0.5). ISABEL_OUTER_MAX_NORMAL_ANGLE sets the face-to-field normal limit in degrees (default 85, strictly between 0 and 90). ISABEL_OUTER_MAX_TRIANGLES sets the per-side budget (default 500000). ISABEL_OUTER_STAGE=surface exports and audits the assembled surface. FIBER_RESULTS_ROOT sets the output parent. PYTHON selects the analysis interpreter.\n' "$0"
}
if (( $# > 1 )); then usage >&2; exit 1; fi
max_trace_steps=10000
case "${1:-outer}" in
    outer)
        stage=${ISABEL_OUTER_STAGE:-all}
        if [[ $stage != all && $stage != surface ]]; then
            printf 'ISABEL_OUTER_STAGE must be all or surface.\n' >&2
            exit 1
        fi
        step=${ISABEL_OUTER_STEP:-0.5}
        if [[ ! $step =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ || ! $step =~ [1-9] ]]; then
            printf 'ISABEL_OUTER_STEP must be a positive decimal step in grid units.\n' >&2
            exit 1
        fi
        max_triangles=${ISABEL_OUTER_MAX_TRIANGLES:-500000}
        if [[ ! $max_triangles =~ ^[1-9][0-9]*$ ]]; then
            printf 'ISABEL_OUTER_MAX_TRIANGLES must be a positive integer per side.\n' >&2
            exit 1
        fi
        name="outer_h${step}_budget${max_triangles}"
        normal_angle=${ISABEL_OUTER_MAX_NORMAL_ANGLE:-85}
        if ! awk -v a="$normal_angle" 'BEGIN { exit !(a ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ && a+0 > 0 && a+0 < 90) }'; then
            printf 'ISABEL_OUTER_MAX_NORMAL_ANGLE must be strictly between 0 and 90 degrees.\n' >&2
            exit 1
        fi
        if [[ $normal_angle != 55 ]]; then name+="_normal${normal_angle}"; fi
        triangle=(--q0 -25,-19 --q1 165,-22 --q2 165,-15)
        ;;
    narrow)
        name=narrow_h1
        triangle=(--q0 -1000,-22 --q1 -900,-22 --q2 -950,-18)
        step=1.0; spatial_tolerance=0.35; max_triangles=50000; max_trace_steps=100000 ;;
    wedge)
        name=wedge_h1
        triangle=(--q0 -1000,-21 --q1 -850,-21 --q2 -850,-15)
        step=1.0; spatial_tolerance=0.20; max_triangles=50000; max_trace_steps=100000 ;;
    loop)
        name=p1200_h1
        triangle=(--q0 -1100,-27 --q1 -1100,-22 --q2 -1200,-22)
        step=1.0; spatial_tolerance=0.08; max_triangles=8000 ;;
    storm)
        name=storm_h5_budget3500
        triangle=(--q0 -1500,-36 --q1 100,-36 --q2 150,10)
        step=5.0; spatial_tolerance=0.4; max_triangles=3500 ;;
    core)
        name=core_h5_tol1p6
        triangle=(--q0 -3300,-35 --q1 -550,-55 --q2 -550,19)
        step=5.0; spatial_tolerance=1.6; max_triangles=15000 ;;
    --help|-h) usage; exit 0 ;;
    *) usage >&2; exit 1 ;;
esac

experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$experiment_dir/../../.." && pwd)
results_root=${FIBER_RESULTS_ROOT:-"$repo_root/build/fiber-volume-results"}
mkdir -p -- "$results_root"
results_root=$(cd -- "$results_root" && pwd)
cd "$repo_root"

model_dir=build/src/hurricane_isabel_3d
source_dir=build/examples/hurricane_isabel
for field in Pf30 TCf30; do
    files=("$model_dir/${field}_3d.mfa")
    if [[ ${stage:-all} != surface ]]; then files+=("$source_dir/${field}_3d.bin"); fi
    for file in "${files[@]}"; do
        [[ -r "$file" ]] || { printf 'Missing Isabel input: %s\n' "$file" >&2; exit 1; }
    done
done
analysis=("${PYTHON:-python3}")
imports='import numpy, vtkmodules.vtkIOLegacy, matplotlib.pyplot'
if [[ ${stage:-all} != surface ]] && ! "${analysis[@]}" -c "$imports" 2>/dev/null; then
    if [[ -n ${PYTHON:-} ]] || ! command -v conda >/dev/null; then
        printf 'Set PYTHON to an interpreter with NumPy, VTK and Matplotlib.\n' >&2
        exit 1
    fi
    analysis=(conda run --no-capture-output -n mfa_env python)
    "${analysis[@]}" -c "$imports"
fi

# A separate directory for each invocation also preserves partial runs.
run_dir=$(mktemp -d "$results_root/hurricane_isabel_t30_${name}_seed40.XXXXXXXX")
printf 'Results: %s\n' "$run_dir"
(
    if [[ $(type -t spack) != function ]]; then
        source "${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}"
    fi
    source ./load-env.sh || exit 1
    printf '1. Build the MFA volume extractor.\n'
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
    targets=(fiber_volume_mfa)
    if [[ ${1:-outer} == outer ]]; then
        targets+=(fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface
                  fiber_stitch_surface fiber_volume_from_patches)
        if [[ $stage == surface ]]; then targets+=(fiber_grid_recovery); fi
    fi
    cmake --build build --target "${targets[@]}" -j 3

    extractor=(build/src/fiber_surface/fiber_volume_mfa
        --first "$model_dir/Pf30_3d.mfa" --second "$model_dir/TCf30_3d.mfa"
        --variable1 0 --variable2 0 --seed-grid 40 --max-seed-starts 1632160)
    printf '2. Inspect the models on the full 202 x 202 x 40 lattice.\n'
    "${extractor[@]}" --inspect 2>&1 | tee "$run_dir/inspect.log"

    if [[ ${1:-outer} == outer ]]; then
        source "$experiment_dir/run_triangle_stages.sh"
        polyline="${triangle[1]};${triangle[3]};${triangle[5]};${triangle[1]}"
        run_triangle_stages "$run_dir" "$stage" "$step" --model mfa "${extractor[@]:1}" \
            --polyline "$polyline" --residual-tolerance 0.005 \
            --max-normal-angle "$normal_angle" --max-triangles "$max_triangles"
    else
        command=("${extractor[@]}" "${triangle[@]}" --step "$step"
            --spatial-tolerance "$spatial_tolerance" --residual-tolerance 0.005
            --max-triangles "$max_triangles" --max-trace-steps "$max_trace_steps"
            --stage all --output "$run_dir/volume")
        printf '%q ' "${command[@]}" > "$run_dir/run.txt"
        printf '\n' >> "$run_dir/run.txt"
        printf '3. Trace the corner fibers, assemble and validate the boundary, then tetrahedralize.\n'
        "${command[@]}" 2>&1 | tee "$run_dir/extract.log" || {
            status=("${PIPESTATUS[@]}")
            # A partial mesh has a diagnostic report and can still be plotted.
            if (( status[0] != 2 || status[1] != 0 )); then exit 1; fi
        }
    fi
)

if [[ ${1:-outer} == outer && $stage == surface ]]; then
    printf 'Assembled boundary and surface audit: %s\n' "$run_dir"
    exit 0
fi
printf '4. Measure the source footprint, render the boundary and compare any validated fill.\n'
"${analysis[@]}" "$experiment_dir/evaluate_hurricane_isabel.py" "$run_dir" "$source_dir"
printf 'Completed: %s\n' "$run_dir"
