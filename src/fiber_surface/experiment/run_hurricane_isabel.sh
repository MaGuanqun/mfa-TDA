#!/usr/bin/env bash
set -eo pipefail

usage() {
    printf 'Usage: %s [narrow|wedge|core|storm|loop] (default: narrow)\nFIBER_RESULTS_ROOT sets the output parent. PYTHON selects the analysis interpreter.\n' "$0"
}
if (( $# > 1 )); then usage >&2; exit 1; fi
max_trace_steps=10000
case "${1:-narrow}" in
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
    for file in "$model_dir/${field}_3d.mfa" "$source_dir/${field}_3d.bin"; do
        [[ -r "$file" ]] || { printf 'Missing Isabel input: %s\n' "$file" >&2; exit 1; }
    done
done
analysis=("${PYTHON:-python3}")
imports='import numpy, vtkmodules.vtkIOLegacy, matplotlib.pyplot'
if ! "${analysis[@]}" -c "$imports" 2>/dev/null; then
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
    cmake --build build --target fiber_volume_mfa -j 3

    extractor=(build/src/fiber_surface/fiber_volume_mfa
        --first "$model_dir/Pf30_3d.mfa" --second "$model_dir/TCf30_3d.mfa"
        --variable1 0 --variable2 0 --seed-grid 40 --max-seed-starts 1632160)
    printf '2. Inspect the models on the full 202 x 202 x 40 lattice.\n'
    "${extractor[@]}" --inspect 2>&1 | tee "$run_dir/inspect.log"

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
)

printf '4. Measure the source footprint, render the boundary and compare any validated fill.\n'
"${analysis[@]}" "$experiment_dir/evaluate_hurricane_isabel.py" "$run_dir" "$source_dir"
printf 'Completed: %s\n' "$run_dir"
