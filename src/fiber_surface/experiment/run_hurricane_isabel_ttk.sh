#!/usr/bin/env bash
set -eo pipefail

usage() {
    cat <<'EOF'
Usage: run_hurricane_isabel_ttk.sh [both | narrow | wedge | MFA_RUN_DIR ...]
Default: compare narrow and wedge on the full 500 x 500 x 100 MFA sample grid.
Named presets reuse the latest completed MFA run, or run the MFA experiment if missing.
TTK_PYTHON selects a Python with topologytoolkit, VTK, NumPy and SciPy.
FIBER_RESULTS_ROOT sets the results parent. TTK_THREADS defaults to 3.
EOF
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if [[ $# == 0 || ( $# == 1 && $1 == both ) ]]; then set -- narrow wedge; fi

experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$experiment_dir/../../.." && pwd)
results_root=${FIBER_RESULTS_ROOT:-"$repo_root/build/fiber-volume-results"}
mkdir -p -- "$results_root"
results_root=$(cd -- "$results_root" && pwd)
threads=${TTK_THREADS:-3}
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo 'TTK_THREADS must be a positive integer.' >&2; exit 1; }
analysis=("${TTK_PYTHON:-python3}")
imports='import vtk, numpy, scipy; from topologytoolkit import ttkFiberSurface'
if ! "${analysis[@]}" -c "$imports" 2>/dev/null; then
    if [[ -n ${TTK_PYTHON:-} ]] || ! command -v conda >/dev/null; then
        echo 'Set TTK_PYTHON to an interpreter with topologytoolkit, VTK, NumPy and SciPy.' >&2
        exit 1
    fi
    analysis=(conda run --no-capture-output -n mfa python)
    "${analysis[@]}" -c "$imports"
fi

latest_mfa_run() {
    local dir newest=
    for dir in "$results_root/hurricane_isabel_t30_${1}_h1_seed40".*; do
        [[ -s "$dir/volume_tet_boundary.vtk" && -f "$dir/volume_report.json" ]] || continue
        if [[ -z $newest || "$dir/volume_report.json" -nt "$newest/volume_report.json" ]]; then
            if python3 -c 'import json, sys
with open(sys.argv[1]) as f:
    report = json.load(f)
sys.exit(not (report.get("complete") and report.get("boundary_valid")))' "$dir/volume_report.json" 2>/dev/null; then
                newest=$dir
            fi
        fi
    done
    printf '%s' "$newest"
}

runs=()
for case_name in "$@"; do
    if [[ $case_name == narrow || $case_name == wedge ]]; then
        latest=$(latest_mfa_run "$case_name")
        if [[ -z $latest ]]; then
            FIBER_RESULTS_ROOT="$results_root" "$experiment_dir/run_hurricane_isabel.sh" "$case_name"
            latest=$(latest_mfa_run "$case_name")
            [[ -n $latest ]] || { echo "No completed MFA run for $case_name" >&2; exit 1; }
        fi
        runs+=("$latest")
    elif [[ -f "$case_name/volume_report.json" && -s "$case_name/volume_tet_boundary.vtk" ]]; then
        runs+=("$(cd -- "$case_name" && pwd)")
    else
        printf 'Unknown preset or missing completed MFA run: %s\n' "$case_name" >&2
        usage >&2
        exit 1
    fi
done

run_dir=$(mktemp -d "$results_root/hurricane_isabel_t30_ttk.XXXXXXXX")
printf 'Results: %s\n' "$run_dir"
cd -- "$repo_root"
(
    if [[ $(type -t spack) != function ]]; then
        source "${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}"
    fi
    source ./load-env.sh || exit 1
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
    cmake --build build --target sample_local_field -j 3
    sampler=build/src/critical_point_tracking/local_reference/sample_local_field
    for field in Pf30 TCf30; do
        "$sampler" --backend mfa --dataset hurricane_isabel_3d \
            --input "build/src/hurricane_isabel_3d/${field}_3d.mfa" \
            --region 0 499 0 499 0 99 --ratio 1 --threads "$threads" --output "$run_dir/$field"
    done
) 2>&1 | tee "$run_dir/sample.log"

"${analysis[@]}" "$experiment_dir/ttk_hurricane_isabel.py" \
    "$run_dir" "${runs[@]}" --threads "$threads" 2>&1 | tee "$run_dir/ttk.log"
printf 'Completed: %s\n' "$run_dir"
