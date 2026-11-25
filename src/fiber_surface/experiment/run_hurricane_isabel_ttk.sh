#!/usr/bin/env bash
set -eo pipefail

usage() {
    cat <<'EOF'
Usage: run_hurricane_isabel_ttk.sh [outer | both | narrow | wedge | MFA_RUN_DIR ...]
Default: compare the broad outer ring on the full 500 x 500 x 100 MFA sample grid.
Use both for the earlier narrow and wedge triangles. Named presets prefer a
completed fill, then a saved partial boundary, and run the MFA code if missing.
ISABEL_OUTER_STEP and ISABEL_OUTER_MAX_NORMAL_ANGLE select the matching outer-ring
MFA run (defaults 0.5 grid units and 85 degrees).
TTK_PYTHON selects a Python with topologytoolkit, VTK, NumPy, SciPy and Matplotlib.
FIBER_RESULTS_ROOT sets the results parent. TTK_THREADS defaults to 3.
EOF
}
if [[ ${1:-} == --help || ${1:-} == -h ]]; then usage; exit 0; fi
if [[ $# == 0 ]]; then set -- outer; fi
if [[ $# == 1 && $1 == both ]]; then set -- narrow wedge; fi

experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$experiment_dir/../../.." && pwd)
results_root=${FIBER_RESULTS_ROOT:-"$repo_root/build/fiber-volume-results"}
mkdir -p -- "$results_root"
results_root=$(cd -- "$results_root" && pwd)
threads=${TTK_THREADS:-3}
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo 'TTK_THREADS must be a positive integer.' >&2; exit 1; }
analysis=("${TTK_PYTHON:-python3}")
imports='import vtk, numpy, scipy, matplotlib.pyplot; from topologytoolkit import ttkFiberSurface'
if ! "${analysis[@]}" -c "$imports" 2>/dev/null; then
    if [[ -n ${TTK_PYTHON:-} ]] || ! command -v conda >/dev/null; then
        echo 'Set TTK_PYTHON to an interpreter with topologytoolkit, VTK, NumPy, SciPy and Matplotlib.' >&2
        exit 1
    fi
    analysis=(conda run --no-capture-output -n mfa python)
    "${analysis[@]}" -c "$imports"
fi

latest_mfa_run() {
    python3 - "$results_root" "$1" "${ISABEL_OUTER_STEP:-0.5}" "${ISABEL_OUTER_MAX_NORMAL_ANGLE:-85}" "${ISABEL_OUTER_MAX_TRIANGLES:-500000}" <<'PYTHON'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
runs = []
preset = sys.argv[2]
step = float(sys.argv[3]) if preset == "outer" else 1
if preset == "outer":
    angle, budget = sys.argv[4:6]
    if not 0 < float(angle) < 90:
        sys.exit("ISABEL_OUTER_MAX_NORMAL_ANGLE must be strictly between 0 and 90 degrees.")
    suffix = "" if angle == "55" else f"_normal{angle}"
    pattern = f"outer_h*_budget{budget}{suffix}"
else:
    pattern = f"{preset}_h1"
for path in root.glob(f"hurricane_isabel_t30_{pattern}_seed40.*/volume_report.json"):
    try:
        report = json.loads(path.read_text())
        if report["step"] != step:
            continue
        if preset == "outer" and (report["spatial_tolerance"] != step or
                                  report.get("fiber_seed_radius", report["spatial_tolerance"]) != step):
            continue
        filled = report["complete"] and report["boundary_valid"]
        boundary = path.parent / ("volume_tet_boundary.vtk" if filled else "volume_boundary.vtk")
        if boundary.stat().st_size:
            runs.append((filled, path.stat().st_mtime_ns, path.parent))
    except (OSError, ValueError, KeyError):
        continue
if runs:
    print(max(runs)[2])
PYTHON
}

runs=()
for case_name in "$@"; do
    if [[ $case_name == outer || $case_name == narrow || $case_name == wedge ]]; then
        latest=$(latest_mfa_run "$case_name")
        if [[ -z $latest ]]; then
            # Exit 2 means MFA wrote a diagnostic boundary. TTK still has a range query.
            FIBER_RESULTS_ROOT="$results_root" "$experiment_dir/run_hurricane_isabel.sh" "$case_name" || {
                status=$?
                (( status == 2 )) || exit "$status"
            }
            latest=$(latest_mfa_run "$case_name")
            [[ -n $latest ]] || { echo "No MFA boundary report for $case_name" >&2; exit 1; }
        fi
        runs+=("$latest")
    elif [[ -f "$case_name/volume_report.json" && ( -s "$case_name/volume_tet_boundary.vtk" || -s "$case_name/volume_boundary.vtk" ) ]]; then
        runs+=("$(cd -- "$case_name" && pwd)")
    else
        printf 'Unknown preset or missing MFA boundary: %s\n' "$case_name" >&2
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

source_grid="$repo_root/build/examples/hurricane_isabel"
source_args=()
if [[ -r "$source_grid/Pf30_3d.bin" && -r "$source_grid/TCf30_3d.bin" ]]; then
    source_args=(--source-grid "$source_grid")
fi
"${analysis[@]}" "$experiment_dir/ttk_hurricane_isabel.py" \
    "$run_dir" "${runs[@]}" --threads "$threads" "${source_args[@]}" 2>&1 | tee "$run_dir/ttk.log"
printf 'Completed: %s\n' "$run_dir"
