#!/usr/bin/env bash
# Analytical quartic_potential_2 -> sampled gradient -> slice seeds -> SFFF
# -> separated branches with the analytical reference's ColorId values.
set -eo pipefail

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$repo_root"
original_args=("$@")
nx=100 ny=100 nt=100
output_dir=""
reference="$repo_root/build/src/quartic_potential_2/quartic_potential_2_original_4_colorid.vtp"
pvpython_bin=${PVPYTHON:-/usr/local/bin/pvpython}
postprocess_python=${POSTPROCESS_PYTHON:-pvpython}
build=1
step=-1
strength=10
max_steps=10000

usage() {
    cat <<'EOF'
Usage: bash feature_flow_fields_analytical.sh [options]
  --grid NX NY NT       Sample counts, including endpoints (default: 100 100 100)
  --output-dir DIR      Default: build/src/quartic_potential_2/sfff_analytical_NXxNYxNT
  --reference-vtp FILE  Reference branches carrying RegionId and ColorId
  --step LENGTH        RK4 arc-length step (default: minimum spatial grid spacing)
  --strength K         Stable FFF convergence strength (default: 10)
  --max-steps N         Integration limit per direction (default: 10000)
  --pvpython EXE       ParaView Python with TTKTopologicalSkeleton and NumPy
                       (default: PVPYTHON or /usr/local/bin/pvpython)
  --postprocess-python EXE  Python with VTK, NumPy and SciPy (default: pvpython
                           on PATH, or POSTPROCESS_PYTHON)
  --no-build           Reuse existing executables
  -h, --help           Show this help

Domain: x,y in [-2,2], t in [0,4]. Samples are exact analytical derivatives.
Tracking uses the sampled field; seeds are extracted from that field by TTK.
Output includes gradient.vff, seeds.csv, tracking.obj/.vtp,
tracking_colorid.vtp, color_matches.csv, branches/*.vtp and branches.csv.
In ParaView, color tracking_colorid.vtp by ColorId using the reference's palette;
RegionId selects individual branches. Geometry coordinates are (x,y,t).
Repeating an invocation overwrites its experiment outputs.
EOF
}
fail() { printf 'Error: %s\n' "$*" >&2; exit 1; }
while (($#)); do
    case "$1" in
        --grid) [[ $# -ge 4 ]] || fail '--grid requires NX NY NT'; nx=$2; ny=$3; nt=$4; shift 4 ;;
        --output-dir|--reference-vtp|--step|--strength|--max-steps|--pvpython|--postprocess-python)
            [[ $# -ge 2 ]] || fail "$1 requires a value"
            case "$1" in
                --output-dir) output_dir=$2 ;;
                --reference-vtp) reference=$2 ;;
                --step) step=$2 ;;
                --strength) strength=$2 ;;
                --max-steps) max_steps=$2 ;;
                --pvpython) pvpython_bin=$2 ;;
                --postprocess-python) postprocess_python=$2 ;;
            esac
            shift 2 ;;
        --no-build) build=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) fail "Unknown option: $1" ;;
    esac
done
for n in "$nx" "$ny" "$nt"; do
    [[ $n =~ ^[0-9]+$ && $n -ge 2 ]] || fail 'Grid dimensions must be integers >= 2'
done
[[ $max_steps =~ ^[0-9]+$ && $max_steps -gt 0 ]] || fail '--max-steps must be positive'
python3 - "$step" "$strength" <<'PY'
import math
import sys
step, strength = map(float, sys.argv[1:])
if not math.isfinite(step) or (step != -1 and step <= 0):
    raise SystemExit('--step must be positive (or -1 for automatic)')
if not math.isfinite(strength) or strength <= 0:
    raise SystemExit('--strength must be finite and positive')
PY
[[ -f $reference ]] || fail "Missing color reference: $reference"
command -v "$pvpython_bin" >/dev/null || fail "Missing pvpython: $pvpython_bin"
command -v "$postprocess_python" >/dev/null || fail "Missing postprocessing Python: $postprocess_python"
output_dir=${output_dir:-"$repo_root/build/src/quartic_potential_2/sfff_analytical_${nx}x${ny}x${nt}"}
mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)

# Activate in the same shell used for CMake and the C++ executables.
if [[ $(type -t spack || true) != function ]]; then
    spack_setup=${SPACK_SETUP_ENV:-"${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"}
    [[ -f $spack_setup ]] || fail 'Initialize Spack shell support or set SPACK_SETUP_ENV'
    source "$spack_setup"
fi
# load-env.sh may deactivate an inactive environment, which returns nonzero.
set +e
source ./load-env.sh
env_status=$?
set -e
[[ $env_status -eq 0 && ${SPACK_ENV##*/} == mfa-TDA-env && -d $MFA_PATH && -d $TORCH_PATH ]] || fail 'Could not activate mfa-TDA-env'

pv=("$pvpython_bin" --no-mpi --force-offscreen-rendering)
post=("$postprocess_python")
if [[ ${postprocess_python##*/} == pvpython* ]]; then
    post+=(--no-mpi --force-offscreen-rendering)
fi
"${pv[@]}" -c 'import numpy, vtk; from paraview.simple import TTKTopologicalSkeleton; print("ParaView/TTK dependencies available")'
"${post[@]}" -c 'import numpy, scipy, vtk; print("VTK/NumPy/SciPy postprocessing dependencies available")'
{
    printf 'cd %q\n' "$repo_root"
    printf 'bash %q ' "$repo_root/feature_flow_fields_analytical.sh"
    printf '%q ' "${original_args[@]}"
    printf '\n'
} > "$output_dir/invocation.sh"
run_logged() {
    local log=$1
    shift
    printf '%q ' "$@" >> "$output_dir/commands.txt"
    printf '\n' >> "$output_dir/commands.txt"
    "$@" 2>&1 | tee "$output_dir/$log.log"
}
: > "$output_dir/commands.txt"
if ((build)); then
    run_logged configure cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-flto=auto -Dmfa_thread=tbb
    run_logged build cmake --build build --target sample_quartic_gradient stable_feature_flow_fields -j 3
fi
sampler="$repo_root/build/src/feature_flow_fields/sample_quartic_gradient"
tracker="$repo_root/build/src/feature_flow_fields/stable_feature_flow_fields"
[[ -x $sampler && -x $tracker ]] || fail 'Missing executables; rerun without --no-build'
run_logged sampling "$sampler" -x "$nx" -y "$ny" -t "$nt" -o "$output_dir/gradient.vff"
run_logged seeds "${pv[@]}" "$repo_root/src/feature_flow_fields/vff_critical_points.py" \
    -i "$output_dir/gradient.vff" -o "$output_dir/seeds.csv"
run_logged tracking "$tracker" -f "$output_dir/gradient.vff" -s "$output_dir/seeds.csv" \
    -b "$output_dir/tracking" -l "$step" -m "$max_steps" -k "$strength"
run_logged conversion "${post[@]}" "$repo_root/src/python/convert_obj_to_vtp.py" \
    -i "$output_dir/tracking.obj" -o "$output_dir/tracking.vtp"
# Infinite threshold copies the closest reference branch's color for every
# branch, as in seperate_ttk_trajectory.py, while also recording match distances.
run_logged coloring "${post[@]}" "$repo_root/src/python/seperate_branch_assign_color2.py" \
    -i "$output_dir/tracking.vtp" -r "$reference" -o "$output_dir/tracking_colorid.vtp" \
    --distance-mode mean --match-distance-threshold inf --seed 0 \
    --matches-csv "$output_dir/color_matches.csv"
run_logged branches "${post[@]}" "$repo_root/src/feature_flow_fields/export_sfff_branches.py" \
    --experiment-dir "$output_dir" --reference-vtp "$reference"
printf '\nColored result: %s/tracking_colorid.vtp\nIndividual branches: %s/branches/\n' "$output_dir" "$output_dir"
