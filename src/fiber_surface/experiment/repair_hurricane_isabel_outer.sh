#!/usr/bin/env bash
# Reproduce the two local refinements for the Isabel t=30 candidate-2 experiment.
set -eo pipefail
if [[ $# != 2 ]]; then
    printf 'Usage: %s FULL_H0.5_NORMAL85_RUN NEW_OUTPUT_DIRECTORY\n' "$0" >&2
    exit 1
fi
full=$(realpath "$1")
out=$(realpath -m "$2")
experiment=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$experiment/../../.."
if [[ $(type -t spack) != function ]]; then
    source "${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}"
fi
source ./load-env.sh
py=${PYTHON:-python3}
"$py" - "$full" <<'PY'
import json, sys
from pathlib import Path
import numpy, scipy
root = Path(sys.argv[1])
report = json.loads((root/'volume_report.json').read_text())
config = (root/'stages/run.conf').read_text().splitlines()
options = next(line.split()[1:] for line in config if line.startswith('options '))
if (report['triangle'] != [[-25,-19],[165,-22],[165,-15]] or report['step'] != .5 or
        report['spatial_tolerance'] != .5 or report['fiber_seed_radius'] != .5 or float(options[7]) != 85 or
        report['box_min'] != [0,0,0] or report['box_max'] != [499,499,99]):
    sys.exit('Expected the full-domain Isabel candidate-2 h=0.5, 85-degree run')
PY
mkdir -p "$(dirname "$out")"
mkdir "$out"
bin=build/src/fiber_surface
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
cmake --build build -j3 --target fiber_repair_boundary fiber_grid_recovery fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface fiber_stitch_surface fiber_volume_from_patches
"$bin/fiber_repair_boundary" "$full/stages" "$full/volume_boundary.ply" "$out/repaired" --repair

mapfile -t models < <("$py" - "$full/stages/run.conf" <<'PY'
import shlex, sys
from pathlib import Path
lines = [shlex.split(line) for line in Path(sys.argv[1]).read_text().splitlines()]
for key in ['first','second']:
    print(next(row[1] for row in lines if row and row[0] == key))
PY
)
partial_stage() {
    local status=0
    "$@" || status=$?
    [[ $status == 0 || $status == 2 ]]
}
crop() {
    local dir="$out/$1" lo=$2 hi=$3
    mkdir -p "$dir"
    "$bin/fiber_find_seeds" --run-dir "$dir/stages" --model mfa --first "${models[0]}" --second "${models[1]}" \
        --variable1 0 --variable2 0 --polyline '-25,-19;165,-22;165,-15;-25,-19' --min "$lo" --max "$hi" \
        --step .125 --fiber-seed-radius .125 --spatial-tolerance .125 --residual-tolerance .005 \
        --max-normal-angle 85 --max-triangles 500000 --max-trace-steps 100000 --seed-grid 7
    partial_stage "$bin/fiber_trace_endpoint_fibers" --run-dir "$dir/stages"
    for side in 0 1 2; do
        partial_stage "$bin/fiber_march_surface" --run-dir "$dir/stages" --segment "$side"
        partial_stage "$bin/fiber_stitch_surface" --run-dir "$dir/stages" --segment "$side"
    done
    "$bin/fiber_volume_from_patches" --run-dir "$dir/stages" --output "$dir/volume"
}
crop isolated 323.5,351.5,42.25 326.5,354.5,43.25
crop main 277.5,348.5,42.5 280.5,351.5,43.5
"$py" "$experiment/join_isabel_refinements.py" "$out/repaired.ply" "$out/isolated/volume_boundary.ply" \
    "$out/main/volume_boundary.ply" "$out/adaptive"
for mesh in adaptive_transition adaptive; do
    "$bin/fiber_grid_recovery" audit "$full/stages" "$out/$mesh.ply" "$out/$mesh" --all-face-samples
done
"$py" - "$out/adaptive.json" <<'PY'
import json, re, sys
from pathlib import Path
report = json.loads(Path(sys.argv[1]).read_text())
topology = dict(re.findall(r'(\w+)=(-?\d+)', report['topology']))
if (any(int(v) for k,v in topology.items() if k != 'Euler') or not report['intersections_checked']
        or report['face_sites_projection_failed'] or report['seam_failures']
        or report['face_sites_max_angle'] > 85 or report['face_sites_max_projection'] > .5):
    sys.exit('The assembled surface failed the topology or fitted-field checks')
PY
printf 'Closed labeled surface saved to %s/adaptive.ply\n' "$out"
