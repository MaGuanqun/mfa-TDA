#!/usr/bin/env bash
# Run after the build in COORDNET_TABLE.md. Leave the smoke reports under build/.
set -eo pipefail
cd "$(dirname -- "${BASH_SOURCE[0]}")/../../.."
if [[ $(type -t spack) != function ]]; then
    source "${SPACK_SETUP_ENV:-${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh}"
fi
source ./load-env.sh
exp=src/fiber_surface/experiment
source "$exp/run_triangle_stages.sh"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run=$(mktemp -d "$PWD/build/coordnet-table-test.XXXXXXXX")
printf 'Smoke test directory: %s\n' "$run"
"${TORCH_PYTHON:-python3}" - "$run" <<'PY'
from pathlib import Path
import sys
import torch

for name, axis in (("x", 2), ("y", 1)):
    model = torch.jit.trace(lambda p: p[:, axis], torch.zeros(3, 3, dtype=torch.float64))
    model.save(str(Path(sys.argv[1]) / (name + ".pt")))
PY
run_triangle_stages "$run" surface 0.2 --model coordnet --first "$run/x.pt" --second "$run/y.pt" \
    --domain-min 0,0,0 --domain-max 1,1,1 --seed-grid 3 --max-triangles 2000 \
    --polyline '-0.42,-0.26;0.69,-0.18;-0.13,0.68;-0.42,-0.26' \
    --residual-tolerance 0.00005 --fiber-field-tolerance 1e-9 --surface-solve-tolerance 1e-10 \
    --derivative-tolerance 1e-10 --max-normal-angle 85
"${TORCH_PYTHON:-python3}" "$exp/sample_coordnet_grid.py" "$run" "$run/grid" --dims 11 13 15
"${METRICS_PYTHON:-python3}" "$exp/evaluate_coordnet_table.py" "$run" "$run/grid" "$run/metrics" \
    --model-name affine-test --ttk-python "${TTK_PYTHON:-${METRICS_PYTHON:-python3}}"
python3 - "$run/metrics/table.csv" <<'PY'
import csv
from pathlib import Path
import sys

with Path(sys.argv[1]).open() as stream:
    ours, ttk = csv.DictReader(stream)
assert [ours["method"], ttk["method"]] == ["Ours", "TTK"]
for row in (ours, ttk):
    assert int(row["triangles"]) > 0 and int(row["projection_samples"]) == 20000
    assert int(row["side_vertices_checked"]) > 0 and float(row["max_vertex_range_residual"]) < 1e-6
    assert 0 <= float(row["max_vertex_overshoot"]) < 1e-6
    assert int(row["failed_projections"]) == 0 and float(row["d95"]) < 1e-6
assert all(int(ours[k]) == 0 for k in ("open_edges", "pinches", "intersections"))
both, missed, excess, neither = [int(ours[k]) for k in ("both", "fitted_only", "mesh_only", "neither")]
assert both + missed == 16 * 15 and both + missed + excess + neither == 11 * 13 * 15
assert excess == 0 and missed <= 16 and float(ours["iou"]) == both / (both + missed + excess)
assert int(ttk["open_edges"]) > 0 and ttk["iou"] == "" and "boundary_edges" in ttk["iou_note"]
print("CoordNet table smoke test passed")
PY
