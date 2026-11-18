"""Exercise volume CLI validation, outputs, and blocked filling of partial meshes."""
import argparse
import json
from pathlib import Path
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--exe', required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)

def run(options, expected):
    result = subprocess.run([a.exe, *options], text=True, capture_output=True)
    assert result.returncode == expected, result.stdout + result.stderr

run(['--q0', '0,0', '--q1', '0,1', '--q2', '1,0'], 1)
run(['--q0', '0,0'], 1)
run(['--q0', '0,0', '--q1', '1,1', '--q2', '2,2'], 1)
run(['--example', 'unknown'], 1)
prefix = str(a.work / 'volume')
run(['--example', 'prism', '--output', prefix], 0)
report = json.loads(Path(prefix + '_report.json').read_text())
assert report['complete'] and abs(report['volume'] - 0.5) < 1e-10
assert report['tetrahedra'] > 0
for suffix in ['_boundary.obj', '_boundary.vtk', '_boundary.poly', '_volume.vtk', '_tet_boundary.vtk']:
    assert Path(prefix + suffix).stat().st_size > 0
# Reusing the prefix must replace a formerly successful volume with an empty
# output when the new boundary is incomplete, avoiding a stale success artifact.
run(['--example', 'prism', '--max-triangles', '4', '--output', prefix], 2)
report = json.loads(Path(prefix + '_report.json').read_text())
assert not report['complete'] and report['tetrahedra'] == 0
assert report['failed_stage'] == 'boundary_validation'
assert 'CELLS 0 0' in Path(prefix + '_volume.vtk').read_text()
print('Volume CLI and partial-result checks passed')
