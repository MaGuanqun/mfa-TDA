"""Exercise volume CLI validation, outputs, and blocked filling of partial meshes."""
import argparse
import json
from pathlib import Path
import struct
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

# Stop before TetGen and inspect the PLY independently of the C++ reader.
run(['--example', 'prism', '--stage', 'boundary', '--output', prefix], 0)
boundary_report = json.loads(Path(prefix + '_report.json').read_text())
assert boundary_report['stage'] == 'boundary' and boundary_report['boundary_valid']
assert boundary_report['tetrahedra'] == 0
ply = Path(prefix + '_boundary.ply')
header, payload = ply.read_bytes().split(b'end_header\n', 1)
assert header.startswith(b'ply\nformat binary_little_endian 1.0\n')
obj = Path(prefix + '_boundary.obj').read_text().splitlines()
points = [tuple(map(float, line.split()[1:])) for line in obj if line.startswith('v ')]
vertices_end = 24 * len(points)
assert list(struct.iter_unpack('<3d', payload[:vertices_end])) == points
faces = []
for line in obj:
    if line.startswith('g patch_'):
        patch = int(line.split('_')[1])
    if line.startswith('f '):
        faces.append((3, *(int(i)-1 for i in line.split()[1:]), patch))
assert list(struct.iter_unpack('<B3iB', payload[vertices_end:])) == faces

# Use a different prefix and a triangle budget too low to re-extract the mesh.
resumed = str(a.work / 'resumed')
run(['--boundary-input', str(ply), '--max-triangles', '4', '--output', resumed], 0)
filled = json.loads(Path(resumed + '_report.json').read_text())
assert filled['stage'] == 'fill' and filled['complete'] and filled['stage_seconds'][0] == 0
for key in ['boundary_triangles', 'corner_fiber_components', 'box_contact_triangles']:
    assert filled[key] == boundary_report[key], key
assert filled['tetrahedra'] > 0 and abs(filled['volume'] - 0.5) < 1e-10
run(['--boundary-input', str(ply), '--stage', 'boundary'], 1)
bad = a.work / 'bad.ply'
bad.write_bytes(header + b'end_header\n' + payload[:-1])
run(['--boundary-input', str(bad), '--output', resumed], 1)
bad.write_bytes(header + b'end_header\n' + payload[:vertices_end+1] +
                struct.pack('<i', -1) + payload[vertices_end+5:])
run(['--boundary-input', str(bad), '--output', resumed], 1)
# A plausible patch label on the wrong face must fail field validation.
offset = vertices_end + 13
bad.write_bytes(header + b'end_header\n' + payload[:offset] +
                bytes([faces[0][-1] % 9 + 1]) + payload[offset+1:])
run(['--boundary-input', str(bad), '--output', resumed], 2)
invalid = json.loads(Path(resumed + '_report.json').read_text())
assert invalid['failed_stage'] == 'boundary_validation' and invalid['tetrahedra'] == 0

# Reusing the prefix must replace a formerly successful volume with an empty
# output when the new boundary is incomplete, avoiding a stale success artifact.
run(['--example', 'prism', '--max-triangles', '4', '--output', prefix], 2)
report = json.loads(Path(prefix + '_report.json').read_text())
assert not report['complete'] and report['tetrahedra'] == 0
assert report['failed_stage'] == 'boundary_validation'
assert 'CELLS 0 0' in Path(prefix + '_volume.vtk').read_text()
run(['--stage', 'fill', '--output', prefix], 2)
print('Volume CLI and partial-result checks passed')
