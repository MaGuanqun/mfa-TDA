"""A missing interior triangle is recoverable while the two open box rims stay intact."""
import json
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile

seed, grid, repair, refine = sys.argv[1:]
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    stages = work / 'stages'
    subprocess.run([seed, '--run-dir', stages, '--model', 'analytical', '--example', 'strip',
                    '--min', '-1,-1,-1', '--max', '1,1,1', '--seed-grid', '3', '--step', '.5',
                    '--polyline', '-.287,-.133;.611,-.227;.094,.539;-.287,-.133'], check=True)
    prism = work / 'prism'
    subprocess.run([grid, 'grid', stages, '.25', prism], check=True)
    header, payload = prism.with_suffix('.ply').read_bytes().split(b'end_header\n', 1)
    vertex_count = int(re.search(rb'element vertex (\d+)', header)[1])
    count = int(re.search(rb'element face (\d+)', header)[1])
    points = list(struct.iter_unpack('<ddd', payload[:24*vertex_count]))
    cells = list(struct.iter_unpack('<BiiiB', payload[24*vertex_count:]))
    # Choose a face whose corners are interior to the side, away from box caps.
    incident_labels = [set() for _ in points]
    for _, a, b, c, label in cells:
        for v in (a, b, c):
            incident_labels[v].add(label)
    missing = next(i for i, cell in enumerate(cells)
                   if all(abs(points[v][2]) < .8 and len(incident_labels[v]) == 1 for v in cell[1:4]))
    damaged = work / 'damaged.ply'
    header = header.replace(f'element face {count}\n'.encode(), f'element face {count-1}\n'.encode())
    offset = vertex_count*24 + missing*14
    damaged.write_bytes(header + b'end_header\n' + payload[:offset] + payload[offset+14:])
    fixed = work / 'fixed'
    subprocess.run([repair, stages, damaged, fixed, '--repair'], check=True)
    subprocess.run([grid, 'audit', stages, fixed.with_suffix('.ply'), fixed], check=True)
    original, result = (json.loads(p.with_suffix('.json').read_text()) for p in (prism, fixed))
    assert result['triangles'] == count
    assert result['topology'] == original['topology']
    assert abs(result['side_area'] - original['side_area']) < 1e-10
    assert result['seam_failures'] == result['failed_projections'] == 0
    fitted = work / 'fitted'
    subprocess.run([refine, stages, fixed.with_suffix('.ply'), fitted.with_suffix('.ply')], check=True)
    subprocess.run([grid, 'audit', stages, fitted.with_suffix('.ply'), fitted, '--all-face-samples'], check=True)
    refined = json.loads(fitted.with_suffix('.json').read_text())
    assert refined['triangles'] > result['triangles']
    assert refined['topology'] == original['topology']
    assert abs(refined['side_area'] - original['side_area']) < 1e-10
    assert refined['seam_chord_sites'] > 0
    assert refined['seam_chord_failures'] == refined['face_sites_projection_failed'] == 0
    assert refined['max_seam_chord_projection'] < 1e-10
