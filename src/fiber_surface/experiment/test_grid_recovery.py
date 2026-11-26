"""An affine bivariate field gives three planar sides with known total area."""
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile

seeds, recover = sys.argv[1:]
corners = [(-.287, -.133), (.611, -.227), (.094, .539)]
polyline = ';'.join(f'{x},{y}' for x, y in corners + corners[:1])
with tempfile.TemporaryDirectory() as tmp:
    stages = Path(tmp) / 'stages'
    subprocess.run([seeds, '--run-dir', stages, '--model', 'analytical', '--example', 'strip',
                    '--min', '-1,-1,-1', '--max', '1,1,1', '--seed-grid', '3',
                    '--polyline', polyline], check=True)
    prefix = Path(tmp) / 'prism_sides'
    cases = [('grid', []), ('grid', ['--refine', '-.4,-.5,-.3', '.5,.61,.28', '.14']),
             ('fit-grid', []), ('carrier-grid', []), ('recover-grid', [])]
    for mode, refinement in cases:
        subprocess.run([recover, mode, stages, '.37', prefix, *refinement], check=True)
        report = json.loads(prefix.with_suffix('.json').read_text())
        topology = {k: int(v) for k, v in re.findall(r'(\w+)=(-?\d+)', report['topology'])}
        perimeter = sum(math.dist(corners[i], corners[(i+1) % 3]) for i in range(3))
        assert math.isclose(report['side_area'], 2 * perimeter, rel_tol=1e-12)
        assert report['intersections_checked'] and topology['intersections'] == 0
        assert topology['nonmanifold_vertices'] == topology['nonmanifold_edges'] == 0
        assert topology['inconsistent_edges'] == topology['duplicate_faces'] == 0
        assert topology['boundary_edges'] > 0 and topology['Euler'] == 0
        assert report['seam_tests'] > 0 and report['seam_failures'] == 0
        assert report['side_vertices_checked'] == report['vertices']
        assert report['max_vertex_range_residual'] < 1e-12
        assert 0 <= report['max_vertex_overshoot'] < 1e-12
        assert report['projection_max'] < 1e-12 and report['reversed_sample_normals'] == 0

    # F(x,y,z)=(x,y). Move the shared A edge 0.2 of the way beyond A on CA.
    # Its AB overshoot is smaller, so checking only the first incident label would fail.
    a, b, c = corners
    xy = [v + .2*(v-w) for v, w in zip(a, c)]
    points = [(*xy, -.5), (*xy, .5), (*b, 0.), (*c, 0.)]
    sample = Path(tmp) / 'shifted_seam.ply'
    with sample.open('wb') as stream:
        stream.write(b'ply\nformat binary_little_endian 1.0\ncomment fiber_volume_boundary 0 0 0 0 0\n'
                     b'element vertex 4\nproperty double x\nproperty double y\nproperty double z\n'
                     b'element face 2\nproperty list uchar int vertex_indices\nproperty uchar patch\nend_header\n')
        for point in points:
            stream.write(struct.pack('<3d', *point))
        for face in [(0, 1, 2, 1), (1, 0, 3, 3)]:
            stream.write(struct.pack('<BiiiB', 3, *face))
    subprocess.run([recover, 'audit', stages, sample, prefix], check=True)
    shifted = json.loads(prefix.with_suffix('.json').read_text())
    assert shifted['side_vertices_checked'] == 4
    assert math.isclose(shifted['max_vertex_overshoot'], .2, abs_tol=1e-14)
    assert math.isclose(shifted['max_vertex_range_residual'], math.dist(xy, a), abs_tol=1e-14)
