"""An affine bivariate field gives three planar sides with known total area."""
import json
import math
from pathlib import Path
import re
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
        assert report['projection_max'] < 1e-12 and report['reversed_sample_normals'] == 0
