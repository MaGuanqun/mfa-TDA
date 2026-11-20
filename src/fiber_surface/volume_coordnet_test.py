"""Check CoordNet volume selection in a cropped box and filling a saved boundary."""
import json
import math
from pathlib import Path
import subprocess
import sys

exe, fixture, work = map(Path, sys.argv[1:])
prefix = work / "volume"
common = [exe, "--first", fixture / "quadratic.pt", "--second", fixture / "affine.pt",
          "--domain-min", "-2,1,10", "--domain-max", "6,3,14",
          "--min", "1.75,1.5,11.75", "--max", "2.25,2.5,12.25",
          "--step", "0.09", "--seed-grid", "3", "--output", prefix]
# Only the bottom edge intersects this box. The retained material has
# 1.75*(x-2) - 3*(y-2) + (t-12) >= 0.3. Integrating over x and t gives
# a mean retained y thickness of 0.4 and a volume of 0.5 * 0.5 * 0.4.
selection = ["--q0", "-100,0.3", "--q1", "100,0.3", "--q2", "0,100"]
report_path = work / "volume_report.json"
subprocess.run([*common, *selection, "--stage", "boundary"], check=True)
boundary = json.loads(report_path.read_text())
assert boundary["boundary_valid"] and boundary["tetrahedra"] == 0, boundary
subprocess.run([*common, *selection, "--stage", "fill"], check=True)
filled = json.loads(report_path.read_text())
assert filled["complete"] and filled["tetrahedra"] > 0, filled
assert math.isclose(filled["volume"], 0.1, rel_tol=1e-9), filled
