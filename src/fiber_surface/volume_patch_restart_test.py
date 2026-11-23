"""Compare saved, independently marched triangle sides with direct volume extraction."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main(binary, work):
    triangle = ("0,0", "1,0", "0,1")
    def run(name, *args, expected=0):
        p = subprocess.run([binary/name, *args], capture_output=True, text=True)
        if p.returncode != expected:
            raise AssertionError(f"{name} returned {p.returncode}: {p.stdout}{p.stderr}")

    for budget, expected in ((30000, 0), (4, 2)):
        folder = work/str(budget)
        folder.mkdir(parents=True, exist_ok=True)
        options = ("--step", "0.2", "--spatial-tolerance", "0.015", "--max-triangles", str(budget))
        run("fiber_find_seeds", "--run-dir", folder, "--model", "analytical", "--example", "strip",
            "--min", "-1,-1,0", "--max", "2,2,1", "--polyline", ";".join((*triangle, triangle[0])), *options)
        run("fiber_trace_endpoint_fibers", "--run-dir", folder)
        run("fiber_march_surface", "--run-dir", folder, "--segment", "3", expected=1)
        for i in (2, 0, 1):
            run("fiber_march_surface", "--run-dir", folder, "--segment", str(i))
            run("fiber_stitch_surface", "--run-dir", folder, "--segment", str(i), expected=expected)
        direct, resumed = folder/"direct", folder/"resumed"
        run("fiber_volume_explicit", "--example", "prism", "--output", direct, *options, expected=expected)
        run("fiber_volume_from_patches", "--run-dir", folder, "--output", resumed, expected=expected)
        for suffix in ("_boundary.ply", "_volume.vtk", "_tet_boundary.vtk"):
            assert Path(str(direct)+suffix).read_bytes() == Path(str(resumed)+suffix).read_bytes(), suffix
        reports = [json.loads(Path(str(prefix)+"_report.json").read_text()) for prefix in (direct, resumed)]
        for r in reports:
            del r["stage_seconds"]
        assert reports[0] == reports[1]


if __name__ == "__main__":
    binary, work = map(Path, sys.argv[1:])
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=work) as run:
        main(binary, Path(run))
