#!/usr/bin/env python3
"""Cross-process restart, export equivalence, and malformed-state regression tests."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess


STAGES = ["fiber_find_seeds", "fiber_trace_endpoint_fibers", "fiber_march_surface",
          "fiber_stitch_surface", "fiber_evaluate_surface"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--fixture", type=Path)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)

    def invoke(name, arguments, expected=0):
        completed = subprocess.run([str(args.bin / name), *map(str, arguments)],
                                   text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if completed.returncode != expected:
            raise AssertionError(f"{name}: expected {expected}, got {completed.returncode}\n{completed.stdout}")
        return completed.stdout

    def stage(number, run, extras=(), expected=0):
        return invoke(STAGES[number - 1], ["--run-dir", run, *extras], expected)

    def pipeline(run, seed_args, partial=False):
        if run.exists():
            shutil.rmtree(run)
        stage(1, run, seed_args)
        stage(2, run)
        stage(3, run)
        march = run / "marching.state"
        before = march.read_bytes()
        stage(4, run, expected=2 if partial else 0)
        stage(5, run, expected=2 if partial else 0)
        mesh = (run / "surface.obj").read_bytes()
        # Restart only the last two executables. Upstream geometry stays fixed.
        stage(4, run, expected=2 if partial else 0)
        stage(5, run, expected=2 if partial else 0)
        assert before == march.read_bytes(), "Stitching modified its input checkpoint"
        assert mesh == (run / "surface.obj").read_bytes(), "Restart changed mesh geometry"
        report = json.loads((run / "report.json").read_text())
        assert report["complete"] is not partial
        artifacts = ["run.conf", "surface_seeds.dat", "endpoint_seeds.dat",
                     "surface_seeds.vtk", "endpoint_seeds.vtk", "seeds_report.json",
                     "fibers.state", "marching.state", "stitched.state",
                     "surface.obj", "surface.vtk", "surface_fibers.vtk",
                     "surface_report.json", "report.json", "evaluation.conf"]
        for prefix in ["endpoint_fibers", "marching", "stitching"]:
            artifacts.extend(f"{prefix}_{name}" for name in
                             ["mesh.vtk", "fibers.vtk", "fronts.vtk", "report.json"])
        for path in artifacts:
            assert (run / path).is_file(), f"Missing stage artifact: {path}"
        assert not any(path.is_dir() for path in run.iterdir()), "Pipeline created output subdirectories"
        return mesh, report

    if args.fixture:
        # Copy only the small synthetic test model; never mutate real datasets.
        model = args.work / "sphere model.mfa"
        shutil.copyfile(args.fixture, model)
        common = ["--first", model, "--variable1", "0", "--variable2", "1",
                  "--a", "1,-2", "--b", "1,0.5", "--min", "-2,-2,-2", "--max", "2,2,2"]
        run = args.work / "mfa sphere"
        mesh, report = pipeline(run, ["--model", "mfa", *common])
        assert report["components"] == 1 and report["fiber_components"] == 1
        prefix = args.work / "mfa_reference"
        invoke("fiber_surface_mfa", [*common, "--output", prefix])
        assert mesh == prefix.with_suffix(".obj").read_bytes(), "Staged MFA differs from in-memory extraction"
        original = model.read_bytes()
        try:
            model.write_bytes(original + b"changed")
            assert "MFA input changed" in stage(5, run, expected=1)
        finally:
            model.write_bytes(original)
    else:
        for name in ["sphere", "closed-sphere", "strip", "disconnected"]:
            run = args.work / name
            mesh, report = pipeline(run, ["--model", "analytical", "--example", name])
            assert report["components"] == (2 if name == "disconnected" else 1)
            prefix = args.work / (name + "_reference")
            invoke("fiber_surface_explicit", ["--example", name, "--output", prefix])
            assert mesh == prefix.with_suffix(".obj").read_bytes(), f"Staged {name} changed geometry"
        run = args.work / "sphere"
        target = run / "marching.state"
        original = target.read_text()
        corrupted = original.splitlines()
        triangle_header = next(i for i, line in enumerate(corrupted) if line.startswith("triangles "))
        corrupted[triangle_header + 1] = "99999999 0 1"
        try:
            for content in [original.replace("FIBER_SURFACE_STATE 1", "FIBER_SURFACE_STATE 99", 1),
                            original[:len(original) // 2], "\n".join(corrupted) + "\n"]:
                target.write_text(content)
                stage(4, run, expected=1)
        finally:
            target.write_text(original)
        seed_file = run / "endpoint_seeds.dat"
        original_seed = seed_file.read_bytes()
        try:
            seed_file.write_bytes(original_seed + b"\n")
            assert "upstream output changed" in stage(4, run, expected=1)
        finally:
            seed_file.write_bytes(original_seed)
        config = run / "run.conf"
        original_config = config.read_text()
        try:
            config.write_text(original_config.replace('example "sphere"', 'example "closed-sphere"'))
            assert "configuration mismatch" in stage(5, run, expected=1)
        finally:
            config.write_text(original_config)
        missing = args.work / "missing"
        stage(3, missing, expected=1)
        pipeline(args.work / "partial", ["--example", "sphere", "--max-triangles", "4"], partial=True)
    print("Separate processes, exact OBJ equivalence, restarts, and invalid-input checks passed")


if __name__ == "__main__":
    main()
