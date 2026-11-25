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
        assert report["step"] == report["fiber_step"]
        seed_report = json.loads((run / "seeds_report.json").read_text())
        shape = seed_report["discovery_grid_shape"]
        assert seed_report["starting_points"] == shape[0] * shape[1] * shape[2]
        assert report["discovery_grid_shape"] == shape
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

    def curve_pipeline(run, program, model, options, closed):
        stage(1, run, ["--model", model, *options])
        for step in [2, 3, 4, 5]:
            stage(step, run)
        mesh = (run / "surface.obj").read_bytes()
        prefix = args.work / (run.name + "_reference")
        invoke(program, [*options, "--output", prefix])
        assert mesh == prefix.with_suffix(".obj").read_bytes()
        assert (run / "surface.vtk").read_bytes() == prefix.with_suffix(".vtk").read_bytes()
        report = json.loads((run / "report.json").read_text())
        assert report["closed"] is closed and report["complete"] and len(report["range_vertices"]) == 3
        stage(4, run)
        stage(5, run)
        assert mesh == (run / "surface.obj").read_bytes()

    if args.fixture:
        # Copy only the small synthetic test model; never mutate real datasets.
        model = args.work / "sphere model.mfa"
        shutil.copyfile(args.fixture, model)
        common = ["--first", model, "--variable1", "0", "--variable2", "1",
                  "--a", "1,-2", "--b", "1,0.5", "--min", "-2,-2,-2", "--max", "2,2,2"]
        run = args.work / "mfa sphere"
        mesh, report = pipeline(run, ["--model", "mfa", *common])
        assert report["components"] == 1 and report["fiber_components"] == 1
        assert report["discovery_grid_shape"] == [7, 7, 7]
        prefix = args.work / "mfa_reference"
        invoke("fiber_surface_mfa", [*common, "--output", prefix])
        assert mesh == prefix.with_suffix(".obj").read_bytes(), "Staged MFA differs from in-memory extraction"
        anisotropic = args.work / "mfa anisotropic"
        stage(1, anisotropic, ["--model", "mfa", "--first", model,
                              "--a", "1,-2", "--b", "1,0.5", "--seed-grid", "5",
                              "--min", "-1,-2,-2", "--max", "1,2,2"])
        seed_report = json.loads((anisotropic / "seeds_report.json").read_text())
        assert seed_report["discovery_grid_shape"] == [5, 10, 10]
        # Old checkpoints can only resume if their traced and marched geometry use the same step.
        config_path = anisotropic / "run.conf"
        original_config = config_path.read_text()
        config_lines = original_config.splitlines()
        for i, line in enumerate(config_lines):
            if line.startswith("options "):
                fields = line.split()
                fields[3] = str(float(fields[1]) / 2)
                config_lines[i] = " ".join(fields)
        try:
            config_path.write_text("\n".join(config_lines) + "\n")
            assert "rerun the seeds stage" in stage(2, anisotropic, expected=1)
        finally:
            config_path.write_text(original_config)
        assert seed_report["starting_points"] == 500
        curve_pipeline(args.work / "mfa elbow", "fiber_surface_mfa", "mfa",
                       ["--first", model, "--variable1", "0", "--variable2", "1",
                        "--polyline", "1,-2;1,0.5;2,0.5", "--min", "-2,-2,-2", "--max", "2,2,2",
                        "--step", "0.1", "--residual-tolerance", "0.005"],
                       closed=False)
        original = model.read_bytes()
        try:
            model.write_bytes(original + b"changed")
            assert "MFA input changed" in stage(5, run, expected=1)
        finally:
            model.write_bytes(original)
    else:
        anisotropic = args.work / "anisotropic"
        stage(1, anisotropic, ["--example", "strip", "--seed-grid", "5",
                               "--min", "-1,-2,-2", "--max", "1,2,2"])
        seed_report = json.loads((anisotropic / "seeds_report.json").read_text())
        assert seed_report["discovery_grid_shape"] == [5, 10, 10]
        # Shape metadata is checked before a later stage accepts the seeds.
        seed_path = anisotropic / "surface_seeds.dat"
        original_seed = seed_path.read_text()
        malformed = original_seed.splitlines()
        header = malformed[0].split()
        header[6] = "6"
        malformed[0] = " ".join(header)
        try:
            seed_path.write_text("\n".join(malformed) + "\n")
            assert "grid shape" in stage(2, anisotropic, expected=1)
        finally:
            seed_path.write_text(original_seed)
        # Existing version-1 discoveries remain readable.
        for name in ["surface_seeds.dat", "endpoint_seeds.dat"]:
            path = anisotropic / name
            lines = path.read_text().splitlines()
            header = lines[0].split()
            header[1] = "1"
            del header[6:9]
            lines[0] = " ".join(header)
            path.write_text("\n".join(lines) + "\n")
        stage(2, anisotropic)
        for name in ["sphere", "closed-sphere", "strip", "disconnected"]:
            run = args.work / name
            options = ["--example", name]
            if name == "strip":
                options += ["--fiber-seed-radius", "0.12"]
            mesh, report = pipeline(run, ["--model", "analytical", *options])
            assert report["fiber_seed_radius"] == (0.12 if name == "strip" else report["spatial_tolerance"])
            assert report["components"] == (2 if name == "disconnected" else 1)
            prefix = args.work / (name + "_reference")
            invoke("fiber_surface_explicit", [*options, "--output", prefix])
            assert mesh == prefix.with_suffix(".obj").read_bytes(), f"Staged {name} changed geometry"
        run = args.work / "triangle loop"
        options = ["--example", "strip", "--polyline", "0,0;1,0;0,1;0,0"]
        curve_pipeline(run, "fiber_surface_explicit", "analytical", options, closed=True)
        # Patch zero must detect changed corner seeds owned by patch one.
        neighbor = run / "segment_1" / "endpoint_seeds.dat"
        original = neighbor.read_bytes()
        try:
            neighbor.write_bytes(original + b"\n")
            error = stage(5, run, expected=1)
            assert "segment_0/fibers.state" in error and "upstream output changed" in error
        finally:
            neighbor.write_bytes(original)
        for path in ["0,0", "0,0;1,0;", "0,0;1,1;0,1;1,0", "0,0;nan,0"]:
            invoke("fiber_surface_explicit", ["--polyline", path], expected=1)
        invoke("fiber_surface_explicit", [*options, "--a", "0,0", "--b", "1,0"], expected=1)
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
