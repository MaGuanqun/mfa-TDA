"""Produce the paper's Ours and TTK rows for one CoordNet range triangle."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy

from join_isabel_refinements import read_ply, write_ply
from sample_coordnet_grid import read_grid, run_parameters

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_isabel_grid import grid_selection

BINARY = Path(__file__).resolve().parents[3] / "build/src/fiber_surface/fiber_grid_recovery"


def mesh_mask(xyz, faces, image):
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(xyz, deep=True))
    cells = vtk.vtkCellArray()
    cells.SetData(3, numpy_to_vtkIdTypeArray(np.ascontiguousarray(faces, dtype=np.int64).ravel(), deep=True))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(cells)
    stencil = vtk.vtkPolyDataToImageStencil()
    stencil.SetInputData(mesh)
    stencil.SetOutputOrigin(image.GetOrigin())
    stencil.SetOutputSpacing(image.GetSpacing())
    stencil.SetOutputWholeExtent(image.GetExtent())
    stencil.SetTolerance(0.)
    mask = vtk.vtkImageStencilToImage()
    mask.SetInputConnection(stencil.GetOutputPort())
    mask.SetInsideValue(1)
    mask.SetOutsideValue(0)
    mask.SetOutputScalarTypeToUnsignedChar()
    mask.Update()
    return vtk_to_numpy(mask.GetOutput().GetPointData().GetScalars()).reshape(image.GetDimensions()[::-1]).astype(bool)


def audit(stages, mesh, prefix):
    subprocess.run([BINARY, "audit", stages, mesh, prefix], check=True)
    report = json.loads(prefix.with_suffix(".json").read_text())
    topology = {k: int(v) for k, v in re.findall(r"(\w+)=(-?\d+)", report["topology"])}
    return report, topology


def measure(method, path, run, dest, grid, image, target):
    xyz, faces, labels, _ = read_ply(path)
    row = dict(method=method, mesh=str(path.resolve()), triangles=len(faces),
               open_edges=0, pinches=0, intersections=0, max_vertex_range_residual=None,
               max_vertex_overshoot=None, d95=None, iou=None,
               side_vertices_checked=0,
               projection_samples=0, failed_projections=0, excluded_degenerate_faces=0,
               both=None, fitted_only=None, mesh_only=None, neither=None, iou_note="")
    if not len(faces):
        row["iou_note"] = "Empty surface"
        return row

    prefix = dest / method.lower()
    report, counts = audit(run / "stages", path, prefix)
    row.update(open_edges=counts["boundary_edges"], pinches=counts["nonmanifold_vertices"],
               max_vertex_range_residual=report["max_vertex_range_residual"],
               max_vertex_overshoot=report["max_vertex_overshoot"],
               side_vertices_checked=report["side_vertices_checked"],
               projection_samples=report["samples"], failed_projections=report["failed_projections"])
    if report["samples"] == 20000 and report["failed_projections"] == 0:
        row["d95"] = report["projection_p95"]

    row["intersections"] = counts["intersections"] if report["intersections_checked"] else None
    if not report["intersections_checked"] and counts["degenerate_faces"]:
        # Same cross-product cutoff as validate_mesh(). Keep the full mesh for every other column.
        a = xyz[faces[:, 1]] - xyz[faces[:, 0]]
        b = xyz[faces[:, 2]] - xyz[faces[:, 0]]
        cutoff = grid["step"] * 1e-8 * np.maximum(np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1))
        keep = np.linalg.norm(np.cross(a, b), axis=1) > cutoff
        excluded = int((~keep).sum())
        if excluded != counts["degenerate_faces"]:
            raise ValueError("Degenerate face count disagrees with the full audit")
        row["excluded_degenerate_faces"] = excluded
        if keep.any():
            subset = dest / (method.lower() + "_nondegenerate")
            write_ply(subset.with_suffix(".ply"), xyz, faces[keep], labels[keep], (0, 0, 0))
            filtered, topology = audit(run / "stages", subset.with_suffix(".ply"), subset)
            if filtered["intersections_checked"]:
                row["intersections"] = topology["intersections"]

    defects = [key for key in ("boundary_edges", "nonmanifold_edges", "nonmanifold_vertices",
                               "degenerate_faces", "duplicate_faces", "invalid_indices") if counts[key]]
    if defects or row["intersections"] != 0:
        if row["intersections"] != 0:
            defects.append("intersections unchecked" if row["intersections"] is None else "intersections")
        row["iou_note"] = ", ".join(defects)
        return row

    inside = mesh_mask(xyz, faces, image)
    both = int(np.count_nonzero(inside & target))
    missed, excess = int(target.sum()) - both, int(inside.sum()) - both
    union = both + missed + excess
    row.update(both=both, fitted_only=missed, mesh_only=excess, neither=int(target.size) - union)
    if union:
        row["iou"] = both / union
    else:
        row["iou_note"] = "Empty union on the evaluation grid"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="native run directory with stages/run.conf")
    parser.add_argument("grid", type=Path, help="output directory from sample_coordnet_grid.py")
    parser.add_argument("output", type=Path, help="new directory for audits and table.csv")
    parser.add_argument("--model-name", required=True, help="model label for both table rows")
    parser.add_argument("--mesh", type=Path, help="native labeled PLY, default RUN/volume_boundary.ply")
    parser.add_argument("--ttk-python", default=sys.executable, help="Python with the ttkFiberSurface bindings")
    parser.add_argument("--threads", type=int, default=3, help="TTK extraction threads")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    grid, fields = read_grid(args.grid)
    if any(grid[k] != val for k, val in run_parameters(args.run).items()):
        parser.error("Grid and stage configuration differ. Sample the selected run again.")
    dims = grid["dimensions"]
    if not hasattr(vtk, "vtkPolyDataToImageStencil"):
        parser.error("Use a Python with VTK ImagingStencil, and select the TTK interpreter with --ttk-python.")
    native = args.mesh or args.run / "volume_boundary.ply"
    if not native.is_file():
        parser.error(f"Missing native boundary: {native}")

    dest = args.output
    dest.mkdir(parents=True, exist_ok=False)
    image = vtk.vtkImageData()
    image.SetDimensions(*dims)
    image.SetOrigin(*grid["box_min"])
    image.SetSpacing(*grid["spacing"])
    query = np.array(grid["triangle"])
    target = np.empty(dims[::-1], dtype=bool)
    for k in range(dims[2]):
        target[k] = grid_selection(fields[0, k], fields[1, k], query)

    reference = dest / "ttk.ply"
    subprocess.run([args.ttk_python, Path(__file__).with_name("ttk_coordnet_grid.py"),
                    args.grid, reference, "--threads", str(args.threads)], check=True)
    ttk_meta = json.loads((dest / "ttk_extraction.json").read_text())
    meta = dict(grid, model=args.model_name, evaluation_points=target.size,
                fitted_selected=int(target.sum()), stencil_vtk_version=vtk.vtkVersion.GetVTKVersion(), **ttk_meta)
    (dest / "evaluation.json").write_text(json.dumps(meta, indent=2) + "\n")
    rows = [dict(model=args.model_name, **measure(method, path, args.run, dest, grid, image, target))
            for method, path in [("Ours", native), ("TTK", reference)]]
    with (dest / "table.csv.tmp").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (dest / "table.csv.tmp").rename(dest / "table.csv")
    print((dest / "table.csv").read_text(), end="")


if __name__ == "__main__":
    main()
