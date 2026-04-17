#!/usr/bin/env python3
"""
Export progressive trajectory snapshots from a VTP polyline file.

For z range [z1, z2] computed from input point coordinates and k steps:
  step_size = (z2 - z1) / k
  thresholds = z1 + i * step_size, i = 1..k

At each threshold z_t, trajectories are clipped to the portion with z <= z_t.
If a segment crosses z_t, an interpolated boundary point is inserted.
Each step is written to a separate .vtp file.
"""

from __future__ import annotations

import argparse
import os
from typing import List, Sequence, Tuple

import vtk

TARGET_COLORIDS = {2, 4, 8}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read trajectory polylines from a VTP and export progressive "
            "partial trajectories for increasing z thresholds."
        )
    )
    parser.add_argument("--input", required=True, help="Input trajectory .vtp file.")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output step .vtp files.",
    )
    parser.add_argument(
        "--prefix",
        default="tracking_step",
        help='Output filename prefix (default: "tracking_step").',
    )
    parser.add_argument(
        "--steps",
        type=int,
        required=True,
        help="Number of progressive steps (k).",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-10,
        help="Numerical tolerance for z comparisons (default: 1e-10).",
    )
    parser.add_argument(
        "--all-trajectories",
        action="store_true",
        help=(
            "Process all trajectories without ColorId pre-extraction. "
            f"By default only ColorId in {sorted(TARGET_COLORIDS)} are processed."
        ),
    )
    return parser.parse_args()


def read_polydata(path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    poly = reader.GetOutput()
    if poly is None:
        raise RuntimeError(f"Failed to read VTP: {path}")
    return poly


def write_polydata(poly: vtk.vtkPolyData, path: str) -> None:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(poly)
    writer.SetDataModeToBinary()
    ok = writer.Write()
    if ok != 1:
        raise RuntimeError(f"Failed to write VTP: {path}")


def compute_z_range(poly: vtk.vtkPolyData) -> Tuple[float, float]:
    pts = poly.GetPoints()
    if pts is None or pts.GetNumberOfPoints() == 0:
        raise ValueError("Input has no points; cannot compute z range.")
    z_min = float("inf")
    z_max = float("-inf")
    for i in range(pts.GetNumberOfPoints()):
        z = float(pts.GetPoint(i)[2])
        if z < z_min:
            z_min = z
        if z > z_max:
            z_max = z
    return z_min, z_max


def get_colorid_cell_array(poly: vtk.vtkPolyData) -> vtk.vtkDataArray | None:
    cell_data = poly.GetCellData()
    if cell_data is None:
        return None
    for name in ("ColorId", "colorId", "color_id"):
        arr = cell_data.GetArray(name)
        if arr is not None:
            return arr
    return None


def colorid_value_as_int(arr: vtk.vtkDataArray | None, cell_id: int) -> int:
    if arr is None:
        return -1
    try:
        return int(arr.GetTuple1(cell_id))
    except Exception:
        return -1


def point_interp(
    p0: Sequence[float], p1: Sequence[float], z_target: float, tol: float
) -> Tuple[float, float, float]:
    z0 = float(p0[2])
    z1 = float(p1[2])
    dz = z1 - z0
    if abs(dz) <= tol:
        # Nearly horizontal in z; fallback to midpoint projected to z_target.
        return (
            0.5 * (float(p0[0]) + float(p1[0])),
            0.5 * (float(p0[1]) + float(p1[1])),
            z_target,
        )
    t = (z_target - z0) / dz
    t = max(0.0, min(1.0, t))
    return (
        float(p0[0]) + t * (float(p1[0]) - float(p0[0])),
        float(p0[1]) + t * (float(p1[1]) - float(p0[1])),
        z_target,
    )


def almost_same_point(a: Sequence[float], b: Sequence[float], tol: float) -> bool:
    return (
        abs(float(a[0]) - float(b[0])) <= tol
        and abs(float(a[1]) - float(b[1])) <= tol
        and abs(float(a[2]) - float(b[2])) <= tol
    )


def append_point_if_new(
    line_pts: List[Tuple[float, float, float]], p: Tuple[float, float, float], tol: float
) -> None:
    if not line_pts or not almost_same_point(line_pts[-1], p, tol):
        line_pts.append(p)


def clip_polyline_below_z(
    pts: Sequence[Tuple[float, float, float]], z_th: float, tol: float
) -> List[List[Tuple[float, float, float]]]:
    """
    Return one or more contiguous sub-polylines that satisfy z <= z_th.
    Crossing points are linearly interpolated.
    """
    if len(pts) < 2:
        return []

    out_lines: List[List[Tuple[float, float, float]]] = []
    current: List[Tuple[float, float, float]] = []

    for i in range(len(pts) - 1):
        p0 = pts[i]
        p1 = pts[i + 1]
        z0 = float(p0[2])
        z1 = float(p1[2])
        below0 = z0 <= z_th + tol
        below1 = z1 <= z_th + tol

        if below0 and below1:
            append_point_if_new(current, p0, tol)
            append_point_if_new(current, p1, tol)
            continue

        if below0 and not below1:
            # Leaving the kept side: close at intersection.
            append_point_if_new(current, p0, tol)
            p_cross = point_interp(p0, p1, z_th, tol)
            append_point_if_new(current, p_cross, tol)
            if len(current) >= 2:
                out_lines.append(current)
            current = []
            continue

        if (not below0) and below1:
            # Entering the kept side: start at intersection.
            p_cross = point_interp(p0, p1, z_th, tol)
            current = []
            append_point_if_new(current, p_cross, tol)
            append_point_if_new(current, p1, tol)
            continue

        # Both points are above z_th: no kept geometry on this segment.
        if len(current) >= 2:
            out_lines.append(current)
        current = []

    if len(current) >= 2:
        out_lines.append(current)
    return out_lines


def polyline_cell_point_ids(poly: vtk.vtkPolyData, cell_id: int) -> List[int]:
    cell = poly.GetCell(cell_id)
    n = cell.GetNumberOfPoints()
    return [cell.GetPointId(i) for i in range(n)]


def extract_lines(
    poly: vtk.vtkPolyData,
    z_th: float,
    tol: float,
    selected_colorids: set[int] | None,
) -> vtk.vtkPolyData:
    in_points = poly.GetPoints()
    if in_points is None:
        out = vtk.vtkPolyData()
        out.SetPoints(vtk.vtkPoints())
        out.SetLines(vtk.vtkCellArray())
        return out

    out_pts = vtk.vtkPoints()
    out_lines = vtk.vtkCellArray()
    out_colorid = vtk.vtkIntArray()
    out_colorid.SetName("ColorId")
    in_colorid = get_colorid_cell_array(poly)

    n_cells = poly.GetNumberOfCells()
    for cid in range(n_cells):
        if poly.GetCellType(cid) != vtk.VTK_POLY_LINE and poly.GetCellType(cid) != vtk.VTK_LINE:
            continue
        cell_colorid = colorid_value_as_int(in_colorid, cid)
        if selected_colorids is not None and cell_colorid not in selected_colorids:
            continue
        ids = polyline_cell_point_ids(poly, cid)
        if len(ids) < 2:
            continue
        trajectory = [tuple(in_points.GetPoint(pid)) for pid in ids]
        clipped_parts = clip_polyline_below_z(trajectory, z_th, tol)

        for part in clipped_parts:
            if len(part) < 2:
                continue
            line_ids = vtk.vtkIdList()
            for p in part:
                new_pid = out_pts.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
                line_ids.InsertNextId(new_pid)
            out_lines.InsertNextCell(line_ids)
            out_colorid.InsertNextValue(cell_colorid)

    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_poly.SetLines(out_lines)
    out_poly.GetCellData().AddArray(out_colorid)
    out_poly.GetCellData().SetActiveScalars("ColorId")
    return out_poly


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input not found: {args.input}")
    if args.steps <= 0:
        raise ValueError("--steps must be positive.")

    os.makedirs(args.output_dir, exist_ok=True)

    poly = read_polydata(args.input)
    z1, z2 = compute_z_range(poly)
    step_size = (z2 - z1) / float(args.steps)

    print(f"Input: {args.input}")
    print(f"Output dir: {args.output_dir}")
    print(f"z-range (from input): [{z1}, {z2}], steps={args.steps}, step_size={step_size}")
    selected_colorids: set[int] | None = None if args.all_trajectories else TARGET_COLORIDS
    if selected_colorids is None:
        print("ColorId: processing all trajectories (no pre-extraction).")
    else:
        print(
            "ColorId: processing only trajectories with ColorId in "
            f"{sorted(selected_colorids)}; values are copied to output cell data."
        )

    for i in range(1, args.steps + 1):
        z_th = z1 + i * step_size
        out_poly = extract_lines(poly, z_th, args.tol, selected_colorids)
        out_name = f"{args.prefix}_{i:04d}.vtp"
        out_path = os.path.join(args.output_dir, out_name)
        write_polydata(out_poly, out_path)
        print(
            f"[{i:04d}/{args.steps:04d}] z_th={z_th:.12g}, "
            f"lines={out_poly.GetNumberOfLines()}, points={out_poly.GetNumberOfPoints()} -> {out_path}"
        )


if __name__ == "__main__":
    main()
