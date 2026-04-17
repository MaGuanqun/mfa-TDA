#!/usr/bin/env python3
"""
Convert a 4D-encoded VTP (x, y, z with point array `t`) into shifted 3D points.

For each point:
    (x, y, z) -> (x + 3*t, y + 3*t, z + 3*t)

The output preserves:
  - topology/cells (including edges/lines)
  - all point data arrays (including `t`)
  - cell and field data
"""

from __future__ import annotations

import argparse
import os

import vtk


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read a VTP with point array 't', shift points by +3*t in x/y/z, "
            "and write a new VTP while preserving arrays and edges."
        )
    )
    parser.add_argument("-i", "--input-vtp", required=True, help="Input .vtp file")
    parser.add_argument("-o", "--output-vtp", required=True, help="Output .vtp file")
    parser.add_argument(
        "--t-array-name",
        default="t",
        help="Point-data array name for 4th dimension (default: t)",
    )
    parser.add_argument(
        "--data-mode",
        choices=["binary", "ascii"],
        default="binary",
        help="VTP writer mode (default: binary)",
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


def write_polydata(poly: vtk.vtkPolyData, path: str, mode: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    if mode == "ascii":
        writer.SetDataModeToAscii()
    else:
        writer.SetDataModeToBinary()
    writer.SetInputData(poly)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write VTP: {path}")


def convert_points_with_t(poly: vtk.vtkPolyData, t_array_name: str) -> vtk.vtkPolyData:
    point_data = poly.GetPointData()
    if point_data is None:
        raise RuntimeError("Input polydata has no point data.")

    t_array = point_data.GetArray(t_array_name)
    if t_array is None:
        available = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]
        raise ValueError(
            f"Point-data array '{t_array_name}' not found. Available arrays: {available}"
        )

    if t_array.GetNumberOfComponents() < 1:
        raise ValueError(f"Array '{t_array_name}' has no components.")

    n_points = poly.GetNumberOfPoints()
    if t_array.GetNumberOfTuples() != n_points:
        raise ValueError(
            f"Tuple count mismatch: points={n_points}, "
            f"{t_array_name} tuples={t_array.GetNumberOfTuples()}"
        )

    # Deep-copy to preserve all existing arrays/cells/topology.
    out_poly = vtk.vtkPolyData()
    out_poly.DeepCopy(poly)

    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(n_points)

    for pid in range(n_points):
        x, y, z = poly.GetPoint(pid)
        t = float(t_array.GetTuple1(pid))
        shift = 3.0 * t
        new_points.SetPoint(pid, x + shift, y + shift, z + shift)

    out_poly.SetPoints(new_points)
    return out_poly


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_vtp):
        raise FileNotFoundError(f"Input file not found: {args.input_vtp}")

    poly = read_polydata(args.input_vtp)
    out_poly = convert_points_with_t(poly, args.t_array_name)
    write_polydata(out_poly, args.output_vtp, args.data_mode)

    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {out_poly.GetNumberOfPoints()}")
    print(f"Cells: {out_poly.GetNumberOfCells()}")


if __name__ == "__main__":
    main()
