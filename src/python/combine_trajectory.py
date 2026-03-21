#!/usr/bin/env pvpython
"""
Combine a range of trajectory .vtp files into a single .vtp.

Each input .vtp is assumed to be a VTK PolyData containing both:
  - points (Verts/Points)
  - edges (Lines/PolyLines)

We preserve all point/cell attribute arrays by concatenating the PolyData
objects without cleaning/merging points.

Example:
  pvpython src/python/combine_trajectory.py \
    --dir /path/to/vtps \
    --start 0 --end 1 \
    --pattern "trajectory_index{i}.vtp" \
    --output /path/to/vtps/trajectory_combined.vtp
"""

from __future__ import annotations

import argparse
import os
from typing import List

import vtk


def build_input_paths(input_dir: str, pattern: str, start: int, end: int) -> List[str]:
    step = 1 if end >= start else -1
    paths: List[str] = []
    for i in range(start, end + step, step):
        fname = pattern.format(i=i)
        paths.append(os.path.join(input_dir, fname))
    return paths


def read_polydata(vtp_path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    return reader.GetOutput()


def combine_polydata(polydata_list: List[vtk.vtkPolyData]) -> vtk.vtkPolyData:
    if not polydata_list:
        raise ValueError("No input PolyData provided.")

    first = polydata_list[0]

    append = vtk.vtkAppendPolyData()
    for poly in polydata_list:
        append.AddInputData(poly)

    append.Update()
    out = append.GetOutput()

    # Preserve active scalar metadata (if present) from the first input.
    first_pd = first.GetPointData()
    if first_pd is not None:
        scalars = first_pd.GetScalars()
        if scalars is not None:
            out.GetPointData().SetActiveScalars(scalars.GetName())

    first_cd = first.GetCellData()
    if first_cd is not None:
        scalars = first_cd.GetScalars()
        if scalars is not None:
            out.GetCellData().SetActiveScalars(scalars.GetName())

    return out


def ensure_edge_values_array(
    out_poly: vtk.vtkPolyData,
    polydata_list: List[vtk.vtkPolyData],
    edge_values_array: str,
) -> None:
    """
    Ensure the combined output contains a per-cell "edge values" array.

    In this repo, `extract_trajectories_from_vtp.py` expects a cell-data array
    (default name: `EdgeValues`) with one tuple per VTK cell (typically VTK_LINE edges).
    """
    out_cd = out_poly.GetCellData()
    if out_cd is None:
        raise RuntimeError("Output polydata has no cell data.")

    expected_total_cells = sum(p.GetNumberOfCells() for p in polydata_list)
    out_arr = out_cd.GetArray(edge_values_array)
    if out_arr is not None and out_arr.GetNumberOfTuples() == expected_total_cells:
        # Already present and sized correctly.
        out_cd.SetActiveScalars(edge_values_array)
        return

    # Find a template array (so we preserve dtype and number of components).
    template_arr = None
    for p in polydata_list:
        a = p.GetCellData().GetArray(edge_values_array)
        if a is not None:
            template_arr = a
            break

    if template_arr is None:
        raise ValueError(
            f"Could not find cell data array '{edge_values_array}' in any input VTP."
        )

    num_components = template_arr.GetNumberOfComponents()
    zeros = [0.0] * num_components

    # Build a new array with correct size.
    new_arr = template_arr.NewInstance()
    new_arr.SetName(edge_values_array)
    new_arr.SetNumberOfComponents(num_components)
    new_arr.SetNumberOfTuples(expected_total_cells)

    offset = 0
    for p in polydata_list:
        cd = p.GetCellData()
        in_arr = cd.GetArray(edge_values_array) if cd is not None else None

        n_cells = p.GetNumberOfCells()
        for j in range(n_cells):
            if in_arr is None:
                new_arr.SetTuple(offset + j, zeros)
            else:
                new_arr.SetTuple(offset + j, in_arr.GetTuple(j))
        offset += n_cells

    # Replace any existing array with the rebuilt one.
    try:
        out_cd.RemoveArray(edge_values_array)
    except Exception:
        # RemoveArray may not exist in some VTK builds; fall back to overwrite via AddArray.
        pass

    out_cd.AddArray(new_arr)
    out_cd.SetActiveScalars(edge_values_array)


def write_vtp(polydata: vtk.vtkPolyData, output_path: str, data_mode: str) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    if data_mode.lower() == "binary":
        writer.SetDataModeToBinary()
    elif data_mode.lower() == "ascii":
        writer.SetDataModeToAscii()
    else:
        raise ValueError("--data-mode must be one of: binary, ascii")

    # Keep arrays exactly as produced by the appender.
    writer.SetInputData(polydata)
    writer.Write()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Combine trajectory_index{i}.vtp files into one VTP.")
    p.add_argument("--dir", required=True, help="Directory containing the trajectory .vtp files.")
    p.add_argument(
        "--start",
        required=True,
        type=int,
        help="Start index (inclusive), e.g. 0 for trajectory_index0.vtp.",
    )
    p.add_argument(
        "--end",
        required=True,
        type=int,
        help="End index (inclusive), e.g. 1 for trajectory_index1.vtp.",
    )
    p.add_argument(
        "--pattern",
        default="trajectory_index{i}.vtp",
        help='Filename pattern with "{i}" placeholder (default: trajectory_index{i}.vtp).',
    )
    p.add_argument(
        "--output",
        required=True,
        help="Output .vtp file path for the combined result.",
    )
    p.add_argument(
        "--data-mode",
        default="binary",
        choices=["binary", "ascii"],
        help="Output VTP data mode (default: binary).",
    )
    p.add_argument(
        "--edge-values-array",
        default="EdgeValues",
        help="Cell-data array name storing per-edge values (default: EdgeValues).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    input_paths = build_input_paths(args.dir, args.pattern, args.start, args.end)
    missing = [p for p in input_paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Missing input VTP file(s): {missing}")

    polydata_list = [read_polydata(p) for p in input_paths]
    out_poly = combine_polydata(polydata_list)
    ensure_edge_values_array(out_poly, polydata_list, args.edge_values_array)
    write_vtp(out_poly, args.output, args.data_mode)

    print(f"Wrote combined VTP: {args.output}")
    print(f"Input files combined: {len(polydata_list)}")


if __name__ == "__main__":
    main()

