#!/usr/bin/env pvpython
"""
Combine trajectory .vtp files into a single .vtp.

Each input .vtp is assumed to be a VTK PolyData containing both:
  - points (Verts/Points)
  - edges (Lines/PolyLines)

We preserve all point/cell attribute arrays by concatenating the PolyData
objects without cleaning/merging points.

Examples:
  pvpython src/python/combine_trajectory.py \
    --dir /path/to/vtps \
    --start 0 --end 1 \
    --pattern "trajectory_index{i}.vtp" \
    --output /path/to/vtps/trajectory_combined.vtp

  pvpython src/python/combine_trajectory.py \
    --dir /path/to/vtps \
    --indices "0,2,5-8,10" \
    --pattern "trajectory_index{i}.vtp" \
    --output /path/to/vtps/trajectory_selected.vtp
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


def parse_indices_string(indices: str) -> List[int]:
    """
    Parse a comma-separated index specification into an ordered list of integers.

    Supported tokens:
      - single index: "4"
      - ascending range: "2-6" -> 2,3,4,5,6
      - descending range: "6-2" -> 6,5,4,3,2
    """
    if not indices or not indices.strip():
        raise ValueError("--indices cannot be empty.")

    parsed: List[int] = []
    for raw_token in indices.split(","):
        token = raw_token.strip()
        if not token:
            raise ValueError(f"Invalid empty token in --indices: {indices!r}")

        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
                raise ValueError(
                    f"Invalid range token {token!r} in --indices. Use forms like 3-7."
                )
            start = int(parts[0].strip())
            end = int(parts[1].strip())
            step = 1 if end >= start else -1
            parsed.extend(range(start, end + step, step))
        else:
            parsed.append(int(token))

    return parsed


def build_input_paths_from_indices(input_dir: str, pattern: str, indices: List[int]) -> List[str]:
    if not indices:
        raise ValueError("No indices provided to build input paths.")
    return [os.path.join(input_dir, pattern.format(i=i)) for i in indices]


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
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--indices",
        type=str,
        help=(
            'Comma-separated indices/ranges, e.g. "0,2,5-8,10". '
            "Use this to directly choose files."
        ),
    )
    mode.add_argument(
        "--start",
        type=int,
        help="Start index (inclusive), e.g. 0 for trajectory_index0.vtp.",
    )
    p.add_argument(
        "--end",
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

    if args.indices is not None:
        indices = parse_indices_string(args.indices)
        input_paths = build_input_paths_from_indices(args.dir, args.pattern, indices)
    else:
        if args.end is None:
            raise ValueError("--end is required when using --start/--end mode.")
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

