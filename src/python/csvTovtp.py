#!/usr/bin/env python3
"""
Convert CSV point data to VTP (VTK PolyData).

The script uses three CSV columns as point coordinates and writes every other
column as point-data arrays in the output .vtp file.
"""

import argparse
import csv
import math
import os
from typing import Dict, List, Sequence, Tuple

import vtk


def detect_position_columns(fieldnames: Sequence[str]) -> Tuple[str, str, str]:
    """Find XYZ columns from common naming conventions."""
    candidates = [
        ("PositionX", "PositionY", "PositionZ"),
        ("x", "y", "z"),
        ("x0", "x1", "x2"),
    ]
    field_set = set(fieldnames)
    for triplet in candidates:
        if all(name in field_set for name in triplet):
            return triplet
    raise ValueError(
        "Could not detect position columns. Expected one of: "
        "['PositionX','PositionY','PositionZ'], ['x','y','z'], ['x0','x1','x2']."
    )


def can_parse_int(values: Sequence[str]) -> bool:
    for value in values:
        text = value.strip()
        if text == "":
            continue
        try:
            int(text)
        except ValueError:
            return False
    return True


def can_parse_float(values: Sequence[str]) -> bool:
    for value in values:
        text = value.strip()
        if text == "":
            continue
        try:
            float(text)
        except ValueError:
            return False
    return True


def infer_column_type(values: Sequence[str]) -> str:
    if can_parse_int(values):
        return "int"
    if can_parse_float(values):
        return "float"
    return "string"


def add_numeric_array(
    point_data: vtk.vtkPointData,
    name: str,
    values: Sequence[str],
    as_int: bool,
) -> None:
    if as_int:
        arr = vtk.vtkIntArray()
        arr.SetName(name)
        for value in values:
            text = value.strip()
            arr.InsertNextValue(int(text) if text != "" else 0)
    else:
        arr = vtk.vtkDoubleArray()
        arr.SetName(name)
        for value in values:
            text = value.strip()
            arr.InsertNextValue(float(text) if text != "" else math.nan)
    point_data.AddArray(arr)


def add_string_array(point_data: vtk.vtkPointData, name: str, values: Sequence[str]) -> None:
    arr = vtk.vtkStringArray()
    arr.SetName(name)
    for value in values:
        arr.InsertNextValue(value)
    point_data.AddArray(arr)


def convert_csv_to_vtp(input_csv: str, output_vtp: str) -> None:
    with open(input_csv, "r", newline="") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header.")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    if not rows:
        raise ValueError("Input CSV has no data rows.")

    x_col, y_col, z_col = detect_position_columns(fieldnames)

    points = vtk.vtkPoints()
    for row in rows:
        try:
            x = float(row[x_col])
            y = float(row[y_col])
            z = float(row[z_col])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid coordinate in row with values "
                f"{x_col}={row.get(x_col)!r}, {y_col}={row.get(y_col)!r}, {z_col}={row.get(z_col)!r}"
            ) from exc
        points.InsertNextPoint(x, y, z)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    verts = vtk.vtkCellArray()
    for idx in range(len(rows)):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, idx)
        verts.InsertNextCell(vertex)
    polydata.SetVerts(verts)

    point_data = polydata.GetPointData()
    data_columns = [name for name in fieldnames if name not in (x_col, y_col, z_col)]
    column_values: Dict[str, List[str]] = {
        name: [row.get(name, "") for row in rows] for name in data_columns
    }

    for name in data_columns:
        values = column_values[name]
        inferred = infer_column_type(values)
        if inferred == "int":
            add_numeric_array(point_data, name, values, as_int=True)
        elif inferred == "float":
            add_numeric_array(point_data, name, values, as_int=False)
        else:
            add_string_array(point_data, name, values)

    out_dir = os.path.dirname(output_vtp)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_vtp)
    writer.SetInputData(polydata)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write output VTP: {output_vtp}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert CSV points + properties to VTP point-data arrays."
    )
    parser.add_argument("-i", "--input-csv", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output-vtp", required=True, help="Output .vtp file path")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert_csv_to_vtp(args.input_csv, args.output_vtp)
    print(f"Wrote: {args.output_vtp}")
