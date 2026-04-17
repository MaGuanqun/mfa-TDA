#!/usr/bin/env python3
"""
Find common points between:
  A) a VTP file of points
  B) a CSV file of points (with CriticalType)

For each point in A, this script finds at most one nearest point in B
within Euclidean distance <= d (no CriticalType constraint for matching).

Matched B points are written to a single CSV file.
"""

import argparse
import csv
import os

import vtk
import numpy as np
from scipy.spatial import cKDTree


def find_csv_coord_columns(fieldnames):
    """Return coordinate column names from a CSV header."""
    # Prefer TTK-style names first.
    ttk_candidates = ("PositionX", "PositionY", "PositionZ")
    if all(name in fieldnames for name in ttk_candidates):
        return ttk_candidates

    # Fallback to x0/x1/x2-style names.
    generic_candidates = ("x0", "x1", "x2")
    if all(name in fieldnames for name in generic_candidates):
        return generic_candidates

    raise ValueError(
        "Could not detect coordinate columns in CSV. "
        "Expected either (PositionX, PositionY, PositionZ) or (x0, x1, x2)."
    )


def read_vtp_points(vtp_path):
    """
    Read VTP points and optional per-point CriticalType.
    Returns list of dicts: {'x','y','z','critical_type','point_id'}
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    poly = reader.GetOutput()

    pts = poly.GetPoints()
    if pts is None:
        return []

    point_data = poly.GetPointData()
    crit_arr = point_data.GetArray("CriticalType")

    out = []
    npts = pts.GetNumberOfPoints()
    for i in range(npts):
        x, y, z = pts.GetPoint(i)
        ctype = None
        if crit_arr is not None:
            ctype = int(crit_arr.GetComponent(i, 0))
        out.append(
            {
                "point_id": i,
                "x": float(x),
                "y": float(y),
                "z": float(z),
                "critical_type": ctype,
            }
        )
    return out


def read_csv_points(csv_path):
    """
    Read CSV points with optional CriticalType.
    Returns:
      - list of dicts with parsed coordinates and original row payload
      - original CSV fieldnames (in source order)
    """
    points = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV file has no header.")
        csv_fieldnames = list(reader.fieldnames)

        x_col, y_col, z_col = find_csv_coord_columns(reader.fieldnames)

        for idx, row in enumerate(reader):
            critical_type = None
            if "CriticalType" in row and row["CriticalType"] not in (None, ""):
                critical_type = int(float(row["CriticalType"]))

            points.append(
                {
                    "row_id": idx,
                    "x": float(row[x_col]),
                    "y": float(row[y_col]),
                    "z": float(row[z_col]),
                    "critical_type": critical_type,
                    "original_row": dict(row),
                }
            )
    return points, csv_fieldnames


def match_points(points_a, points_b, threshold, z_k=None, z_eps=1e-6):
    """
    For each point in B, find whether it is within `threshold` of any point in A.
    Matching is distance-only; no CriticalType filter is applied.
    If `z_k` is provided, points in B are additionally filtered to those whose
    `z` coordinate is within `z_eps` of `z_k`.

    Output includes all matched B points. Since we iterate B exactly once, each
    B row_id appears at most once in the output.
    """
    if not points_a or not points_b:
        return []

    a_coords = np.array([[pa["x"], pa["y"], pa["z"]] for pa in points_a], dtype=float)
    tree = cKDTree(a_coords)
    matches = []

    for pb in points_b:
        if z_k is not None:
            if abs(pb["z"] - z_k) > z_eps:
                continue

        dist, idx = tree.query(
            [pb["x"], pb["y"], pb["z"]], distance_upper_bound=threshold
        )
        # `cKDTree.query(..., distance_upper_bound=...)` returns `idx` even when
        # out of bounds; we explicitly require a finite distance.
        if not np.isfinite(dist) or idx >= len(points_a):
            continue

        matches.append({"RowId": pb["row_id"], **pb["original_row"]})

    return matches


def write_single_csv(rows, output_csv, csv_fieldnames):
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["RowId"] + [name for name in csv_fieldnames if name != "RowId"],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Find points in VTP close to points in CSV."
    )
    parser.add_argument("-a", "--input-vtp", required=True, help="Input VTP file A")
    parser.add_argument("-b", "--input-csv", required=True, help="Input CSV file B")
    parser.add_argument(
        "-d", "--distance", type=float, required=True, help="Distance threshold"
    )
    parser.add_argument(
        "--k",
        type=float,
        default=None,
        help="Optional Z (PositionZ) value; only keep CSV points with z ~= k",
    )
    parser.add_argument(
        "--k-eps",
        type=float,
        default=1e-3,
        help="Epsilon for Z filtering when --k is set (absolute tolerance)",
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        required=True,
        help="Output CSV file path",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.distance < 0:
        raise ValueError("--distance must be non-negative.")
    if args.k_eps < 0:
        raise ValueError("--k-eps must be non-negative.")

    points_a = read_vtp_points(args.input_vtp)
    points_b, csv_fieldnames = read_csv_points(args.input_csv)
    matched_rows = match_points(
        points_a, points_b, args.distance, z_k=args.k, z_eps=args.k_eps
    )
    write_single_csv(matched_rows, args.output_csv, csv_fieldnames)

    total = len(matched_rows)
    print(f"A points: {len(points_a)}")
    print(f"B points: {len(points_b)}")
    print(f"Matched points: {total}")
    print(f"Output CSV: {args.output_csv}")


if __name__ == "__main__":
    main()