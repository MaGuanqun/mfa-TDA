#!/usr/bin/env python3
"""
Assign RegionId from a VTP to each CSV row by nearest 3D Euclidean distance.

- VTP point coordinates: (x, y, z) — z is the third dimension (e.g. time).
- CSV columns: PositionX, PositionY, PositionZ (or detected equivalents).

For every CSV point, find the closest VTP point in 3D and copy its RegionId.
No separate time vs spatial rules — one 3D distance only.
"""

from __future__ import annotations

import argparse
import os
from typing import Sequence, Tuple

import numpy as np
import pandas as pd
import vtk
from scipy.spatial import cKDTree
from vtk.util.numpy_support import vtk_to_numpy


def _require_columns(columns: Sequence[str], required: Sequence[str]) -> None:
    missing = [c for c in required if c not in columns]
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")


def _pick_y_column(columns: Sequence[str]) -> str:
    for cand in ("PositionY", "Position", "Y", "y"):
        if cand in columns:
            return cand
    raise ValueError(
        "Could not detect Y column. Expected one of: PositionY, Position, Y, y"
    )


def _pick_x_column(columns: Sequence[str]) -> str:
    for cand in ("PositionX", "X", "x"):
        if cand in columns:
            return cand
    raise ValueError("Could not detect X column. Expected PositionX/X/x.")


def _pick_z_column(columns: Sequence[str]) -> str:
    for cand in ("PositionZ", "t", "T", "time", "Time", "z", "Z"):
        if cand in columns:
            return cand
    raise ValueError(
        "Could not detect Z (3rd dim) column. Expected PositionZ or similar."
    )


def read_vtp_xyz_region(
    vtp_path: str, region_array_name: str = "RegionId"
) -> Tuple[np.ndarray, np.ndarray]:
    """Returns vtp_xyz (N,3) float64, region (N,) int64."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    poly = reader.GetOutput()

    vtk_pts = poly.GetPoints()
    if vtk_pts is None:
        raise ValueError(f"VTP has no points: {vtp_path}")

    pts = vtk_to_numpy(vtk_pts.GetData()).astype(np.float64, copy=False)
    if pts.ndim != 2 or pts.shape[1] < 3:
        raise ValueError(
            f"VTP points must be 3D. Got shape {pts.shape} from: {vtp_path}"
        )
    pts = pts[:, :3]

    reg_arr = poly.GetPointData().GetArray(region_array_name)
    if reg_arr is None:
        raise ValueError(
            f"VTP point-data array '{region_array_name}' not found in: {vtp_path}"
        )
    region = vtk_to_numpy(reg_arr).astype(np.int64, copy=False)
    if region.shape[0] != pts.shape[0]:
        raise ValueError(
            f"Region array length ({region.shape[0]}) != #points ({pts.shape[0]})."
        )
    return pts, region


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assign RegionId by nearest 3D point in VTP (Euclidean distance)."
    )
    p.add_argument("--input-vtp", "-i", required=True, help="Input VTP with RegionId.")
    p.add_argument("--input-csv", "-c", required=True, help="Input CSV points.")
    p.add_argument("--output-csv", "-o", required=True, help="Output CSV with RegionId.")
    p.add_argument(
        "--region-array",
        default="RegionId",
        help="RegionId array name in VTP point data (default: RegionId).",
    )
    p.add_argument(
        "--max-3d-distance",
        type=float,
        default=None,
        help="If set, rows with nearest 3D distance above this get RegionId=-1.",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=-1,
        help="Workers for scipy KD-tree query (-1 = all CPUs, 1 = single-thread).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_3d_distance is not None and args.max_3d_distance < 0:
        raise ValueError("--max-3d-distance must be >= 0 if set.")

    vtp_xyz, vtp_region = read_vtp_xyz_region(
        args.input_vtp, region_array_name=args.region_array
    )
    tree = cKDTree(vtp_xyz)

    df = pd.read_csv(args.input_csv)
    cols = list(df.columns)
    x_col = _pick_x_column(cols)
    y_col = _pick_y_column(cols)
    z_col = _pick_z_column(cols)
    _require_columns(cols, [x_col, y_col, z_col])

    csv_xyz = df[[x_col, y_col, z_col]].to_numpy(dtype=np.float64, copy=False)

    # scipy>=1.6: query(..., workers=-1) uses all CPUs; older scipy has no workers kw
    try:
        w = args.workers if args.workers is not None else -1
        dists, nn_idx = tree.query(csv_xyz, k=1, workers=w)
    except TypeError:
        dists, nn_idx = tree.query(csv_xyz, k=1)

    nn_idx = np.asarray(nn_idx, dtype=np.int64).reshape(-1)
    dists = np.asarray(dists, dtype=np.float64).reshape(-1)
    assigned_region = vtp_region[nn_idx].astype(np.int64, copy=False)

    if args.max_3d_distance is not None:
        far = dists > args.max_3d_distance
        if np.any(far):
            assigned_region = assigned_region.copy()
            assigned_region[far] = -1

    out_df = df.copy()
    out_df["RegionId"] = assigned_region
    out_df["NearestDist3D"] = dists

    out_dir = os.path.dirname(args.output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)

    unmatched_mask = out_df["RegionId"].to_numpy(dtype=np.int64, copy=False) == -1
    unmatched_count = int(np.sum(unmatched_mask))
    matched = int(np.sum(assigned_region >= 0))

    print(f"VTP points:     {vtp_xyz.shape[0]}")
    print(f"Input CSV rows: {df.shape[0]}")
    print(f"Matched rows:   {matched}")
    print(f"Unmatched rows: {df.shape[0] - matched}")
    print(f"Has RegionId=-1: {unmatched_count > 0}")
    if unmatched_count > 0:
        print("Unmatched rows detail (RegionId == -1):")
        print(out_df.loc[unmatched_mask].to_string(index=False))
    print(f"Output CSV:     {args.output_csv}")


if __name__ == "__main__":
    main()
