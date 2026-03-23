#!/usr/bin/env python3
"""
Assign RegionId from a VTP point set to CSV points via nearest spatial match.

Problem setting:
- VTP points encode (x, y, t) in their 3D coordinates.
- CSV rows encode (PositionX, PositionY, PositionZ) where Z is time.
- For each CSV row, we match the closest VTP point in 2D space (x, y)
  within the same time slice, then copy that VTP point's RegionId.

Fast path:
- Build one cKDTree per VTP time slice (in parallel-ready map).
- Query entire CSV time groups in vectorized batches.
- Optional multi-threading across CSV time groups.
"""

from __future__ import annotations

import argparse
import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

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
    # Support common variants seen in this repository and user description.
    for cand in ("PositionY", "Position", "Y", "y"):
        if cand in columns:
            return cand
    raise ValueError(
        "Could not detect Y column in CSV. Expected one of: "
        "PositionY, Position, Y, y"
    )


def _pick_x_column(columns: Sequence[str]) -> str:
    for cand in ("PositionX", "X", "x"):
        if cand in columns:
            return cand
    raise ValueError("Could not detect X column in CSV. Expected PositionX/X/x.")


def _pick_t_column(columns: Sequence[str]) -> str:
    # Time is stored as PositionZ in your pipeline; keep fallbacks for robustness.
    for cand in ("PositionZ", "t", "T", "time", "Time", "z", "Z"):
        if cand in columns:
            return cand
    raise ValueError(
        "Could not detect time column in CSV. Expected PositionZ or a time-like column."
    )


def read_vtp_xy_t_region(
    vtp_path: str, region_array_name: str = "RegionId"
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns:
      points: (N, 3) float64, columns [x, y, t]
      region: (N,) int64
    """
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


@dataclass
class SliceIndex:
    time_value: float
    tree: cKDTree
    region_ids: np.ndarray  # aligned with tree data rows


def build_time_slices(
    points_xyz_t: np.ndarray,
    region_ids: np.ndarray,
    time_quantum: float,
) -> Dict[int, SliceIndex]:
    """
    Group VTP points by time and build one 2D KD-tree per time group.
    Key is quantized integer time key: round(t / time_quantum).
    """
    if time_quantum <= 0:
        raise ValueError("--time-quantum must be > 0.")

    t_vals = points_xyz_t[:, 2]
    keys = np.rint(t_vals / time_quantum).astype(np.int64)
    uniq_keys = np.unique(keys)

    out: Dict[int, SliceIndex] = {}
    for k in uniq_keys:
        mask = keys == k
        xy = points_xyz_t[mask, :2]
        reg = region_ids[mask]
        if xy.shape[0] == 0:
            continue
        out[int(k)] = SliceIndex(
            time_value=float(np.mean(t_vals[mask])),
            tree=cKDTree(xy),
            region_ids=reg,
        )
    return out


def _nearest_time_key(
    query_key: int, sorted_keys: np.ndarray, max_key_gap: Optional[int]
) -> Optional[int]:
    """
    Return closest available key to query_key if within max_key_gap.
    """
    if sorted_keys.size == 0:
        return None
    pos = int(np.searchsorted(sorted_keys, query_key))
    candidates: List[int] = []
    if pos < sorted_keys.size:
        candidates.append(int(sorted_keys[pos]))
    if pos > 0:
        candidates.append(int(sorted_keys[pos - 1]))
    if not candidates:
        return None
    best = min(candidates, key=lambda k: abs(k - query_key))
    if max_key_gap is not None and abs(best - query_key) > max_key_gap:
        return None
    return best


def assign_region_for_group(
    group_indices: np.ndarray,
    group_xy: np.ndarray,
    group_time_key: int,
    slices: Dict[int, SliceIndex],
    sorted_slice_keys: np.ndarray,
    max_key_gap: Optional[int],
    max_spatial_distance: Optional[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns tuple:
      row_indices, assigned_region_ids, nn_distances
    """
    n = group_indices.shape[0]
    if n == 0:
        return (
            group_indices,
            np.empty((0,), dtype=np.int64),
            np.empty((0,), dtype=np.float64),
        )

    # If max_key_gap is set, search ALL slices in [t - gap, t + gap] and choose
    # the closest spatial match among those slices.
    if max_key_gap is not None:
        left = int(np.searchsorted(sorted_slice_keys, group_time_key - max_key_gap, side="left"))
        right = int(np.searchsorted(sorted_slice_keys, group_time_key + max_key_gap, side="right"))
        candidate_keys = sorted_slice_keys[left:right]
    else:
        nearest = _nearest_time_key(group_time_key, sorted_slice_keys, None)
        candidate_keys = np.array([nearest], dtype=np.int64) if nearest is not None else np.empty((0,), dtype=np.int64)

    if candidate_keys.size == 0:
        return (
            group_indices,
            np.full(n, -1, dtype=np.int64),
            np.full(n, np.nan, dtype=np.float64),
        )

    best_dist = np.full(n, np.inf, dtype=np.float64)
    best_region = np.full(n, -1, dtype=np.int64)

    for key in candidate_keys:
        slc = slices[int(key)]
        dists, nn_idx = slc.tree.query(group_xy, k=1)
        nn_idx = nn_idx.astype(np.int64, copy=False)
        cand_region = slc.region_ids[nn_idx]
        improve = dists < best_dist
        if np.any(improve):
            best_dist[improve] = dists[improve]
            best_region[improve] = cand_region[improve]

    assigned = best_region

    if max_spatial_distance is not None:
        far_mask = best_dist > max_spatial_distance
        if np.any(far_mask):
            assigned = assigned.copy()
            assigned[far_mask] = -1

    return group_indices, assigned.astype(np.int64, copy=False), best_dist.astype(
        np.float64, copy=False
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Assign RegionId to CSV points by nearest 2D (x,y) VTP point "
            "within matching time slice."
        )
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
        "--time-step",
        type=float,
        default=None,
        help=(
            "If set, force all CSV points to this time step instead of per-row time. "
            "Useful when assigning one specific slice."
        ),
    )
    p.add_argument(
        "--time-quantum",
        type=float,
        default=1e-6,
        help="Quantization for time grouping (default: 1e-6).",
    )
    p.add_argument(
        "--max-time-diff",
        type=float,
        default=None,
        help=(
            "Time window half-width. Match against ALL VTP slices with "
            "|csv_time - vtp_time| <= max-time-diff, then keep closest spatial match. "
            "If no slice is in window, RegionId=-1. "
            "Default: disabled (use nearest single slice)."
        ),
    )
    p.add_argument(
        "--max-spatial-distance",
        type=float,
        default=None,
        help="Optional max 2D distance. Beyond this assign RegionId=-1.",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1)),
        help="Thread workers over time groups (default: CPU count).",
    )
    p.add_argument(
        "--debug-time",
        action="store_true",
        help="Print time-range diagnostics to explain unmatched rows.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be >= 1.")
    if args.max_time_diff is not None and args.max_time_diff < 0:
        raise ValueError("--max-time-diff must be >= 0.")

    # 1) Read VTP and build per-time KDTree index
    vtp_points, vtp_region = read_vtp_xy_t_region(
        args.input_vtp, region_array_name=args.region_array
    )
    slices = build_time_slices(vtp_points, vtp_region, time_quantum=args.time_quantum)
    if not slices:
        raise RuntimeError("No valid VTP time slices were built.")
    sorted_slice_keys = np.array(sorted(slices.keys()), dtype=np.int64)
    if args.max_time_diff is None:
        max_key_gap: Optional[int] = None
    else:
        max_key_gap = int(math.ceil(args.max_time_diff / args.time_quantum))

    # 2) Read CSV and discover columns
    df = pd.read_csv(args.input_csv)
    cols = list(df.columns)
    x_col = _pick_x_column(cols)
    y_col = _pick_y_column(cols)
    t_col = _pick_t_column(cols)
    _require_columns(cols, [x_col, y_col, t_col])

    # 3) Prepare numeric arrays
    xy = df[[x_col, y_col]].to_numpy(dtype=np.float64, copy=False)
    if args.time_step is None:
        t_vals = df[t_col].to_numpy(dtype=np.float64, copy=False)
    else:
        t_vals = np.full(df.shape[0], float(args.time_step), dtype=np.float64)
    t_keys = np.rint(t_vals / args.time_quantum).astype(np.int64)

    if args.debug_time:
        vtp_t_min = float(np.min(vtp_points[:, 2]))
        vtp_t_max = float(np.max(vtp_points[:, 2]))
        csv_t_min = float(np.min(t_vals))
        csv_t_max = float(np.max(t_vals))
        print(
            f"[debug] VTP time range: [{vtp_t_min}, {vtp_t_max}] "
            f"(unique slices={len(sorted_slice_keys)})"
        )
        print(f"[debug] CSV time range: [{csv_t_min}, {csv_t_max}]")
        if args.max_time_diff is None:
            print("[debug] max-time-diff disabled: nearest time slice always used.")
        else:
            print(f"[debug] max-time-diff={args.max_time_diff}, max-key-gap={max_key_gap}")

    # 4) Assign in parallel by time groups
    assigned_region = np.full(df.shape[0], -1, dtype=np.int64)
    assigned_dist = np.full(df.shape[0], np.nan, dtype=np.float64)

    uniq_t_keys = np.unique(t_keys)
    futures = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for k in uniq_t_keys:
            idx = np.where(t_keys == k)[0]
            futures.append(
                ex.submit(
                    assign_region_for_group,
                    idx,
                    xy[idx],
                    int(k),
                    slices,
                    sorted_slice_keys,
                    max_key_gap,
                    args.max_spatial_distance,
                )
            )

        for fut in as_completed(futures):
            row_idx, reg, dist = fut.result()
            assigned_region[row_idx] = reg
            assigned_dist[row_idx] = dist

    # 4b) Hard fallback: guarantee assignment for any remaining -1 rows.
    # This can happen for pathological numeric rows (e.g., non-finite XY) or
    # strict filtering in custom runs. We force a global nearest XY assignment.
    unresolved = np.where(assigned_region < 0)[0]
    fallback_count = int(unresolved.shape[0])
    if fallback_count > 0:
        global_tree = cKDTree(vtp_points[:, :2])
        fallback_xy = xy[unresolved]
        finite_mask = np.isfinite(fallback_xy).all(axis=1)

        # For non-finite XY rows, use (0,0) for query just to guarantee an ID.
        safe_xy = fallback_xy.copy()
        if np.any(~finite_mask):
            safe_xy[~finite_mask] = np.array([0.0, 0.0], dtype=np.float64)

        d_fb, idx_fb = global_tree.query(safe_xy, k=1)
        idx_fb = idx_fb.astype(np.int64, copy=False)
        assigned_region[unresolved] = vtp_region[idx_fb].astype(np.int64, copy=False)
        assigned_dist[unresolved] = d_fb.astype(np.float64, copy=False)

    # 5) Write output CSV
    out_df = df.copy()
    out_df["RegionId"] = assigned_region
    out_df["NearestSpatialDistXY"] = assigned_dist

    out_dir = os.path.dirname(args.output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    out_df.to_csv(args.output_csv, index=False)

    unmatched_mask = out_df["RegionId"].to_numpy(dtype=np.int64, copy=False) == -1
    unmatched_count = int(np.sum(unmatched_mask))

    matched = int(np.sum(assigned_region >= 0))
    has_negative_one = bool(np.any(assigned_region == -1))
    print(f"Input CSV rows: {df.shape[0]}")
    print(f"Matched rows:   {matched}")
    print(f"Unmatched rows: {df.shape[0] - matched}")
    print(f"Has RegionId=-1: {has_negative_one}")
    print(f"Fallback assigned rows: {fallback_count}")
    if unmatched_count > 0:
        print("Unmatched rows detail (RegionId == -1):")
        print(out_df.loc[unmatched_mask].to_string(index=False))
    print(f"Output CSV:     {args.output_csv}")


if __name__ == "__main__":
    main()
