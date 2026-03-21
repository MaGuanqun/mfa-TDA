#!/usr/bin/env python3
"""
For each point P in a VTP file (with point-data array "t"):
1) Find the two neighboring regular-grid indices A/B on [grid_min, grid_max]
   with grid_count points.
2) Load CSV files for indices A and B (3D points + CriticalType).
3) In each CSV, find the nearest point to P within epsilon.
4) If both exist, keep the closer one.

Output one CSV containing extracted matched points with columns:
  PositionX, PositionY, PositionZ, t, CriticalType, RegionId

Memory-aware behavior:
- CSV/KDTree data are loaded on demand and cached with LRU eviction.
- At most --max-cached-indices index datasets are kept in memory.
"""

import argparse
import csv
import math
import os
from collections import OrderedDict

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


def read_vtp_points(vtp_path, t_array_name="t", region_array_name="RegionId"):
    """
    Read VTP points and per-point t value.
    Returns list of dicts: {'x','y','z','t','region_id','point_id'}
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    poly = reader.GetOutput()

    pts = poly.GetPoints()
    if pts is None:
        return []

    point_data = poly.GetPointData()
    t_arr = point_data.GetArray(t_array_name)
    if t_arr is None:
        raise ValueError(f"VTP point-data array '{t_array_name}' not found.")

    region_arr = point_data.GetArray(region_array_name)
    if region_arr is None:
        raise ValueError(
            f"VTP point-data array '{region_array_name}' not found. "
            f"Compute it first (e.g. via 'assign_component_id_to_point.py')."
        )

    out = []
    npts = pts.GetNumberOfPoints()
    for i in range(npts):
        x, y, z = pts.GetPoint(i)
        t_val = float(t_arr.GetComponent(i, 0))
        region_id = int(region_arr.GetComponent(i, 0))
        out.append(
            {
                "point_id": i,
                "x": float(x),
                "y": float(y),
                "z": float(z),
                "t": t_val,
                "region_id": region_id,
            }
        )
    return out


def load_csv_index_data(csv_path):
    """
    Load one CSV file and build KDTree.
    Returns dict with keys: coords, crit, tree
    """
    coords = []
    crit_types = []

    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV file has no header.")

        if "CriticalType" not in reader.fieldnames:
            raise ValueError("CSV must contain a 'CriticalType' column.")

        x_col, y_col, z_col = find_csv_coord_columns(reader.fieldnames)

        for row in reader:
            coords.append([float(row[x_col]), float(row[y_col]), float(row[z_col])])
            crit_types.append(int(float(row["CriticalType"])))

    if not coords:
        return {"coords": np.empty((0, 3), dtype=float), "crit": np.empty((0,), dtype=int), "tree": None}

    coord_arr = np.asarray(coords, dtype=float)
    crit_arr = np.asarray(crit_types, dtype=int)
    tree = cKDTree(coord_arr)
    return {"coords": coord_arr, "crit": crit_arr, "tree": tree}


def grid_step(grid_min, grid_max, grid_count):
    if grid_count < 2:
        raise ValueError("--grid-count must be >= 2.")
    return (grid_max - grid_min) / float(grid_count - 1)


def bracket_indices(t_val, gmin, gmax, gstep, gcount):
    """Return neighboring grid indices (A, B) with A <= B."""
    if t_val <= gmin:
        return 0, 0
    if t_val >= gmax:
        last = gcount - 1
        return last, last

    pos = (t_val - gmin) / gstep
    left = int(math.floor(pos))
    right = int(math.ceil(pos))
    left = max(0, min(gcount - 1, left))
    right = max(0, min(gcount - 1, right))
    return left, right


class IndexDataCache:
    """LRU cache for per-index CSV KDTree data."""

    def __init__(self, csv_template, max_cached_indices):
        self.csv_template = csv_template
        self.max_cached_indices = max_cached_indices
        self.cache = OrderedDict()
        self.missing = set()

    def get(self, idx):
        if idx in self.cache:
            self.cache.move_to_end(idx)
            return self.cache[idx]
        if idx in self.missing:
            return None

        csv_path = self.csv_template.format(idx=idx)
        if not os.path.exists(csv_path):
            self.missing.add(idx)
            return None

        data = load_csv_index_data(csv_path)
        self.cache[idx] = data
        self.cache.move_to_end(idx)

        while len(self.cache) > self.max_cached_indices:
            self.cache.popitem(last=False)
        return data


def nearest_within_eps(index_data, point_xyz, epsilon):
    """Return (distance, matched_xyz, matched_critical_type) or None."""
    if index_data is None or index_data["tree"] is None:
        return None

    dist, idx = index_data["tree"].query(point_xyz, distance_upper_bound=epsilon)
    if not np.isfinite(dist):
        return None
    if idx >= index_data["coords"].shape[0]:
        return None

    matched_xyz = index_data["coords"][int(idx)]
    matched_type = int(index_data["crit"][int(idx)])
    return float(dist), matched_xyz, matched_type


def match_points(points_vtp, cache, epsilon, gmin, gmax, gcount, log_matches=True):
    """
    For each VTP point, query bracketing grid indices A/B and keep the closer
    valid nearest match (within epsilon). Returns output rows.
    """
    rows = []
    gstep = grid_step(gmin, gmax, gcount)

    for p in points_vtp:
        pxyz = [p["x"], p["y"], p["z"]]
        a_idx, b_idx = bracket_indices(p["t"], gmin, gmax, gstep, gcount)

        a_data = cache.get(a_idx)
        a_match = nearest_within_eps(a_data, pxyz, epsilon)

        if b_idx == a_idx:
            b_match = None
        else:
            b_data = cache.get(b_idx)
            b_match = nearest_within_eps(b_data, pxyz, epsilon)

        chosen = None
        chosen_idx = None
        if a_match is not None and b_match is not None:
            if a_match[0] <= b_match[0]:
                chosen = a_match
                chosen_idx = a_idx
            else:
                chosen = b_match
                chosen_idx = b_idx
        elif a_match is not None:
            chosen = a_match
            chosen_idx = a_idx
        elif b_match is not None:
            chosen = b_match
            chosen_idx = b_idx

        if chosen is None:
            continue

        _, matched_xyz, matched_type = chosen
        rows.append(
            {
                "PositionX": float(matched_xyz[0]),
                "PositionY": float(matched_xyz[1]),
                "PositionZ": float(matched_xyz[2]),
                "t": float(gmin + chosen_idx * gstep),
                "CriticalType": matched_type,
                # RegionId comes from the VTP point, not from the matched CSV point.
                "RegionId": int(p["region_id"]),
            }
        )
        if log_matches:
            print(
                "Matched "
                f"vtp_id={p['point_id']} "
                f"grid_idx={chosen_idx} "
                f"pos=({matched_xyz[0]:.6f},{matched_xyz[1]:.6f},{matched_xyz[2]:.6f}) "
                f"t={gmin + chosen_idx * gstep:.6f} "
                f"CriticalType={matched_type}",
                flush=True,
            )

    return rows


def write_single_csv(rows, output_csv):
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "PositionX",
                "PositionY",
                "PositionZ",
                "t",
                "CriticalType",
                "RegionId",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="4D matching: for each VTP point, query bracketing grid-index CSVs and save closest match."
    )
    parser.add_argument("-i", "--input-vtp", required=True, help="Input VTP file")
    parser.add_argument(
        "-c",
        "--csv-template",
        required=True,
        help="CSV path template containing {idx}, e.g. /path/mfa_8_{idx}.csv",
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        required=True,
        help="Output CSV file path",
    )
    parser.add_argument(
        "-e",
        "--epsilon",
        type=float,
        required=True,
        help="Distance threshold for nearest neighbor acceptance",
    )
    parser.add_argument("--t-array", default="t", help="VTP point-data t array name")
    parser.add_argument(
        "--region-array",
        default="RegionId",
        help="VTP point-data RegionId array name (default: RegionId).",
    )
    parser.add_argument("--grid-min", type=float, default=0.0, help="Grid min value")
    parser.add_argument("--grid-max", type=float, default=89.0, help="Grid max value")
    parser.add_argument("--grid-count", type=int, default=81, help="Number of grid points")
    parser.add_argument(
        "--max-cached-indices",
        type=int,
        default=4,
        help="Max number of index CSV datasets cached in memory",
    )
    parser.add_argument(
        "--no-match-log",
        action="store_true",
        help="Disable printing each matched point during runtime",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.epsilon < 0:
        raise ValueError("--epsilon must be non-negative.")
    if args.max_cached_indices < 1:
        raise ValueError("--max-cached-indices must be >= 1.")

    points_vtp = read_vtp_points(
        args.input_vtp, t_array_name=args.t_array, region_array_name=args.region_array
    )
    cache = IndexDataCache(args.csv_template, args.max_cached_indices)
    matched_rows = match_points(
        points_vtp,
        cache,
        args.epsilon,
        args.grid_min,
        args.grid_max,
        args.grid_count,
        log_matches=(not args.no_match_log),
    )
    write_single_csv(matched_rows, args.output_csv)

    total = len(matched_rows)
    print(f"VTP points: {len(points_vtp)}")
    print(f"Matched points: {total}")
    print(f"Output CSV: {args.output_csv}")


if __name__ == "__main__":
    main()