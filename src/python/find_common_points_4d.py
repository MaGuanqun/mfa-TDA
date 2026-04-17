#!/usr/bin/env python3
"""
For each point P in a VTP file (with point-data array "t"):
1) Choose candidate slice indices from the regular grid on [grid_min, grid_max]
   with grid_count points: by default the two bracketing neighbors of P's `t`; if
   --t-slice-tol is set, every slice whose grid time is within that tolerance of
   P's `t` is included (useful when many time levels may contain the same feature).
2) Load CSV files for each candidate index (3D points + CriticalType).
3) In each CSV, find CSV point(s) in **3D** (PositionX/Y/Z vs VTP x,y,z) within distance
   <= epsilon (same role as --distance in find_common_points.py).
4) Default --match-mode nearest: among all candidate slices, output the single closest
   CSV point. With --match-mode all_within_epsilon: output every CSV row within epsilon
   on every candidate slice (multiple rows per VTP point when many neighbors exist).

The "4D" workflow is **time via slice choice**, not via a 4D distance: epsilon is
**spatial only**; time enters only through which slice CSV(s) are searched (--t-slice-tol
or bracketing).

Output one CSV containing extracted matched points with columns:
  VtpPointId, SliceIndex, PositionX, PositionY, PositionZ, t, CriticalType, ColorId
  RegionId is copied from the VTP only if point-data --region-array exists; otherwise
  the column is omitted. Use --without-region-id to omit it even when present.

ColorId is taken from the VTP point closest in 3D to each output row's position
(same nearest-neighbor idea as assigning VTP fields to CSV points elsewhere).

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


def read_vtp_points(
    vtp_path,
    t_array_name="t",
    region_array_name="RegionId",
    color_array_name="ColorId",
    force_omit_region_id=False,
):
    """
    Read VTP points and per-point t value.

    Returns (points, copy_region_to_csv) where points are dicts with optional 'region_id'
    if the VTP has region_array_name (unless force_omit_region_id is True).
    copy_region_to_csv is True iff output rows should include RegionId.
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    poly = reader.GetOutput()

    pts = poly.GetPoints()
    if pts is None:
        return [], False

    point_data = poly.GetPointData()
    t_arr = point_data.GetArray(t_array_name)
    if t_arr is None:
        raise ValueError(f"VTP point-data array '{t_array_name}' not found.")

    region_arr = None
    if not force_omit_region_id:
        region_arr = point_data.GetArray(region_array_name)

    color_arr = point_data.GetArray(color_array_name)
    if color_arr is None:
        raise ValueError(
            f"VTP point-data array '{color_array_name}' not found. "
            f"Compute it first (e.g. branch coloring / 'seperate_diffferent_branches.py')."
        )

    out = []
    npts = pts.GetNumberOfPoints()
    for i in range(npts):
        x, y, z = pts.GetPoint(i)
        t_val = float(t_arr.GetComponent(i, 0))
        color_id = int(color_arr.GetComponent(i, 0))
        rec = {
            "point_id": i,
            "x": float(x),
            "y": float(y),
            "z": float(z),
            "t": t_val,
            "color_id": color_id,
        }
        if region_arr is not None:
            rec["region_id"] = int(region_arr.GetComponent(i, 0))
        out.append(rec)

    copy_region = region_arr is not None
    return out, copy_region


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


def grid_time_for_index(idx, gmin, gstep):
    return float(gmin + float(idx) * gstep)


def candidate_slice_indices(t_val, gmin, gmax, gstep, gcount, t_slice_tol=None):
    """
    Indices into the CSV template to query for a VTP time `t_val`.

    Default (t_slice_tol is None or negative): bracketing pair (or one index if
    at endpoint / degenerate).

    With t_slice_tol >= 0: all k with |grid_time(k) - t_val| <= t_slice_tol.
    If that set is empty, fall back to bracketing (same as default).
    """
    a, b = bracket_indices(t_val, gmin, gmax, gstep, gcount)
    if t_slice_tol is None or t_slice_tol < 0:
        if a == b:
            return [a]
        return sorted({a, b})

    in_tol = [
        k
        for k in range(gcount)
        if abs(grid_time_for_index(k, gmin, gstep) - float(t_val)) <= t_slice_tol
    ]
    if not in_tol:
        if a == b:
            return [a]
        return sorted({a, b})
    return sorted(set(in_tol))


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
    """
    Return (distance, matched_xyz, matched_critical_type) or None.
    `epsilon` is the cKDTree `distance_upper_bound`: max **3D** Euclidean distance
    from `point_xyz` to CSV coordinates (same idea as find_common_points.py --distance).
    """
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


def all_within_eps(index_data, point_xyz, epsilon):
    """
    Return a list of (distance, matched_xyz, matched_critical_type) for every CSV
    point within **3D** Euclidean distance <= epsilon (empty list if none).
    """
    if index_data is None or index_data["tree"] is None:
        return []
    tree = index_data["tree"]
    coords = index_data["coords"]
    crit = index_data["crit"]
    q = np.asarray(point_xyz, dtype=float)
    idxs = tree.query_ball_point(q, r=float(epsilon), p=2)
    idxs = np.atleast_1d(np.asarray(idxs, dtype=np.intp))
    out = []
    eps = float(epsilon)
    for j in idxs:
        j = int(j)
        if j < 0 or j >= coords.shape[0]:
            continue
        d = float(np.linalg.norm(coords[j] - q))
        if d <= eps + 1e-12:
            out.append((d, coords[j].copy(), int(crit[j])))
    out.sort(key=lambda x: x[0])
    return out


def match_points(
    points_vtp,
    cache,
    epsilon,
    gmin,
    gmax,
    gcount,
    log_matches=True,
    include_region_id=True,
    t_slice_tol=None,
    match_mode="nearest",
):
    """
    For each VTP point, query one or more slice CSVs (bracketing and/or --t-slice-tol).

    match_mode 'nearest': one output row — smallest 3D distance among all slices.
    match_mode 'all_within_epsilon': one row per CSV point within epsilon on each slice.

    Each row's ColorId is the ColorId of the VTP vertex nearest in 3D to the
    matched (PositionX, PositionY, PositionZ); RegionId comes from the driving VTP point.
    """
    rows = []
    gstep = grid_step(gmin, gmax, gcount)
    all_mode = match_mode == "all_within_epsilon"

    if points_vtp:
        vtp_xyz = np.array([[p["x"], p["y"], p["z"]] for p in points_vtp], dtype=float)
        vtp_color = np.array([p["color_id"] for p in points_vtp], dtype=int)
        vtp_tree = cKDTree(vtp_xyz)
    else:
        vtp_tree = None

    def _emit_row(matched_xyz, matched_type, slice_idx, vtp_id, p):
        _, nn_vtp = vtp_tree.query(matched_xyz)
        color_id_out = int(vtp_color[int(nn_vtp)])
        row = {
            "VtpPointId": int(vtp_id),
            "SliceIndex": int(slice_idx),
            "PositionX": float(matched_xyz[0]),
            "PositionY": float(matched_xyz[1]),
            "PositionZ": float(matched_xyz[2]),
            "t": float(grid_time_for_index(slice_idx, gmin, gstep)),
            "CriticalType": matched_type,
            "ColorId": color_id_out,
        }
        if include_region_id:
            row["RegionId"] = int(p["region_id"])
        rows.append(row)
        if log_matches:
            print(
                "Matched "
                f"vtp_id={vtp_id} "
                f"slice={slice_idx} "
                f"pos=({matched_xyz[0]:.6f},{matched_xyz[1]:.6f},{matched_xyz[2]:.6f}) "
                f"t={grid_time_for_index(slice_idx, gmin, gstep):.6f} "
                f"CriticalType={matched_type} "
                f"ColorId={color_id_out}",
                flush=True,
            )

    for p in points_vtp:
        pxyz = [p["x"], p["y"], p["z"]]
        cand = candidate_slice_indices(p["t"], gmin, gmax, gstep, gcount, t_slice_tol)
        vtp_id = p["point_id"]

        if all_mode:
            for k in cand:
                data = cache.get(k)
                for dist, matched_xyz, matched_type in all_within_eps(data, pxyz, epsilon):
                    _emit_row(matched_xyz, matched_type, k, vtp_id, p)
            continue

        chosen = None
        chosen_idx = None
        for k in cand:
            data = cache.get(k)
            m = nearest_within_eps(data, pxyz, epsilon)
            if m is None:
                continue
            if chosen is None or m[0] < chosen[0]:
                chosen = m
                chosen_idx = k

        if chosen is None:
            continue

        _, matched_xyz, matched_type = chosen
        _emit_row(matched_xyz, matched_type, chosen_idx, vtp_id, p)

    return rows


def write_single_csv(rows, output_csv, include_region_id=True):
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fieldnames = [
        "VtpPointId",
        "SliceIndex",
        "PositionX",
        "PositionY",
        "PositionZ",
        "t",
        "CriticalType",
    ]
    if include_region_id:
        fieldnames.append("RegionId")
    fieldnames.append("ColorId")

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "4D pipeline matching: VTP point `t` selects time-slice CSVs; within each slice, "
            "find CSV point(s) within 3D distance <= epsilon. Default: single closest match "
            "across slices; optional: emit every neighbor within epsilon per slice."
        )
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
        help=(
            "Max 3D Euclidean distance (VTP x,y,z vs CSV PositionX/Y/Z) for a match on each "
            "candidate time slice; same meaning as --distance in find_common_points.py. "
            "With --match-mode all_within_epsilon, every CSV point within this radius is output."
        ),
    )
    parser.add_argument(
        "--match-mode",
        choices=("nearest", "all_within_epsilon"),
        default="nearest",
        help=(
            "nearest: one CSV row per VTP point (closest match over candidate slices). "
            "all_within_epsilon: all CSV rows within --epsilon on each candidate slice "
            "(multiple rows per VTP point possible)."
        ),
    )
    parser.add_argument("--t-array", default="t", help="VTP point-data t array name")
    parser.add_argument(
        "--region-array",
        default="RegionId",
        help="VTP point-data RegionId array name when present (default: RegionId).",
    )
    parser.add_argument(
        "--without-region-id",
        action="store_true",
        help="Omit RegionId from the output even if the VTP has --region-array.",
    )
    parser.add_argument(
        "--color-array",
        default="ColorId",
        help="VTP point-data ColorId array name (default: ColorId).",
    )
    parser.add_argument("--grid-min", type=float, default=0.0, help="Grid min value")
    parser.add_argument("--grid-max", type=float, default=89.0, help="Grid max value")
    parser.add_argument("--grid-count", type=int, default=81, help="Number of grid points")
    parser.add_argument(
        "--t-slice-tol",
        type=float,
        default=None,
        help=(
            "If set (>=0), query every slice index whose grid time is within this "
            "distance of the VTP point's t (same units as grid_min/max). "
            "With nearest match-mode, the smallest 3D distance among those slices wins; "
            "with all_within_epsilon, every in-radius CSV point on every such slice is emitted. "
            "If no slice falls in the window, fall back to bracketing neighbors only."
        ),
    )
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
    if args.t_slice_tol is not None and args.t_slice_tol < 0:
        raise ValueError("--t-slice-tol must be non-negative if set.")

    match_mode = args.match_mode

    points_vtp, include_region = read_vtp_points(
        args.input_vtp,
        t_array_name=args.t_array,
        region_array_name=args.region_array,
        color_array_name=args.color_array,
        force_omit_region_id=args.without_region_id,
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
        include_region_id=include_region,
        t_slice_tol=args.t_slice_tol,
        match_mode=match_mode,
    )
    write_single_csv(matched_rows, args.output_csv, include_region_id=include_region)

    total = len(matched_rows)
    print(f"VTP points: {len(points_vtp)}")
    print(f"Matched points: {total}")
    print(f"Output CSV: {args.output_csv}")


if __name__ == "__main__":
    main()