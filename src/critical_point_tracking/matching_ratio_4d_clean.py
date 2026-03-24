import argparse
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def load_csv_points(csv_path, has_header=True):
    skip = 1 if has_header else 0
    pts = np.loadtxt(csv_path, delimiter=",", skiprows=skip, usecols=(0, 1, 2))
    pts = np.atleast_2d(pts)
    return pts.astype(float, copy=False)


def load_vtp_points_and_t(vtp_path, t_array_name="t"):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()

    polydata = reader.GetOutput()
    vtk_points = polydata.GetPoints()
    if vtk_points is None:
        raise ValueError(f"No points found in VTP: {vtp_path}")

    xyz = vtk_to_numpy(vtk_points.GetData()).astype(float, copy=False)
    if xyz.shape[1] < 3:
        raise ValueError(f"VTP points do not contain x,y,z columns: {vtp_path}")
    xyz = xyz[:, :3]

    point_data = polydata.GetPointData()
    t_arr = point_data.GetArray(t_array_name)
    if t_arr is None:
        for i in range(point_data.GetNumberOfArrays()):
            name = point_data.GetArrayName(i)
            if name is not None and name.lower() == "t":
                t_arr = point_data.GetArray(i)
                break
    if t_arr is None:
        available = []
        for i in range(point_data.GetNumberOfArrays()):
            name = point_data.GetArrayName(i)
            if name is not None:
                available.append(name)
        raise ValueError(
            f"Cannot find point-data array '{t_array_name}' in {vtp_path}. "
            f"Available arrays: {available}"
        )

    t = vtk_to_numpy(t_arr).astype(float, copy=False).reshape(-1)
    if t.shape[0] != xyz.shape[0]:
        raise ValueError(
            f"Point count mismatch in {vtp_path}: xyz={xyz.shape[0]}, t={t.shape[0]}"
        )
    return xyz, t


def build_indexed_csv_resolver(csv_pattern):
    if "{index}" in csv_pattern:
        return lambda idx: csv_pattern.format(index=idx)

    stem, ext = os.path.splitext(csv_pattern)
    if ext.lower() != ".csv":
        raise ValueError("csv_pattern must end with .csv or contain '{index}'.")

    m = re.match(r"^(.*)_(\d+)$", stem)
    prefix = m.group(1) if m else stem
    return lambda idx: f"{prefix}_{idx}.csv"


def index_to_time(index, t_min, t_max, num_knots):
    return t_min + (t_max - t_min) * (index / (num_knots - 1))


def build_tree(points):
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("scipy is required for matching. Please install scipy.") from exc
    return cKDTree(points) if points.size else None


def match_with_time_and_space(
    src_xyz, src_t, tgt_xyz, tgt_t, spatial_threshold, time_threshold, time_eps
):
    n_src = src_xyz.shape[0]
    if n_src == 0 or tgt_xyz.size == 0:
        return np.zeros(n_src, dtype=bool)

    tree = build_tree(tgt_xyz)
    neighbors = tree.query_ball_point(src_xyz, r=spatial_threshold)
    matched = np.zeros(n_src, dtype=bool)

    for i, ids in enumerate(neighbors):
        if not ids:
            continue
        if np.any(np.abs(tgt_t[ids] - src_t[i]) <= (time_threshold + time_eps)):
            matched[i] = True
    return matched


def load_all_csv_points(csv_pattern, has_header, t_min, t_max, num_knots):
    csv_path_for_idx = build_indexed_csv_resolver(csv_pattern)
    all_xyz = []
    all_t = []
    loaded = 0
    missing = 0

    for idx in range(num_knots):
        path = csv_path_for_idx(idx)
        if not os.path.exists(path):
            missing += 1
            continue
        pts = load_csv_points(path, has_header=has_header)
        if pts.size == 0:
            continue
        t_val = index_to_time(idx, t_min=t_min, t_max=t_max, num_knots=num_knots)
        all_xyz.append(pts)
        all_t.append(np.full(pts.shape[0], t_val, dtype=float))
        loaded += 1

    if not all_xyz:
        raise FileNotFoundError("No CSV slices were loaded from the given csv_pattern.")

    xyz = np.vstack(all_xyz)
    t = np.concatenate(all_t)
    print(f"Loaded CSV slices: {loaded}, missing: {missing}, total CSV points: {xyz.shape[0]}")
    return xyz, t


def compute_match_ratio_4d(
    vtp_path,
    csv_pattern,
    time_threshold,
    spatial_threshold,
    csv_has_header,
    t_min,
    t_max,
    num_knots,
    t_array_name,
    time_eps,
):
    if time_threshold < 0 or spatial_threshold < 0:
        raise ValueError("time_threshold and spatial_threshold must be non-negative.")
    if num_knots < 2 or t_max <= t_min:
        raise ValueError("Invalid knot/time configuration.")

    vtp_xyz, vtp_t = load_vtp_points_and_t(vtp_path, t_array_name=t_array_name)
    csv_xyz, csv_t = load_all_csv_points(
        csv_pattern=csv_pattern,
        has_header=csv_has_header,
        t_min=t_min,
        t_max=t_max,
        num_knots=num_knots,
    )

    m_vtp = match_with_time_and_space(
        vtp_xyz, vtp_t, csv_xyz, csv_t, spatial_threshold, time_threshold, time_eps
    )
    matched_vtp = int(np.count_nonzero(m_vtp))
    total_vtp = int(vtp_xyz.shape[0])
    ratio_vtp = (matched_vtp / total_vtp) if total_vtp else 0.0

    m_csv = match_with_time_and_space(
        csv_xyz, csv_t, vtp_xyz, vtp_t, spatial_threshold, time_threshold, time_eps
    )
    matched_csv = int(np.count_nonzero(m_csv))
    total_csv = int(csv_xyz.shape[0])
    ratio_csv = (matched_csv / total_csv) if total_csv else 0.0

    print("VTP(second) result on CSV(first) result:")
    print(f"Matched: {matched_vtp} / {total_vtp}")
    print(f"Ratio:   {ratio_vtp:.6f}")
    print("CSV(first) result on VTP(second) result:")
    print(f"Matched: {matched_csv} / {total_csv}")
    print(f"Ratio:   {ratio_csv:.6f}")
    return ratio_vtp, ratio_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="4D matching: |t1-t2|<s1 and ||xyz1-xyz2||<s2."
    )
    parser.add_argument("vtp_file", type=str, help="Path to VTP file with xyz+t.")
    parser.add_argument(
        "csv_pattern",
        type=str,
        help="CSV pattern: name_{index}.csv, name_58.csv, or name.csv.",
    )
    parser.add_argument("time_threshold", type=float, help="s1: |dt| threshold.")
    parser.add_argument("spatial_threshold", type=float, help="s2: xyz distance threshold.")
    parser.add_argument(
        "--time-eps",
        type=float,
        default=1e-12,
        help="Tolerance for time compare: |dt| <= s1 + eps.",
    )
    parser.add_argument("--header", action="store_true", help="CSV files have header row.")
    parser.add_argument("--t-min", type=float, default=0.0, help="Minimum t value.")
    parser.add_argument("--t-max", type=float, default=89.0, help="Maximum t value.")
    parser.add_argument("--num-knots", type=int, default=161, help="Number of time knots.")
    parser.add_argument(
        "--t-array-name",
        type=str,
        default="t",
        help="VTP point-data array name storing t.",
    )
    args = parser.parse_args()

    compute_match_ratio_4d(
        vtp_path=args.vtp_file,
        csv_pattern=args.csv_pattern,
        time_threshold=args.time_threshold,
        spatial_threshold=args.spatial_threshold,
        csv_has_header=args.header,
        t_min=args.t_min,
        t_max=args.t_max,
        num_knots=args.num_knots,
        t_array_name=args.t_array_name,
        time_eps=args.time_eps,
    )
