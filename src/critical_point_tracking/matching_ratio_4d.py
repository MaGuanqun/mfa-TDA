import argparse
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def load_csv_points(csv_path, has_header=True):
    """Load (x, y, z) points from a CSV file."""
    skip = 1 if has_header else 0
    pts = np.loadtxt(csv_path, delimiter=",", skiprows=skip, usecols=(0, 1, 2))
    pts = np.atleast_2d(pts)
    return pts.astype(float, copy=False)


def t_to_knot_position(t_values, t_min, t_max, num_knots):
    """Convert t to floating knot position in [0, num_knots - 1]."""
    if num_knots < 2:
        raise ValueError("num_knots must be >= 2.")
    if t_max <= t_min:
        raise ValueError("t_max must be > t_min.")

    scale = (num_knots - 1) / (t_max - t_min)
    pos = (np.asarray(t_values, dtype=float) - t_min) * scale
    return np.clip(pos, 0.0, float(num_knots - 1))


def get_left_right_indices(knot_positions, num_knots):
    """For each knot position, return (left_idx, right_idx)."""
    left = np.floor(knot_positions).astype(int)
    right = np.ceil(knot_positions).astype(int)
    left = np.clip(left, 0, num_knots - 1)
    right = np.clip(right, 0, num_knots - 1)
    return left, right


def load_vtp_points_and_t(vtp_path, t_array_name="t"):
    """Load xyz coordinates from VTP and t values from point-data array."""
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
    """
    Build resolver for indexed CSV path.
    Supported forms:
      1) '/path/name_{index}.csv'
      2) '/path/name_58.csv' -> '/path/name_{idx}.csv'
      3) '/path/name.csv'    -> '/path/name_{idx}.csv'
    """
    if "{index}" in csv_pattern:
        return lambda idx: csv_pattern.format(index=idx)

    stem, ext = os.path.splitext(csv_pattern)
    if ext.lower() != ".csv":
        raise ValueError("csv_pattern must end with .csv or contain '{index}'.")

    m = re.match(r"^(.*)_(\d+)$", stem)
    prefix = m.group(1) if m else stem
    return lambda idx: f"{prefix}_{idx}.csv"


def build_xy_tree(points):
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("scipy is required for fast matching. Please install scipy.") from exc
    return cKDTree(points[:, :3]) if points.size else None


def batch_match_with_tree(source_pts, target_pts, target_tree, xy_threshold, z_threshold):
    """
    Batched match:
    - find XY neighbors by KD-tree radius query
    - check any neighbor also matches Z threshold
    Returns boolean array per source point.
    """
    n_src = source_pts.shape[0]
    matched = np.zeros(n_src, dtype=bool)
    if n_src == 0 or target_pts.size == 0 or target_tree is None:
        return matched

    neighbors_list = target_tree.query_ball_point(source_pts[:, :2], r=xy_threshold)
    for i, neighbors in enumerate(neighbors_list):
        if not neighbors:
            continue
        if np.any(np.abs(target_pts[neighbors, 2] - source_pts[i, 2]) < z_threshold):
            matched[i] = True
    return matched


def compute_match_ratio_4d(
    vtp_path,
    csv_pattern,
    xy_threshold,
    z_threshold,
    csv_has_header,
    t_min,
    t_max,
    num_knots,
    t_array_name,
):
    vtp_xyz, vtp_t = load_vtp_points_and_t(vtp_path, t_array_name=t_array_name)
    knot_pos = t_to_knot_position(vtp_t, t_min=t_min, t_max=t_max, num_knots=num_knots)
    left_idx, right_idx = get_left_right_indices(knot_pos, num_knots=num_knots)
    csv_path_for_idx = build_indexed_csv_resolver(csv_pattern)

    needed_indices = np.unique(np.concatenate((left_idx, right_idx)))

    csv_cache = {}
    csv_tree = {}
    for idx in needed_indices:
        path = csv_path_for_idx(int(idx))
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing CSV for time-slice index {int(idx)}: {path}")
        pts = load_csv_points(path, has_header=csv_has_header)
        csv_cache[int(idx)] = pts
        csv_tree[int(idx)] = build_xy_tree(pts)

    # Direction 1 (fast): each VTP point checks left/right CSV slices
    pair_keys = np.stack((left_idx, right_idx), axis=1)
    unique_pairs = np.unique(pair_keys, axis=0)
    matched_vtp_mask = np.zeros(vtp_xyz.shape[0], dtype=bool)

    for li, ri in unique_pairs:
        li = int(li)
        ri = int(ri)
        mask = (left_idx == li) & (right_idx == ri)
        src_pts = vtp_xyz[mask]
        if src_pts.size == 0:
            continue

        m_left = batch_match_with_tree(
            src_pts, csv_cache[li], csv_tree[li], xy_threshold, z_threshold
        )
        if ri == li:
            m_any = m_left
        else:
            m_right = batch_match_with_tree(
                src_pts, csv_cache[ri], csv_tree[ri], xy_threshold, z_threshold
            )
            m_any = m_left | m_right
        matched_vtp_mask[mask] = m_any

    matched_vtp = int(np.count_nonzero(matched_vtp_mask))
    total_vtp = int(vtp_xyz.shape[0])
    ratio_vtp = (matched_vtp / total_vtp) if total_vtp else 0.0

    # Direction 2 (fast): each CSV slice checks relevant VTP subset
    matched_csv = 0
    total_csv = 0
    for idx in needed_indices:
        idx = int(idx)
        csv_pts = csv_cache[idx]
        total_csv += int(csv_pts.shape[0])

        candidate_mask = (left_idx == idx) | (right_idx == idx)
        vtp_candidates = vtp_xyz[candidate_mask]
        if vtp_candidates.size == 0 or csv_pts.size == 0:
            continue

        vtp_tree = build_xy_tree(vtp_candidates)
        m = batch_match_with_tree(
            csv_pts, vtp_candidates, vtp_tree, xy_threshold, z_threshold
        )
        matched_csv += int(np.count_nonzero(m))

    ratio_csv = (matched_csv / total_csv) if total_csv else 0.0

    print("VTP result on left/right sliced CSV result:")
    print(f"Matched: {matched_vtp} / {total_vtp}")
    print(f"Ratio:   {ratio_vtp:.6f}")
    print("Left/right sliced CSV result on VTP result:")
    print(f"Matched: {matched_csv} / {total_csv}")
    print(f"Ratio:   {ratio_csv:.6f}")

    return ratio_vtp, ratio_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="4D matching: each point checks left/right CSV slices by t."
    )
    parser.add_argument("vtp_file", type=str, help="Path to VTP file with xyz + t.")
    parser.add_argument(
        "csv_pattern",
        type=str,
        help="CSV path pattern: name_{index}.csv, name_58.csv, or name.csv.",
    )
    parser.add_argument("xy_threshold", type=float, help="XY-plane distance threshold.")
    parser.add_argument("z_threshold", type=float, help="Z distance threshold.")
    parser.add_argument("--header", action="store_true", help="CSV files have header row.")
    parser.add_argument("--t-min", type=float, default=0.0, help="Minimum t value.")
    parser.add_argument("--t-max", type=float, default=89.0, help="Maximum t value.")
    parser.add_argument(
        "--num-knots",
        type=int,
        default=161,
        help="Number of uniform knots in [t-min, t-max].",
    )
    parser.add_argument(
        "--t-array-name",
        type=str,
        default="t",
        help="Point-data array name in VTP that stores t.",
    )
    args = parser.parse_args()

    compute_match_ratio_4d(
        vtp_path=args.vtp_file,
        csv_pattern=args.csv_pattern,
        xy_threshold=args.xy_threshold,
        z_threshold=args.z_threshold,
        csv_has_header=args.header,
        t_min=args.t_min,
        t_max=args.t_max,
        num_knots=args.num_knots,
        t_array_name=args.t_array_name,
    )
import argparse
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def load_csv_points(csv_path, has_header=True):
    """Load (x, y, z) points from a CSV file."""
    skip = 1 if has_header else 0
    pts = np.loadtxt(csv_path, delimiter=",", skiprows=skip, usecols=(0, 1, 2))
    pts = np.atleast_2d(pts)
    return pts.astype(float, copy=False)


def t_to_knot_position(t_values, t_min, t_max, num_knots):
    """Convert t to floating knot position in [0, num_knots - 1]."""
    if num_knots < 2:
        raise ValueError("num_knots must be >= 2.")
    if t_max <= t_min:
        raise ValueError("t_max must be > t_min.")

    scale = (num_knots - 1) / (t_max - t_min)
    pos = (np.asarray(t_values, dtype=float) - t_min) * scale
    return np.clip(pos, 0.0, float(num_knots - 1))


def get_left_right_indices(knot_positions, num_knots):
    """For each knot position, return (left_idx, right_idx)."""
    left = np.floor(knot_positions).astype(int)
    right = np.ceil(knot_positions).astype(int)
    left = np.clip(left, 0, num_knots - 1)
    right = np.clip(right, 0, num_knots - 1)
    return left, right


def load_vtp_points_and_t(vtp_path, t_array_name="t"):
    """Load xyz coordinates from VTP and t values from point-data array."""
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
    """
    Build resolver for indexed CSV path.
    Supported forms:
      1) '/path/name_{index}.csv'
      2) '/path/name_58.csv' -> '/path/name_{idx}.csv'
      3) '/path/name.csv'    -> '/path/name_{idx}.csv'
    """
    if "{index}" in csv_pattern:
        return lambda idx: csv_pattern.format(index=idx)

    stem, ext = os.path.splitext(csv_pattern)
    if ext.lower() != ".csv":
        raise ValueError("csv_pattern must end with .csv or contain '{index}'.")

    m = re.match(r"^(.*)_(\d+)$", stem)
    prefix = m.group(1) if m else stem
    return lambda idx: f"{prefix}_{idx}.csv"


def point_matches_any_slice(point_xyz, candidate_csv_slices, xy_threshold, z_threshold):
    """Check if one point matches any of the provided candidate slice point sets."""
    x, y, z = point_xyz
    xy = np.array([x, y], dtype=float)

    for csv_pts in candidate_csv_slices:
        if csv_pts.size == 0:
            continue
        diff_xy = csv_pts[:, :2] - xy
        dist_xy = np.sqrt(np.sum(diff_xy**2, axis=1))
        mask = dist_xy < xy_threshold
        if np.any(np.abs(csv_pts[mask, 2] - z) < z_threshold):
            return True
    return False


def compute_match_ratio_4d(
    vtp_path,
    csv_pattern,
    xy_threshold,
    z_threshold,
    csv_has_header,
    t_min,
    t_max,
    num_knots,
    t_array_name,
):
    vtp_xyz, vtp_t = load_vtp_points_and_t(vtp_path, t_array_name=t_array_name)
    knot_pos = t_to_knot_position(vtp_t, t_min=t_min, t_max=t_max, num_knots=num_knots)
    left_idx, right_idx = get_left_right_indices(knot_pos, num_knots=num_knots)
    csv_path_for_idx = build_indexed_csv_resolver(csv_pattern)

    needed_indices = np.unique(np.concatenate((left_idx, right_idx)))
    csv_cache = {}
    for idx in needed_indices:
        path = csv_path_for_idx(int(idx))
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing CSV for time-slice index {int(idx)}: {path}")
        csv_cache[int(idx)] = load_csv_points(path, has_header=csv_has_header)

    # Direction 1: each VTP point checks both left/right CSV slices
    matched_vtp = 0
    total_vtp = vtp_xyz.shape[0]
    for i in range(total_vtp):
        li = int(left_idx[i])
        ri = int(right_idx[i])
        candidates = [csv_cache[li]]
        if ri != li:
            candidates.append(csv_cache[ri])
        if point_matches_any_slice(vtp_xyz[i], candidates, xy_threshold, z_threshold):
            matched_vtp += 1
    ratio_vtp = (matched_vtp / total_vtp) if total_vtp else 0.0

    # Direction 2: each CSV point in a used slice checks VTP points from neighboring slices
    idx_to_vtp_mask = {}
    for idx in needed_indices:
        idx_to_vtp_mask[int(idx)] = (left_idx == idx) | (right_idx == idx)

    matched_csv = 0
    total_csv = 0
    for idx in needed_indices:
        csv_pts = csv_cache[int(idx)]
        vtp_candidates = vtp_xyz[idx_to_vtp_mask[int(idx)]]
        total_csv += csv_pts.shape[0]
        for p in csv_pts:
            if point_matches_any_slice(
                p,
                [vtp_candidates],
                xy_threshold=xy_threshold,
                z_threshold=z_threshold,
            ):
                matched_csv += 1
    ratio_csv = (matched_csv / total_csv) if total_csv else 0.0

    print("VTP result on left/right sliced CSV result:")
    print(f"Matched: {matched_vtp} / {total_vtp}")
    print(f"Ratio:   {ratio_vtp:.6f}")
    print("Left/right sliced CSV result on VTP result:")
    print(f"Matched: {matched_csv} / {total_csv}")
    print(f"Ratio:   {ratio_csv:.6f}")

    return ratio_vtp, ratio_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="4D matching: each point checks left/right CSV slices by t."
    )
    parser.add_argument("vtp_file", type=str, help="Path to VTP file with xyz + t.")
    parser.add_argument(
        "csv_pattern",
        type=str,
        help="CSV path pattern: name_{index}.csv, name_58.csv, or name.csv.",
    )
    parser.add_argument("xy_threshold", type=float, help="XY-plane distance threshold.")
    parser.add_argument("z_threshold", type=float, help="Z distance threshold.")
    parser.add_argument("--header", action="store_true", help="CSV files have header row.")
    parser.add_argument("--t-min", type=float, default=0.0, help="Minimum t value.")
    parser.add_argument("--t-max", type=float, default=89.0, help="Maximum t value.")
    parser.add_argument(
        "--num-knots",
        type=int,
        default=161,
        help="Number of uniform knots in [t-min, t-max].",
    )
    parser.add_argument(
        "--t-array-name",
        type=str,
        default="t",
        help="Point-data array name in VTP that stores t.",
    )
    args = parser.parse_args()

    compute_match_ratio_4d(
        vtp_path=args.vtp_file,
        csv_pattern=args.csv_pattern,
        xy_threshold=args.xy_threshold,
        z_threshold=args.z_threshold,
        csv_has_header=args.header,
        t_min=args.t_min,
        t_max=args.t_max,
        num_knots=args.num_knots,
        t_array_name=args.t_array_name,
    )
import argparse
import os
import re

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def load_csv_points(csv_path, has_header=True):
    """Load (x, y, z) points from a CSV file."""
    skip = 1 if has_header else 0
    pts = np.loadtxt(csv_path, delimiter=",", skiprows=skip, usecols=(0, 1, 2))
    pts = np.atleast_2d(pts)
    return pts.astype(float, copy=False)


def map_t_to_index(t_values, t_min, t_max, num_knots):
    """Map t values to uniform knot indices in [0, num_knots-1]."""
    if num_knots < 2:
        raise ValueError("num_knots must be >= 2.")
    denom = (t_max - t_min)
    if denom <= 0:
        raise ValueError("t_max must be > t_min.")

    idx_float = (np.asarray(t_values, dtype=float) - t_min) / denom * (num_knots - 1)
    idx = np.rint(idx_float).astype(int)
    return np.clip(idx, 0, num_knots - 1)


def load_vtp_points_and_t(vtp_path, t_array_name="t"):
    """
    Load xyz coordinates from VTP and t values from point-data array.
    Returns:
      xyz: (N, 3)
      t:   (N,)
    """
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
        n_arrays = point_data.GetNumberOfArrays()
        for i in range(n_arrays):
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


def match_ratio(source_pts, target_pts, xy_threshold, z_threshold):
    """For each source point, check if any target point matches in xy/z thresholds."""
    if source_pts.size == 0:
        return 0, 0, 0.0
    if target_pts.size == 0:
        return 0, source_pts.shape[0], 0.0

    try:
        from scipy.spatial import cKDTree

        tree = cKDTree(target_pts[:, :2])
        neighbors_list = tree.query_ball_point(source_pts[:, :2], r=xy_threshold)
        matched = 0
        for i, neighbors in enumerate(neighbors_list):
            if not neighbors:
                continue
            z_src = source_pts[i, 2]
            if np.any(np.abs(target_pts[neighbors, 2] - z_src) < z_threshold):
                matched += 1
    except ImportError:
        target_xy = target_pts[:, :2]
        target_z = target_pts[:, 2]
        matched = 0
        for x, y, z in source_pts:
            diff_xy = target_xy - np.array([x, y])
            dist_xy = np.sqrt(np.sum(diff_xy**2, axis=1))
            candidates = dist_xy < xy_threshold
            if np.any(np.abs(target_z[candidates] - z) < z_threshold):
                matched += 1

    total = source_pts.shape[0]
    return matched, total, matched / total


def build_indexed_csv_resolver(csv_pattern):
    """
    Build resolver for indexed CSV path.
    Supported forms:
      1) '/path/name_{index}.csv'
      2) '/path/name_58.csv' -> '/path/name_{idx}.csv'
      3) '/path/name.csv'    -> '/path/name_{idx}.csv'
    """
    if "{index}" in csv_pattern:
        return lambda idx: csv_pattern.format(index=idx)

    stem, ext = os.path.splitext(csv_pattern)
    if ext.lower() != ".csv":
        raise ValueError("csv_pattern must end with .csv or contain '{index}'.")

    m = re.match(r"^(.*)_(\d+)$", stem)
    if m:
        prefix = m.group(1)
    else:
        prefix = stem

    return lambda idx: f"{prefix}_{idx}.csv"


def compute_match_ratio_4d(
    vtp_path,
    csv_pattern,
    xy_threshold,
    z_threshold,
    csv_has_header,
    t_min,
    t_max,
    num_knots,
    t_array_name,
):
    vtp_xyz, vtp_t = load_vtp_points_and_t(vtp_path, t_array_name=t_array_name)
    vtp_idx = map_t_to_index(vtp_t, t_min=t_min, t_max=t_max, num_knots=num_knots)
    csv_path_for_idx = build_indexed_csv_resolver(csv_pattern)

    unique_idx = np.unique(vtp_idx)
    csv_cache = {}
    for idx in unique_idx:
        path = csv_path_for_idx(int(idx))
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"Missing CSV for time-slice index {int(idx)}: {path}"
            )
        csv_cache[int(idx)] = load_csv_points(path, has_header=csv_has_header)

    # Direction 1: VTP points -> their corresponding CSV slice
    matched_vtp = 0
    total_vtp = 0
    for idx in unique_idx:
        mask = vtp_idx == idx
        src = vtp_xyz[mask]
        tgt = csv_cache[int(idx)]
        m, n, _ = match_ratio(src, tgt, xy_threshold, z_threshold)
        matched_vtp += m
        total_vtp += n
    ratio_vtp = (matched_vtp / total_vtp) if total_vtp else 0.0

    # Direction 2: CSV slice points -> corresponding VTP points from same time index
    matched_csv = 0
    total_csv = 0
    for idx in unique_idx:
        mask = vtp_idx == idx
        tgt = vtp_xyz[mask]
        src = csv_cache[int(idx)]
        m, n, _ = match_ratio(src, tgt, xy_threshold, z_threshold)
        matched_csv += m
        total_csv += n
    ratio_csv = (matched_csv / total_csv) if total_csv else 0.0

    print("VTP result on sliced CSV result:")
    print(f"Matched: {matched_vtp} / {total_vtp}")
    print(f"Ratio:   {ratio_vtp:.6f}")
    print("Sliced CSV result on VTP result:")
    print(f"Matched: {matched_csv} / {total_csv}")
    print(f"Ratio:   {ratio_csv:.6f}")

    return ratio_vtp, ratio_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="4D matching: use VTP point t to match against indexed CSV slices."
    )
    parser.add_argument("vtp_file", type=str, help="Path to VTP file with xyz + t.")
    parser.add_argument(
        "csv_pattern",
        type=str,
        help="CSV path pattern. Supports name_{index}.csv, name_58.csv, or name.csv.",
    )
    parser.add_argument("xy_threshold", type=float, help="XY-plane distance threshold.")
    parser.add_argument("z_threshold", type=float, help="Z distance threshold.")
    parser.add_argument(
        "--header", action="store_true", help="Indicate that CSV files have a header row."
    )
    parser.add_argument("--t-min", type=float, default=0.0, help="Minimum t value.")
    parser.add_argument("--t-max", type=float, default=89.0, help="Maximum t value.")
    parser.add_argument(
        "--num-knots",
        type=int,
        default=161,
        help="Number of uniform knots in [t-min, t-max].",
    )
    parser.add_argument(
        "--t-array-name",
        type=str,
        default="t",
        help="Point-data array name in VTP that stores t.",
    )

    args = parser.parse_args()

    compute_match_ratio_4d(
        vtp_path=args.vtp_file,
        csv_pattern=args.csv_pattern,
        xy_threshold=args.xy_threshold,
        z_threshold=args.z_threshold,
        csv_has_header=args.header,
        t_min=args.t_min,
        t_max=args.t_max,
        num_knots=args.num_knots,
        t_array_name=args.t_array_name,
    )
import argparse
import numpy as np
import vtk
import os

def load_obj_vertices(obj_path):
    """Load vertex positions (x, y, z) from an .obj file."""
    vertices = []
    with open(obj_path, 'r') as f:
        for line in f:
            if line.startswith('v '):
                parts = line.split()
                if len(parts) >= 4:
                    x, y, z = map(float, parts[1:4])
                    vertices.append((x, y, z))
    return np.array(vertices, dtype=float)


def load_csv_points(csv_path, has_header=True):
    """Load (x, y, z) points from a CSV file."""
    skip = 1 if has_header else 0
    return np.loadtxt(csv_path, delimiter=',', skiprows=skip, usecols=(0, 1, 2))


def load_vtp_points(vtp_path):
    """Load (x, y, z) point coordinates from a .vtp file."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()

    polydata = reader.GetOutput()
    vtk_points = polydata.GetPoints()

    if vtk_points is None:
        raise ValueError(f"No points found in VTP: {vtp_path}")

    n = vtk_points.GetNumberOfPoints()
    pts = np.zeros((n, 3), dtype=float)

    for i in range(n):
        pts[i] = vtk_points.GetPoint(i)

    return pts

def compute_match_ratio(obj_pts, csv_pts, xy_threshold, z_threshold, csv_has_header):
    if obj_pts.size == 0:
        raise ValueError("No vertices found in the OBJ file.")
    if csv_pts.size == 0:
        raise ValueError("No points found in the CSV file.")

    # Try KD-tree acceleration
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(obj_pts[:, :2])
        use_kdtree = True
        print("Using KD-tree for nearest neighbor search.")
    except ImportError:
        use_kdtree = False
        print("scipy not available; using slower method for nearest neighbor search.")

    matched = 0
    n_csv = csv_pts.shape[0]

    if use_kdtree:
        neighbors_list = tree.query_ball_point(csv_pts[:, :2], r=xy_threshold)

        for i, neighbors in enumerate(neighbors_list):
            if not neighbors:
                continue
            z_csv = csv_pts[i, 2]
            if np.any(np.abs(obj_pts[neighbors, 2] - z_csv) < z_threshold):
                matched += 1
    else:
        # Slower fallback
        obj_xy = obj_pts[:, :2]
        obj_z = obj_pts[:, 2]

        for x, y, z in csv_pts:
            diff_xy = obj_xy - np.array([x, y])
            dist_xy = np.sqrt(np.sum(diff_xy**2, axis=1))
            candidates = dist_xy < xy_threshold
            if np.any(np.abs(obj_z[candidates] - z) < z_threshold):
                matched += 1

    ratio = matched / n_csv
    print(f"Matched: {matched} / {n_csv}")
    print(f"Ratio:   {ratio:.6f}")

    return ratio


def compute_match_ratio_different(obj_path, csv_path, xy_threshold, z_threshold, csv_has_header):
    # obj_pts = load_obj_vertices(obj_path)
    obj_pts = load_vtp_points(obj_path)
    
    ext = os.path.splitext(csv_path)[1].lower()
    if ext == ".csv":
        csv_pts= load_csv_points(csv_path, has_header=csv_has_header)
    elif ext == ".vtp":
        csv_pts =load_vtp_points(csv_path)
    # csv_pts = load_csv_points(csv_path, has_header=csv_has_header)
    print("discrete(second) result on our(first) result:")
    compute_match_ratio(obj_pts, csv_pts, xy_threshold, z_threshold, csv_has_header)
    print("our(first) result on discrete(second) result:")
    compute_match_ratio(csv_pts, obj_pts, xy_threshold, z_threshold, csv_has_header)

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare OBJ vertices and CSV points.")
    parser.add_argument("obj_file", type=str, help="Path to OBJ file.")
    parser.add_argument("csv_file", type=str, help="Path to CSV file.")
    parser.add_argument("xy_threshold", type=float, help="XY-plane distance threshold.")
    parser.add_argument("z_threshold", type=float, help="Z distance threshold.")
    parser.add_argument("--header", action="store_true",
                        help="Indicate that CSV has a header row.")

    args = parser.parse_args()

    print("Comparing points between OBJ and CSV files:")

    compute_match_ratio_different(
        obj_path=args.obj_file,
        csv_path=args.csv_file,
        xy_threshold=args.xy_threshold,
        z_threshold=args.z_threshold,
        csv_has_header=args.header,
    )