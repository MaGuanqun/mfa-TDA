#!/usr/bin/env python3
"""
Separate branches in an input VTP, then assign ColorId per branch by matching each
input branch to the closest branch in a reference VTP (reference must already carry
ColorId, e.g. from ``seperate_diffferent_branches.py``).

Pipeline for the input VTP (same as ``seperate_diffferent_branches.py``):
  junction split -> break loops -> enforce monotonic t -> connected components (RegionId)

For each input RegionId, find the closest reference branch (by RegionId, not
ColorId — ColorId can repeat when there are many branches) and copy that branch's
ColorId onto the input branch.

Default matching uses subsampled per-point mean distance with KD-tree candidate
pruning (accurate and fast for 1k+ branches). Use ``--distance-mode centroid`` for
a faster but coarser match.

Usage:
  pvpython src/python/seperate_ttk_trajectory.py \\
    --input-vtp in.vtp \\
    --reference-vtp ref_colorid.vtp \\
    --output-vtp out.vtp
"""

from __future__ import annotations

import argparse
import os
import random
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import vtk
from scipy.spatial import cKDTree
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from seperate_diffferent_branches import (  # noqa: E402
    break_all_loops,
    compute_connectivity_regions,
    enforce_monopoly_in_trajectories,
    split_high_degree_points,
    vtk_read_polydata,
    vtk_write_polydata,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Separate branches in an input VTP and assign ColorId by matching each branch "
            "to the closest branch in a reference VTP (reference ColorId required)."
        )
    )
    p.add_argument("--input-vtp", "-i", required=True, help="Input .vtp to separate and color.")
    p.add_argument(
        "--reference-vtp",
        "-r",
        required=True,
        help="Reference .vtp with ColorId per branch (from seperate_diffferent_branches.py).",
    )
    p.add_argument("--output-vtp", "-o", required=True, help="Output .vtp path.")
    p.add_argument("--data-mode", default="binary", choices=["binary", "ascii"], help="Output VTP data mode.")
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed for branch-splitting randomness (loop breaks; default: 0).",
    )
    p.add_argument(
        "--t-array-name",
        default="t",
        help="Point-data array name for time/scalar (default: t).",
    )
    p.add_argument(
        "--region-array-name",
        default="RegionId",
        help="Connectivity region array name on the input (default: RegionId).",
    )
    p.add_argument(
        "--color-array-name",
        default="ColorId",
        help="Color array name on reference and output (default: ColorId).",
    )
    p.add_argument(
        "--ref-region-array-name",
        default="RegionId",
        help="Reference branch id array (default: RegionId). Branches are matched by this, not ColorId.",
    )
    p.add_argument(
        "--distance-dim",
        choices=["3d", "4d"],
        default="3d",
        help="Match branches in 3D (x,y,z) or 4D (x,y,z,t). Default: 3d.",
    )
    p.add_argument(
        "--distance-mode",
        choices=["mean", "max", "centroid"],
        default="mean",
        help=(
            "Branch distance metric. "
            '"mean" (default) uses subsampled per-point nearest distances; '
            '"max" is directed Hausdorff-like; '
            '"centroid" compares branch centroids only (fastest, least accurate).'
        ),
    )
    p.add_argument(
        "--sample-points",
        type=int,
        default=48,
        help="Max query points sampled per input branch for matching (default: 48).",
    )
    p.add_argument(
        "--candidate-k",
        type=int,
        default=32,
        help=(
            "Nearest reference points queried per sample to build branch candidates "
            "(default: 32). Only used for mean/max modes."
        ),
    )
    p.add_argument(
        "--verbose-matches",
        action="store_true",
        help="Print every input-branch match (default: summary only when >50 branches).",
    )
    return p.parse_args()


def get_points_xyz(poly: vtk.vtkPolyData) -> np.ndarray:
    npts = poly.GetNumberOfPoints()
    if npts == 0:
        return np.zeros((0, 3), dtype=np.float64)

    xyz = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64, copy=False)
    if xyz.shape[1] != 3:
        raise RuntimeError(f"Expected 3D point coordinates, got shape {xyz.shape}")
    return xyz


def get_match_points(
    poly: vtk.vtkPolyData,
    t_array_name: str,
    distance_dim: str,
) -> np.ndarray:
    xyz = get_points_xyz(poly)
    if distance_dim == "3d":
        return xyz

    npts = poly.GetNumberOfPoints()
    t_arr = poly.GetPointData().GetArray(t_array_name)
    if t_arr is None:
        available = [
            poly.GetPointData().GetArrayName(i)
            for i in range(poly.GetPointData().GetNumberOfArrays())
        ]
        raise RuntimeError(
            f'Point-data array "{t_array_name}" not found (required for --distance-dim 4d). '
            f"Available arrays: {available}"
        )

    t_np = vtk_to_numpy(t_arr).astype(np.float64, copy=False)
    if t_np.ndim == 2:
        if t_np.shape[1] != 1:
            raise RuntimeError(
                f'T-array "{t_array_name}" must be scalar per point, got shape {t_np.shape}'
            )
        t_np = t_np[:, 0]
    if t_np.shape[0] != npts:
        raise RuntimeError(
            f'T-array length mismatch: points={npts}, "{t_array_name}"={t_np.shape[0]}'
        )
    return np.column_stack((xyz, t_np))


def group_points_by_scalar_array(
    poly: vtk.vtkPolyData,
    array_name: str,
) -> Dict[int, np.ndarray]:
    arr = poly.GetPointData().GetArray(array_name)
    if arr is None:
        raise RuntimeError(
            f'Point-data array "{array_name}" not found on {poly.GetNumberOfPoints()} points.'
        )

    npts = poly.GetNumberOfPoints()
    if arr.GetNumberOfTuples() != npts:
        raise RuntimeError(
            f'Array "{array_name}" length mismatch: points={npts}, array={arr.GetNumberOfTuples()}'
        )

    groups: Dict[int, List[int]] = defaultdict(list)
    for pid in range(npts):
        groups[int(arr.GetTuple1(pid))].append(pid)
    return {key: np.asarray(pids, dtype=np.int64) for key, pids in groups.items()}


@dataclass
class ReferenceBranchIndex:
    """Precomputed reference branches for fast nearest-branch lookup."""

    region_ids: np.ndarray
    region_to_color: Dict[int, int]
    centroids: np.ndarray
    centroid_tree: cKDTree
    branch_trees: Dict[int, cKDTree]
    global_tree: cKDTree
    global_point_regions: np.ndarray


def _region_to_color_map(
    poly: vtk.vtkPolyData,
    branches: Dict[int, np.ndarray],
    color_array_name: str,
) -> Dict[int, int]:
    color_arr = poly.GetPointData().GetArray(color_array_name)
    if color_arr is None:
        raise RuntimeError(f'Reference missing PointData["{color_array_name}"].')
    mapping: Dict[int, int] = {}
    for region_id, pids in branches.items():
        mapping[int(region_id)] = int(color_arr.GetTuple1(int(pids[0])))
    return mapping


def prepare_reference_branches(
    ref_poly: vtk.vtkPolyData,
    region_array_name: str,
    color_array_name: str,
) -> Tuple[Dict[int, np.ndarray], Dict[int, int], vtk.vtkPolyData]:
    """
    Group reference points by branch id (RegionId).

    ColorId may repeat across branches when there are many components; matching
    must use RegionId and only then copy the matched branch's ColorId.
    """
    region_arr = ref_poly.GetPointData().GetArray(region_array_name)
    if region_arr is None:
        ref_poly = compute_connectivity_regions(ref_poly, region_array_name=region_array_name)

    branches = group_points_by_scalar_array(ref_poly, region_array_name)
    region_to_color = _region_to_color_map(ref_poly, branches, color_array_name)
    return branches, region_to_color, ref_poly


def build_reference_branch_index(
    ref_points: np.ndarray,
    ref_branches: Dict[int, np.ndarray],
) -> ReferenceBranchIndex:
    region_ids = np.asarray(sorted(ref_branches.keys()), dtype=np.int64)
    centroids = np.empty((len(region_ids), ref_points.shape[1]), dtype=np.float64)
    branch_trees: Dict[int, cKDTree] = {}
    global_chunks: List[np.ndarray] = []
    global_region_chunks: List[np.ndarray] = []

    for i, region_id in enumerate(region_ids):
        branch_pids = ref_branches[int(region_id)]
        branch_pts = ref_points[branch_pids]
        centroids[i] = branch_pts.mean(axis=0)
        branch_trees[int(region_id)] = cKDTree(branch_pts)
        global_chunks.append(branch_pts)
        global_region_chunks.append(np.full(branch_pts.shape[0], int(region_id), dtype=np.int64))

    global_points = np.vstack(global_chunks) if global_chunks else np.empty((0, ref_points.shape[1]))
    global_point_regions = (
        np.concatenate(global_region_chunks) if global_region_chunks else np.empty(0, dtype=np.int64)
    )

    return ReferenceBranchIndex(
        region_ids=region_ids,
        region_to_color={},
        centroids=centroids,
        centroid_tree=cKDTree(centroids) if len(region_ids) else cKDTree(np.empty((0, ref_points.shape[1]))),
        branch_trees=branch_trees,
        global_tree=cKDTree(global_points) if global_points.shape[0] else cKDTree(np.empty((0, ref_points.shape[1]))),
        global_point_regions=global_point_regions,
    )


def subsample_branch_points(
    query_point_ids: np.ndarray,
    query_points: np.ndarray,
    max_points: int,
) -> np.ndarray:
    if max_points <= 0 or query_point_ids.size <= max_points:
        return query_points[query_point_ids]

    idx = np.linspace(0, query_point_ids.size - 1, max_points, dtype=np.int64)
    return query_points[query_point_ids[idx]]


def candidate_region_ids(
    query_pts: np.ndarray,
    ref_index: ReferenceBranchIndex,
    candidate_k: int,
) -> np.ndarray:
    if query_pts.shape[0] == 0 or ref_index.global_point_regions.size == 0:
        return ref_index.region_ids

    k = min(max(candidate_k, 1), ref_index.global_point_regions.size)
    try:
        _, nn_idx = ref_index.global_tree.query(query_pts, k=k, workers=-1)
    except TypeError:
        _, nn_idx = ref_index.global_tree.query(query_pts, k=k)

    nn_idx = np.atleast_2d(np.asarray(nn_idx, dtype=np.int64))
    hit_regions = ref_index.global_point_regions[nn_idx.reshape(-1)]
    candidates = np.unique(hit_regions)
    if candidates.size == 0:
        return ref_index.region_ids
    return candidates


def _query_knn_dists(tree: cKDTree, query_pts: np.ndarray) -> np.ndarray:
    if query_pts.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    try:
        dists, _ = tree.query(query_pts, k=1, workers=-1)
    except TypeError:
        dists, _ = tree.query(query_pts, k=1)
    return np.atleast_1d(np.asarray(dists, dtype=np.float64)).reshape(-1)


def score_branch_distance(
    query_pts: np.ndarray,
    ref_tree: cKDTree,
    distance_mode: str,
) -> float:
    dists = _query_knn_dists(ref_tree, query_pts)
    if dists.size == 0:
        return float("inf")
    if distance_mode == "max":
        return float(np.max(dists))
    return float(np.mean(dists))


def find_closest_reference_branch(
    query_point_ids: np.ndarray,
    query_points: np.ndarray,
    ref_index: ReferenceBranchIndex,
    distance_mode: str,
    sample_points: int,
    candidate_k: int,
) -> Tuple[int, int, float]:
    """
    Returns (region_id, color_id, distance).
    """
    if ref_index.region_ids.size == 0:
        raise RuntimeError("No reference branches available for matching.")

    if query_point_ids.size == 0:
        raise RuntimeError("Input branch has no points.")

    query_pts = subsample_branch_points(query_point_ids, query_points, sample_points)

    if distance_mode == "centroid":
        centroid = query_pts.mean(axis=0)
        dist, idx = ref_index.centroid_tree.query(centroid, k=1)
        region_id = int(ref_index.region_ids[int(idx)])
        color_id = int(ref_index.region_to_color[region_id])
        return region_id, color_id, float(dist)

    candidates = candidate_region_ids(query_pts, ref_index, candidate_k)
    centroid = query_pts.mean(axis=0)
    _, centroid_order = ref_index.centroid_tree.query(
        centroid,
        k=min(len(ref_index.region_ids), max(len(candidates), 1)),
    )
    ordered = [int(ref_index.region_ids[int(i)]) for i in np.atleast_1d(centroid_order)]
    seen: set[int] = set()
    search_order: List[int] = []
    for region_id in list(candidates) + ordered:
        rid = int(region_id)
        if rid not in seen and rid in ref_index.branch_trees:
            seen.add(rid)
            search_order.append(rid)

    best_region = -1
    best_dist = float("inf")
    for region_id in search_order:
        dist = score_branch_distance(
            query_pts=query_pts,
            ref_tree=ref_index.branch_trees[region_id],
            distance_mode=distance_mode,
        )
        if dist < best_dist or (
            dist == best_dist and (best_region < 0 or region_id < best_region)
        ):
            best_dist = dist
            best_region = int(region_id)

    if best_region < 0:
        raise RuntimeError("Failed to match input branch to any reference branch.")

    color_id = int(ref_index.region_to_color[best_region])
    return best_region, color_id, best_dist


def separate_branches(
    poly: vtk.vtkPolyData,
    rng: random.Random,
    t_array_name: str,
    region_array_name: str,
) -> vtk.vtkPolyData:
    split = split_high_degree_points(poly, degree_threshold=2)
    no_loops = break_all_loops(split, rng, t_array_name)
    monopoly = enforce_monopoly_in_trajectories(no_loops, t_array_name=t_array_name)
    return compute_connectivity_regions(monopoly, region_array_name=region_array_name)


def assign_color_ids_from_reference(
    poly: vtk.vtkPolyData,
    ref_poly: vtk.vtkPolyData,
    region_array_name: str,
    ref_region_array_name: str,
    color_array_name: str,
    t_array_name: str,
    distance_dim: str,
    distance_mode: str,
    sample_points: int,
    candidate_k: int,
) -> Tuple[Dict[int, Tuple[int, float]], int]:
    ref_branches, region_to_color, ref_poly = prepare_reference_branches(
        ref_poly,
        region_array_name=ref_region_array_name,
        color_array_name=color_array_name,
    )
    input_regions = group_points_by_scalar_array(poly, region_array_name)

    input_points = get_match_points(poly, t_array_name, distance_dim)
    ref_points = get_match_points(ref_poly, t_array_name, distance_dim)
    ref_index = build_reference_branch_index(ref_points, ref_branches)
    ref_index.region_to_color = region_to_color

    region_to_color_out: Dict[int, int] = {}
    match_info: Dict[int, Tuple[int, float]] = {}

    for region_id in sorted(input_regions.keys()):
        query_pids = input_regions[region_id]
        _, color_id, dist = find_closest_reference_branch(
            query_point_ids=query_pids,
            query_points=input_points,
            ref_index=ref_index,
            distance_mode=distance_mode,
            sample_points=sample_points,
            candidate_k=candidate_k,
        )
        region_to_color_out[region_id] = color_id
        match_info[region_id] = (color_id, dist)

    pd = poly.GetPointData()
    cd = poly.GetCellData()
    region_points = pd.GetArray(region_array_name)
    region_cells = cd.GetArray(region_array_name)
    if region_points is None and region_cells is None:
        raise RuntimeError(f"Missing '{region_array_name}' on input polydata.")

    if pd.GetArray(color_array_name) is not None:
        pd.RemoveArray(color_array_name)
    if cd.GetArray(color_array_name) is not None:
        cd.RemoveArray(color_array_name)

    max_region = max(region_to_color_out.keys(), default=0)
    lookup = np.zeros(max_region + 1, dtype=np.int32)
    for region_id, color_id in region_to_color_out.items():
        if region_id >= 0:
            lookup[region_id] = int(color_id)

    if region_points is not None and region_points.GetNumberOfTuples() == poly.GetNumberOfPoints():
        region_np = vtk_to_numpy(region_points).astype(np.int64, copy=False).reshape(-1)
        point_colors = np.zeros(region_np.shape[0], dtype=np.int32)
        valid = (region_np >= 0) & (region_np < lookup.size)
        point_colors[valid] = lookup[region_np[valid]]
    else:
        point_colors = np.zeros(poly.GetNumberOfPoints(), dtype=np.int32)
        poly.BuildLinks()
        cell_ids = vtk.vtkIdList()
        for pid in range(poly.GetNumberOfPoints()):
            poly.GetPointCells(pid, cell_ids)
            if cell_ids.GetNumberOfIds() > 0 and region_cells is not None:
                rid = int(region_cells.GetTuple1(int(cell_ids.GetId(0))))
                if 0 <= rid < lookup.size:
                    point_colors[pid] = lookup[rid]

    if region_cells is not None and region_cells.GetNumberOfTuples() == poly.GetNumberOfCells():
        region_cell_np = vtk_to_numpy(region_cells).astype(np.int64, copy=False).reshape(-1)
        cell_colors = np.zeros(region_cell_np.shape[0], dtype=np.int32)
        valid_cells = (region_cell_np >= 0) & (region_cell_np < lookup.size)
        cell_colors[valid_cells] = lookup[region_cell_np[valid_cells]]
    else:
        cell_colors = np.zeros(poly.GetNumberOfCells(), dtype=np.int32)

    col_p = numpy_to_vtk(point_colors, deep=True)
    col_p.SetName(color_array_name)
    col_p.SetNumberOfComponents(1)

    col_c = numpy_to_vtk(cell_colors, deep=True)
    col_c.SetName(color_array_name)
    col_c.SetNumberOfComponents(1)

    pd.AddArray(col_p)
    cd.AddArray(col_c)
    pd.SetActiveScalars(color_array_name)
    return match_info, len(ref_branches)


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_vtp):
        raise FileNotFoundError(f"Input VTP not found: {args.input_vtp}")
    if not os.path.exists(args.reference_vtp):
        raise FileNotFoundError(f"Reference VTP not found: {args.reference_vtp}")
    if int(args.sample_points) <= 0:
        raise ValueError("--sample-points must be > 0.")
    if int(args.candidate_k) <= 0:
        raise ValueError("--candidate-k must be > 0.")

    rng = random.Random(int(args.seed))
    random.seed(int(args.seed))

    ref = vtk_read_polydata(args.reference_vtp)
    if ref.GetPointData().GetArray(args.color_array_name) is None:
        raise RuntimeError(
            f"Reference VTP missing PointData['{args.color_array_name}']. "
            "Run seperate_diffferent_branches.py on the reference first."
        )

    src = vtk_read_polydata(args.input_vtp)
    labeled = separate_branches(
        src,
        rng=rng,
        t_array_name=args.t_array_name,
        region_array_name=args.region_array_name,
    )

    match_info, ref_branch_count = assign_color_ids_from_reference(
        labeled,
        ref_poly=ref,
        region_array_name=args.region_array_name,
        ref_region_array_name=args.ref_region_array_name,
        color_array_name=args.color_array_name,
        t_array_name=args.t_array_name,
        distance_dim=args.distance_dim,
        distance_mode=args.distance_mode,
        sample_points=int(args.sample_points),
        candidate_k=int(args.candidate_k),
    )

    vtk_write_polydata(labeled, args.output_vtp, args.data_mode)

    dist_label = "distance_3d" if args.distance_dim == "3d" else "distance_4d"
    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {labeled.GetNumberOfPoints()}  Cells: {labeled.GetNumberOfCells()}")
    print(f"Input branches: {len(match_info)}  Reference branches: {ref_branch_count}")
    print(
        f"Distance dim: {args.distance_dim}  Distance mode: {args.distance_mode}  "
        f"sample_points: {args.sample_points}  candidate_k: {args.candidate_k}"
    )
    print_matches = args.verbose_matches or len(match_info) <= 50
    if print_matches:
        print(f"Matches (input RegionId -> reference ColorId, {dist_label}):")
        for region_id in sorted(match_info.keys()):
            color_id, dist = match_info[region_id]
            print(f"  RegionId={region_id} -> ColorId={color_id}, {dist_label}={dist:.8g}")
    else:
        print(f"Matches: {len(match_info)} branches assigned (use --verbose-matches to print all).")


if __name__ == "__main__":
    main()
