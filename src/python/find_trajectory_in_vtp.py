#!/usr/bin/env python3
"""
Find the closest trajectory in B for each trajectory/component in A.

Distance is computed in 4D: (x, y, z, t), where `t` is from point-data array.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read two VTP trajectory files A and B, then for each connected "
            "component in A, find the closest trajectory in B under 4D "
            "distance (x,y,z,t)."
        )
    )
    parser.add_argument("--input-a", required=True, help="Smaller VTP A.")
    parser.add_argument("--input-b", required=True, help="Larger VTP B.")
    parser.add_argument(
        "--t-array-name",
        default="t",
        help='Point-data array name for time dimension in both files (default: "t").',
    )
    parser.add_argument(
        "--region-array-name",
        default="RegionId",
        help='Point-data array name in A for connected component id (default: "RegionId").',
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="XYZ search radius for candidate B points (default: 1e-6).",
    )
    parser.add_argument(
        "--nearest-k",
        type=int,
        default=8,
        help=(
            "For each A point, query this many nearest B points to build "
            "candidate trajectories before exact 4D scoring (default: 8)."
        ),
    )
    parser.add_argument(
        "--distance-mode",
        choices=["mean", "max"],
        default="mean",
        help=(
            'How to score A->B trajectory distance from per-point nearest 4D '
            'distances. "mean" is average; "max" is directed Hausdorff-like.'
        ),
    )
    return parser.parse_args()


def read_polydata(vtp_path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    poly = reader.GetOutput()
    if poly is None:
        raise RuntimeError(f"Failed to read VTP: {vtp_path}")
    return poly


def get_points_4d(poly: vtk.vtkPolyData, t_array_name: str) -> np.ndarray:
    npts = poly.GetNumberOfPoints()
    if npts == 0:
        return np.zeros((0, 4), dtype=np.float64)

    xyz = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64, copy=False)
    if xyz.shape[1] != 3:
        raise RuntimeError(f"Expected 3D point coordinates, got shape {xyz.shape}")

    t_arr = poly.GetPointData().GetArray(t_array_name)
    if t_arr is None:
        available = [
            poly.GetPointData().GetArrayName(i)
            for i in range(poly.GetPointData().GetNumberOfArrays())
        ]
        raise RuntimeError(
            f'Point-data array "{t_array_name}" not found. Available arrays: {available}'
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


def compute_region_ids(poly: vtk.vtkPolyData) -> np.ndarray:
    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputData(poly)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.ColorRegionsOn()
    connectivity.Update()
    out = connectivity.GetOutput()
    region = out.GetPointData().GetArray("RegionId")
    if region is None:
        raise RuntimeError("Connectivity filter did not produce RegionId.")
    region_np = vtk_to_numpy(region)
    if region_np.shape[0] != poly.GetNumberOfPoints():
        raise RuntimeError("RegionId length mismatch after connectivity filter.")
    return region_np.astype(np.int64, copy=False)


def get_region_ids_for_a(poly_a: vtk.vtkPolyData, region_array_name: str) -> np.ndarray:
    arr = poly_a.GetPointData().GetArray(region_array_name)
    if arr is not None:
        region_np = vtk_to_numpy(arr)
        if region_np.ndim == 2 and region_np.shape[1] == 1:
            region_np = region_np[:, 0]
        if region_np.shape[0] != poly_a.GetNumberOfPoints():
            raise RuntimeError(
                f'Region array "{region_array_name}" length mismatch with points in A.'
            )
        return region_np.astype(np.int64, copy=False)
    return compute_region_ids(poly_a)


def build_b_point_to_trajectories(
    poly_b: vtk.vtkPolyData,
) -> tuple[list[set[int]], dict[int, np.ndarray]]:
    npts = poly_b.GetNumberOfPoints()
    point_to_traj = [set() for _ in range(npts)]
    traj_to_points: dict[int, list[int]] = defaultdict(list)
    lines = poly_b.GetLines()

    if lines is None or lines.GetNumberOfCells() == 0:
        return point_to_traj, {}

    id_list = vtk.vtkIdList()
    lines.InitTraversal()
    cell_id = 0
    while lines.GetNextCell(id_list):
        for j in range(id_list.GetNumberOfIds()):
            pid = id_list.GetId(j)
            point_to_traj[pid].add(cell_id)
            traj_to_points[cell_id].append(pid)
        cell_id += 1

    traj_to_points_np = {
        tid: np.asarray(pids, dtype=np.int64) for tid, pids in traj_to_points.items()
    }
    return point_to_traj, traj_to_points_np


def build_candidate_trajectories_for_component(
    comp_point_ids: np.ndarray,
    a_points_4d: np.ndarray,
    b_locator: vtk.vtkStaticPointLocator,
    b_point_to_traj: list[set[int]],
    nearest_k: int,
) -> set[int]:
    candidates: set[int] = set()
    nearest_ids = vtk.vtkIdList()

    for pid in comp_point_ids:
        p4 = a_points_4d[pid]
        b_locator.FindClosestNPoints(nearest_k, p4[:3], nearest_ids)
        for i in range(nearest_ids.GetNumberOfIds()):
            bpid = nearest_ids.GetId(i)
            candidates.update(b_point_to_traj[bpid])
    return candidates


def directed_component_to_trajectory_distance(
    comp_point_ids: np.ndarray,
    a_points_4d: np.ndarray,
    traj_point_ids: np.ndarray,
    b_points_4d: np.ndarray,
    distance_mode: str,
) -> float:
    traj_pts = b_points_4d[traj_point_ids]
    per_point_min_dists = []

    for pid in comp_point_ids:
        p4 = a_points_4d[pid]
        diffs = traj_pts - p4
        dist2 = np.einsum("ij,ij->i", diffs, diffs)
        per_point_min_dists.append(float(np.sqrt(np.min(dist2))))

    if distance_mode == "max":
        return float(np.max(per_point_min_dists))
    return float(np.mean(per_point_min_dists))


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_a):
        raise FileNotFoundError(f"A not found: {args.input_a}")
    if not os.path.exists(args.input_b):
        raise FileNotFoundError(f"B not found: {args.input_b}")
    if args.tolerance < 0:
        raise ValueError("--tolerance must be non-negative.")
    if args.nearest_k <= 0:
        raise ValueError("--nearest-k must be positive.")

    poly_a = read_polydata(args.input_a)
    poly_b = read_polydata(args.input_b)

    a_points_4d = get_points_4d(poly_a, args.t_array_name)
    b_points_4d = get_points_4d(poly_b, args.t_array_name)

    region_ids = get_region_ids_for_a(poly_a, args.region_array_name)
    region_to_points: dict[int, list[int]] = defaultdict(list)
    for point_id, region_id in enumerate(region_ids):
        region_to_points[int(region_id)].append(point_id)

    b_locator = vtk.vtkStaticPointLocator()
    b_locator.SetDataSet(poly_b)
    b_locator.BuildLocator()
    b_point_to_traj, b_traj_to_points = build_b_point_to_trajectories(poly_b)

    print(f"Loaded A: points={poly_a.GetNumberOfPoints()}, regions={len(region_to_points)}")
    print(
        f"Loaded B: points={poly_b.GetNumberOfPoints()}, trajectories(lines)={poly_b.GetLines().GetNumberOfCells() if poly_b.GetLines() else 0}"
    )
    print(f"Candidate search XYZ radius: {args.tolerance}")
    print(f"Nearest-k candidate search: {args.nearest_k}")
    print(f"Distance mode: {args.distance_mode}")
    print("Matches (A_region_id -> closest B_trajectory_id, distance):")

    matched_count = 0
    for region_id in sorted(region_to_points.keys()):
        comp_point_ids = np.asarray(region_to_points[region_id], dtype=np.int64)
        candidate_ids = build_candidate_trajectories_for_component(
            comp_point_ids=comp_point_ids,
            a_points_4d=a_points_4d,
            b_locator=b_locator,
            b_point_to_traj=b_point_to_traj,
            nearest_k=args.nearest_k,
        )
        if not candidate_ids:
            candidate_ids = set(b_traj_to_points.keys())

        best_b_id = None
        best_dist = float("inf")
        for b_tid in candidate_ids:
            traj_point_ids = b_traj_to_points.get(b_tid)
            if traj_point_ids is None or traj_point_ids.size == 0:
                continue
            d = directed_component_to_trajectory_distance(
                comp_point_ids=comp_point_ids,
                a_points_4d=a_points_4d,
                traj_point_ids=traj_point_ids,
                b_points_4d=b_points_4d,
                distance_mode=args.distance_mode,
            )
            if (d < best_dist) or (d == best_dist and (best_b_id is None or b_tid < best_b_id)):
                best_dist = d
                best_b_id = b_tid

        if best_b_id is not None:
            matched_count += 1
            print(
                f"A_region_id={region_id}, A_id={region_id}, B_trajectory_id={best_b_id}, distance_4d={best_dist:.8g}"
            )

    print(f"Matched components in A: {matched_count}/{len(region_to_points)}")


if __name__ == "__main__":
    main()
