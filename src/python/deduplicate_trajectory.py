#!/usr/bin/env pvpython
"""
Remove redundant trajectories from a combined .vtp (VTK PolyData with LINE / POLY_LINE).

1) Optional: merge vertices that lie within a spatial tolerance (same physical node, different
   point ids) and merge duplicate edges.

2) Drop shorter trajectories that lie entirely on a longer polyline (within tolerance, monotone
   along arc length).

3) Trim middle / prefix / suffix overlap: if only part of a shorter curve lies on an already
   kept polyline, remove that overlapping run and keep the remaining one or two polylines.
   (Without this, two curves that share only a middle segment both stay — neither is fully
   subsumed by the other.)

Typical use: after combine_trajectory.py, several traces overlap; keep one copy of the
shared geometry.

After ``seperate_diffferent_branches.py``, each branch has its own ``RegionId`` but
junctions may share the same ``(x,y,z)``. Default vertex merge unions those points in
**xyz only**, re-gluing the LINE graph so **one connected component ≠ one branch**.
Use ``--merge-same-region-only`` (or ``--no-merge-nodes``) if you want one component
per trajectory/branch.

``--recompute-region-id`` rebuilds ``RegionId`` from the current LINE graph (and turns
on region-aware merging when vertex merge is enabled). You get one id per **graph**
component—still only ~20 if the mesh was already merged into ~20 components; run it on
output right after branch splitting, before any step that welds junctions in xyz only.

  pvpython src/python/deduplicate_trajectory.py -i all.vtp -o dedup.vtp --tolerance 1e-4
  pvpython src/python/deduplicate_trajectory.py -i all.vtp -o dedup.vtp --t-array t --tolerance 1e-3

With ``--t-array``, distances default to **4D** Euclidean on ``(x,y,z,time_scale*t)``; vertex merge
stays **3D**. Use ``--spatial-only-distance`` to keep 3D distances while still ordering by ``t``.
  pvpython src/python/deduplicate_trajectory.py -i all.vtp -o dedup.vtp --no-merge-nodes

Pairwise rule (default): drop the shorter trajectory if more than half its arc length lies
within ``--tolerance`` of a longer trajectory (``--overlap-length-frac 0.5``). Set to ``0``
to disable.

Performance: geometry uses NumPy-vectorized segment sweeps (not Python nested loops). Use
``--no-trim-middle`` if you only need whole-trajectory removal (fewer iterations). Vertex
merge cost scales with point count; loosen ``--node-tolerance`` only if needed.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def read_polydata(path: str) -> vtk.vtkPolyData:
    abs_path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(abs_path):
        raise FileNotFoundError(
            f"VTP input does not exist: {abs_path}\n"
            f"  (resolved from {path!r}, cwd={os.getcwd()})"
        )
    if not os.path.isfile(abs_path):
        raise ValueError(f"VTP input is not a regular file: {abs_path}")
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(abs_path)
    reader.Update()
    out = reader.GetOutput()
    if out is None:
        raise RuntimeError(f"VTK reader produced no output for: {abs_path}")
    if out.GetNumberOfPoints() == 0 and out.GetNumberOfCells() == 0:
        raise RuntimeError(
            f"VTK read returned empty PolyData (0 points, 0 cells): {abs_path}. "
            "File may be missing, corrupt, or not a valid .vtp."
        )
    return out


def write_polydata(poly: vtk.vtkPolyData, path: str, data_mode: str) -> None:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(poly)
    if data_mode == "binary":
        writer.SetDataModeToBinary()
    elif data_mode == "ascii":
        writer.SetDataModeToAscii()
    else:
        writer.SetDataModeToAppended()
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write: {path}")


def _is_line_like(ctype: int) -> bool:
    return ctype in (vtk.VTK_LINE, vtk.VTK_POLY_LINE)


def expand_line_segments(poly: vtk.vtkPolyData) -> List[Tuple[int, int, int]]:
    """Each consecutive pair in a LINE / POLY_LINE becomes one undirected segment (cid, p, q)."""
    segs: List[Tuple[int, int, int]] = []
    idlist = vtk.vtkIdList()
    for cid in range(poly.GetNumberOfCells()):
        ctype = poly.GetCellType(cid)
        if not _is_line_like(ctype):
            continue
        poly.GetCellPoints(cid, idlist)
        n = idlist.GetNumberOfIds()
        if n < 2:
            continue
        for k in range(n - 1):
            p = int(idlist.GetId(k))
            q = int(idlist.GetId(k + 1))
            if p != q:
                segs.append((cid, p, q))
    return segs


class UnionFind:
    def __init__(self, n: int) -> None:
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            self.p[ra] = rb
        elif self.r[ra] > self.r[rb]:
            self.p[rb] = ra
        else:
            self.p[rb] = ra
            self.r[ra] += 1


def _all_points_xyz_array(poly: vtk.vtkPolyData) -> np.ndarray:
    pts = poly.GetPoints()
    if pts is None or pts.GetNumberOfPoints() == 0:
        return np.zeros((0, 3), dtype=np.float64)
    return vtk_to_numpy(pts.GetData()).astype(np.float64, copy=False)


def _union_indices_within_tol(
    xyz: np.ndarray,
    tol: float,
    region_id: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    For each point index i, return the smallest index j in its cluster such that all
    pairwise distances in the cluster are connected via edges of length <= tol
    (transitive closure of the "within tol" relation).

    If ``region_id`` is given (length n), only merge pairs with the same region id.
    Use this after ``seperate_diffferent_branches.py``, which duplicates junction points at
    the same (x,y,z); blind spatial merge would glue branches back together.
    """
    n = len(xyz)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    uf = UnionFind(n)
    if tol <= 0:
        return np.arange(n, dtype=np.int64)

    inv = 1.0 / tol
    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    for i in range(n):
        ix = int(np.floor(xyz[i, 0] * inv))
        iy = int(np.floor(xyz[i, 1] * inv))
        iz = int(np.floor(xyz[i, 2] * inv))
        grid[(ix, iy, iz)].append(i)

    tol2 = tol * tol
    for i in range(n):
        ix = int(np.floor(xyz[i, 0] * inv))
        iy = int(np.floor(xyz[i, 1] * inv))
        iz = int(np.floor(xyz[i, 2] * inv))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((ix + dx, iy + dy, iz + dz), []):
                        if j <= i:
                            continue
                        d2 = float(np.sum((xyz[i] - xyz[j]) ** 2))
                        if d2 <= tol2:
                            if region_id is not None and int(region_id[i]) != int(
                                region_id[j]
                            ):
                                continue
                            uf.union(i, j)

    root_min: Dict[int, int] = {}
    for i in range(n):
        r = uf.find(i)
        prev = root_min.get(r, n)
        if i < prev:
            root_min[r] = i
    canon = np.empty(n, dtype=np.int64)
    for i in range(n):
        r = uf.find(i)
        canon[i] = root_min[r]
    return canon


def merge_graph_by_proximity(
    poly: vtk.vtkPolyData,
    node_tol: float,
    region_ids: Optional[np.ndarray] = None,
) -> vtk.vtkPolyData:
    """
    Merge vertices closer than `node_tol` (same cluster gets one point, position averaged).
    Collapse duplicate edges (same undirected merged endpoints). Keeps first cell's cell data.

    ``region_ids`` (optional): only merge point pairs that share the same id (e.g. RegionId).
    """
    segs = expand_line_segments(poly)
    if not segs:
        return poly

    npts = poly.GetNumberOfPoints()
    xyz = _all_points_xyz_array(poly)
    if len(xyz) != npts:
        raise RuntimeError("Point count mismatch.")

    canon = _union_indices_within_tol(xyz, node_tol, region_ids)
    unique_old = sorted(set(int(canon[i]) for i in range(npts)))
    old_rep_to_new = {old: j for j, old in enumerate(unique_old)}
    n_new = len(unique_old)

    sums = np.zeros((n_new, 3), dtype=np.float64)
    counts = np.zeros(n_new, dtype=np.int64)
    for i in range(npts):
        j = old_rep_to_new[int(canon[i])]
        sums[j] += xyz[i]
        counts[j] += 1
    means = sums / np.maximum(counts, 1).reshape(-1, 1)

    out_pts = vtk.vtkPoints()
    out_pts.SetData(numpy_to_vtk(means.astype(np.float64)))

    old_pd = poly.GetPointData()
    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_pd = out_poly.GetPointData()

    for ai in range(old_pd.GetNumberOfArrays()):
        src = old_pd.GetArray(ai)
        if src is None:
            continue
        dst = src.NewInstance()
        dst.DeepCopy(src)
        dst.SetNumberOfTuples(n_new)
        for new_j, old_pid in enumerate(unique_old):
            dst.SetTuple(new_j, old_pid, src)
        out_pd.AddArray(dst)
    if old_pd.GetScalars() is not None and old_pd.GetScalars().GetName():
        out_pd.SetActiveScalars(old_pd.GetScalars().GetName())

    old_cd = poly.GetCellData()
    edge_seen: Set[Tuple[int, int]] = set()
    out_lines = vtk.vtkCellArray()
    line = vtk.vtkLine()
    kept_cell_ids: List[int] = []

    for cid, p, q in segs:
        rp = old_rep_to_new[int(canon[p])]
        rq = old_rep_to_new[int(canon[q])]
        if rp == rq:
            continue
        a, b = (rp, rq) if rp < rq else (rq, rp)
        if (a, b) in edge_seen:
            continue
        edge_seen.add((a, b))
        line.GetPointIds().SetId(0, rp)
        line.GetPointIds().SetId(1, rq)
        out_lines.InsertNextCell(line)
        kept_cell_ids.append(cid)

    out_poly.SetLines(out_lines)
    out_cd = out_poly.GetCellData()
    n_cells = len(kept_cell_ids)
    for ai in range(old_cd.GetNumberOfArrays()):
        src = old_cd.GetArray(ai)
        if src is None:
            continue
        dst = src.NewInstance()
        dst.DeepCopy(src)
        dst.SetNumberOfTuples(n_cells)
        for new_cid, old_cid in enumerate(kept_cell_ids):
            dst.SetTuple(new_cid, old_cid, src)
        out_cd.AddArray(dst)
    if old_cd.GetScalars() is not None and old_cd.GetScalars().GetName():
        out_cd.SetActiveScalars(old_cd.GetScalars().GetName())

    out_poly.BuildCells()
    out_poly.BuildLinks()
    return out_poly


def build_adjacency(npts: int, segs: Sequence[Tuple[int, int, int]]) -> List[Set[int]]:
    adj: List[Set[int]] = [set() for _ in range(npts)]
    for _cid, p, q in segs:
        adj[p].add(q)
        adj[q].add(p)
    return adj


def component_segments(
    segs: Sequence[Tuple[int, int, int]], npts: int
) -> Dict[int, List[Tuple[int, int, int]]]:
    uf = UnionFind(npts)
    for _cid, p, q in segs:
        uf.union(p, q)
    by_root: Dict[int, List[Tuple[int, int, int]]] = defaultdict(list)
    for item in segs:
        cid, p, q = item
        by_root[uf.find(p)].append(item)
    return dict(by_root)


def assign_region_ids_from_line_graph(
    poly: vtk.vtkPolyData,
    array_name: str,
    propagate_to_line_cells: bool = True,
) -> Tuple[int, int]:
    """
    Replace point-data ``array_name`` with dense labels 0 .. K-1 from VTK_LINE /
    VTK_POLY_LINE connectivity (one id per connected component). Points not used by any
    line cell each get their own id.

    Optionally writes the same ``array_name`` on **cell data** for line-like cells
    (value = RegionId of the cell's first point).

    Returns ``(n_line_graph_components, n_distinct_region_ids)``.
    """
    npts = poly.GetNumberOfPoints()
    segs = expand_line_segments(poly)
    by_comp = component_segments(segs, npts) if segs else {}

    region = np.full(npts, -1, dtype=np.int32)
    next_rid = 0
    for _root in sorted(by_comp.keys()):
        comp_segs = by_comp[_root]
        for _c, p, q in comp_segs:
            region[p] = next_rid
            region[q] = next_rid
        next_rid += 1

    n_line = len(by_comp)
    for i in range(npts):
        if int(region[i]) < 0:
            region[i] = next_rid
            next_rid += 1

    n_distinct = int(next_rid)

    pd = poly.GetPointData()
    if pd.GetArray(array_name) is not None:
        pd.RemoveArray(array_name)
    pa = numpy_to_vtk(region.ravel(), deep=1, array_type=vtk.VTK_INT)
    pa.SetName(array_name)
    pd.AddArray(pa)

    if propagate_to_line_cells:
        idlist = vtk.vtkIdList()
        cd = poly.GetCellData()
        if cd.GetArray(array_name) is not None:
            cd.RemoveArray(array_name)
        ca = vtk.vtkIntArray()
        ca.SetName(array_name)
        n_cells = poly.GetNumberOfCells()
        ca.SetNumberOfTuples(n_cells)
        for cid in range(n_cells):
            ctype = poly.GetCellType(cid)
            if not _is_line_like(ctype):
                ca.SetTuple1(cid, -1)
                continue
            poly.GetCellPoints(cid, idlist)
            if idlist.GetNumberOfIds() < 1:
                ca.SetTuple1(cid, -1)
                continue
            p0 = int(idlist.GetId(0))
            ca.SetTuple1(cid, int(region[p0]))
        cd.AddArray(ca)

    return n_line, n_distinct


def _edge_set_from_segs(comp_segs: Sequence[Tuple[int, int, int]]) -> Set[Tuple[int, int]]:
    es: Set[Tuple[int, int]] = set()
    for _c, p, q in comp_segs:
        es.add((p, q))
        es.add((q, p))
    return es


def _try_order_by_increasing_t(
    comp_nodes: Set[int], t_vals: np.ndarray, edge_set: Set[Tuple[int, int]]
) -> Optional[List[int]]:
    """If sorting by t yields consecutive graph edges, return that order."""
    nodes = sorted(comp_nodes, key=lambda pid: float(t_vals[pid]))
    for i in range(len(nodes) - 1):
        if (nodes[i], nodes[i + 1]) not in edge_set:
            return None
    return nodes


def _walk_simple_path(
    start: int,
    comp_nodes: Set[int],
    adj: List[Set[int]],
    max_steps: int,
) -> List[int]:
    """Walk degree-1 chain from `start` until branch or end."""
    order: List[int] = []
    prev = -1
    curr = start
    visited: Set[int] = set()
    for _ in range(max_steps):
        if curr in visited:
            break
        order.append(curr)
        visited.add(curr)
        nxts = [v for v in adj[curr] if v != prev and v in comp_nodes]
        if not nxts:
            break
        if len(nxts) > 1:
            break
        prev, curr = curr, nxts[0]
    return order


def order_path_points(
    comp_segs: Sequence[Tuple[int, int, int]],
    adj: List[Set[int]],
    comp_nodes: Set[int],
    t_vals: Optional[np.ndarray],
) -> Optional[List[int]]:
    """Return an ordered list of point ids along the polyline (only consecutive pairs are edges)."""
    edge_set = _edge_set_from_segs(comp_segs)
    deg = {pid: len(adj[pid] & comp_nodes) for pid in comp_nodes}
    endpoints = [pid for pid in comp_nodes if deg.get(pid, 0) == 1]

    if t_vals is not None:
        tord = _try_order_by_increasing_t(comp_nodes, t_vals, edge_set)
        if tord is not None and len(tord) == len(comp_nodes):
            return tord

    # Path with two endpoints: walk from one end; reverse if t decreases along the walk.
    if len(endpoints) == 2:
        w = _walk_simple_path(endpoints[0], comp_nodes, adj, len(comp_nodes) + 2)
        if len(w) == len(comp_nodes):
            if t_vals is not None and len(w) >= 2:
                ts = t_vals[np.array(w, dtype=np.int64)]
                if float(ts[-1]) < float(ts[0]):
                    w = list(reversed(w))
            return w

    # One endpoint (dead end + branch) or tree-like: try walking from each degree-1 node.
    if len(endpoints) >= 1:
        best: Optional[List[int]] = None
        for ep in endpoints:
            w = _walk_simple_path(ep, comp_nodes, adj, len(comp_nodes) + 2)
            if len(w) == len(comp_nodes):
                best = w
                break
        if best is not None:
            if t_vals is not None and len(best) >= 2:
                ts = t_vals[np.array(best, dtype=np.int64)]
                if float(ts[-1]) < float(ts[0]):
                    best = list(reversed(best))
            return best

    # Cycle: every vertex has degree 2 in the component.
    if not endpoints and comp_nodes:
        n = len(comp_nodes)
        if n >= 3 and all(deg.get(pid, 0) == 2 for pid in comp_nodes):
            start = min(comp_nodes)
            neighs = sorted([v for v in adj[start] if v in comp_nodes])
            if len(neighs) == 2:
                for first in neighs:
                    order_c: List[int] = [start]
                    prev = start
                    curr = first
                    ok = True
                    for _ in range(n - 1):
                        order_c.append(curr)
                        nxts = [v for v in adj[curr] if v != prev and v in comp_nodes]
                        if len(nxts) != 1:
                            ok = False
                            break
                        prev, curr = curr, nxts[0]
                    if ok and len(order_c) == n and start in adj[order_c[-1]]:
                        w = order_c
                        if t_vals is not None and len(w) >= 2:
                            ts = t_vals[np.array(w, dtype=np.int64)]
                            if float(ts[-1]) < float(ts[0]):
                                w = list(reversed(w))
                        return w

    # Branching graph: cannot build one edge-consistent polyline without extra rules.
    return None


def points_xyz(poly: vtk.vtkPolyData, pids: Sequence[int]) -> np.ndarray:
    pts = poly.GetPoints()
    out = np.zeros((len(pids), 3), dtype=np.float64)
    for i, pid in enumerate(pids):
        out[i] = pts.GetPoint(pid)
    return out


def points_coords(
    poly: vtk.vtkPolyData,
    pids: Sequence[int],
    t_vals: np.ndarray,
    time_scale: float = 1.0,
) -> np.ndarray:
    """
    Build (N, 4) coordinates [x, y, z, time_scale * t] for distance in 4D Euclidean space.
    """
    pts = poly.GetPoints()
    n = len(pids)
    out = np.zeros((n, 4), dtype=np.float64)
    for i, pid in enumerate(pids):
        out[i, 0], out[i, 1], out[i, 2] = pts.GetPoint(pid)
        out[i, 3] = time_scale * float(t_vals[pid])
    return out


def trajectory_coords(
    poly: vtk.vtkPolyData,
    order: Sequence[int],
    t_vals: Optional[np.ndarray],
    time_scale: float,
    use_4d_distance: bool,
) -> np.ndarray:
    """(N,3) spatial-only or (N,4) with t as 4th coordinate for distances."""
    if use_4d_distance and t_vals is not None:
        return points_coords(poly, order, t_vals, time_scale)
    return points_xyz(poly, order)


def polyline_arc_length(poly_xyz: np.ndarray) -> np.ndarray:
    """Cumulative arc length at each vertex; length len(poly_xyz)."""
    if len(poly_xyz) < 2:
        return np.zeros(len(poly_xyz), dtype=np.float64)
    seg = np.linalg.norm(np.diff(poly_xyz, axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


def _point_at_arc_length(poly_xyz: np.ndarray, cum: np.ndarray, s: float) -> np.ndarray:
    """Point at arc length s along the polyline (piecewise linear)."""
    if len(poly_xyz) < 2:
        return poly_xyz[0].copy()
    L = float(cum[-1])
    if s <= 0:
        return poly_xyz[0].copy()
    if s >= L:
        return poly_xyz[-1].copy()
    idx = int(np.searchsorted(cum, s, side="right") - 1)
    idx = max(0, min(idx, len(poly_xyz) - 2))
    seg = float(cum[idx + 1] - cum[idx])
    if seg <= 1e-15:
        return poly_xyz[idx].copy()
    t = (s - float(cum[idx])) / seg
    return poly_xyz[idx] + t * (poly_xyz[idx + 1] - poly_xyz[idx])


def overlap_arc_length_fraction(
    shorter_xyz: np.ndarray, longer_xyz: np.ndarray, tol: float
) -> float:
    """
    Fraction of the shorter polyline's total arc length that lies within `tol` of the
    longer polyline (uniform sampling along arc length of the shorter curve).
    """
    if len(shorter_xyz) < 2 or len(longer_xyz) < 2:
        return 0.0
    cum = polyline_arc_length(shorter_xyz)
    L = float(cum[-1])
    if L <= 1e-15:
        return 0.0
    if _aabbs_min_separation_squared(shorter_xyz, longer_xyz) > tol * tol:
        return 0.0
    n_seg = max(24, min(500, int(np.ceil(L / max(tol * 0.25, 1e-12)))))
    mids = (np.arange(n_seg, dtype=np.float64) + 0.5) * (L / n_seg)
    ddim = int(shorter_xyz.shape[1])
    pts = np.empty((n_seg, ddim), dtype=np.float64)
    for k in range(n_seg):
        pts[k] = _point_at_arc_length(shorter_xyz, cum, mids[k])
    d = distances_points_to_polyline(pts, longer_xyz)
    seg_len = L / n_seg
    overlap_len = float(np.sum(seg_len * (d <= tol)))
    return overlap_len / L


def _sep_axis_aligned_boxes(
    q_lo: np.ndarray,
    q_hi: np.ndarray,
    p_lo: np.ndarray,
    p_hi: np.ndarray,
) -> float:
    """Minimum Euclidean distance between two AABBs (squared). Dimension-agnostic (3D or 4D)."""
    d2 = 0.0
    for ax in range(int(q_lo.shape[0])):
        if q_hi[ax] < p_lo[ax]:
            t = float(p_lo[ax] - q_hi[ax])
            d2 += t * t
        elif p_hi[ax] < q_lo[ax]:
            t = float(q_lo[ax] - p_hi[ax])
            d2 += t * t
    return d2


def _aabbs_min_separation_squared(q_xyz: np.ndarray, p_xyz: np.ndarray) -> float:
    """Minimum Euclidean distance between two point clouds' axis-aligned bounding boxes (squared)."""
    if len(q_xyz) == 0 or len(p_xyz) == 0:
        return np.inf
    q_lo = np.min(q_xyz, axis=0)
    q_hi = np.max(q_xyz, axis=0)
    p_lo = np.min(p_xyz, axis=0)
    p_hi = np.max(p_xyz, axis=0)
    return _sep_axis_aligned_boxes(q_lo, q_hi, p_lo, p_hi)


def distances_and_projection_to_polyline(
    query_xyz: np.ndarray, poly_xyz: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Per query point: minimum distance to the polyline and arc-length parameter of the
    closest point on that polyline. Works in R^D for D=3 (xyz) or D=4 (xyz + t).
    Fully vectorized over queries (one segment loop).
    """
    nq = len(query_xyz)
    if len(poly_xyz) == 0:
        return np.full(nq, np.inf, dtype=np.float64), np.zeros(nq, dtype=np.float64)
    if len(poly_xyz) == 1:
        d = np.linalg.norm(query_xyz - poly_xyz[0], axis=1)
        return d, np.zeros(nq, dtype=np.float64)

    best_d = np.full(nq, np.inf, dtype=np.float64)
    best_s = np.zeros(nq, dtype=np.float64)
    cum = polyline_arc_length(poly_xyz)
    Q = query_xyz
    for i in range(len(poly_xyz) - 1):
        a = poly_xyz[i]
        b = poly_xyz[i + 1]
        ab = b - a
        denom = float(np.dot(ab, ab))
        if denom <= 0:
            d = np.linalg.norm(Q - a, axis=1)
            s_seg = np.full(nq, cum[i], dtype=np.float64)
        else:
            t = np.dot(Q - a, ab) / denom
            t = np.clip(t, 0.0, 1.0)
            proj = a + t[:, None] * ab
            d = np.linalg.norm(Q - proj, axis=1)
            s_seg = cum[i] + t * (cum[i + 1] - cum[i])
        mask = d < best_d
        best_d[mask] = d[mask]
        best_s[mask] = s_seg[mask]
    return best_d, best_s


def distances_points_to_polyline(query_xyz: np.ndarray, poly_xyz: np.ndarray) -> np.ndarray:
    """Minimum distance from each query point to the polyline (union of segments)."""
    d, _ = distances_and_projection_to_polyline(query_xyz, poly_xyz)
    return d


def projection_arc_lengths(query_xyz: np.ndarray, poly_xyz: np.ndarray) -> np.ndarray:
    """Arc length along `poly_xyz` of the closest point to each query (for monotonicity)."""
    _, s = distances_and_projection_to_polyline(query_xyz, poly_xyz)
    return s


def is_subsumed_by(
    short_xyz: np.ndarray,
    long_xyz: np.ndarray,
    tol: float,
    mono_rel: float,
    overlap_fraction: float = 1.0,
    long_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> bool:
    """
    True if enough of `short_xyz` lies on the polyline `long_xyz` (within tol) and
    projections along arc length are monotone on the close points.
    overlap_fraction=1.0 requires every sample within tol (strict).
    """
    if len(short_xyz) < 2 or len(long_xyz) < 2:
        return False
    q_lo = np.min(short_xyz, axis=0)
    q_hi = np.max(short_xyz, axis=0)
    if long_bounds is not None:
        p_lo, p_hi = long_bounds
        if _sep_axis_aligned_boxes(q_lo, q_hi, p_lo, p_hi) > tol * tol:
            return False
    elif _aabbs_min_separation_squared(short_xyz, long_xyz) > tol * tol:
        return False
    d, s = distances_and_projection_to_polyline(short_xyz, long_xyz)
    close = d <= tol
    if overlap_fraction >= 1.0 - 1e-12:
        if float(np.max(d)) > tol:
            return False
        mask = np.ones(len(short_xyz), dtype=bool)
    else:
        if float(np.mean(close)) < overlap_fraction:
            return False
        mask = close
        if int(np.sum(mask)) < 2:
            return False
    s = s[mask]
    if len(s) < 2:
        return False
    ds = np.diff(s)
    inc = bool(np.all(ds >= -mono_rel))
    dec = bool(np.all(ds <= mono_rel))
    return inc or dec


def _monotonic_1d(s: np.ndarray, mono_rel: float) -> bool:
    if len(s) < 2:
        return True
    ds = np.diff(s)
    return bool(np.all(ds >= -mono_rel) or np.all(ds <= mono_rel))


def split_once_against_reference(
    order: List[int],
    xyz: np.ndarray,
    ref_xyz: np.ndarray,
    tol: float,
    mono_rel: float,
    ref_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> Optional[List[Tuple[List[int], np.ndarray]]]:
    """
    If a contiguous run of samples lies on `ref_xyz` (close + monotonic along ref), remove
    that run and return one or two polylines. Returns None if there is nothing to remove.
    """
    n = len(xyz)
    if n < 2 or len(ref_xyz) < 2:
        return None
    q_lo = np.min(xyz, axis=0)
    q_hi = np.max(xyz, axis=0)
    if ref_bounds is not None:
        p_lo, p_hi = ref_bounds
        if _sep_axis_aligned_boxes(q_lo, q_hi, p_lo, p_hi) > tol * tol:
            return None
    elif _aabbs_min_separation_squared(xyz, ref_xyz) > tol * tol:
        return None

    d, s = distances_and_projection_to_polyline(xyz, ref_xyz)

    runs: List[Tuple[int, int]] = []
    k = 0
    while k < n:
        if d[k] > tol:
            k += 1
            continue
        k0 = k
        while k + 1 < n and d[k + 1] <= tol:
            k += 1
        ss = s[k0 : k + 1]
        if not _monotonic_1d(ss, mono_rel):
            k += 1
            continue
        runs.append((k0, k))
        k += 1

    if not runs:
        return None

    i, j = max(runs, key=lambda t: (t[1] - t[0], t[1]))

    if i == 0 and j == n - 1:
        return []

    out: List[Tuple[List[int], np.ndarray]] = []

    if i > 0 and j < n - 1:
        if i >= 2:
            out.append((order[:i], xyz[:i].copy()))
        if n - j - 1 >= 2:
            out.append((order[j + 1 :], xyz[j + 1 :].copy()))
        return out if out else None

    if i == 0 and j < n - 1:
        if n - j - 1 >= 2:
            return [(order[j + 1 :], xyz[j + 1 :].copy())]
        return None

    if i > 0 and j == n - 1:
        if i >= 2:
            return [(order[:i], xyz[:i].copy())]
        return None

    return None


def trim_polyline_against_reference(
    order: List[int],
    xyz: np.ndarray,
    ref_xyz: np.ndarray,
    tol: float,
    mono_rel: float,
    overlap_fraction: float,
    ref_bounds: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> List[Tuple[List[int], np.ndarray]]:
    """
    Remove all geometry of (order, xyz) that duplicates `ref_xyz` (fully subsumed, or
    iterative prefix / interior / suffix overlap removal).
    """
    if len(xyz) < 2:
        return []

    if is_subsumed_by(
        xyz,
        ref_xyz,
        tol,
        mono_rel,
        overlap_fraction=overlap_fraction,
        long_bounds=ref_bounds,
    ):
        return []

    work: List[Tuple[List[int], np.ndarray]] = [(list(order), np.array(xyz, copy=True))]
    # Each successful split shortens the polyline or shortens overlap; cap passes tightly.
    max_passes = min(128, max(24, len(xyz) + 16))
    for _ in range(max_passes):
        new_work: List[Tuple[List[int], np.ndarray]] = []
        changed = False
        for o, x in work:
            if len(x) < 2:
                continue
            if is_subsumed_by(
                x,
                ref_xyz,
                tol,
                mono_rel,
                overlap_fraction=overlap_fraction,
                long_bounds=ref_bounds,
            ):
                changed = True
                continue
            nxt = split_once_against_reference(
                o, x, ref_xyz, tol, mono_rel, ref_bounds=ref_bounds
            )
            if nxt is None:
                new_work.append((o, x))
                continue
            changed = True
            for piece in nxt:
                if len(piece[1]) >= 2:
                    new_work.append(piece)
        if not changed:
            break
        work = new_work
        if not work:
            return []

    return [p for p in work if len(p[1]) >= 2]


def collect_edges_along_order(
    order: Sequence[int], edge_map: Dict[Tuple[int, int], List[int]]
) -> List[int]:
    """Map consecutive point pairs to original cell ids (first match)."""
    cids: List[int] = []
    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        key = (a, b) if a < b else (b, a)
        lst = edge_map.get(key)
        if not lst:
            raise RuntimeError(f"Missing edge ({a},{b}) in segment map.")
        cids.append(lst[0])
    return cids


def deduplicate_trajectories(
    poly: vtk.vtkPolyData,
    tolerance: float,
    t_array_name: Optional[str],
    mono_rel: float,
    merge_close_nodes: bool = True,
    node_tolerance: Optional[float] = None,
    overlap_fraction: float = 1.0,
    trim_middle_overlap: bool = True,
    overlap_length_frac_min: float = 0.5,
    time_scale: float = 1.0,
    spatial_only_distance: bool = False,
    merge_same_region_only: bool = False,
    region_array_name: str = "RegionId",
    recompute_region_id: bool = False,
    propagate_cell_region_id: bool = True,
) -> Tuple[vtk.vtkPolyData, Dict[str, int]]:
    nt = node_tolerance if node_tolerance is not None else tolerance

    region_id_recomputed = 0
    line_graph_components_after_recompute = 0
    distinct_region_ids_after_recompute = 0
    if recompute_region_id:
        line_graph_components_after_recompute, distinct_region_ids_after_recompute = (
            assign_region_ids_from_line_graph(
                poly,
                region_array_name,
                propagate_to_line_cells=propagate_cell_region_id,
            )
        )
        region_id_recomputed = 1

    point_data_has_region_id = (
        poly.GetPointData().GetArray(region_array_name) is not None
    )

    n_pts_before = poly.GetNumberOfPoints()
    segs0 = expand_line_segments(poly)
    n_seg_before = len(segs0)
    graph_components_before_merge = (
        len(component_segments(segs0, n_pts_before)) if segs0 and n_pts_before > 0 else 0
    )

    region_for_merge: Optional[np.ndarray] = None
    if merge_close_nodes and nt > 0:
        if merge_same_region_only:
            rarr = poly.GetPointData().GetArray(region_array_name)
            if rarr is None:
                raise ValueError(
                    f"--merge-same-region-only requires point-data array {region_array_name!r}."
                )
            region_for_merge = vtk_to_numpy(rarr).astype(np.int64).reshape(-1)
            if len(region_for_merge) != poly.GetNumberOfPoints():
                raise ValueError(
                    f"{region_array_name!r} length {len(region_for_merge)} != "
                    f"number of points {poly.GetNumberOfPoints()}."
                )
        poly = merge_graph_by_proximity(poly, nt, region_for_merge)

    n_pts_after = poly.GetNumberOfPoints()
    segs = expand_line_segments(poly)
    n_seg_after = len(segs)
    if not segs:
        n_pts_now = poly.GetNumberOfPoints()
        n_cells_now = poly.GetNumberOfCells()
        raise ValueError(
            "No VTK_LINE / VTK_POLY_LINE segments in input. "
            f"The mesh has {n_pts_now} points and {n_cells_now} cells but no line/polylines "
            "this script can use. Check that the .vtp loaded correctly and contains trajectory "
            "edges (not only vertices or triangles)."
        )

    npts = poly.GetNumberOfPoints()
    adj = build_adjacency(npts, segs)
    by_comp = component_segments(segs, npts)
    graph_components_after_merge = len(by_comp)

    edge_map: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for cid, p, q in segs:
        key = (p, q) if p < q else (q, p)
        edge_map[key].append(cid)

    t_vals: Optional[np.ndarray] = None
    if t_array_name:
        arr = poly.GetPointData().GetArray(t_array_name)
        if arr is None:
            raise ValueError(f"Point-data array {t_array_name!r} not found.")
        t_vals = vtk_to_numpy(arr).astype(np.float64, copy=False).reshape(-1)

    use_4d_distance = t_vals is not None and not spatial_only_distance

    # Per-component ordered points: (N,3) or (N,4) with t as 4th dim for distances.
    traj: List[Tuple[float, List[int], np.ndarray]] = []
    for _root, comp_segs in by_comp.items():
        nodes = set()
        for _c, p, q in comp_segs:
            nodes.add(p)
            nodes.add(q)
        if len(nodes) < 2:
            continue
        order = order_path_points(comp_segs, adj, nodes, t_vals)
        if order is None or len(order) < 2:
            continue
        coords = trajectory_coords(poly, order, t_vals, time_scale, use_4d_distance)
        arc = float(polyline_arc_length(coords)[-1])
        traj.append((arc, order, coords))

    if not traj:
        raise ValueError("Could not extract any trajectory with >= 2 points.")

    traj.sort(key=lambda x: -x[0])
    trajectories_extracted = len(traj)
    removed_pairwise_overlap = 0
    if overlap_length_frac_min > 0.0:
        n = len(traj)
        mark = [False] * n
        for i in range(n):
            _ai, _oi, xyz_i = traj[i]
            for j in range(i + 1, n):
                if mark[j]:
                    continue
                _aj, _oj, xyz_j = traj[j]
                ov = overlap_arc_length_fraction(xyz_j, xyz_i, tolerance)
                if ov >= overlap_length_frac_min:
                    mark[j] = True
        removed_pairwise_overlap = int(sum(mark))
        traj = [traj[k] for k in range(n) if not mark[k]]

    if not traj:
        raise ValueError(
            "All trajectories were removed by pairwise overlap-length rule; "
            "lower --overlap-length-frac or increase --tolerance."
        )

    # Longest arc first. For each trajectory, either drop it (fully on a kept polyline),
    # trim overlapping runs against kept polylines, or keep it as new.
    kept: List[Tuple[List[int], np.ndarray]] = []
    kept_bounds: List[Tuple[np.ndarray, np.ndarray]] = []
    removed_entirely = 0
    split_into_multiple = 0
    for _arc, order, xyz in traj:
        parts: List[Tuple[List[int], np.ndarray]] = [(order, xyz)]
        if trim_middle_overlap:
            for (_ref_o, ref_x), rb in zip(kept, kept_bounds):
                new_parts: List[Tuple[List[int], np.ndarray]] = []
                for o, x in parts:
                    new_parts.extend(
                        trim_polyline_against_reference(
                            o,
                            x,
                            ref_x,
                            tolerance,
                            mono_rel,
                            overlap_fraction,
                            ref_bounds=rb,
                        )
                    )
                parts = new_parts
        else:
            redundant = False
            for (_ref_o, ref_x), rb in zip(kept, kept_bounds):
                if is_subsumed_by(
                    xyz,
                    ref_x,
                    tolerance,
                    mono_rel,
                    overlap_fraction=overlap_fraction,
                    long_bounds=rb,
                ):
                    redundant = True
                    break
            if redundant:
                parts = []
            else:
                parts = [(order, xyz)]
        n_kept_from_this = 0
        for p in parts:
            if len(p[1]) >= 2:
                xz = p[1]
                kept.append(p)
                kept_bounds.append((np.min(xz, axis=0), np.max(xz, axis=0)))
                n_kept_from_this += 1
        if n_kept_from_this == 0:
            removed_entirely += 1
        elif n_kept_from_this > 1:
            split_into_multiple += 1

    stats: Dict[str, int] = {
        "input_points_before_merge": n_pts_before,
        "output_points_after_merge": n_pts_after,
        "points_merged_down": max(0, n_pts_before - n_pts_after),
        "input_segments_before_merge": n_seg_before,
        "output_segments_after_merge": n_seg_after,
        "duplicate_segments_collapsed": max(0, n_seg_before - n_seg_after),
        "graph_components_before_merge": graph_components_before_merge,
        "graph_components_after_merge": graph_components_after_merge,
        "merge_same_region_used": 1 if (merge_same_region_only and region_for_merge is not None) else 0,
        "trajectories_extracted": trajectories_extracted,
        "removed_pairwise_overlap": removed_pairwise_overlap,
        "input_trajectories": len(traj),
        "output_polylines": len(kept),
        "removed_entirely": removed_entirely,
        "split_into_multiple": split_into_multiple,
        "distance_uses_4d": 1 if use_4d_distance else 0,
        "point_data_has_region_id": 1 if point_data_has_region_id else 0,
        "region_id_recomputed": region_id_recomputed,
        "line_graph_components_after_recompute": line_graph_components_after_recompute,
        "distinct_region_ids_after_recompute": distinct_region_ids_after_recompute,
    }

    # Rebuild polydata from kept trajectories (VTK_LINE per edge, preserve cell data).
    old_pts = poly.GetPoints()
    old_pd = poly.GetPointData()
    old_cd = poly.GetCellData()

    point_ids: List[int] = []
    seen: Set[int] = set()
    for order, _ in kept:
        for pid in order:
            if pid not in seen:
                seen.add(pid)
                point_ids.append(pid)

    old_to_new = {pid: i for i, pid in enumerate(point_ids)}
    out_pts = vtk.vtkPoints()
    out_pts.SetNumberOfPoints(len(point_ids))
    for old_pid, new_pid in old_to_new.items():
        out_pts.SetPoint(new_pid, old_pts.GetPoint(old_pid))

    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_pd = out_poly.GetPointData()
    for ai in range(old_pd.GetNumberOfArrays()):
        src = old_pd.GetArray(ai)
        if src is None:
            continue
        dst = src.NewInstance()
        dst.DeepCopy(src)
        dst.SetNumberOfTuples(len(point_ids))
        for old_pid, new_pid in old_to_new.items():
            dst.SetTuple(new_pid, old_pid, src)
        out_pd.AddArray(dst)
    if old_pd.GetScalars() is not None and old_pd.GetScalars().GetName():
        out_pd.SetActiveScalars(old_pd.GetScalars().GetName())

    out_lines = vtk.vtkCellArray()
    line = vtk.vtkLine()
    kept_cell_ids: List[int] = []
    for order, _xyz in kept:
        cids = collect_edges_along_order(order, edge_map)
        for k in range(len(order) - 1):
            line.GetPointIds().SetId(0, old_to_new[order[k]])
            line.GetPointIds().SetId(1, old_to_new[order[k + 1]])
            out_lines.InsertNextCell(line)
            kept_cell_ids.append(cids[k])

    out_poly.SetLines(out_lines)
    out_cd = out_poly.GetCellData()
    n_out_cells = len(kept_cell_ids)
    for ai in range(old_cd.GetNumberOfArrays()):
        src = old_cd.GetArray(ai)
        if src is None:
            continue
        dst = src.NewInstance()
        dst.DeepCopy(src)
        dst.SetNumberOfTuples(n_out_cells)
        for new_cid, old_cid in enumerate(kept_cell_ids):
            dst.SetTuple(new_cid, old_cid, src)
        out_cd.AddArray(dst)
    if old_cd.GetScalars() is not None and old_cd.GetScalars().GetName():
        out_cd.SetActiveScalars(old_cd.GetScalars().GetName())

    out_poly.BuildCells()
    out_poly.BuildLinks()
    return out_poly, stats


def main() -> None:
    p = argparse.ArgumentParser(
        description="Remove trajectory polylines that duplicate part of a longer trajectory (spatial tolerance)."
    )
    p.add_argument("-i", "--input-vtp", required=True, help="Input .vtp with trajectories.")
    p.add_argument("-o", "--output-vtp", required=True, help="Output .vtp path.")
    p.add_argument(
        "--tolerance",
        type=float,
        default=1e-4,
        help=(
            "Threshold distance for overlap / subsumption. With --t-array (and without "
            "--spatial-only-distance), this is Euclidean distance in 4D (x,y,z,time_scale*t)."
        ),
    )
    p.add_argument(
        "--t-array",
        type=str,
        default=None,
        help=(
            "Point-data array name for time t. Used to order vertices along each component, "
            "and (unless --spatial-only-distance) as the 4th dimension for all distance tests."
        ),
    )
    p.add_argument(
        "--time-scale",
        type=float,
        default=1.0,
        help="Multiply t by this before 4D distance when --t-array is set (align units with xyz).",
    )
    p.add_argument(
        "--spatial-only-distance",
        action="store_true",
        help="With --t-array, use only (x,y,z) for distances; t is still used for vertex ordering.",
    )
    p.add_argument(
        "--mono-rel",
        type=float,
        default=1e-9,
        help="Relative tolerance for monotonicity of arc-length projections.",
    )
    p.add_argument(
        "--no-merge-nodes",
        action="store_true",
        help="Do not merge vertices that are closer than the node tolerance (only polyline subsumption).",
    )
    p.add_argument(
        "--merge-same-region-only",
        action="store_true",
        help=(
            "When merging close vertices, only merge pairs with the same point-data id "
            "(see --region-array). Use after seperate_diffferent_branches.py so junctions "
            "at identical (x,y,z) on different branches stay split."
        ),
    )
    p.add_argument(
        "--region-array",
        type=str,
        default="RegionId",
        help="Point-data array name for --merge-same-region-only (default: RegionId).",
    )
    p.add_argument(
        "--recompute-region-id",
        action="store_true",
        help=(
            "Before vertex merge: replace RegionId (--region-array) from LINE/POLY_LINE "
            "connected components (0..K-1) and on line cells if present. With merging on, "
            "region-aware merge is enabled automatically so thousands of branches are not "
            "glued at coincident junctions."
        ),
    )
    p.add_argument(
        "--no-cell-region-id",
        action="store_true",
        help="With --recompute-region-id, only overwrite point-data RegionId (not cell data).",
    )
    p.add_argument(
        "--node-tolerance",
        type=float,
        default=None,
        help="Vertex merge distance; defaults to --tolerance if omitted.",
    )
    p.add_argument(
        "--overlap-fraction",
        type=float,
        default=1.0,
        help=(
            "Fraction of shorter polyline samples that must lie within --tolerance of the "
            "longer polyline (1.0 = all samples; use e.g. 0.85 for noisy / partial overlap)."
        ),
    )
    p.add_argument(
        "--no-trim-middle",
        action="store_true",
        help="Only remove whole trajectories fully contained in a longer one (no middle/prefix/suffix trim).",
    )
    p.add_argument(
        "--overlap-length-frac",
        type=float,
        default=0.5,
        help=(
            "If > 0: before the main pass, drop the shorter trajectory in each pair when the "
            "fraction of the shorter's arc length within --tolerance of the longer polyline "
            "is >= this value (default 0.5). Set to 0 to disable."
        ),
    )
    p.add_argument(
        "--data-mode",
        choices=("ascii", "binary", "appended"),
        default="binary",
        help="VTK XML writer data mode.",
    )
    p.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print trajectory removal summary.",
    )
    args = p.parse_args()

    merge_same_region_only = args.merge_same_region_only or (
        args.recompute_region_id and not args.no_merge_nodes
    )

    poly = read_polydata(args.input_vtp)
    out, stats = deduplicate_trajectories(
        poly,
        tolerance=args.tolerance,
        t_array_name=args.t_array,
        mono_rel=args.mono_rel,
        merge_close_nodes=not args.no_merge_nodes,
        node_tolerance=args.node_tolerance,
        overlap_fraction=args.overlap_fraction,
        trim_middle_overlap=not args.no_trim_middle,
        overlap_length_frac_min=args.overlap_length_frac,
        time_scale=args.time_scale,
        spatial_only_distance=args.spatial_only_distance,
        merge_same_region_only=merge_same_region_only,
        region_array_name=args.region_array,
        recompute_region_id=args.recompute_region_id,
        propagate_cell_region_id=not args.no_cell_region_id,
    )
    os.makedirs(os.path.dirname(os.path.abspath(args.output_vtp)) or ".", exist_ok=True)
    write_polydata(out, args.output_vtp, args.data_mode)
    print(f"Wrote {args.output_vtp}")
    if not args.quiet:
        if stats.get("distance_uses_4d"):
            print(
                "Distance metric: 4D Euclidean on (x, y, z, "
                f"{args.time_scale} * t) using point-data array {args.t_array!r}. "
                "Vertex/edge merge still uses xyz only."
            )
        if stats.get("region_id_recomputed"):
            print(
                "Recomputed RegionId from line graph: "
                f"{stats['line_graph_components_after_recompute']} LINE component(s), "
                f"{stats['distinct_region_ids_after_recompute']} distinct point ids "
                f"(includes isolated points)."
            )
        print(
            "Vertex/edge merge: "
            f"{stats['input_points_before_merge']} → {stats['output_points_after_merge']} points "
            f"({stats['points_merged_down']} merged into another); "
            f"{stats['input_segments_before_merge']} → {stats['output_segments_after_merge']} "
            f"line segments ({stats['duplicate_segments_collapsed']} duplicate segments dropped)."
        )
        if stats.get("removed_pairwise_overlap", 0) or stats["trajectories_extracted"] != stats[
            "input_trajectories"
        ]:
            print(
                "Pairwise overlap-length rule: "
                f"{stats['trajectories_extracted']} trajectories → {stats['input_trajectories']} "
                f"({stats['removed_pairwise_overlap']} shorter curves dropped, "
                "shorter arc within --tolerance of a longer one over the length fraction)."
            )
        print(
            "Trajectory pass (after merge): "
            f"{stats['input_trajectories']} line-graph component(s) (ordered polylines) → "
            f"{stats['output_polylines']} polyline(s); "
            f"{stats['removed_entirely']} removed entirely as a duplicate of a longer kept curve; "
            f"{stats['split_into_multiple']} split into multiple polylines (overlap trim)."
        )
        if (
            stats.get("point_data_has_region_id")
            and not stats.get("merge_same_region_used")
            and not args.no_merge_nodes
            and stats["graph_components_before_merge"]
            > max(1, 2 * stats["graph_components_after_merge"])
        ):
            print(
                "Hint: Input has point-data "
                f"{args.region_array!r} but vertex merge ignored it; many branches can share one "
                "(x,y,z), so the LINE graph was reconnected into fewer components. "
                "For one component per branch, add --merge-same-region-only "
                f"(or keep --region-array {args.region_array!r} if renamed), or use --no-merge-nodes."
            )
        if stats["removed_entirely"] == 0:
            print(
                "Why 0 'removed entirely'? Duplicates are often eliminated earlier (vertex/edge merge, "
                "or two copies become one connected component). Pass 2 only counts whole curves "
                "dropped because they lie on a longer kept polyline. Try larger --tolerance if "
                "merge shows no duplicate segments collapsed either."
            )


if __name__ == "__main__":
    main()
