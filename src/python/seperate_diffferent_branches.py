#!/usr/bin/env python3
"""
Split branches in a VTP (points + edges) into separate connected components,
then assign ColorId per component: if at most 20 branches, a random
permutation of 0..n-1; if more than 20, random integers in 0..59 (may repeat).

Definition used here:
- An "edge" is a line-like cell (VTK_LINE / VTK_POLY_LINE).
- A "branch" is a connected component after we *split* any junction point that
  is incident to more than 2 edges (degree > 2). Splitting is done by replacing
  that shared point with per-edge duplicated points (same coordinates and point
  data), so that the junction no longer glues multiple branches together.

Pipeline order: junction split → explode ``POLY_LINE`` to ``LINE`` → break all
**loops (cycles)** in the line graph until none remain: for each cycle, pick two
anchor vertices by **min** and **max** of point array ``t`` (name configurable)
when that array exists, otherwise by **min** and **max z**. Duplicate those
vertices so the two arcs between the anchors become two disjoint open paths.
Tie-breaking for ties on min/max uses the same RNG as ``ColorId`` (``--seed``).

Output:
- Preserves all existing point-data and cell-data arrays.
- Adds:
    - CellData:  ColorId (int; 0..n-1 if n<=20 components, else 0..59 by default)
    - PointData: ColorId (same)
  based on connected components computed on the post-split geometry.

Usage:
  python3 seperate_diffferent_branches.py --input-vtp in.vtp --output-vtp out.vtp
"""

from __future__ import annotations

import argparse
import os
import random
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import vtk

# Above this many connected components, ColorId is drawn randomly in [0, color_mod-1] (default 0..59).
_PERMUTATION_MAX_COMPONENTS = 20


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Split junction points (degree>2) so each branch becomes its own connected component, "
            "then assign ColorId per branch (<=20: random permutation 0..n-1; >20: random 0..59) and write a new VTP."
        )
    )
    p.add_argument("--input-vtp", "-i", required=True, help="Input .vtp (VTK PolyData) containing points and edges.")
    p.add_argument("--output-vtp", "-o", required=True, help="Output .vtp path.")
    p.add_argument("--data-mode", default="binary", choices=["binary", "ascii"], help="Output VTP data mode.")
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Seed for all randomness (ColorId, tie-breaks when breaking loops; default: 0).",
    )
    p.add_argument(
        "--t-array-name",
        default="t",
        help="Point-data array name for time/scalar when breaking loops (default: t). If missing, use z coordinate.",
    )
    p.add_argument(
        "--color-mod",
        type=int,
        default=48,
        help="ColorId range is 0..color_mod-1 (default: 48).",
    )
    p.add_argument(
        "--region-array-name",
        default="RegionId",
        help="Connectivity region array name to use/produce (default: RegionId).",
    )
    p.add_argument(
        "--color-array-name",
        default="ColorId",
        help="Output array name for random component color ids (default: ColorId).",
    )
    return p.parse_args()


def vtk_read_polydata(vtp_path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    out = vtk.vtkPolyData()
    out.ShallowCopy(reader.GetOutput())
    return out


def vtk_write_polydata(polydata: vtk.vtkPolyData, output_path: str, data_mode: str) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    if data_mode.lower() == "binary":
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAscii()
    writer.SetInputData(polydata)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write VTP: {output_path}")


def _is_edge_cell_type(cell_type: int) -> bool:
    return cell_type in (vtk.VTK_LINE, vtk.VTK_POLY_LINE)


def compute_point_edge_degree(poly: vtk.vtkPolyData) -> List[int]:
    """
    Degree = number of incident edge-cells (line / polyline) for each point.
    """
    npts = poly.GetNumberOfPoints()
    deg = [0] * npts
    ncells = poly.GetNumberOfCells()
    idlist = vtk.vtkIdList()
    for cid in range(ncells):
        ctype = poly.GetCellType(cid)
        if not _is_edge_cell_type(ctype):
            continue
        poly.GetCellPoints(cid, idlist)
        seen: set[int] = set()
        for j in range(idlist.GetNumberOfIds()):
            pid = int(idlist.GetId(j))
            # Avoid double-counting if a polyline repeats a point id (rare).
            if pid not in seen:
                deg[pid] += 1
                seen.add(pid)
    return deg


def _new_like_array(arr: vtk.vtkDataArray) -> vtk.vtkDataArray:
    new_arr = arr.NewInstance()
    new_arr.SetName(arr.GetName())
    new_arr.SetNumberOfComponents(arr.GetNumberOfComponents())
    return new_arr


def _copy_point_tuple(src_pd: vtk.vtkPointData, dst_pd: vtk.vtkPointData, src_pid: int) -> None:
    na = src_pd.GetNumberOfArrays()
    for ai in range(na):
        src_arr = src_pd.GetArray(ai)
        if src_arr is None:
            continue
        dst_arr = dst_pd.GetArray(src_arr.GetName())
        if dst_arr is None:
            raise RuntimeError(f"Missing destination point array '{src_arr.GetName()}'.")
        ncomp = src_arr.GetNumberOfComponents()
        tup = [0.0] * ncomp
        src_arr.GetTuple(src_pid, tup)
        dst_arr.InsertNextTuple(tup)


def _copy_cell_tuple(src_cd: vtk.vtkCellData, dst_cd: vtk.vtkCellData, src_cid: int) -> None:
    na = src_cd.GetNumberOfArrays()
    for ai in range(na):
        src_arr = src_cd.GetArray(ai)
        if src_arr is None:
            continue
        dst_arr = dst_cd.GetArray(src_arr.GetName())
        if dst_arr is None:
            raise RuntimeError(f"Missing destination cell array '{src_arr.GetName()}'.")
        ncomp = src_arr.GetNumberOfComponents()
        tup = [0.0] * ncomp
        src_arr.GetTuple(src_cid, tup)
        dst_arr.InsertNextTuple(tup)


@dataclass(frozen=True)
class _PointOcc:
    cell_id: int
    local_index: int  # index within the cell's point-id list


def split_high_degree_points(poly: vtk.vtkPolyData, degree_threshold: int = 2) -> vtk.vtkPolyData:
    """
    For any point with edge-degree > degree_threshold, duplicate that point per
    incident cell occurrence and rewrite cell connectivity to use the duplicates.

    This removes junction sharing, so branches become separate components.
    """
    npts = poly.GetNumberOfPoints()
    if npts == 0:
        out = vtk.vtkPolyData()
        out.ShallowCopy(poly)
        return out

    deg = compute_point_edge_degree(poly)
    to_split = [d > degree_threshold for d in deg]
    if not any(to_split):
        out = vtk.vtkPolyData()
        out.DeepCopy(poly)
        return out

    src_pts = poly.GetPoints()
    src_pd = poly.GetPointData()
    src_cd = poly.GetCellData()

    out_pts = vtk.vtkPoints()
    out_pts.SetDataType(src_pts.GetDataType())

    out_pd = vtk.vtkPointData()
    out_cd = vtk.vtkCellData()

    # Prepare destination point arrays (same schema).
    for ai in range(src_pd.GetNumberOfArrays()):
        a = src_pd.GetArray(ai)
        if a is None or a.GetName() is None:
            continue
        out_pd.AddArray(_new_like_array(a))

    # Prepare destination cell arrays (same schema).
    for ai in range(src_cd.GetNumberOfArrays()):
        a = src_cd.GetArray(ai)
        if a is None or a.GetName() is None:
            continue
        out_cd.AddArray(_new_like_array(a))

    old_to_new_shared: Dict[int, int] = {}

    def add_point_from_old(old_pid: int) -> int:
        x, y, z = src_pts.GetPoint(old_pid)
        new_pid = out_pts.InsertNextPoint(float(x), float(y), float(z))
        _copy_point_tuple(src_pd, out_pd, old_pid)
        return int(new_pid)

    # Rebuild ALL cell arrays while rewriting point ids as needed.
    out_verts = vtk.vtkCellArray()
    out_lines = vtk.vtkCellArray()
    out_polys = vtk.vtkCellArray()
    out_strips = vtk.vtkCellArray()

    idlist = vtk.vtkIdList()
    ncells = poly.GetNumberOfCells()

    for cid in range(ncells):
        poly.GetCellPoints(cid, idlist)
        nids = idlist.GetNumberOfIds()
        if nids <= 0:
            continue

        new_ids = [0] * nids
        for j in range(nids):
            old_pid = int(idlist.GetId(j))
            if 0 <= old_pid < npts and to_split[old_pid]:
                # Split point: create a fresh duplicate for this occurrence.
                new_ids[j] = add_point_from_old(old_pid)
            else:
                # Not split: reuse one shared copy to preserve connectivity.
                mapped = old_to_new_shared.get(old_pid)
                if mapped is None:
                    mapped = add_point_from_old(old_pid)
                    old_to_new_shared[old_pid] = mapped
                new_ids[j] = mapped

        ctype = poly.GetCellType(cid)
        if ctype == vtk.VTK_VERTEX:
            out_verts.InsertNextCell(1)
            out_verts.InsertCellPoint(new_ids[0])
        elif ctype == vtk.VTK_POLY_VERTEX:
            out_verts.InsertNextCell(nids)
            for pid in new_ids:
                out_verts.InsertCellPoint(pid)
        elif ctype in (vtk.VTK_LINE, vtk.VTK_POLY_LINE):
            out_lines.InsertNextCell(nids)
            for pid in new_ids:
                out_lines.InsertCellPoint(pid)
        elif ctype == vtk.VTK_TRIANGLE_STRIP:
            out_strips.InsertNextCell(nids)
            for pid in new_ids:
                out_strips.InsertCellPoint(pid)
        else:
            # Fallback: keep as polys (works for triangles/quads/polygons).
            out_polys.InsertNextCell(nids)
            for pid in new_ids:
                out_polys.InsertCellPoint(pid)

        # Preserve cell data tuple aligned with the output cell ordering.
        _copy_cell_tuple(src_cd, out_cd, cid)

    out = vtk.vtkPolyData()
    out.SetPoints(out_pts)
    out.SetVerts(out_verts)
    out.SetLines(out_lines)
    out.SetPolys(out_polys)
    out.SetStrips(out_strips)
    out.GetPointData().ShallowCopy(out_pd)
    out.GetCellData().ShallowCopy(out_cd)

    # Ensure internal links are consistent (needed for some downstream filters).
    out.BuildCells()
    out.BuildLinks()
    return out


def _expand_line_segments(poly: vtk.vtkPolyData) -> List[Tuple[int, int, int]]:
    """
    Each consecutive pair in a LINE / POLY_LINE cell becomes one undirected segment.
    Returns list of (cell_id, p, q) with p != q.
    """
    segs: List[Tuple[int, int, int]] = []
    idlist = vtk.vtkIdList()
    ncells = poly.GetNumberOfCells()
    for cid in range(ncells):
        ctype = poly.GetCellType(cid)
        if not _is_edge_cell_type(ctype):
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


def _edge_key(p: int, q: int) -> Tuple[int, int]:
    return (p, q) if p < q else (q, p)


def _seg_sig(cid: int, p: int, q: int) -> Tuple[int, int, int]:
    """Canonical segment id for membership in sets (cell id + sorted endpoints)."""
    if p <= q:
        return (cid, p, q)
    return (cid, q, p)


def _append_vtk_cell(ctype: int, ids: List[int], va: vtk.vtkCellArray, la: vtk.vtkCellArray, pa: vtk.vtkCellArray, sa: vtk.vtkCellArray) -> None:
    if ctype == vtk.VTK_VERTEX:
        va.InsertNextCell(1)
        va.InsertCellPoint(ids[0])
    elif ctype == vtk.VTK_POLY_VERTEX:
        va.InsertNextCell(len(ids))
        for x in ids:
            va.InsertCellPoint(x)
    elif ctype == vtk.VTK_LINE:
        if len(ids) != 2:
            raise ValueError("LINE cell must have exactly 2 point ids")
        la.InsertNextCell(2)
        la.InsertCellPoint(ids[0])
        la.InsertCellPoint(ids[1])
    elif ctype == vtk.VTK_TRIANGLE:
        pa.InsertNextCell(3)
        for x in ids[:3]:
            pa.InsertCellPoint(x)
    elif ctype == vtk.VTK_QUAD:
        pa.InsertNextCell(4)
        for x in ids[:4]:
            pa.InsertCellPoint(x)
    elif ctype == vtk.VTK_POLYGON:
        pa.InsertNextCell(len(ids))
        for x in ids:
            pa.InsertCellPoint(x)
    elif ctype == vtk.VTK_TRIANGLE_STRIP:
        sa.InsertNextCell(len(ids))
        for x in ids:
            sa.InsertCellPoint(x)
    else:
        pa.InsertNextCell(len(ids))
        for x in ids:
            pa.InsertCellPoint(x)


def rebuild_explode_polylines(poly: vtk.vtkPolyData) -> vtk.vtkPolyData:
    """
    Turn each POLY_LINE into consecutive VTK_LINE cells (cell data duplicated).
    Other cell types are preserved. Needed so cycle cuts never assign two
    different copies of the same junction to one polyline cell.
    """
    idlist = vtk.vtkIdList()
    ncells = poly.GetNumberOfCells()
    rows: List[Tuple[int, List[int], int]] = []

    for cid in range(ncells):
        ctype = poly.GetCellType(cid)
        poly.GetCellPoints(cid, idlist)
        n = idlist.GetNumberOfIds()
        ids = [int(idlist.GetId(i)) for i in range(n)]
        if ctype == vtk.VTK_POLY_LINE and n >= 3:
            for k in range(n - 1):
                rows.append((vtk.VTK_LINE, [ids[k], ids[k + 1]], cid))
        elif ctype == vtk.VTK_POLY_LINE and n == 2:
            rows.append((vtk.VTK_LINE, ids, cid))
        else:
            rows.append((ctype, ids, cid))

    src_cd = poly.GetCellData()
    out = vtk.vtkPolyData()
    out.SetPoints(poly.GetPoints())
    out.GetPointData().ShallowCopy(poly.GetPointData())

    va = vtk.vtkCellArray()
    la = vtk.vtkCellArray()
    pa = vtk.vtkCellArray()
    sa = vtk.vtkCellArray()
    out_cd = vtk.vtkCellData()
    for ai in range(src_cd.GetNumberOfArrays()):
        a = src_cd.GetArray(ai)
        if a is None or a.GetName() is None:
            continue
        out_cd.AddArray(_new_like_array(a))

    for ct, ids, src_cid in rows:
        _append_vtk_cell(ct, ids, va, la, pa, sa)
        _copy_cell_tuple(src_cd, out_cd, src_cid)

    out.SetVerts(va)
    out.SetLines(la)
    out.SetPolys(pa)
    out.SetStrips(sa)
    out.GetCellData().ShallowCopy(out_cd)
    out.BuildCells()
    out.BuildLinks()
    return out


def _find_parallel_pair(segments: List[Tuple[int, int, int]]) -> Optional[Tuple[int, int]]:
    """If two or more segments share the same undirected edge, return (u, v) with u < v."""
    cnt: Dict[Tuple[int, int], int] = defaultdict(int)
    for _, p, q in segments:
        cnt[_edge_key(p, q)] += 1
    for (u, v), c in cnt.items():
        if c >= 2:
            return (u, v)
    return None


def _find_cycle_dfs(
    n: int,
    adj: List[List[int]],
) -> Optional[List[int]]:
    """
    Find one simple cycle; return vertex list [v0, v1, ..., vm] such that
    consecutive pairs and (vm, v0) are edges (unique-edge graph).
    """
    color = [0] * n  # 0 white, 1 gray, 2 black
    parent = [-1] * n

    def dfs(u: int, pu: int) -> Optional[List[int]]:
        color[u] = 1
        for v in adj[u]:
            if v == pu:
                continue
            if color[v] == 0:
                parent[v] = u
                cyc = dfs(v, u)
                if cyc is not None:
                    return cyc
            elif color[v] == 1:
                # Back edge (u, v); v is ancestor -> cycle along parent chain u -> ... -> v
                path_uv: List[int] = []
                x = u
                while x != v and x != -1:
                    path_uv.append(x)
                    x = parent[x]
                if x != v:
                    continue
                path_uv.append(v)
                # Cyclic order: v ... u, closing edge (u, v)
                return list(reversed(path_uv))
        color[u] = 2
        return None

    for s in range(n):
        if color[s] != 0:
            continue
        parent[s] = -1
        cyc = dfs(s, -1)
        if cyc is not None:
            return cyc
    return None


def _cycle_edges_from_vertices(cycle_verts: List[int]) -> Set[Tuple[int, int]]:
    """Undirected edge keys for a closed walk given ordered unique vertices on the cycle."""
    L = len(cycle_verts)
    if L < 2:
        return set()
    keys: Set[Tuple[int, int]] = set()
    for i in range(L):
        p = cycle_verts[i]
        q = cycle_verts[(i + 1) % L]
        keys.add(_edge_key(p, q))
    return keys


def _arc_edge_keys(cycle_verts: List[int], ia: int, ib: int) -> Set[Tuple[int, int]]:
    """Edges on the forward walk from cycle_verts[ia] to cycle_verts[ib] (exclusive of wrapping past ib)."""
    L = len(cycle_verts)
    keys: Set[Tuple[int, int]] = set()
    i = ia
    while i != ib:
        j = (i + 1) % L
        keys.add(_edge_key(cycle_verts[i], cycle_verts[j]))
        i = j
    return keys


def _pick_cycle_anchors(
    cycle_verts: List[int],
    poly: vtk.vtkPolyData,
    t_arr: Optional[vtk.vtkDataArray],
    rng: random.Random,
) -> Tuple[int, int]:
    """Distinct vertices a, b on the cycle for min/max scalar (t or z)."""
    pts = poly.GetPoints()
    uniq = list(dict.fromkeys(cycle_verts))  # preserve order, unique

    def scalar_at(pid: int) -> float:
        if t_arr is not None and 0 <= pid < t_arr.GetNumberOfTuples():
            return float(t_arr.GetTuple1(pid))
        x, y, z = pts.GetPoint(pid)
        return float(z)

    vals = {pid: scalar_at(pid) for pid in uniq}
    mn = min(vals.values())
    mx = max(vals.values())
    cand_min = [p for p, v in vals.items() if v == mn]
    cand_max = [p for p, v in vals.items() if v == mx]
    a = int(rng.choice(cand_min))
    b_choices = [p for p in cand_max if p != a]
    if not b_choices:
        b_choices = cand_max
    b = int(rng.choice(b_choices))
    if a == b and len(uniq) >= 2:
        others = [p for p in uniq if p != a]
        b = int(rng.choice(others))
    return a, b


def _duplicate_point(
    poly: vtk.vtkPolyData,
    old_pid: int,
) -> int:
    """Append a copy of point old_pid (coordinates + all point data); return new index."""
    pts = poly.GetPoints()
    pd = poly.GetPointData()
    x, y, z = pts.GetPoint(old_pid)
    new_id = pts.InsertNextPoint(float(x), float(y), float(z))
    for ai in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(ai)
        if arr is None or arr.GetName() is None:
            continue
        ncomp = arr.GetNumberOfComponents()
        tup = [0.0] * ncomp
        arr.GetTuple(old_pid, tup)
        arr.InsertNextTuple(tup)
    return int(new_id)


def _remap_segment_endpoints(
    p: int,
    q: int,
    a: int,
    b: int,
    a0: int,
    a1: int,
    b0: int,
    b1: int,
    in_arc0: bool,
    in_arc1: bool,
) -> Tuple[int, int]:
    def map_one(x: int, arc0: bool, arc1: bool) -> int:
        if x == a:
            if arc0:
                return a0
            if arc1:
                return a1
            return a0
        if x == b:
            if arc0:
                return b0
            if arc1:
                return b1
            return b0
        return x

    return map_one(p, in_arc0, in_arc1), map_one(q, in_arc0, in_arc1)


def _replace_poly_cell_point_ids(poly: vtk.vtkPolyData, cell_id: int, new_ids: List[int]) -> None:
    idl = vtk.vtkIdList()
    idl.SetNumberOfIds(len(new_ids))
    for i, pid in enumerate(new_ids):
        idl.SetId(i, int(pid))
    if hasattr(poly, "ReplaceCell"):
        poly.ReplaceCell(cell_id, idl)
        return
    raise RuntimeError(
        "This VTK build has no vtkPolyData.ReplaceCell; upgrade VTK or rebuild cells manually."
    )


def break_one_cycle(poly: vtk.vtkPolyData, rng: random.Random, t_array_name: str) -> bool:
    """
    If the line graph contains a cycle, break one cycle into two paths by
    duplicating the min/max (t or z) vertices. Returns True if a cycle was broken.
    """
    n = poly.GetNumberOfPoints()
    if n < 2:
        return False

    segments = _expand_line_segments(poly)
    if not segments:
        return False

    par = _find_parallel_pair(segments)
    cycle_verts: Optional[List[int]] = None
    arc0_seg: Set[Tuple[int, int, int]] = set()
    arc1_seg: Set[Tuple[int, int, int]] = set()
    arc0_edges: Set[Tuple[int, int]] = set()
    arc1_edges: Set[Tuple[int, int]] = set()

    if par is not None:
        u, v = par
        cycle_verts = [u, v]
        matching = sorted(
            _seg_sig(cid, p, q) for cid, p, q in segments if _edge_key(p, q) == (u, v)
        )
        if len(matching) < 2:
            return False
        mid = len(matching) // 2
        arc0_seg = set(matching[:mid])
        arc1_seg = set(matching[mid:])
    else:
        adj_unique: List[List[int]] = [[] for _ in range(n)]
        seen_e: Set[Tuple[int, int]] = set()
        for _, p, q in segments:
            ek = _edge_key(p, q)
            if ek in seen_e:
                continue
            seen_e.add(ek)
            adj_unique[p].append(q)
            adj_unique[q].append(p)
        cycle_verts = _find_cycle_dfs(n, adj_unique)

    if cycle_verts is None or len(cycle_verts) < 2:
        return False

    L = len(cycle_verts)
    if cycle_verts[0] == cycle_verts[-1]:
        cycle_verts = cycle_verts[:-1]
        L = len(cycle_verts)
    if L < 2:
        return False

    pd = poly.GetPointData()
    t_arr = pd.GetArray(t_array_name)
    a, b = _pick_cycle_anchors(cycle_verts, poly, t_arr, rng)

    if not arc0_seg and not arc1_seg:
        ia = cycle_verts.index(a)
        ib = cycle_verts.index(b)
        arc0_edges = _arc_edge_keys(cycle_verts, ia, ib)
        arc1_edges = _cycle_edges_from_vertices(cycle_verts) - arc0_edges

    idlist = vtk.vtkIdList()
    ncells = poly.GetNumberOfCells()
    to_update: List[Tuple[int, int, int, bool, bool]] = []
    for cid in range(ncells):
        if poly.GetCellType(cid) != vtk.VTK_LINE:
            continue
        poly.GetCellPoints(cid, idlist)
        if idlist.GetNumberOfIds() != 2:
            continue
        p = int(idlist.GetId(0))
        q = int(idlist.GetId(1))
        sig = _seg_sig(cid, p, q)
        if arc0_seg or arc1_seg:
            in0 = sig in arc0_seg
            in1 = sig in arc1_seg
        else:
            ek = _edge_key(p, q)
            in0 = ek in arc0_edges
            in1 = ek in arc1_edges
        if in0 or in1:
            to_update.append((cid, p, q, in0, in1))

    if not to_update:
        return False

    a0 = _duplicate_point(poly, a)
    a1 = _duplicate_point(poly, a)
    b0 = _duplicate_point(poly, b)
    b1 = _duplicate_point(poly, b)

    for cid, p, q, in0, in1 in to_update:
        p2, q2 = _remap_segment_endpoints(p, q, a, b, a0, a1, b0, b1, in0, in1)
        _replace_poly_cell_point_ids(poly, cid, [p2, q2])
    poly.Modified()
    poly.BuildCells()
    poly.BuildLinks()
    return True


def break_all_loops(poly: vtk.vtkPolyData, rng: random.Random, t_array_name: str) -> vtk.vtkPolyData:
    """Repeatedly break cycles until the line graph is acyclic (no parallel 2-cycles, no graph cycles)."""
    out = vtk.vtkPolyData()
    out.DeepCopy(poly)
    out = rebuild_explode_polylines(out)
    safety = 0
    max_iter = max(1000, out.GetNumberOfCells() * 10)
    while safety < max_iter:
        safety += 1
        if not break_one_cycle(out, rng, t_array_name):
            break
    return out


def compute_connectivity_regions(poly: vtk.vtkPolyData, region_array_name: str) -> vtk.vtkPolyData:
    """
    Compute connected components (regions) on the polydata. Produces a RegionId
    array (int) in both CellData and PointData.
    """
    conn = vtk.vtkConnectivityFilter()
    conn.SetInputData(poly)
    conn.SetExtractionModeToAllRegions()
    conn.ColorRegionsOn()
    conn.Update()

    out = vtk.vtkPolyData()
    out.ShallowCopy(conn.GetOutput())

    # vtkConnectivityFilter uses "RegionId" as the default name.
    # If user asked for a different name, rename in both point/cell data.
    if region_array_name != "RegionId":
        for data_obj in (out.GetPointData(), out.GetCellData()):
            arr = data_obj.GetArray("RegionId")
            if arr is not None:
                arr.SetName(region_array_name)
    return out


def add_random_color_ids(
    poly: vtk.vtkPolyData,
    region_array_name: str,
    color_array_name: str,
    rng: random.Random,
    color_mod: int,
) -> None:
    """
    Add ColorId arrays (point and cell) based on connected components.

    Policy:
    - If number of components K <= 20, assign ColorId as a random permutation of
      0..K-1 (each branch gets a unique index; order depends on --seed).
    - If K > 20, assign random ColorId in [0, color_mod-1] per component (default 0..59).
    """
    if color_mod <= 0:
        raise ValueError("--color-mod must be > 0.")

    pd = poly.GetPointData()
    cd = poly.GetCellData()

    region_cells = cd.GetArray(region_array_name)
    region_points = pd.GetArray(region_array_name)
    if region_cells is None and region_points is None:
        raise RuntimeError(f"Missing '{region_array_name}' on both points and cells.")

    # Build mapping from observed region ids.
    region_ids: set[int] = set()
    if region_cells is not None:
        for i in range(region_cells.GetNumberOfTuples()):
            region_ids.add(int(region_cells.GetTuple1(i)))
    if region_points is not None:
        for i in range(region_points.GetNumberOfTuples()):
            region_ids.add(int(region_points.GetTuple1(i)))

    reg_to_color: Dict[int, int] = {}
    sorted_regions = sorted(region_ids)
    k = len(sorted_regions)
    if k <= _PERMUTATION_MAX_COMPONENTS:
        perm = list(range(k))
        rng.shuffle(perm)
        for idx, rid in enumerate(sorted_regions):
            reg_to_color[rid] = int(perm[idx])
    else:
        for rid in sorted_regions:
            reg_to_color[rid] = rng.randrange(color_mod)

    # Remove any existing ColorId arrays to avoid duplicates.
    if pd.GetArray(color_array_name) is not None:
        pd.RemoveArray(color_array_name)
    if cd.GetArray(color_array_name) is not None:
        cd.RemoveArray(color_array_name)

    col_p = vtk.vtkIntArray()
    col_p.SetName(color_array_name)
    col_p.SetNumberOfComponents(1)
    col_p.SetNumberOfTuples(poly.GetNumberOfPoints())

    col_c = vtk.vtkIntArray()
    col_c.SetName(color_array_name)
    col_c.SetNumberOfComponents(1)
    col_c.SetNumberOfTuples(poly.GetNumberOfCells())

    # Cell colors
    if region_cells is not None and region_cells.GetNumberOfTuples() == poly.GetNumberOfCells():
        for i in range(poly.GetNumberOfCells()):
            rid = int(region_cells.GetTuple1(i))
            col_c.SetValue(i, int(reg_to_color.get(rid, 0)))
    else:
        # Fallback: if cell RegionId missing, set 0
        for i in range(poly.GetNumberOfCells()):
            col_c.SetValue(i, 0)

    # Point colors
    if region_points is not None and region_points.GetNumberOfTuples() == poly.GetNumberOfPoints():
        for i in range(poly.GetNumberOfPoints()):
            rid = int(region_points.GetTuple1(i))
            col_p.SetValue(i, int(reg_to_color.get(rid, 0)))
    else:
        # Fallback: derive point color from any incident cell's RegionId.
        poly.BuildLinks()
        cell_ids = vtk.vtkIdList()
        for pid in range(poly.GetNumberOfPoints()):
            poly.GetPointCells(pid, cell_ids)
            if cell_ids.GetNumberOfIds() > 0 and region_cells is not None:
                cid = int(cell_ids.GetId(0))
                rid = int(region_cells.GetTuple1(cid))
                col_p.SetValue(pid, int(reg_to_color.get(rid, 0)))
            else:
                col_p.SetValue(pid, 0)

    pd.AddArray(col_p)
    cd.AddArray(col_c)
    pd.SetActiveScalars(color_array_name)


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_vtp):
        raise FileNotFoundError(f"Input VTP not found: {args.input_vtp}")

    rng = random.Random(int(args.seed))
    random.seed(int(args.seed))

    src = vtk_read_polydata(args.input_vtp)

    # 1) Split junction points (degree > 2 on edge cells).
    split = split_high_degree_points(src, degree_threshold=2)

    # 2) Break all graph loops into two open trajectories (min/max t or z anchors).
    no_loops = break_all_loops(split, rng, args.t_array_name)

    # 3) Connected components on the split geometry.
    labeled = compute_connectivity_regions(no_loops, region_array_name=args.region_array_name)

    # 4) ColorId per component (both points and cells); all randomness uses rng / random.seed above.
    add_random_color_ids(
        labeled,
        region_array_name=args.region_array_name,
        color_array_name=args.color_array_name,
        rng=rng,
        color_mod=int(args.color_mod),
    )

    # 5) Write.
    vtk_write_polydata(labeled, args.output_vtp, args.data_mode)

    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {labeled.GetNumberOfPoints()}  Cells: {labeled.GetNumberOfCells()}")
    print(f"Added arrays: PointData['{args.color_array_name}'], CellData['{args.color_array_name}']")


if __name__ == "__main__":
    main()

