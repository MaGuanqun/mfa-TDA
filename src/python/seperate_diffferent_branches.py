#!/usr/bin/env python3
"""
Split branches in a VTP (points + edges) into separate connected components,
then assign a random ColorId (0..59) per component.

Definition used here:
- An "edge" is a line-like cell (VTK_LINE / VTK_POLY_LINE).
- A "branch" is a connected component after we *split* any junction point that
  is incident to more than 2 edges (degree > 2). Splitting is done by replacing
  that shared point with per-edge duplicated points (same coordinates and point
  data), so that the junction no longer glues multiple branches together.

Output:
- Preserves all existing point-data and cell-data arrays.
- Adds:
    - CellData:  ColorId (int, 0..59) on all cells
    - PointData: ColorId (int, 0..59) on all points
  based on connected components computed on the post-split geometry.

Usage:
  python3 seperate_diffferent_branches.py --input-vtp in.vtp --output-vtp out.vtp
"""

from __future__ import annotations

import argparse
import os
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import vtk


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Split junction points (degree>2) so each branch becomes its own connected component, "
            "then assign a random ColorId (0..59) per component and write a new VTP."
        )
    )
    p.add_argument("--input-vtp", "-i", required=True, help="Input .vtp (VTK PolyData) containing points and edges.")
    p.add_argument("--output-vtp", "-o", required=True, help="Output .vtp path.")
    p.add_argument("--data-mode", default="binary", choices=["binary", "ascii"], help="Output VTP data mode.")
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for ColorId mapping (default: 0, deterministic).",
    )
    p.add_argument(
        "--color-mod",
        type=int,
        default=60,
        help="ColorId range is 0..color_mod-1 (default: 60).",
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
    seed: int,
    color_mod: int,
) -> None:
    """
    Add ColorId arrays (point and cell) by mapping each connected component id
    (RegionId) -> random integer in [0, color_mod-1].
    """
    if color_mod <= 0:
        raise ValueError("--color-mod must be > 0.")

    rng = random.Random(seed)

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
    for rid in sorted(region_ids):
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

    src = vtk_read_polydata(args.input_vtp)

    # 1) Split junction points (degree > 2 on edge cells).
    split = split_high_degree_points(src, degree_threshold=2)

    # 2) Connected components on the split geometry.
    labeled = compute_connectivity_regions(split, region_array_name=args.region_array_name)

    # 3) Random ColorId per component (both points and cells).
    add_random_color_ids(
        labeled,
        region_array_name=args.region_array_name,
        color_array_name=args.color_array_name,
        seed=int(args.seed),
        color_mod=int(args.color_mod),
    )

    # 4) Write.
    vtk_write_polydata(labeled, args.output_vtp, args.data_mode)

    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {labeled.GetNumberOfPoints()}  Cells: {labeled.GetNumberOfCells()}")
    print(f"Added arrays: PointData['{args.color_array_name}'], CellData['{args.color_array_name}']")


if __name__ == "__main__":
    main()

