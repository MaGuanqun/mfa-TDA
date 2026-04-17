#!/usr/bin/env pvpython
"""
Assign a connected-component index to every point in a VTP.

This script uses ParaView's `Connectivity` filter (same basic pattern as
`src/contour/count_betti_num.py`) to compute connected components, reads the
resulting `RegionId` array, and writes a new VTP that contains:
  - all original point/cell data arrays
  - plus a per-point array `RegionId` (connected component id).

Expected input VTP:
  - VTK PolyData with `Points` and connectivity given via `Lines` /
    `PolyLines` cells (edges).
"""

from __future__ import annotations

import argparse
import os
from collections import deque

import vtk
from paraview.simple import Connectivity, XMLPolyDataReader, servermanager


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute connected components and write RegionId to each point.")
    p.add_argument("--input-vtp", required=True, help="Input VTP PolyData containing points and edges (Lines/PolyLines).")
    p.add_argument("--output-vtp", required=True, help="Output VTP path with RegionId added to point data.")
    p.add_argument("--region-array-name", default="RegionId", help="Name of RegionId array (default: RegionId).")
    p.add_argument("--data-mode", default="binary", choices=["binary", "ascii"], help="Output VTP data mode (default: binary).")
    return p.parse_args()


def vtk_read_polydata(vtp_path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    return reader.GetOutput()


def vtk_write_polydata(polydata: vtk.vtkPolyData, output_path: str, data_mode: str) -> None:
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_path)
    if data_mode.lower() == "binary":
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAscii()
    writer.SetInputData(polydata)
    writer.Write()


def _compute_region_ids_from_polydata_graph(poly: vtk.vtkPolyData, array_name: str) -> vtk.vtkIntArray:
    """
    Robust fallback that labels connected components directly from input polydata.
    This guarantees one tuple per input point, including isolated points.
    """
    npts = poly.GetNumberOfPoints()
    adjacency = [[] for _ in range(npts)]
    touched = [False] * npts

    id_list = vtk.vtkIdList()
    ncells = poly.GetNumberOfCells()
    for cid in range(ncells):
        ctype = poly.GetCellType(cid)
        # Support both LINE and POLY_LINE connectivity.
        if ctype not in (vtk.VTK_LINE, vtk.VTK_POLY_LINE):
            continue
        poly.GetCellPoints(cid, id_list)
        m = id_list.GetNumberOfIds()
        if m <= 0:
            continue
        # Any point that appears in a line/polyline is part of graph connectivity.
        for i in range(m):
            pid = int(id_list.GetId(i))
            if 0 <= pid < npts:
                touched[pid] = True
        # Consecutive points along a polyline are adjacent.
        for i in range(m - 1):
            p0 = int(id_list.GetId(i))
            p1 = int(id_list.GetId(i + 1))
            if p0 == p1 or p0 < 0 or p1 < 0 or p0 >= npts or p1 >= npts:
                continue
            adjacency[p0].append(p1)
            adjacency[p1].append(p0)

    region = [-1] * npts
    next_region_id = 0

    # Connected components for points used by graph cells.
    for seed in range(npts):
        if not touched[seed] or region[seed] != -1:
            continue
        q = deque([seed])
        region[seed] = next_region_id
        while q:
            u = q.popleft()
            for v in adjacency[u]:
                if region[v] == -1:
                    region[v] = next_region_id
                    q.append(v)
        next_region_id += 1

    # Isolated points become singleton components.
    for pid in range(npts):
        if region[pid] == -1:
            region[pid] = next_region_id
            next_region_id += 1

    arr = vtk.vtkIntArray()
    arr.SetName(array_name)
    arr.SetNumberOfTuples(npts)
    for pid, rid in enumerate(region):
        arr.SetTuple1(pid, rid)
    return arr


def main() -> None:
    args = parse_args()

    if not os.path.exists(args.input_vtp):
        raise FileNotFoundError(f"Input VTP not found: {args.input_vtp}")

    # 1) Use ParaView to compute connectivity / region ids.
    vtp_reader = XMLPolyDataReader(registrationName="Input", FileName=args.input_vtp)
    connectivity = Connectivity(registrationName="Connectivity1", Input=vtp_reader)
    connectivity.UpdatePipeline()

    connectivity_data = servermanager.Fetch(connectivity)
    region_arr = connectivity_data.GetPointData().GetArray(args.region_array_name)
    if region_arr is None:
        raise RuntimeError(
            f"Connectivity filter did not produce point array '{args.region_array_name}'. "
            f"Available point arrays: "
            f"{[connectivity_data.GetPointData().GetArrayName(i) for i in range(connectivity_data.GetPointData().GetNumberOfArrays())]}"
        )

    # 2) Read original PolyData with VTK and attach RegionId (preserves original arrays).
    poly = vtk_read_polydata(args.input_vtp)
    npts = poly.GetNumberOfPoints()

    pd = poly.GetPointData()
    existing = pd.GetArray(args.region_array_name)
    if existing is not None:
        pd.RemoveArray(args.region_array_name)

    tuple_count = region_arr.GetNumberOfTuples()
    if tuple_count == npts:
        # Copy to detach from the ParaView output object.
        new_region = region_arr.NewInstance()
        new_region.DeepCopy(region_arr)
        new_region.SetName(args.region_array_name)
        print(
            f"Using ParaView Connectivity output for '{args.region_array_name}' "
            f"({tuple_count} tuples)."
        )
    else:
        print(
            f"Warning: Connectivity tuple count mismatch ({tuple_count} vs {npts}). "
            "Falling back to direct graph-based component labeling on input polydata."
        )
        new_region = _compute_region_ids_from_polydata_graph(poly, args.region_array_name)

    pd.AddArray(new_region)
    pd.SetActiveScalars(args.region_array_name)

    vtk_write_polydata(poly, args.output_vtp, args.data_mode)

    rid_min, rid_max = new_region.GetRange()
    # region_arr.GetRange() returns (min, max) as floats for some array types.
    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {npts}")
    print(f"RegionId range: {rid_min} .. {rid_max} (components = {int(rid_max) + 1})")


if __name__ == "__main__":
    main()

