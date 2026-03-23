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
    if region_arr.GetNumberOfTuples() != npts:
        raise RuntimeError(
            f"RegionId tuple count mismatch: RegionId has {region_arr.GetNumberOfTuples()} tuples, "
            f"but input PolyData has {npts} points."
        )

    pd = poly.GetPointData()
    existing = pd.GetArray(args.region_array_name)
    if existing is not None:
        pd.RemoveArray(args.region_array_name)

    # Copy to detach from the ParaView output object.
    new_region = region_arr.NewInstance()
    new_region.DeepCopy(region_arr)
    new_region.SetName(args.region_array_name)
    pd.AddArray(new_region)
    pd.SetActiveScalars(args.region_array_name)

    vtk_write_polydata(poly, args.output_vtp, args.data_mode)

    rid_min, rid_max = region_arr.GetRange()
    # region_arr.GetRange() returns (min, max) as floats for some array types.
    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {npts}")
    print(f"RegionId range: {rid_min} .. {rid_max} (components = {int(rid_max) + 1})")


if __name__ == "__main__":
    main()

