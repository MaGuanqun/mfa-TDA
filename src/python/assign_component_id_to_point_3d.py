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
from vtk.util.numpy_support import vtk_to_numpy
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

    # 1) Read original PolyData with VTK (we preserve this as output geometry/data).
    poly = vtk_read_polydata(args.input_vtp)
    npts = poly.GetNumberOfPoints()

    # Add an explicit point-id array so we can map connectivity output back to original points.
    id_filter = vtk.vtkIdFilter()
    id_filter.SetInputData(poly)
    id_filter.SetPointIds(True)
    id_filter.SetCellIds(False)
    try:
        # VTK 9+
        id_filter.SetPointIdsArrayName("OriginalPointId")
    except AttributeError:
        # Older VTK
        id_filter.SetIdsArrayName("OriginalPointId")
    id_filter.Update()

    tmp_with_ids = vtk.vtkPolyData()
    tmp_with_ids.ShallowCopy(id_filter.GetOutput())

    # Write a temporary input for ParaView so XMLPolyDataReader can consume it.
    tmp_input_path = args.output_vtp + ".tmp_input_with_ids.vtp"
    vtk_write_polydata(tmp_with_ids, tmp_input_path, "binary")

    # 2) Use ParaView to compute connectivity / region ids.
    vtp_reader = XMLPolyDataReader(registrationName="Input", FileName=args.input_vtp)
    vtp_reader.FileName = [tmp_input_path]
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

    orig_pid_arr = connectivity_data.GetPointData().GetArray("OriginalPointId")
    if orig_pid_arr is None:
        raise RuntimeError("Connectivity output is missing 'OriginalPointId' point array.")

    # Build full-length RegionId on original point indexing.
    # We first mark unconnected/orphan points with -1, then remap them to a
    # dedicated positive RegionId after connected components are assigned.
    full_region = vtk.vtkIntArray()
    full_region.SetName(args.region_array_name)
    full_region.SetNumberOfTuples(npts)
    full_region.Fill(-1)

    region_np = vtk_to_numpy(region_arr)
    orig_pid_np = vtk_to_numpy(orig_pid_arr)
    if region_np.shape[0] != orig_pid_np.shape[0]:
        raise RuntimeError(
            f"Connectivity arrays length mismatch: "
            f"{args.region_array_name} has {region_np.shape[0]} tuples but OriginalPointId has {orig_pid_np.shape[0]}."
        )

    for i in range(region_np.shape[0]):
        pid = int(orig_pid_np[i])
        if 0 <= pid < npts:
            full_region.SetValue(pid, int(region_np[i]))

    pd = poly.GetPointData()
    existing = pd.GetArray(args.region_array_name)
    if existing is not None:
        pd.RemoveArray(args.region_array_name)

    # Copy to detach and preserve the original input point indexing.
    new_region = full_region.NewInstance()
    new_region.DeepCopy(full_region)
    new_region.SetName(args.region_array_name)

    orphan_count = 0
    for i in range(npts):
        if int(new_region.GetValue(i)) < 0:
            orphan_count += 1

    orphan_region_id = None
    if orphan_count > 0:
        connected_max = int(region_arr.GetRange()[1]) if region_arr.GetNumberOfTuples() > 0 else -1
        orphan_region_id = connected_max + 1
        for i in range(npts):
            if int(new_region.GetValue(i)) < 0:
                new_region.SetValue(i, orphan_region_id)

    # Final verification: output must not contain RegionId == -1.
    remaining_negative = 0
    for i in range(npts):
        if int(new_region.GetValue(i)) < 0:
            remaining_negative += 1
    if remaining_negative > 0:
        raise RuntimeError(
            f"Verification failed: found {remaining_negative} points with RegionId < 0 "
            f"after orphan remapping."
        )

    pd.AddArray(new_region)
    pd.SetActiveScalars(args.region_array_name)

    vtk_write_polydata(poly, args.output_vtp, args.data_mode)
    try:
        os.remove(tmp_input_path)
    except OSError:
        pass

    rid_min, rid_max = new_region.GetRange()

    # region_arr.GetRange() returns (min, max) as floats for some array types.
    print(f"Wrote: {args.output_vtp}")
    print(f"Points: {npts}")
    print(f"RegionId range: {rid_min} .. {rid_max}")
    print(f"Orphan points (unconnected): {orphan_count}")
    print("Verification: RegionId == -1 count: 0")
    if orphan_region_id is not None:
        print(f"Assigned orphan RegionId: {orphan_region_id}")
    if rid_max >= 0:
        print(f"Total components (including orphan group if present): {int(rid_max) + 1}")


if __name__ == "__main__":
    main()

