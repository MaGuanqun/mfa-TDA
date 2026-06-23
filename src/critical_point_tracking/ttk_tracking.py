#!/usr/bin/env python
from paraview.simple import *
from paraview import servermanager as sm
import vtk
import vtk.util.numpy_support as VN
import argparse
import numpy as np
import sys
import os


# def plugin_log(is_server):
#     if is_server == 0:
#         LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals())
#     else:
#         LoadPlugin("/home/u1435513-gma/apps/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals())


def spatial_bounds_xy(image, tol=1e-8):
    """Return (xmin, xmax, ymin, ymax) for the spatial domain of a vtkImageData."""
    ox, oy, _ = image.GetOrigin()
    sx, sy, _ = image.GetSpacing()
    nx, ny, _ = image.GetDimensions()
    xmin = ox
    xmax = ox + (nx - 1) * sx
    ymin = oy
    ymax = oy + (ny - 1) * sy
    return xmin, xmax, ymin, ymax, tol


def boundary_point_mask(polydata, xmin, xmax, ymin, ymax, tol):
    """True for points on the spatial domain boundary (x/y only; z is time)."""
    bnd_arr = polydata.GetPointData().GetArray("IsOnBoundary")
    if bnd_arr is not None:
        return VN.vtk_to_numpy(bnd_arr).ravel() == 1

    coords = VN.vtk_to_numpy(polydata.GetPoints().GetData())
    return (
        (np.abs(coords[:, 0] - xmin) <= tol)
        | (np.abs(coords[:, 0] - xmax) <= tol)
        | (np.abs(coords[:, 1] - ymin) <= tol)
        | (np.abs(coords[:, 1] - ymax) <= tol)
    )


def remove_boundary_points(polydata, domain_image):
    """Drop boundary points and line cells that touch them."""
    xmin, xmax, ymin, ymax, tol = spatial_bounds_xy(domain_image)
    on_boundary = boundary_point_mask(polydata, xmin, xmax, ymin, ymax, tol)
    n_pts = polydata.GetNumberOfPoints()
    n_removed = int(on_boundary.sum())
    if n_removed == 0:
        return polydata

    kept_map = [-1] * n_pts
    kept_pts = vtk.vtkPoints()
    for pid in range(n_pts):
        if not on_boundary[pid]:
            kept_map[pid] = kept_pts.InsertNextPoint(polydata.GetPoint(pid))

    point_data = polydata.GetPointData()
    kept_point_arrays = []
    for ai in range(point_data.GetNumberOfArrays()):
        arr = point_data.GetArray(ai)
        kept_arr = arr.NewInstance()
        kept_arr.SetName(arr.GetName())
        kept_arr.SetNumberOfComponents(arr.GetNumberOfComponents())
        kept_point_arrays.append(kept_arr)

    for pid in range(n_pts):
        new_pid = kept_map[pid]
        if new_pid == -1:
            continue
        for ai in range(point_data.GetNumberOfArrays()):
            kept_point_arrays[ai].InsertNextTuple(point_data.GetArray(ai).GetTuple(pid))

    cell_data = polydata.GetCellData()
    original_cell_arrays = [cell_data.GetArray(ai) for ai in range(cell_data.GetNumberOfArrays())]
    kept_cell_arrays = []
    for arr in original_cell_arrays:
        kept_arr = arr.NewInstance()
        kept_arr.SetName(arr.GetName())
        kept_arr.SetNumberOfComponents(arr.GetNumberOfComponents())
        kept_cell_arrays.append(kept_arr)

    kept_lines = vtk.vtkCellArray()
    id_list = vtk.vtkIdList()
    n_cells_removed = 0

    for cid in range(polydata.GetNumberOfCells()):
        cell_type = polydata.GetCellType(cid)
        if cell_type not in (vtk.VTK_LINE, vtk.VTK_POLY_LINE):
            continue

        polydata.GetCellPoints(cid, id_list)
        pt_ids = [id_list.GetId(i) for i in range(id_list.GetNumberOfIds())]
        if any(on_boundary[pid] for pid in pt_ids):
            n_cells_removed += 1
            continue

        new_ids = [kept_map[pid] for pid in pt_ids]
        if cell_type == vtk.VTK_LINE:
            kept_lines.InsertNextCell(2)
            kept_lines.InsertCellPoint(new_ids[0])
            kept_lines.InsertCellPoint(new_ids[1])
        else:
            kept_lines.InsertNextCell(len(new_ids))
            for new_pid in new_ids:
                kept_lines.InsertCellPoint(new_pid)

        for ai, arr in enumerate(original_cell_arrays):
            kept_cell_arrays[ai].InsertNextTuple(arr.GetTuple(cid))

    filtered = vtk.vtkPolyData()
    filtered.SetPoints(kept_pts)
    filtered.SetLines(kept_lines)

    filtered_pd = filtered.GetPointData()
    for kept_arr in kept_point_arrays:
        filtered_pd.AddArray(kept_arr)

    filtered_cd = filtered.GetCellData()
    for kept_arr in kept_cell_arrays:
        filtered_cd.AddArray(kept_arr)

    print(
        f"Removed {n_removed} boundary point(s) and {n_cells_removed} "
        f"line cell(s) touching the domain boundary"
    )
    return filtered


def compute_tracking(input_file, output_file, z_translation=1.0):
    
    print(f"Input file: {input_file}")
    
    # Load the VTI file and select all point arrays (each array = one time step)
    timeTrackingvti = XMLImageDataReader(FileName=[input_file])



    print("read the input file")
    all_point_arrays = list(timeTrackingvti.PointData.keys())
    timeTrackingvti.PointArrayStatus = all_point_arrays
    timeTrackingvti.TimeArray = 'None'
    timeTrackingvti.UpdatePipeline()

    print(f"Point arrays found: {all_point_arrays}")

    # Determine z-translation from FieldData scalar values (same approach as
    # extract_all_critical_points.py which sets z = field_data_array[0] per step).
    # If FieldData carries the actual time values, derive the per-step spacing;
    # otherwise fall back to the user-supplied z_translation argument.
    fetched_raw = sm.Fetch(timeTrackingvti)
    dims = fetched_raw.GetDimensions()  # (nx, ny, nz); nz==1 for 2-D data
    relative_destruction_cost = 2.0 / min(dims[0], dims[1])
    print(f"Grid dimensions: {dims}, Relativedestructioncost: {relative_destruction_cost}")

    if len(all_point_arrays) >= 2:
        fd0 = fetched_raw.GetFieldData().GetArray(all_point_arrays[0])
        fd1 = fetched_raw.GetFieldData().GetArray(all_point_arrays[1])
        if fd0 is not None and fd1 is not None:
            v0 = VN.vtk_to_numpy(fd0)[0]
            v1 = VN.vtk_to_numpy(fd1)[0]
            derived = abs(float(v1) - float(v0))
            if derived > 0.0:
                z_translation = derived
                print(f"Z translation derived from FieldData: {z_translation}")

    # Apply array preconditioning (required for TTK filters)
    precond = TTKArrayPreconditioning(Input=timeTrackingvti)
    precond.UpdatePipeline()
    
    print( "z_translation", z_translation)

    # Tetrahedralize so TTK treats the 2-D image grid as an unstructured mesh.
    # This matches the same step used in extract_all_critical_points.py.
    tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=precond)
    tetrahedralize1.UpdatePipeline()

    # Compute tracking trajectories across all time-step arrays.
    # ForceZtranslation lifts each time step along Z (since the dataset is 2-D
    # and all z-coordinates are 0), using the same spacing convention as
    # extract_all_critical_points.py (z = field_data scalar value per step).
    tTKTrackingFromFields1 = TTKTrackingFromFields(Input=tetrahedralize1)
    tTKTrackingFromFields1.Set(
    Relativedestructioncost=relative_destruction_cost,
    Persistencethreshold = 0,
    ForceZtranslation=1,
    ZTranslation=z_translation,
)
    
    tTKTrackingFromFields1.UpdatePipeline()

    # Extract surface to obtain a clean poly-data output for all critical types
    extractSurface1 = ExtractSurface(Input=tTKTrackingFromFields1)
    extractSurface1.UpdatePipeline()

    tracking_surface = sm.Fetch(extractSurface1)
    tracking_filtered = remove_boundary_points(tracking_surface, fetched_raw)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(tracking_filtered)
    writer.Write()
    print(f"Tracking result saved to: {output_file}")


parser = argparse.ArgumentParser(description='TTK Tracking From Fields.')
parser.add_argument('-i', '--input_name',   type=str,   default='file_name.vti', help='input VTI file containing time-step arrays')
parser.add_argument('-o', '--output_name',  type=str,   default='tracking.vtp',  help='output VTP file for tracking trajectories')
parser.add_argument('-z', '--z_translation', type=float, default=1.0,            help='Z spacing between consecutive time steps (for 2D data)')
parser.add_argument('-s', '--server',       type=int,   default=0,               help='0 = local PC, 1 = remote server')

args = parser.parse_args()

# plugin_log(is_server=args.server)

compute_tracking(args.input_name, args.output_name, args.z_translation)
