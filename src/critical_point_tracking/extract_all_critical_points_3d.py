# from vtkmodules.vtkCommonCore import vtkLogger
# vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_ERROR)

from paraview.simple import *
from paraview import servermanager as sm
import vtk.util.numpy_support as VN
import argparse
import numpy as np
import csv
import sys
import os



def plugin_log(is_server):
    if(is_server==0):
        LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals()) 
    else:
        LoadPlugin("/home/u1435513-gma/apps/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",remote=False,ns=globals())


def compute_critical_points_3d(
    input,
    output,
    scalar=None,
    association="POINTS",
    drop_boundary=True,
    write_type=True,
    write_boundary=True,
):
    vti = XMLImageDataReader(FileName=[input])

    if association.upper() == "CELLS":
        available_arrays = list(vti.CellData.keys())
        vti.CellArrayStatus = available_arrays
        assoc_key = "CELLS"
    else:
        available_arrays = list(vti.PointData.keys())
        vti.PointArrayStatus = available_arrays
        assoc_key = "POINTS"

    vti.TimeArray = 'None'
    vti.UpdatePipeline()

    if not available_arrays:
        raise RuntimeError(f"No {assoc_key.lower()} data arrays found in input: {input}")

    # Your generator (`binary_time_data_convert_3d.py`) writes the scalar as PointData scalars named "scalars".
    if scalar:
        chosen = scalar
    elif "scalars" in available_arrays:
        chosen = "scalars"
    else:
        chosen = available_arrays[0]
    if chosen not in available_arrays:
        raise RuntimeError(
            f"Requested scalar '{chosen}' not found in {assoc_key.lower()} data arrays: {available_arrays}"
        )

    precond = TTKArrayPreconditioning(Input=vti)
    precond.UpdatePipeline()

    tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=precond)

    ttk_cp = TTKScalarFieldCriticalPoints(registrationName='TTKScalarFieldCriticalPoints', Input=tetrahedralize1)
    ttk_cp.ScalarField = [assoc_key, chosen]
    ttk_cp.InputOffsetField = [assoc_key, chosen]
    ttk_cp.UpdatePipeline()

    fetched = sm.Fetch(ttk_cp)
    pts_vtk = fetched.GetPoints()
    if pts_vtk is None or pts_vtk.GetData() is None:
        critical_points_pos = np.empty((0, 3))
    else:
        critical_points_pos = VN.vtk_to_numpy(pts_vtk.GetData())

    criTypeArray_vtk = fetched.GetPointData().GetArray("CriticalType")
    criBoundaryArray_vtk = fetched.GetPointData().GetArray("IsOnBoundary")

    critical_points_type = VN.vtk_to_numpy(criTypeArray_vtk) if criTypeArray_vtk is not None else None
    critical_points_boundary = VN.vtk_to_numpy(criBoundaryArray_vtk) if criBoundaryArray_vtk is not None else None

    if drop_boundary and critical_points_boundary is not None:
        mask = critical_points_boundary != 1
        critical_points_pos = critical_points_pos[mask]
        if critical_points_type is not None:
            critical_points_type = critical_points_type[mask]
        if critical_points_boundary is not None:
        # keep boundary array aligned with filtered points
            critical_points_boundary = critical_points_boundary[mask]

    with open(output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        header = ["PositionX", "PositionY", "PositionZ"]
        if write_type:
            header.append("CriticalType")
        if write_boundary:
            header.append("IsOnBoundary")
        csvwriter.writerow(header)

        for i in range(critical_points_pos.shape[0]):
            row = [
                critical_points_pos[i, 0],
                critical_points_pos[i, 1],
                critical_points_pos[i, 2],
            ]
            if write_type:
                row.append(int(critical_points_type[i]) if critical_points_type is not None else "")
            if write_boundary:
                row.append(int(critical_points_boundary[i]) if critical_points_boundary is not None else "")
            csvwriter.writerow(row)

    
   
    
    
parser = argparse.ArgumentParser(description='TTK-critical points.')


parser.add_argument('-i', '--input_name', type=str, default='file_name.vti', help='input file to compute critical points tracking')
parser.add_argument('-o', '--output_name', type=str, default='file_name.csv', help='output CSV with critical point positions')
parser.add_argument('-s', '--server', type=int, default=0, help='0 for local pc, 1 for remote server')
parser.add_argument('--scalar', type=str, default=None, help='scalar array name to use (default: first array found)')
parser.add_argument('--association', type=str, default="POINTS", choices=["POINTS", "CELLS"], help='where the scalar lives')
parser.add_argument('--keep_boundary', action='store_true', help='do not remove boundary critical points')
parser.add_argument('--positions_only', action='store_true', help='write only XYZ columns (no type/boundary columns)')


args = parser.parse_args()
plugin_log(is_server=args.server)

input_file=args.input_name
output_file=args.output_name

print(output_file)

compute_critical_points_3d(
    input_file,
    output_file,
    scalar=args.scalar,
    association=args.association,
    drop_boundary=(not args.keep_boundary),
    write_type=(not args.positions_only),
    write_boundary=(not args.positions_only),
)