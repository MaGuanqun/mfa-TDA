#!/usr/bin/env pvpython
import os
import argparse
import numpy as np
import csv
from typing import Optional

import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from paraview.simple import *  # noqa: F401,F403
from paraview import servermanager as sm


# ---------------------------------------------------------------------
# TTK plugin loading
# ---------------------------------------------------------------------
def plugin_log(is_server: int):
    """
    Load TTK plugin on local machine (0) or remote server (1).
    Adjust paths to your ParaView installations.
    """
    if is_server == 0:
        LoadPlugin(
            "/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64"
            "/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",
            remote=False,
            ns=globals(),
        )
    else:
        LoadPlugin(
            "/home/u1435513-gma/apps/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64"
            "/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",
            remote=False,
            ns=globals(),
        )


# ---------------------------------------------------------------------
# Dataset geometry & domain (match binary_time_data_convert.py)
# dim = [Nx, Ny, Nz]
# ---------------------------------------------------------------------
def dataset_span_and_domain(name: str, up_sample_ratio: int):
    if name == "vortex_street_3d":
        base_dim = np.array([80, 10, 15], dtype=int)
        dom_min = np.array([-0.5, -0.5, 13.5], dtype=float)
        dom_max = np.array([ 7.5,  0.5, 15.0], dtype=float)
    elif name == "boussinesq_3d":
        base_dim = np.array([10, 30, 15], dtype=int)
        dom_min = np.array([-0.5, -0.5, 0.0], dtype=float)
        dom_max = np.array([ 0.5,  2.5, 1.5], dtype=float)
    elif name == "fluid":
        base_dim = np.array([10, 10, 10], dtype=int)
        dom_min = np.array([0.0, 0.0, 0.0], dtype=float)
        dom_max = np.array([1.0, 1.0, 1.0], dtype=float)
    elif name == "cylinder":
        base_dim = np.array([40, 10, 10], dtype=int)
        dom_min = np.array([1.5, 0.5, 0.0], dtype=float)
        dom_max = np.array([5.5, 1.5, 1.0], dtype=float)
    elif name == "cylinder2":
        base_dim = np.array([23, 10, 10], dtype=int)
        dom_min = np.array([3.2, 0.5, 0.0], dtype=float)
        dom_max = np.array([5.5, 1.5, 1.0], dtype=float)
    elif name == "cylinder3":
        base_dim = np.array([20, 10, 10], dtype=int)
        dom_min = np.array([3.5, 0.5, 0.0], dtype=float)
        dom_max = np.array([5.5, 1.5, 5.0], dtype=float)
    else:
        raise ValueError(f"Unknown dataset name: {name}")

    dim = base_dim * int(up_sample_ratio)  # [Nx, Ny, Nz]
    Nx, Ny, Nz = map(int, dim)
    return np.array([Nx, Ny, Nz], dtype=int), dom_min, dom_max


# ---------------------------------------------------------------------
# Bin → vtkImageData (same layout as binary_time_data_convert.py)
# ---------------------------------------------------------------------
def load_bin_to_vtkImageData(
    input_bin: str,
    global_dims: np.ndarray,
    dom_min: np.ndarray,
    dom_max: np.ndarray,
    float_type: str,
    z_start_idx: int,
    z_count: Optional[int],
) -> vtk.vtkImageData:
    print(f"[load_bin_to_vtkImageData] Reading bin: {input_bin}")
    dtype = np.float32 if float_type == "float32" else np.float64
    data = np.fromfile(input_bin, dtype=dtype)

    Nx_g, Ny_g, Nz_global = map(int, global_dims)
    if z_count is None:
        if data.size % (Nx_g * Ny_g) != 0:
            raise ValueError(
                f"Cannot infer z_count. data.size={data.size}, Nx*Ny={Nx_g * Ny_g}"
            )
        z_count = data.size // (Nx_g * Ny_g)

    expected = Nx_g * Ny_g * z_count
    if data.size != expected:
        raise ValueError(
            f"Bin size mismatch: data.size={data.size}, expected={expected} "
            f"= Nx({Nx_g}) * Ny({Ny_g}) * z_count({z_count})."
        )

    print(
        f"[load_bin_to_vtkImageData] Chunk dims: Nx={Nx_g}, Ny={Ny_g}, "
        f"z_start={z_start_idx}, z_count={z_count}"
    )

    # Reshape like binary_time_data_convert: (Nz, Ny, Nx) = (T, Ny, Nx)
    arr = data.reshape((z_count, Ny_g, Nx_g))

    distance = dom_max - dom_min
    dx = distance[0] / (Nx_g - 1) if Nx_g > 1 else 1.0
    dy = distance[1] / (Ny_g - 1) if Ny_g > 1 else 1.0
    dz = distance[2] / (Nz_global - 1) if Nz_global > 1 else 1.0

    img = vtk.vtkImageData()
    img.SetDimensions((Nx_g, Ny_g, 1))
    img.SetSpacing((dx, dy, 1.0))
    img.SetOrigin(dom_min[0], dom_min[1], dom_min[2])

    for local_z in range(z_count):
        global_z = z_start_idx + local_z
        arr_name = f"{global_z:04d}"

        # Exactly as in binary_time_data_convert.py:
        # slice_2d = np_array[z, :, :]  # (Ny, Nx)
        # flat = slice_2d.T.ravel(order='F')
        slice_2d = arr[local_z]              # (Ny, Nx)
        flat = slice_2d.T.ravel(order="F")   # (Nx*Ny, )

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(arr_name)
        img.GetPointData().AddArray(vtk_arr)

        # Encode physical z (time) in field data
        z_array = vtk.vtkDoubleArray()
        z_array.SetName(arr_name)
        z_array.SetNumberOfComponents(1)
        z_coord = dom_min[2] + dz * global_z
        z_array.InsertNextValue(z_coord)
        img.GetFieldData().AddArray(z_array)

    return img


# ---------------------------------------------------------------------
# TTK: critical points for one chunk
# ---------------------------------------------------------------------
def critical_points_from_vtkImageData_chunk(
    img: vtk.vtkImageData,
    csvwriter: csv.writer,
) -> None:
    print("[critical_chunk] Creating TrivialProducer from vtkImageData")
    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(img)

    precond = TTKArrayPreconditioning(Input=producer)
    precond.UpdatePipeline()
    tet = Tetrahedralize(Input=precond)

    point_data = img.GetPointData()
    arrays = [
        point_data.GetArrayName(i)
        for i in range(point_data.GetNumberOfArrays())
        if point_data.GetArrayName(i) is not None
    ]
    print("[critical_chunk] Scalar fields in this chunk:", arrays)

    for arr_name in arrays:
        print(f"[critical_chunk] Processing scalar field: {arr_name}")
        crit = TTKScalarFieldCriticalPoints(Input=tet)
        crit.ScalarField = ["POINTS", arr_name]
        crit.InputOffsetField = ["POINTS", arr_name]
        crit.UpdatePipeline()

        out = sm.Fetch(crit)

        pts = vtk_to_numpy(out.GetPoints().GetData())
        zvals = vtk_to_numpy(out.GetFieldData().GetArray(arr_name))
        if zvals.size > 0:
            pts[:, 2] = zvals[0]  # overwrite z with physical time

        bndry = vtk_to_numpy(out.GetPointData().GetArray("IsOnBoundary"))
        mask = bndry != 1
        pts = pts[mask]

        for p in pts:
            # csvwriter.writerow([p[0], p[1], p[2]])
            csvwriter.writerow([f"{p[0]:.16g}",f"{p[1]:.16}",f"{p[2]:.16g}"])


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Read INR-sampled bin and compute critical points with TTK."
    )
    ap.add_argument("--input_bin", type=str, required=True,
                    help="Input binary chunk file produced by INR sampler.")
    ap.add_argument("--dataset", type=str, required=True,
                    help="Dataset name (e.g., vortex_street_3d).")
    ap.add_argument("--up_sample_ratio", type=int, default=2,
                    help="Upsampling ratio used in sampling.")
    ap.add_argument("--z_start_idx", type=int, required=True,
                    help="Global time index of first slice in this bin.")
    ap.add_argument("--z_count", type=int, required=True,
                    help="Number of slices stored in this bin.")
    ap.add_argument("--float_type", type=str, default="float64",
                    choices=["float32", "float64"],
                    help="Bin float type.")
    ap.add_argument("--output_csv", type=str, required=True,
                    help="CSV to write/append critical points.")
    ap.add_argument("--append", action="store_true",
                    help="Append to CSV if it exists (no header).")
    ap.add_argument("--server", type=int, default=0,
                    help="0 for local ParaView, 1 for remote server build.")
    args = ap.parse_args()

    plugin_log(args.server)

    global_dims, dom_min, dom_max = dataset_span_and_domain(
        args.dataset, args.up_sample_ratio
    )
    print(
        f"[main] dataset={args.dataset}, global_dims={global_dims}, "
        f"dom_min={dom_min}, dom_max={dom_max}"
    )

    img = load_bin_to_vtkImageData(
        input_bin=args.input_bin,
        global_dims=global_dims,
        dom_min=dom_min,
        dom_max=dom_max,
        float_type=args.float_type,
        z_start_idx=args.z_start_idx,
        z_count=args.z_count,
    )

    mode = "a" if args.append and os.path.exists(args.output_csv) else "w"
    write_header = (mode == "w")

    out_dir = os.path.dirname(args.output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output_csv, mode, newline="") as f_csv:
        writer = csv.writer(f_csv)
        if write_header:
            writer.writerow(["PositionX", "PositionY", "PositionZ"])
        critical_points_from_vtkImageData_chunk(img, writer)

    print(
        f"[main] Done. Critical points written to {args.output_csv} "
        f"(append={args.append})."
    )


if __name__ == "__main__":
    main()