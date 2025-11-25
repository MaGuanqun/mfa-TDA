#!/usr/bin/env pvpython
"""
Chunked MFA sampling → TTK critical points (NO VTI WRITTEN).

Workflow (chunked along z):
    For each z-chunk of size `chunk_size` (default 50):
        1. Call C++ sampler:
               write_vtk <bin_path> <z_start_idx> <z_end_idx> <step_size>
        2. Load that bin (only those slices) into vtkImageData in memory.
        3. Run TTK ScalarFieldCriticalPoints for each slice in that chunk.
        4. Append critical points (non-boundary) to a single CSV.

At the end:
    - You get ONE CSV with all critical points over the full volume.
    - Only one bin file exists on disk at any moment (reused & removed).
"""

import os
import argparse
import csv
import subprocess
import tempfile
from typing import Optional

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from paraview.simple import *  # noqa: F401,F403
from paraview import servermanager as sm

import gc  # <<< for explicit garbage collection


# ---------------------------------------------------------------------------
# TTK plugin loading
# ---------------------------------------------------------------------------

def plugin_log(is_server: int) -> None:
    """Load TTK plugin either on local machine (0) or remote server (1)."""
    if is_server == 0:
        # local workstation
        LoadPlugin(
            "/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64"
            "/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",
            remote=False,
            ns=globals(),
        )
    else:
        # SCI cluster
        LoadPlugin(
            "/home/u1435513-gma/apps/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64"
            "/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",
            remote=False,
            ns=globals(),
        )


# ---------------------------------------------------------------------------
# C++ MFA sampler wrapper
# ---------------------------------------------------------------------------

def run_cpp_sampler(
    exe_path: str,
    out_bin: str,
    z_start: int,
    z_end: int,
    step_size: int,
    mfa_file: str = "mfa_data.mfa",
) -> None:
    """
    Run the C++ MFA sampler executable and produce a bin file for a chunk of slices.

    Assumed C++ signature:
        write_vtk <out_bin> <z_start_idx> <z_end_idx> <step_size>
    """
    cmd = [
        exe_path,
        "-t", out_bin,
        "-m", str(3),
        "-d", str(4),
        "-g", str(0),
        "-z", str(0),
        "-b", str(1),
        "-u", str(step_size)+"-"+str(step_size)+"-"+str(step_size),
        "-x", str(z_start),
        "-y", str(z_end),
        "-f", mfa_file
    ]
    print("[run_cpp_sampler] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"[run_cpp_sampler] Finished, bin written to: {out_bin}")


# ---------------------------------------------------------------------------
# Bin → vtkImageData (chunk, in memory)
# ---------------------------------------------------------------------------

def load_bin_to_vtkImageData(
    input_bin: str,
    global_dims: np.ndarray,   # [Nx_global, Ny_global, Nz_global]
    dom_min: np.ndarray,
    dom_max: np.ndarray,
    float_type: str,
    z_start_idx: int,
    z_count: Optional[int],
) -> vtk.vtkImageData:
    """
    Convert a binary file of a *chunk* of slices into vtkImageData in memory.

    Assumptions:
      - The chunk contains `z_count` slices: z_start_idx .. z_start_idx + z_count - 1.
      - Each slice is Ny_global x Nx_global.
      - Data layout in the bin is (z_count, Ny, Nx) in C-order.
      - global_dims gives the *full* volume dims, so we can compute dz correctly.
    """
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

    # reshape to (z_count, Ny, Nx)
    arr = data.reshape((z_count, Ny_g, Nx_g))

    # full-domain spacings (use Nz_global here!)
    distance = dom_max - dom_min
    dx = distance[0] / (Nx_g - 1) if Nx_g > 1 else 1.0
    dy = distance[1] / (Ny_g - 1) if Ny_g > 1 else 1.0
    dz = distance[2] / (Nz_global - 1) if Nz_global > 1 else 1.0

    # image is XY plane; each slice becomes a scalar array
    img = vtk.vtkImageData()
    img.SetDimensions((Nx_g, Ny_g, 1))
    img.SetSpacing((dx, dy, 1.0))
    img.SetOrigin(dom_min)

    for local_z in range(z_count):
        global_z = z_start_idx + local_z
        arr_name = f"{global_z:04d}"

        # Flatten XY with x fastest: transpose to (x,y) then Fortran-flatten
        flat = arr[local_z].T.ravel(order="F")

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(arr_name)
        img.GetPointData().AddArray(vtk_arr)

        # store physical z coordinate in field data
        z_array = vtk.vtkDoubleArray()
        z_array.SetName(arr_name)
        z_array.SetNumberOfComponents(1)
        z_coord = dom_min[2] + dz * global_z
        z_array.InsertNextValue(z_coord)
        img.GetFieldData().AddArray(z_array)

        # free per-slice flat array
        del flat, vtk_arr, z_array

    # ### MEMORY CLEANUP (numpy side)
    del arr, data
    gc.collect()

    return img


# ---------------------------------------------------------------------------
# TTK critical point extraction for a chunk (append to CSV)
# ---------------------------------------------------------------------------

def critical_points_from_vtkImageData_chunk(
    img: vtk.vtkImageData,
    csvwriter: csv.writer,
) -> None:
    """
    Run TTK on a vtkImageData chunk and append critical points to CSV via csvwriter.
    """

    # 1) Wrap raw vtkImageData in a ParaView pipeline source
    print("[critical_chunk] Creating TrivialProducer from vtkImageData")
    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(img)

    # 2) Preconditioning + tetrahedralization
    precond = TTKArrayPreconditioning(Input=producer)
    precond.UpdatePipeline()
    tet = Tetrahedralize(Input=precond)

    # 3) Get scalar array names from the image's point data
    point_data = img.GetPointData()
    arrays = [
        point_data.GetArrayName(i)
        for i in range(point_data.GetNumberOfArrays())
        if point_data.GetArrayName(i) is not None
    ]
    print("[critical_chunk] Scalar fields in this chunk:", arrays)

    # 4) Run TTK critical points per scalar field
    for arr_name in arrays:
        print(f"[critical_chunk] Processing scalar field: {arr_name}")
        crit = TTKScalarFieldCriticalPoints(Input=tet)
        crit.ScalarField = ["POINTS", arr_name]
        crit.InputOffsetField = ["POINTS", arr_name]
        crit.UpdatePipeline()

        out = sm.Fetch(crit)

        # Positions
        pts = vtk_to_numpy(out.GetPoints().GetData())

        # z/t coordinate stored in field data with same array name
        z_array = out.GetFieldData().GetArray(arr_name)
        if z_array is not None:
            zvals = vtk_to_numpy(z_array)
            if zvals.size > 0:
                pts[:, 2] = zvals[0]
        else:
            zvals = np.array([])

        # Remove boundary critical points
        bndry_arr = out.GetPointData().GetArray("IsOnBoundary")
        if bndry_arr is not None:
            bndry = vtk_to_numpy(bndry_arr)
            mask = bndry != 1
            pts = pts[mask]

        # Append to CSV
        for p in pts:
            csvwriter.writerow([p[0], p[1], p[2]])

        # ### MEMORY CLEANUP for this scalar field
        del pts, zvals
        if bndry_arr is not None:
            del bndry, bndry_arr
        del out
        Delete(crit)
        del crit
        gc.collect()

    # ### MEMORY CLEANUP for this chunk's pipeline
    Delete(tet)
    Delete(precond)
    Delete(producer)
    del tet, precond, producer, point_data, arrays
    gc.collect()


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Chunked MFA bin → TTK critical points (no VTI; 50-slice C++ chunks)."
        )
    )

    parser.add_argument(
        "--cpp_exe",
        type=str,
        required=True,
        help="Path to the C++ MFA sampler executable (e.g., write_vtk).",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=50,
        help="Number of z-slices per chunk (default: 50).",
    )
    parser.add_argument(
        "--step_size",
        type=int,
        default=1,
        help="Sampling step size used in both C++ and Python.",
    )
    parser.add_argument(
        "--float_type",
        type=str,
        default="float32",
        choices=["float32", "float64"],
        help="Binary float type written by the C++ sampler.",
    )
    parser.add_argument(
        "--function",
        type=str,
        default="vortex_street_3d",
        help="Function / dataset name (controls base dims & domain extents).",
    )
    parser.add_argument(
        "--server",
        type=int,
        default=0,
        help="0 for local machine, 1 for server (controls TTK plugin path).",
    )

    parser.add_argument(
        "-o",
        "--output_csv",
        type=str,
        default="critical_points_all.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--mfa_file",
        type=str,
        default="mfa_data.mfa",
        help="Input MFA file for the C++ sampler.",
    )

    args = parser.parse_args()

    # Load TTK plugin
    plugin_log(args.server)

    # Domain & base dims
    if args.function == "rotating_gaussian":
        dom_min = np.array([-2.0, -2.0, 0.0])
        dom_max = np.array([2.0, 2.0, 4.0])
        dim = np.array([100, 100, 100])
    elif args.function in ("quartic_potential", "quartic_potential_2"):
        dom_min = np.array([-2.0, -2.0, 0.0])
        dom_max = np.array([2.0, 2.0, 4.0])
        dim = np.array([100, 100, 100])
    elif args.function == "vortex_street":
        dim = np.array([100, 80, 50])
        dom_min = np.array([0.0, 0.0, 0.0])
        dom_max = np.array([99.0, 79.0, 49.0])
    elif args.function =='vortex_street_3d':
        dim = np.array([80,10,18])
        dom_min = np.array([-0.5, -0.5, 13.5])
        dom_max = np.array([7.5, 0.5, 15.0])
    elif args.function =='boussinesq_3d':
        dim = np.array([7,27,7])
        dom_min = np.array([-0.5, -0.5, 0.0])
        dom_max = np.array([0.5, 2.5, 1.5])
    else:
        raise ValueError(f"Unknown function name: {args.function}")

    # Global dims after step_size
    global_dims = args.step_size * dim
    Nx_g, Ny_g, Nz_g = map(int, global_dims)
    print(
        f"[main] function={args.function}, base dim={dim}, "
        f"global_dims={global_dims} (Nx, Ny, Nz)"
    )

    # Temp bin file reused for all chunks
    bin_path = os.path.join(tempfile.gettempdir(), "mfa_tmp_chunk.bin")
    print(f"[main] Using temp bin file: {bin_path}")

    # Open output CSV once, append per-chunk
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["PositionX", "PositionY", "PositionZ"])

        z = 0
        while z < Nz_g:
            z_start = z
            z_end = min(z_start + args.chunk_size, Nz_g)
            z_count = z_end - z_start

            print(
                f"[main] Processing chunk: z_start={z_start}, z_end={z_end}, "
                f"z_count={z_count}"
            )

            # 1) C++ sampler writes only this chunk
            run_cpp_sampler(
                exe_path=args.cpp_exe,
                out_bin=bin_path,
                z_start=z_start,
                z_end=z_end,
                step_size=args.step_size,
                mfa_file=args.mfa_file
            )

            # 2) Bin → vtkImageData for this chunk
            img = load_bin_to_vtkImageData(
                input_bin=bin_path,
                global_dims=global_dims,
                dom_min=dom_min,
                dom_max=dom_max,
                float_type=args.float_type,
                z_start_idx=z_start,
                z_count=z_count,
            )

            # 3) TTK → append critical points
            critical_points_from_vtkImageData_chunk(img, writer)

            # ### MEMORY CLEANUP: drop img and force GC per chunk
            del img
            gc.collect()

            z = z_end

    # Clean up temp bin
    try:
        os.remove(bin_path)
        print("[main] Removed temp bin:", bin_path)
    except OSError as e:
        print("[main] Warning: could not remove temp bin:", e)

    print(f"[main] Done. All critical points saved to: {args.output_csv}")


if __name__ == "__main__":
    main()
