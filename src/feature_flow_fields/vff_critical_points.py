#!/usr/bin/env python
"""
Compute critical points of every 2-D slice of a 3-D vector field stored in a
``.vff`` file (produced by ``src/convert/write_gradient.cpp::save_vff``), using
TTK's TopologicalSkeleton (discrete vector field topology), and collect all
slice critical points into a single CSV.

Pipeline (per slice along the last/"z" axis):
  1. Extract the 2-D (x, y) slice of the vector field. We keep only the in-plane
     gradient components (vx, vy) and zero the third one, so TTK sees a planar
     vector field whose zeros are the spatial critical points at that z.
  2. Build a 2-D vtkImageData placed at z = 0 in memory (fed to TTK directly via
     a TrivialProducer, so no per-slice .vti is written/read).
  3. Run TTKTopologicalSkeleton on the vector field; output port 0 is the set of
     critical points (see the TTK example:
     https://topology-tool-kit.github.io/examples/discreteVectorFieldTopology/).
  4. The critical points come out with z = 0; set z to the real coordinate of
     this slice so the points are lifted back into 3-D.

All slices' critical points are concatenated and written to one CSV.

.vff layout (little-endian):
  "VFF1", uint32 dtype(0=f64), uint32 D, uint32 C, uint32 n[D],
  float64 dmin[D], float64 dmax[D], then C*prod(n) float64 values,
  AoS (components innermost): offset = (i0 + n0*(i1 + n1*(...)))*C + c.

Run with pvpython (ParaView's python, which provides paraview.simple + TTK):
  pvpython src/feature_flow_fields/vff_critical_points.py \
      -i tracking_result_16.vff \
      -o tracking_result_16_cpt.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import struct

import numpy as np
import vtk
import vtk.util.numpy_support as VN
from paraview import servermanager as sm
from paraview.simple import (
    Delete,
    OutputPort,
    TrivialProducer,
    TTKTopologicalSkeleton,
)


# def plugin_log(is_server: int) -> None:
#     """Load the TTK ParaView plugin (mirrors extract_all_critical_points.py)."""
#     from paraview.simple import LoadPlugin

#     if is_server == 0:
#         plugin = (
#             "/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/"
#             "lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so"
#         )
#     else:
#         plugin = (
#             "/home/u1435513-gma/apps/"
#             "ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64/"
#             "lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so"
#         )
#     try:
#         LoadPlugin(plugin, remote=False, ns=globals())
#     except Exception as exc:  # noqa: BLE001 - plugin may already be available
#         print(f"Warning: could not load TTK plugin ({plugin}): {exc}")


def read_vff(path: str):
    """Read a .vff file. Returns (data, n, dmin, dmax) with data shaped
    (..., n[1], n[0], C) in C-order (x fastest, components innermost)."""
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic != b"VFF1":
            raise ValueError(f"{path}: not a .vff file (bad magic {magic!r})")
        dtype, D, C = struct.unpack("<III", f.read(12))
        if dtype != 0:
            raise ValueError(f"{path}: unsupported dtype {dtype} (only 0=float64)")
        n = list(struct.unpack("<" + "I" * D, f.read(4 * D)))
        dmin = list(struct.unpack("<" + "d" * D, f.read(8 * D)))
        dmax = list(struct.unpack("<" + "d" * D, f.read(8 * D)))
        count = C * int(np.prod(n))
        data = np.fromfile(f, dtype="<f8", count=count)
    if data.size != count:
        raise ValueError(
            f"{path}: expected {count} float64 values, read {data.size}"
        )
    # AoS, x fastest, components innermost -> reshape to (n[D-1], ..., n[1], n[0], C)
    shape = list(reversed(n)) + [C]
    data = data.reshape(shape)
    return data, n, dmin, dmax


def axis_spacing(dmin: float, dmax: float, n: int) -> float:
    return (dmax - dmin) / (n - 1) if n > 1 else 0.0


def build_slice_image(slice_data: np.ndarray, nx: int, ny: int,
                      xmin: float, ymin: float, dx: float, dy: float) -> vtk.vtkImageData:
    """Build a 2-D vtkImageData (z = 0) with a 3-component "VectorField" point
    array holding (vx, vy, 0). slice_data is shaped (ny, nx, C)."""
    image = vtk.vtkImageData()
    image.SetDimensions(nx, ny, 1)
    image.SetOrigin(xmin, ymin, 0.0)
    image.SetSpacing(dx if dx != 0.0 else 1.0, dy if dy != 0.0 else 1.0, 1.0)

    # (ny, nx, C) -> (ny*nx, C) keeps x fastest, matching VTK image point order.
    flat = slice_data.reshape(ny * nx, slice_data.shape[-1])
    vecs = np.zeros((ny * nx, 3), dtype=np.float64)
    vecs[:, 0] = flat[:, 0]
    vecs[:, 1] = flat[:, 1]
    # third (out-of-plane / temporal) component zeroed: planar field for TTK.

    arr = VN.numpy_to_vtk(vecs, deep=True)
    arr.SetName("VectorField")
    image.GetPointData().AddArray(arr)
    image.GetPointData().SetActiveVectors("VectorField")
    return image


def critical_points_for_image(image: vtk.vtkImageData, simplification_threshold: float):
    """Run TTKTopologicalSkeleton on an in-memory vtkImageData (fed via a
    TrivialProducer, no disk I/O) and return (positions Nx3, types N) for the
    critical points (output port 0)."""
    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(image)
    producer.UpdatePipeline()

    skeleton = TTKTopologicalSkeleton(Input=producer)
    skeleton.VectorField = ["POINTS", "VectorField"]
    if simplification_threshold > 0.0:
        skeleton.RunSimplification = 1
        skeleton.SimplificationThreshold = simplification_threshold
    else:
        skeleton.RunSimplification = 0
    skeleton.UpdatePipeline()

    cp = sm.Fetch(OutputPort(skeleton, 0))
    if cp is None or cp.GetNumberOfPoints() == 0:
        positions, types = np.empty((0, 3)), np.empty((0,))
    else:
        positions = VN.vtk_to_numpy(cp.GetPoints().GetData()).astype(np.float64)
        type_arr = cp.GetPointData().GetArray("CriticalType")
        if type_arr is not None:
            types = VN.vtk_to_numpy(type_arr).astype(np.float64).ravel()
        else:
            types = np.full(positions.shape[0], np.nan)

        keep = ~boundary_mask(cp, positions, image)
        positions, types = positions[keep], types[keep]

    # Drop the per-slice pipeline objects so they don't accumulate over slices.
    Delete(skeleton)
    Delete(producer)
    return positions, types


def boundary_mask(cp, positions: np.ndarray, image: vtk.vtkImageData,
                  tol: float = 1e-8) -> np.ndarray:
    """True for critical points on the slice boundary. Prefer TTK's
    "IsOnBoundary" array; fall back to the (x, y) extent of the slice."""
    bnd = cp.GetPointData().GetArray("IsOnBoundary")
    if bnd is not None:
        return VN.vtk_to_numpy(bnd).ravel() == 1

    xmin, xmax, ymin, ymax, _, _ = image.GetBounds()
    return (
        (np.abs(positions[:, 0] - xmin) <= tol)
        | (np.abs(positions[:, 0] - xmax) <= tol)
        | (np.abs(positions[:, 1] - ymin) <= tol)
        | (np.abs(positions[:, 1] - ymax) <= tol)
    )


def compute_all_slices(vff_path: str, output_csv: str, t_stride: int,
                       simplification_threshold: float) -> None:
    data, n, dmin, dmax = read_vff(vff_path)
    D = len(n)
    if D != 3:
        raise ValueError(
            f"{vff_path}: expected a 3-D vector field (D=3), got D={D}. "
            "Slicing is defined along the last (z) axis of a 3-D field."
        )
    C = data.shape[-1]
    if C < 2:
        raise ValueError(f"{vff_path}: need >=2 vector components, got C={C}")

    nx, ny, nz = n[0], n[1], n[2]
    dx = axis_spacing(dmin[0], dmax[0], nx)
    dy = axis_spacing(dmin[1], dmax[1], ny)
    dz = axis_spacing(dmin[2], dmax[2], nz)

    print(f"Loaded {vff_path}: grid n={n}, C={C}")
    print(f"  dmin={dmin}  dmax={dmax}")
    print(f"  spacing dx={dx} dy={dy} dz={dz}  slicing along z ({nz} slices)")

    all_pos = np.empty((0, 3))
    all_types = np.empty((0,))
    all_slice = np.empty((0,))

    slice_indices = range(0, nz, max(1, t_stride))
    for k in slice_indices:
        slice_z = dmin[2] + k * dz
        slice_data = data[k]  # (ny, nx, C)

        image = build_slice_image(slice_data, nx, ny, dmin[0], dmin[1], dx, dy)
        positions, types = critical_points_for_image(image, simplification_threshold)
        if positions.shape[0] > 0:
            positions[:, 2] = slice_z  # lift critical points back to 3-D
            all_pos = np.concatenate((all_pos, positions), axis=0)
            all_types = np.concatenate((all_types, types), axis=0)
            all_slice = np.concatenate(
                (all_slice, np.full(positions.shape[0], k)), axis=0
            )

        print(f"  slice {k:5d} (z={slice_z:.6g}): {positions.shape[0]} critical points")

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PositionX", "PositionY", "PositionZ"])#, "CriticalType", "SliceIndex"
        for i in range(all_pos.shape[0]):
            ctype = all_types[i]
            writer.writerow([
                all_pos[i, 0],
                all_pos[i, 1],
                all_pos[i, 2]
            ]) #                "" if np.isnan(ctype) else int(ctype),
                #int(all_slice[i]),

    print(f"Wrote {all_pos.shape[0]} critical points from "
          f"{len(list(slice_indices))} slices to {output_csv}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Per-slice critical points of a 3-D .vff vector field via TTK."
    )
    p.add_argument("-i", "--input_name", required=True, help="input .vff file")
    p.add_argument("-o", "--output_name", required=True, help="output .csv file")
    p.add_argument("--t-stride", type=int, default=1,
                   help="process every Nth slice along z (default: 1 = all).")
    p.add_argument("--simplification-threshold", type=float, default=0.0,
                   help="TTK VectorSimplification threshold; 0 disables (default).")
    # p.add_argument("-s", "--server", type=int, default=0,
    #                help="0 = local PC, 1 = remote server (for TTK plugin path).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not os.path.exists(args.input_name):
        raise FileNotFoundError(f"Input .vff not found: {args.input_name}")
    # plugin_log(is_server=args.server)
    compute_all_slices(
        vff_path=args.input_name,
        output_csv=args.output_name,
        t_stride=int(args.t_stride),
        simplification_threshold=float(args.simplification_threshold),
    )


if __name__ == "__main__":
    main()
