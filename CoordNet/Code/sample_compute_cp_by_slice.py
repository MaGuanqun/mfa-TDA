#!/usr/bin/env pvpython
"""
inr_bin_critical_point.py

Chunked critical point extraction directly from a float64 TorchScript
INR model (CoordNet, application='super-spatial-temporal').

- No big .dat file, no big full 3D volume.
- We sample only a few t-slices at a time, push them into vtkImageData,
  and run TTK on each chunk, appending critical points to one CSV.

Assumptions:
  * Dataset is one of: vortex_street_3d, boussinesq_3d, fluid, cylinder, cylinder2, cylinder3.
  * Model is a TorchScript file taking 3D coords (t,y,x) in [-1,1]^3.
  * You run with pvpython so ParaView + TTK are available.
"""

import os
import csv
import argparse

import numpy as np
import torch
from torch.utils.data import DataLoader

import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from paraview.simple import *          # noqa: F401,F403
from paraview import servermanager as sm


# ---------------------------------------------------------------------------
# TTK plugin loading (same idea as in mfa_bin_critical_point.py)
# ---------------------------------------------------------------------------

def plugin_log(is_server: int):
    """Load TTK plugin either on local machine (0) or remote server (1)."""
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


# ---------------------------------------------------------------------------
# Dataset geometry & sampling config
# ---------------------------------------------------------------------------

def dataset_span_and_domain(name: str, up_sample_ratio: int):
    """
    Mimic ScalarDataSet.span_num() and domain ranges you used
    in the MFA pipeline, then apply up_sample_ratio.
    """
    if name == "vortex_street_3d":
        base_span = np.array([80, 10, 15])      # [Ny_blocks, Nx_blocks, Nt_blocks]
        dom_min = np.array([-0.5, -0.5, 13.5])  # [y_min, x_min, t_min]
        dom_max = np.array([ 7.5,  0.5, 15.0])  # [y_max, x_max, t_max]
    elif name == "boussinesq_3d":
        base_span = np.array([10, 30, 15])
        dom_min = np.array([-0.5, -0.5, 0.0])
        dom_max = np.array([ 0.5,  2.5, 1.5])
    elif name == "fluid":
        base_span = np.array([10, 10, 10])
        dom_min = np.array([0.0, 0.0, 0.0])
        dom_max = np.array([1.0, 1.0, 1.0])
    elif name == "cylinder":
        base_span = np.array([40, 10, 10])
        dom_min = np.array([1.5, 0.5, 0.0])
        dom_max = np.array([5.5, 1.5, 1.0])
    elif name == "cylinder2":
        base_span = np.array([23, 10, 10])
        dom_min = np.array([3.2, 0.5, 0.0])
        dom_max = np.array([5.5, 1.5, 1.0])
    elif name == "cylinder3":
        base_span = np.array([20, 10, 10])
        dom_min = np.array([3.5, 0.5, 0.0])
        dom_max = np.array([5.5, 1.5, 5.0])
    else:
        raise ValueError(f"Unsupported dataset: {name}")

    # Up-sampled sampling grid (super-spatial-temporal)
    sample_size = up_sample_ratio * base_span
    # Interpret sample_size as [H, W, T] (rows, cols, time)
    H, W, T = map(int, sample_size)
    return H, W, T, dom_min, dom_max


# ---------------------------------------------------------------------------
# Build coords for a chunk of time steps
# ---------------------------------------------------------------------------

def build_coords_chunk(H, W, T, z_start, z_end, device, dtype=torch.float64):
    """
    Build normalized (t,y,x) coords in [-1,1]^3 for time indices [z_start, z_end).
    Loops in a simple order: for t, for y, for x. That order is *only* an indexing
    choice; the model sees explicit coords so the order doesn't matter.
    """
    z_count = z_end - z_start
    total = z_count * H * W

    # Preallocate
    coords = np.empty((total, 3), dtype=np.float64)

    # Normalization factors (match get_mgrid scaling)
    inv_T = 1.0 / max(T - 1, 1)
    inv_H = 1.0 / max(H - 1, 1)
    inv_W = 1.0 / max(W - 1, 1)

    idx = 0
    for t in range(z_start, z_end):
        t_norm = (t * inv_T - 0.5) * 2.0
        for y in range(H):
            y_norm = (y * inv_H - 0.5) * 2.0
            for x in range(W):
                x_norm = (x * inv_W - 0.5) * 2.0
                coords[idx, 0] = t_norm
                coords[idx, 1] = y_norm
                coords[idx, 2] = x_norm
                idx += 1

    coords_tensor = torch.from_numpy(coords).to(device=device, dtype=dtype)
    return coords_tensor  # shape [z_count * H * W, 3]


# ---------------------------------------------------------------------------
# Evaluate INR model on a chunk
# ---------------------------------------------------------------------------

def eval_inr_chunk(model, H, W, T, z_start, z_end, batch_size, device, preferred_dtype):
    """
    Evaluate model on time steps [z_start, z_end) and return values_chunk
    of shape [z_count, H, W] in float64.
    """
    coords_tensor = build_coords_chunk(H, W, T, z_start, z_end, device, dtype=preferred_dtype)
    z_count = z_end - z_start

    loader = DataLoader(
        dataset=coords_tensor,
        batch_size=batch_size,
        shuffle=False,
        pin_memory=(device.type == "cuda"),
    )

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device, non_blocking=(device.type == "cuda")).to(preferred_dtype)
            out = model(batch)  # [B,1] or [B]
            preds.append(out.view(-1).detach().cpu().to(torch.float64).numpy())
        if device.type == "cuda":
            torch.cuda.synchronize()

    preds = np.concatenate(preds, axis=0)  # length = z_count * H * W

    # Rebuild [z_count, H, W]
    values_chunk = np.empty((z_count, H, W), dtype=np.float64)
    per_slice = H * W
    for i in range(z_count):
        v = preds[i * per_slice:(i + 1) * per_slice]
        values_chunk[i] = v.reshape(H, W)   # rows=y, cols=x

    return values_chunk  # [z_count, H, W]


# ---------------------------------------------------------------------------
# Build vtkImageData from INR chunk (no .vti on disk)
# ---------------------------------------------------------------------------

def build_vtk_image_from_chunk(values_chunk, H, W, T,
                               dom_min, dom_max, z_start_idx):
    """
    values_chunk: [z_count, H, W] numpy float64
    H, W, T: global dimensions (rows, cols, total time)
    dom_min, dom_max: physical domain [y_min, x_min, t_min] / [y_max, x_max, t_max]
    z_start_idx: global time index of first slice in this chunk
    Returns vtk.vtkImageData with:
      - dimensions (W, H, 1)  (x,y,1)
      - one scalar array per slice, named by global z index ("%04d")
      - field data arrays storing physical t for each slice
    """
    z_count = values_chunk.shape[0]

    # Spacing using full global dims
    y_dist = dom_max[0] - dom_min[0]
    x_dist = dom_max[1] - dom_min[1]
    t_dist = dom_max[2] - dom_min[2]

    dy = y_dist / (H - 1) if H > 1 else 1.0
    dx = x_dist / (W - 1) if W > 1 else 1.0
    dz = t_dist / (T - 1) if T > 1 else 1.0

    image_data = vtk.vtkImageData()
    image_data.SetDimensions((W, H, 1))               # (Nx, Ny, 1)
    image_data.SetSpacing((dx, dy, 1.0))
    # Origin at (x_min, y_min, t_min); z handled via field data
    image_data.SetOrigin(dom_min[1], dom_min[0], dom_min[2])

    for local_z in range(z_count):
        slice_2d = values_chunk[local_z]             # [H, W]
        # Flatten so x is fastest; Fortran-like over transpose
        flat = slice_2d.T.ravel(order="F")

        global_z_idx = z_start_idx + local_z
        arr_name = f"{global_z_idx:04d}"

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(arr_name)
        image_data.GetPointData().AddArray(vtk_arr)

        # Store physical t coordinate in field data under same name
        t_coord = dom_min[2] + dz * global_z_idx
        z_array = vtk.vtkDoubleArray()
        z_array.SetName(arr_name)
        z_array.SetNumberOfComponents(1)
        z_array.InsertNextValue(t_coord)
        image_data.GetFieldData().AddArray(z_array)

    return image_data


# ---------------------------------------------------------------------------
# Critical point extraction on a vtkImageData chunk
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

        pts = vtk_to_numpy(out.GetPoints().GetData())
        zvals = vtk_to_numpy(out.GetFieldData().GetArray(arr_name))
        if zvals.size > 0:
            pts[:, 2] = zvals[0]   # overwrite z with physical t

        bndry = vtk_to_numpy(out.GetPointData().GetArray("IsOnBoundary"))
        mask = bndry != 1
        pts = pts[mask]

        for p in pts:
            csvwriter.writerow([p[0], p[1], p[2]])


# ---------------------------------------------------------------------------
# TorchScript loading
# ---------------------------------------------------------------------------

def load_ts_model(ts_path: str, device: torch.device, in_dim: int = 3):
    print(f"[inr] Loading TorchScript model: {ts_path}")
    model = torch.jit.load(ts_path, map_location=device)
    model.eval()
    try:
        model = model.to(device)
    except Exception:
        pass

    # Check if it accepts float64; if not, fall back to float32
    preferred_dtype = None
    with torch.no_grad():
        try:
            probe64 = torch.zeros(1, in_dim, dtype=torch.float64, device=device)
            _ = model(probe64)
            preferred_dtype = torch.float64
            print("[inr] model accepted float64 inputs")
        except Exception as e64:
            print(f"[inr] float64 probe failed: {e64}")
            probe32 = torch.zeros(1, in_dim, dtype=torch.float32, device=device)
            _ = model(probe32)
            preferred_dtype = torch.float32
            print("[inr] model accepted float32 inputs")

    return model, preferred_dtype


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Chunked critical point extraction from INR TorchScript model."
    )

    ap.add_argument("--ts_model", type=str, required=True,
                    help="Path to float64 TorchScript model (.pt).")
    ap.add_argument("--dataset", type=str, required=True,
                    help="Dataset name (vortex_street_3d, boussinesq_3d, etc.).")
    ap.add_argument("--output_csv", type=str, required=True,
                    help="Output CSV path for critical points.")
    ap.add_argument("--chunk_size", type=int, default=5,
                    help="Number of time slices per chunk.")
    ap.add_argument("--batch_size", type=int, default=100000,
                    help="Batch size for INR evaluation.")
    ap.add_argument("--up_sample_ratio", type=int, default=2,
                    help="Upsampling ratio (same as in ScalarDataSet.GetTestingData).")
    ap.add_argument("--server", type=int, default=0,
                    help="0 for local PC, 1 for remote ParaView build.")
    ap.add_argument("--no_cuda", action="store_true",
                    help="Force CPU even if CUDA is available.")

    args = ap.parse_args()

    use_cuda = (not args.no_cuda) and torch.cuda.is_available()
    device = torch.device("cuda") if use_cuda else torch.device("cpu")
    print(f"[main] device = {device}")

    # 1) Load TTK plugin
    plugin_log(is_server=args.server)

    # 2) Dataset geometry
    H, W, T, dom_min, dom_max = dataset_span_and_domain(
        args.dataset, args.up_sample_ratio
    )
    print(f"[main] dataset={args.dataset}, H={H}, W={W}, T={T}")
    print(f"[main] dom_min={dom_min}, dom_max={dom_max}")

    # 3) Load TorchScript model
    model, preferred_dtype = load_ts_model(args.ts_model, device, in_dim=3)

    # 4) Prepare CSV
    out_dir = os.path.dirname(args.output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output_csv, "w", newline="") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["PositionX", "PositionY", "PositionZ"])

        # 5) Process time in chunks
        z = 0
        while z < T:
            z_start = z
            z_end = min(z_start + args.chunk_size, T)
            z_count = z_end - z_start
            print(
                f"[main] Processing chunk: z_start={z_start}, "
                f"z_end={z_end}, z_count={z_count}"
            )

            # 5a) Evaluate INR for this chunk
            values_chunk = eval_inr_chunk(
                model=model,
                H=H,
                W=W,
                T=T,
                z_start=z_start,
                z_end=z_end,
                batch_size=args.batch_size,
                device=device,
                preferred_dtype=preferred_dtype,
            )

            # 5b) Wrap into vtkImageData
            img = build_vtk_image_from_chunk(
                values_chunk=values_chunk,
                H=H,
                W=W,
                T=T,
                dom_min=dom_min,
                dom_max=dom_max,
                z_start_idx=z_start,
            )

            # 5c) Run TTK critical points and append to CSV
            critical_points_from_vtkImageData_chunk(img, writer)

            # Next chunk
            z = z_end

    print(f"[main] Done. Critical points saved to: {args.output_csv}")


if __name__ == "__main__":
    main()
