#!/usr/bin/env python
import os
import argparse
import subprocess
import numpy as np
import torch
from torch.utils.data import DataLoader


# ---------------------------------------------------------------------
# Dataset geometry: span_num + domain (H, W, T)
# These H, W, T match ScalarDataSet.span_num() in dataio.py:
#   vortex_street_3d: [80, 10, 15]
#   boussinesq_3d:   [10, 30, 15]
#   fluid:           [10, 10, 10]
#   cylinder*:       [...]
# ---------------------------------------------------------------------
def dataset_span_and_domain(name: str, up_sample_ratio: int):
    if name == "vortex_street_3d":
        base_span = np.array([80, 10, 15])   # [H, W, T]
        dom_min  = np.array([-0.5, -0.5, 13.5])
        dom_max  = np.array([ 7.5,  0.5, 15.0])
    elif name == "boussinesq_3d":
        base_span = np.array([10, 30, 15])
        dom_min  = np.array([-0.5, -0.5, 0.0])
        dom_max  = np.array([ 0.5,  2.5, 1.5])
    elif name == "fluid":
        base_span = np.array([10, 10, 10])
        dom_min  = np.array([0.0, 0.0, 0.0])
        dom_max  = np.array([1.0, 1.0, 1.0])
    elif name == "cylinder":
        base_span = np.array([40, 10, 10])
        dom_min  = np.array([1.5, 0.5, 0.0])
        dom_max  = np.array([5.5, 1.5, 1.0])
    elif name == "cylinder2":
        base_span = np.array([23, 10, 10])
        dom_min  = np.array([3.2, 0.5, 0.0])
        dom_max  = np.array([5.5, 1.5, 1.0])
    elif name == "cylinder3":
        base_span = np.array([20, 10, 10])
        dom_min  = np.array([3.5, 0.5, 0.0])
        dom_max  = np.array([5.5, 1.5, 5.0])
    else:
        raise ValueError(f"Unknown dataset {name}")

    # upsample like in GetTestingData: span_num() * up_sample_ratio
    H0, W0, T0 = base_span
    H, W, T = (base_span * int(up_sample_ratio)).astype(int)
    return H, W, T, dom_min, dom_max


# ---------------------------------------------------------------------
# Build coords for a chunk [z_start, z_end), matching get_mgrid EXACTLY
#
# Original GetTestingData() does:
#   sample_size = span_num() * up_sample_ratio = [H,W,T]
#   get_mgrid([T, H, W], dim=3)
#
# get_mgrid([T,H,W]) internally uses np.mgrid[:T, :W, :H] and returns a
# flattened array in the order that we reproduce here analytically.
#
# Critically:
#   - coords are in [-1,1]^3
#   - flattening order is consistent with original inf() + binary_time_data_convert
# ---------------------------------------------------------------------
def build_coords_chunk(H, W, T, z_start, z_end, torch_dtype):
    """
    Return a CPU tensor with coords for t in [z_start, z_end), in the
    SAME order and values as original get_mgrid([T,H,W], dim=3).
    """
    z_count = z_end - z_start
    total = z_count * H * W

    coords = np.empty((total, 3), dtype=np.float64)

    invT = 1.0 / max(T - 1, 1)
    invH = 1.0 / max(H - 1, 1)
    invW = 1.0 / max(W - 1, 1)

    idx = 0
    for t in range(z_start, z_end):
        # time index normalized [0,1] → [-1,1]
        t_norm = (t * invT - 0.5) * 2.0

        # We want to reproduce the flattened ordering that comes from:
        # np.mgrid[:T, :W, :H] → shape (T,W,H) → ravel('C')
        # Then original inf() remaps these to a [T,H,W] array using a
        # specific permutation. The net effect is:
        #   global linear index j = t*(H*W) + y*W + x
        #   coords[j] must equal get_mgrid(T,H,W)[j]
        #
        # For each (y,x), define:
        #   j_slice = y*W + x
        #   w_idx   = j_slice // H  in [0, W-1]
        #   h_idx   = j_slice %  H  in [0, H-1]
        #
        # These (w_idx,h_idx) are exactly the (W,H) indices used in the
        # original get_mgrid grid. So we compute coords from them:
        for y in range(H):
            for x in range(W):
                j_slice = y * W + x
                w_idx = j_slice // H   # 0..W-1
                h_idx = j_slice %  H   # 0..H-1

                x_norm = (w_idx * invW - 0.5) * 2.0  # second coord
                y_norm = (h_idx * invH - 0.5) * 2.0  # third coord

                coords[idx, 0] = t_norm
                coords[idx, 1] = x_norm
                coords[idx, 2] = y_norm
                idx += 1

    np_dtype = np.float64 if torch_dtype == torch.float64 else np.float32
    return torch.from_numpy(coords.astype(np_dtype))


# ---------------------------------------------------------------------
# Load TorchScript INR (float64 if possible)
# ---------------------------------------------------------------------
def load_ts_model(path, device):
    print(f"[pipeline] Loading TorchScript model: {path}")
    model = torch.jit.load(path, map_location=device)
    model.eval()

    # Try float64 first
    try:
        test = torch.zeros(1, 3, dtype=torch.float64, device=device)
        _ = model(test)
        dtype = torch.float64
        print("[pipeline] model accepted float64 inputs")
    except Exception as e64:
        print(f"[pipeline] float64 probe failed: {e64}")
        test = torch.zeros(1, 3, dtype=torch.float32, device=device)
        _ = model(test)
        dtype = torch.float32
        print("[pipeline] model accepted float32 inputs")

    return model, dtype


# ---------------------------------------------------------------------
# Evaluate INR on one chunk: coords_cpu → values[z_count,H,W]
# ---------------------------------------------------------------------
def eval_inr_chunk(model, coords_cpu, H, W, batch_size, device, torch_dtype):
    loader = DataLoader(
        coords_cpu,
        batch_size=batch_size,
        shuffle=False,
        pin_memory=(device.type == "cuda"),
    )

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device, non_blocking=(device.type == "cuda")).to(torch_dtype)
            out = model(batch)
            preds.append(out.view(-1).detach().cpu().numpy())

        if device.type == "cuda":
            torch.cuda.synchronize()

    preds = np.concatenate(preds)
    # We have total = z_count * H * W entries, in exactly the same
    # (t,y,x) linear order as the original inf() logic reconstructs.
    values = preds.reshape(-1, H, W)  # (z_count, H, W)
    return values


# ---------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Sample INR in chunks and call pvpython to compute critical points."
    )

    ap.add_argument("--ts_model", type=str, required=True,
                    help="Path to TorchScript INR model (.pt).")
    ap.add_argument("--dataset", type=str, required=True,
                    help="Dataset name, e.g. vortex_street_3d.")
    ap.add_argument("--up_sample_ratio", type=int, default=2,
                    help="Upsampling ratio for sampling grid.")
    ap.add_argument("--chunk_size", type=int, default=5,
                    help="Number of time slices per chunk.")
    ap.add_argument("--batch_size", type=int, default=100000,
                    help="Batch size for INR evaluation.")
    ap.add_argument("--output_dir", type=str, required=True,
                    help="Directory to store final CSV (and temporary bins).")
    ap.add_argument("--output_csv_name", type=str, default="critical_points_inr.csv",
                    help="Name of the final critical-point CSV.")
    ap.add_argument("--pvpython", type=str, required=True,
                    help="Path to pvpython executable.")
    ap.add_argument("--pv_script", type=str, required=True,
                    help="Path to compute_critical_points_from_bin.py (pvpython script).")
    ap.add_argument("--server", type=int, default=0,
                    help="0 for local ParaView, 1 for remote server build.")
    ap.add_argument("--no_cuda", action="store_true",
                    help="Force CPU even if CUDA is available.")

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    final_csv = os.path.join(args.output_dir, args.output_csv_name)

    use_cuda = (not args.no_cuda) and torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    print("[pipeline] DEVICE:", device)

    H, W, T, dom_min, dom_max = dataset_span_and_domain(
        args.dataset, args.up_sample_ratio
    )
    print(f"[pipeline] dataset={args.dataset}, H={H}, W={W}, T={T}")
    print(f"[pipeline] dom_min={dom_min}, dom_max={dom_max}")

    model, torch_dtype = load_ts_model(args.ts_model, device)
    print(f"[pipeline] Using dtype {torch_dtype} for model inputs")

    first_chunk = True
    z = 0

    while z < T:
        z_start = z
        z_end = min(z_start + args.chunk_size, T)
        z_count = z_end - z_start

        print(f"\n[pipeline] === CHUNK {z_start} → {z_end - 1} (z_count={z_count}) ===")

        # 1) Build coords on CPU, matching GetTestingData()/get_mgrid
        coords_cpu = build_coords_chunk(H, W, T, z_start, z_end, torch_dtype)

        # 2) Evaluate INR on this chunk
        values_chunk = eval_inr_chunk(
            model=model,
            coords_cpu=coords_cpu,
            H=H,
            W=W,
            batch_size=args.batch_size,
            device=device,
            torch_dtype=torch_dtype,
        )

        # 3) Save temporary bin in float64 (for TTK)
        bin_name = f"inr_chunk_{z_start:04d}_{z_end - 1:04d}.bin"
        bin_path = os.path.join(args.output_dir, bin_name)
        print(f"[pipeline] Writing bin: {bin_path}")
        values_chunk.astype("<f8").ravel(order="C").tofile(bin_path)

        # 4) Call pvpython to compute critical points for this chunk
        cmd = [
            args.pvpython,
            args.pv_script,
            "--input_bin", bin_path,
            "--dataset", args.dataset,
            "--up_sample_ratio", str(args.up_sample_ratio),
            "--z_start_idx", str(z_start),
            "--z_count", str(z_count),
            "--float_type", "float64",
            "--output_csv", final_csv,
            "--server", str(args.server),
        ]
        if not first_chunk:
            cmd.append("--append")

        print("[pipeline] Running pvpython:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        first_chunk = False

        # 5) Delete bin after use to save disk
        print("[pipeline] Deleting bin:", bin_path)
        try:
            os.remove(bin_path)
        except OSError:
            pass

        z = z_end

    print("\n[pipeline] All done!")
    print("[pipeline] Final CSV:", final_csv)


if __name__ == "__main__":
    main()