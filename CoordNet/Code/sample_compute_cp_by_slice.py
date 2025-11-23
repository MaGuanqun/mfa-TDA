#!/usr/bin/env python
import os
import argparse
import subprocess
import numpy as np
import torch
from torch.utils.data import DataLoader


# ---------------------------------------------------------------------
# Dataset geometry
# ---------------------------------------------------------------------
def dataset_span_and_domain(name: str, up_sample_ratio: int):
    if name == "vortex_street_3d":
        base_span = np.array([80, 10, 15])
        dom_min  = np.array([-0.5, -0.5, 13.5])
        dom_max  = np.array([ 7.5,  0.5, 15.0])
    elif name == "boussinesq_3d":
        base_span = np.array([10, 30, 15])
        dom_min  = np.array([-0.5, -0.5, 0.0])
        dom_max  = np.array([ 0.5,  2.5, 1.5])
    else:
        raise ValueError(f"Unknown dataset {name}")

    H, W, T = (base_span * up_sample_ratio).astype(int)
    return H, W, T, dom_min, dom_max


# ---------------------------------------------------------------------
# Build coords on CPU
# ---------------------------------------------------------------------
def build_coords_chunk(H, W, T, z_start, z_end, dtype):
    z_count = z_end - z_start
    total = z_count * H * W

    coords = np.empty((total, 3), dtype=np.float64)

    invT = 1.0 / max(T - 1, 1)
    invH = 1.0 / max(H - 1, 1)
    invW = 1.0 / max(W - 1, 1)

    idx = 0
    for t in range(z_start, z_end):
        tn = (t * invT - 0.5) * 2
        for y in range(H):
            yn = (y * invH - 0.5) * 2
            for x in range(W):
                xn = (x * invW - 0.5) * 2
                coords[idx] = (tn, yn, xn)
                idx += 1

    return torch.from_numpy(coords).to(dtype=dtype)   # CPU tensor


# ---------------------------------------------------------------------
# Load TorchScript INR
# ---------------------------------------------------------------------
def load_ts_model(path, device):
    model = torch.jit.load(path, map_location=device)
    model.eval()

    # test float64
    try:
        test = torch.zeros(1, 3, dtype=torch.float64, device=device)
        _ = model(test)
        dtype = torch.float64
    except:
        dtype = torch.float32

    return model, dtype


# ---------------------------------------------------------------------
# Evaluate chunk
# ---------------------------------------------------------------------
def eval_inr_chunk(model, coords_cpu, H, W, batch_size, device, dtype):
    loader = DataLoader(coords_cpu,
                        batch_size=batch_size,
                        shuffle=False,
                        pin_memory=(device.type == "cuda"))

    preds = []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device, non_blocking=True).to(dtype)
            out = model(batch)
            preds.append(out.view(-1).detach().cpu().numpy())

    preds = np.concatenate(preds)
    values = preds.reshape(-1, H, W)
    return values


# ---------------------------------------------------------------------
# Main Orchestrator
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ts_model", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--up_sample_ratio", type=int, default=2)
    ap.add_argument("--chunk_size", type=int, default=5)
    ap.add_argument("--batch_size", type=int, default=100000)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--pvpython", required=True)
    ap.add_argument("--pv_script", required=True)
    ap.add_argument("--output_csv_name", default="critical_points_inr.csv")
    ap.add_argument("--server", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    final_csv = os.path.join(args.output_dir, args.output_csv_name)

    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    print("DEVICE:", device)

    H, W, T, _, _ = dataset_span_and_domain(args.dataset, args.up_sample_ratio)

    model, dtype = load_ts_model(args.ts_model, device)
    print("Model dtype:", dtype)

    first = True
    z = 0

    while z < T:
        z_start = z
        z_end   = min(z + args.chunk_size, T)
        z_count = z_end - z_start

        print(f"\n=== CHUNK {z_start} → {z_end-1} ===")

        # 1. Build coords (CPU)
        coords_cpu = build_coords_chunk(H, W, T, z_start, z_end, dtype)

        # 2. Evaluate model
        values = eval_inr_chunk(model, coords_cpu, H, W,
                                args.batch_size, device, dtype)

        # 3. Save bin temporarily
        bin_path = os.path.join(
            args.output_dir,
            f"inr_chunk_{z_start:04d}_{z_end-1:04d}.bin"
        )
        values.astype("<f8").ravel(order="C").tofile(bin_path)

        # 4. Call pvpython to extract CPs
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
        if not first:
            cmd.append("--append")
        first = False

        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

        # 5. DELETE BIN FILE AFTER USE
        print("Deleting bin:", bin_path)
        try:
            os.remove(bin_path)
        except:
            pass

        z = z_end

    print("\nAll done!")
    print("Final CSV:", final_csv)


if __name__ == "__main__":
    main()