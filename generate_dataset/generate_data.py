#!/usr/bin/env python3

import argparse
import numpy as np
import torch

import function  # this is your /mnt/data/function.py


def tensor_to_raw(v: torch.Tensor, filename: str) -> None:
    """
    Save a tensor to a .raw file as little-endian float64.

    The tensor will be flattened in row-major order. For a tensor
    with shape (nt, ny, nx), this means x is fastest, then y, then t.
    """
    if isinstance(v, torch.Tensor):
        v = v.detach().cpu().contiguous().view(-1).numpy()
    else:
        v = np.asarray(v).ravel()

    # ensure little-endian float32
    v = v.astype('<f8', copy=False)

    with open(filename, 'wb') as f:
        v.tofile(f)


def generate_quartic_dataset(args):
    
    sx, sy, st = function.sample_size(args.function)  
    
    # Get domain from your function.py
    x_min, x_max, y_min, y_max, t_min, t_max = function.range(args.function)

    # 1D coordinates
    xs = torch.linspace(x_min, x_max, steps=sx)
    ys = torch.linspace(y_min, y_max, steps=sy)
    ts = torch.linspace(t_min, t_max, steps=st)

    # Make 3D grids:
    # T, Y, X all have shape (nt, ny, nx)
    T, Y, X = torch.meshgrid(ts, ys, xs, indexing='ij')

    # Compute scalar field
    if(args.function == 'quartic_potential_2'):
        values = function.compute_quartic_potential_2(X, Y, T)  # shape: (nt, ny, nx)

    return values


def main():
    parser = argparse.ArgumentParser(
        description="Generate a quartic_potential_2 3D dataset and save to .raw"
    )
    parser.add_argument(
        "--function", type=str, default=None,
        help="Number of samples along x (default: from sample_size in function.py)"
    )
    parser.add_argument(
        "--output", type=str, default="quartic_potential_2.raw",
        help="Output .raw filename (default: quartic_potential_2.raw)"
    )

    args = parser.parse_args()

    values = generate_quartic_dataset(args)

    print(f"Value tensor shape: {tuple(values.shape)}")
    
    tensor_to_raw(values, args.output)
    print(f"Saved raw data to {args.output}")
    print("Flattening order: x fastest, then y, then t (overall layout: (t, y, x) in memory).")


if __name__ == "__main__":
    main()