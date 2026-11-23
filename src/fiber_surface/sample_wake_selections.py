"""Estimate the paper's two wake selections directly from the full CoordNet field."""

import json
from pathlib import Path

import numpy as np
import torch


def main():
    root = Path(__file__).resolve().parents[2] / "CoordNet"
    model_dir = root / "Exp" / "vortex_street_3d_2_variables"
    models = []
    for component in ("u", "v"):
        path = model_dir / f"{component}-super-spatial-temporal-64-5-float64.pt"
        model = torch.jit.load(str(path), map_location="cpu").eval()
        data = np.memmap(root / "Data" / f"vortex_street_3d_{component}.bin", dtype="<f4", mode="r")
        models.append((model, float(data.min()), float(data.max())))

    lo = torch.tensor([-0.5, -0.5, 13.5], dtype=torch.float64)
    hi = torch.tensor([7.5, 0.5, 15.0], dtype=torch.float64)
    n = 100_000
    seed = 2026
    rng = torch.Generator().manual_seed(seed)
    results = []
    for time in (None, 13.8, 14.1, 14.4):
        xyz = lo + (hi - lo) * torch.rand((n, 3), dtype=torch.float64, generator=rng)
        if time is not None:
            xyz[:, 2] = time
        inputs = (2 * (xyz - lo) / (hi - lo) - 1).flip((1,))
        up = down = 0
        for batch in inputs.split(1024):
            velocity = []
            for model, minimum, maximum in models:
                output = model(batch).reshape(-1)
                if len(output) != len(batch) or not torch.isfinite(output).all():
                    raise ValueError("Expected one finite velocity per sample")
                velocity.append(minimum + (output + 1) * ((maximum - minimum) / 2))
            u, v = velocity
            # At fixed u the triangle height interpolates from 0 to 0.35 to 0.
            height = (0.35 / 0.15) * torch.minimum(u - 0.60, 0.90 - u)
            up += int(((0 <= v) & (v <= height)).sum())
            down += int(((0 <= -v) & (-v <= height)).sum())
        measure = "volume" if time is None else "area"
        extent = float((hi - lo).prod() if time is None else (hi - lo)[:2].prod())
        results.append({"time": time, "up_hits": up, "down_hits": down,
                        f"up_{measure}": extent * up / n, f"down_{measure}": extent * down / n})
    print(json.dumps({"domain_min": lo.tolist(), "domain_max": hi.tolist(),
                      "seed": seed, "samples_per_estimate": n, "results": results}, indent=2))


if __name__ == "__main__":
    torch.set_num_threads(2)
    with torch.inference_mode():
        main()
