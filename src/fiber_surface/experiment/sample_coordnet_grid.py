"""Sample the two saved CoordNet fields at the nodes used for TTK and mesh IoU."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex

import numpy as np


def run_parameters(run):
    config = (run / "stages/run.conf").read_bytes()
    lines = [shlex.split(line) for line in config.decode().splitlines()]
    fields = {row[0]: row[1:] for row in lines if row}
    if (fields["model"] != ["coordnet"] or fields.get("polyline") != ["4"]
            or fields["FIBER_SURFACE_CONFIG"] != ["3"]):
        raise ValueError("Expected a staged CoordNet run with a three-sided query")
    start = lines.index(["polyline", "4"]) + 1
    polygon = np.array(lines[start:start + 4], dtype=float)
    if polygon.shape != (4, 2) or not np.array_equal(polygon[0], polygon[3]):
        raise ValueError("Expected a closed triangle")
    a, b = polygon[1:3] - polygon[0]
    if a[0]*b[1] <= a[1]*b[0]:
        raise ValueError("Expected counterclockwise corners in model-output units")
    return dict(config_sha256=hashlib.sha256(config).hexdigest(),
                models=[fields[k][0] for k in ("first", "second")],
                model_digests=[fields[k][1] for k in ("first", "second")],
                triangle=polygon[:3].tolist(), step=float(fields["options"][0]),
                **{k: list(map(float, fields[k]))
                   for k in ("box_min", "box_max", "domain_min", "domain_max")})


def read_grid(folder):
    grid = json.loads((folder / "grid.json").read_text())
    fields = np.load(folder / "fields.npy", mmap_mode="r")
    if fields.shape != (2, *grid["dimensions"][::-1]) or fields.dtype != np.dtype("<f8"):
        raise ValueError("Sample layout differs from grid.json")
    return grid, fields


def sample(run, dest, dims):
    import torch

    meta = run_parameters(run)
    if min(dims) < 2:
        raise ValueError("Each axis needs at least two grid nodes")
    models = []
    for name, digest in zip(meta["models"], meta["model_digests"]):
        # Same FNV-1a digest as file_digest() in pipeline_io.cpp.
        h = 14695981039346656037
        for byte in Path(name).read_bytes():
            h = ((h ^ byte) * 1099511628211) & 0xffffffffffffffff
        if f"{h:016x}" != digest:
            raise ValueError(f"Model changed since seed discovery: {name}")
        models.append(torch.jit.load(name, map_location="cpu").eval())

    dest.mkdir(parents=True, exist_ok=False)
    nx, ny, nt = dims
    axes = [np.linspace(a, b, n) for a, b, n in zip(meta["box_min"], meta["box_max"], dims)]
    scale = 2 / (np.array(meta["domain_max"]) - meta["domain_min"])
    samples = np.lib.format.open_memmap(dest / "fields.npy", mode="w+", dtype="<f8",
                                        shape=(2, nt, ny, nx))
    flat = samples.reshape(2, -1)
    torch.set_num_threads(1)
    with torch.inference_mode():
        for start in range(0, math.prod(dims), 2048):
            ids = np.arange(start, min(start + 2048, flat.shape[1]))
            xyz = np.column_stack([axes[0][ids % nx], axes[1][ids // nx % ny], axes[2][ids // (nx*ny)]])
            normalized = (xyz - meta["domain_min"]) * scale - 1
            inputs = torch.frombuffer(normalized[:, ::-1].copy(), dtype=torch.float64).reshape(-1, 3)
            for i, model in enumerate(models):
                values = model(inputs).reshape(-1)
                if (values.numel() != len(ids) or values.dtype != torch.float64
                        or not torch.isfinite(values).all()):
                    raise ValueError(f"Expected one finite float64 scalar per point from {meta['models'][i]}")
                flat[i, start:start + len(ids)] = values.tolist()
    samples.flush()
    spacing = [(b - a) / (n - 1) for a, b, n in zip(meta["box_min"], meta["box_max"], dims)]
    meta.update(dimensions=dims, spacing=spacing, torch_version=torch.__version__)
    (dest / "grid.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="directory containing stages/run.conf")
    parser.add_argument("output", type=Path, help="new directory for fields.npy and grid.json")
    parser.add_argument("--dims", type=int, nargs=3, required=True, metavar=("NX", "NY", "NT"))
    args = parser.parse_args()
    sample(args.run, args.output, args.dims)
