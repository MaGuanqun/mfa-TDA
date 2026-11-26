"""Export The Well's Gray–Scott Bubbles test trajectory for the CoordNet experiment."""

import argparse
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np


FILENAME = "gray_scott_reaction_diffusion_bubbles_F_0.098_k_0.057.hdf5"
SOURCE = "https://huggingface.co/datasets/polymathic-ai/gray_scott_reaction_diffusion/resolve/main/data/test/" + FILENAME
SHA256 = "db6bac94402d5c7a46e8c335d80e77c1f293fabda0d6521ba1dfbbe0c65c2ba1"
DATA = Path(__file__).resolve().parents[1] / "Data"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--frames", type=int, choices=(64, 128), default=128)
    args = parser.parse_args()
    digest = hashlib.sha256()
    with args.archive.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != SHA256:
        raise ValueError("Unexpected Bubbles test archive checksum")
    DATA.mkdir(exist_ok=True)
    trajectory, first = 1, 2
    stop = first + args.frames
    shape = (128, 128, args.frames)
    prefix = "gray_scott_bubbles" if args.frames == 128 else "gray_scott_bubbles_first64"
    with h5py.File(args.archive) as f, TemporaryDirectory(dir=DATA) as tmp:
        dims = f["dimensions"]
        axes = [dims["x"][:], dims["y"][:], dims["time"][first:stop]]
        if any(len(v) != n or not np.allclose(np.diff(v), v[1]-v[0], rtol=1e-5) for v, n in zip(axes, shape)):
            raise ValueError(f"Expected uniform x, y, t axes of lengths {shape}")
        meta = dict(source=SOURCE, source_sha256=SHA256, split="test", trajectory=trajectory,
                    frame_start=first, frame_stop_exclusive=stop, F=float(f.attrs["F"]),
                    k=float(f.attrs["k"]), binary_order=["time", "y", "x"],
                    dimensions_xyz=list(shape), origin=[float(v[0]) for v in axes],
                    upper=[float(v[-1]) for v in axes], fields={})
        for name in ("A", "B"):
            field = f["t0_fields"][name][trajectory, first:stop].transpose(0, 2, 1)
            if field.shape != shape[::-1] or not np.isfinite(field).all() or field.max() <= field.min():
                raise ValueError(f"Invalid field {name}")
            path = Path(tmp) / f"{prefix}_{name}.bin"
            field.astype("<f4", copy=False).tofile(path)
            meta["fields"][name] = dict(file=path.name, minimum=float(field.min()), maximum=float(field.max()))
        (Path(tmp) / f"{prefix}.json").write_text(json.dumps(meta, indent=2) + "\n")
        for p in Path(tmp).iterdir():
            p.replace(DATA / p.name)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
