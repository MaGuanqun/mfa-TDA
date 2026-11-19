"""Export the signed u and v fields used by the vortex street CoordNet fits."""

from pathlib import Path
from tempfile import TemporaryDirectory

from netCDF4 import Dataset
import numpy as np


def main():
    data_dir = Path(__file__).resolve().parents[1] / "Data"
    with Dataset(data_dir / "cylinder2d.nc") as source, TemporaryDirectory(dir=data_dir) as staging:
        for component in ("u", "v"):
            variable = source.variables[component]
            if variable.dimensions != ("tdim", "ydim", "xdim"):
                raise ValueError(f"Unexpected {component} axes: {variable.dimensions}")
            field = np.ma.filled(variable[1351:1501, :, :], -2.0).astype("<f4", copy=False)
            if field.shape != (150, 80, 640):
                raise ValueError(f"Unexpected {component} grid: {field.shape}")
            path = Path(staging) / f"vortex_street_3d_{component}.bin"
            field.tofile(path)
            print(f"{path.name}: {field.shape}, min={field.min()}, max={field.max()}")
        for path in Path(staging).iterdir():
            path.replace(data_dir / path.name)


if __name__ == "__main__":
    main()
