"""Build the upward-flow reference surface with ttkFiberSurface on the sampled INR."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from experiment.evaluate_curved_prism import save_polydata
from experiment.ttk_geometry import extract_triangle


ROOT = Path(__file__).resolve().parents[2]
WAKE = ROOT / "build/fiber-volume-results/vortex_street_up"
LO = np.array([-.5, -.5, 13.5])
HI = np.array([7.5, .5, 15.])
PHYSICAL_QUERY = np.array([[.60, 0.], [.90, 0.], [.75, .35]])


class WakeField:
    def __init__(self):
        model_dir = ROOT / "CoordNet/Exp/vortex_street_3d_2_variables"
        self.models = []
        extrema = []
        for component in ("u", "v"):
            name = f"{component}-super-spatial-temporal-64-5-float64.pt"
            self.models.append(torch.jit.load(str(model_dir / name), map_location="cpu").eval())
            path = ROOT / f"CoordNet/Data/vortex_street_3d_{component}.bin"
            values = np.memmap(path, dtype="<f4", mode="r")
            extrema.append([float(values.min()), float(values.max())])
        lower, upper = np.array(extrema).T
        self.query = 2*(PHYSICAL_QUERY-lower)/(upper-lower)-1

    @torch.inference_mode()
    def __call__(self, xyz):
        result = np.empty((len(xyz), 2))
        for start in range(0, len(xyz), 2048):
            batch = slice(start, start+2048)
            scaled = (2*(xyz[batch]-LO)/(HI-LO)-1)[:, ::-1].copy()
            tensor = torch.from_numpy(scaled)
            for i, model in enumerate(self.models):
                result[batch, i] = model(tensor).reshape(-1).numpy()
        if not np.isfinite(result).all():
            raise ValueError("Nonfinite INR output")
        return result

    def selected(self, uv):
        inside = np.ones(len(uv), dtype=bool)
        for a, b in zip(self.query, np.roll(self.query, -1, axis=0)):
            inside &= (uv-a) @ np.array([a[1]-b[1], b[0]-a[0]]) >= 0
        return inside


def scalars(dataset, uv):
    for i, name in enumerate(("u", "v")):
        a = numpy_to_vtk(uv[:, i].copy(), deep=True)
        a.SetName(name)
        dataset.GetPointData().AddArray(a)


def model_grid(field, h):
    dims = np.rint((HI-LO)/h).astype(int)+1
    axes = [np.linspace(a, b, n) for a, b, n in zip(LO, HI, dims)]
    t, y, x = np.meshgrid(axes[2], axes[1], axes[0], indexing="ij")
    return field(np.column_stack((x.ravel(), y.ravel(), t.ravel()))).reshape(*dims[::-1], 2)


def extract_reference(folder, grid, query):
    dims = np.array(grid.shape[:3][::-1])
    domain = vtk.vtkImageData()
    domain.SetOrigin(*LO)
    domain.SetSpacing(*((HI-LO)/(dims-1)))
    domain.SetDimensions(*dims)
    scalars(domain, grid.reshape(-1, 2))
    xyz, faces, patch, edges, labels, collapsed = extract_triangle(domain, ("u", "v"), query, 2, True)
    if len(np.unique(labels)) != 3:
        raise ValueError("Expected reference edges for each wake corner")
    folder.mkdir(exist_ok=True, parents=True)
    save_polydata(folder / "sides.vtp", xyz, faces, patch)
    save_polydata(folder / "fibers.vtp", xyz, edges, labels, lines=True)
    metadata = dict(dimensions=dims.tolist(), spacing=list(domain.GetSpacing()),
                    implicit_tetrahedra=int(6*np.prod(dims-1)), triangle_model=query.tolist(),
                    triangle_velocity=PHYSICAL_QUERY.tolist(), vtk=vtk.vtkVersion.GetVTKVersion(),
                    comparison="ttkFiberSurface on the image's implicit tetrahedralization",
                    merge="identical coordinates at TTK output precision", collapsed_triangles=collapsed,
                    side_triangles=np.bincount(patch, minlength=4)[1:].tolist(),
                    fiber_segments=np.bincount(labels).tolist())
    (folder / "reference.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(folder, metadata, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=Path, default=WAKE / "ttk")
    folder = parser.parse_args().folder
    torch.set_num_threads(2)
    field = WakeField()
    uv = model_grid(field, .025)
    folder.mkdir(exist_ok=True, parents=True)
    np.save(folder / "grid_values.npy", uv)
    for name, stride in (("grid005", 2), ("grid0025", 1)):
        extract_reference(folder / name, uv[::stride, ::stride, ::stride], field.query)
