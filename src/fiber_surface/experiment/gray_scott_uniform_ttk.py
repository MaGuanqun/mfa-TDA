"""Evaluate Gray–Scott CoordNets on a uniform grid and write the TTK output port."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from gray_scott_bubbles import DATA, NAMES, QUERY, write_json


def bounds(meta, time_max):
    return np.array([*meta["origin"][:2], 0.]), np.array([*meta["upper"][:2], time_max])


def sample_component(folder, name, dims, device, meta, time_max):
    import torch

    lo, hi = bounds(meta, time_max)
    x, y, t = [np.linspace(a, b, n) for a, b, n in zip(lo, hi, dims)]
    dataset = Path(meta["fields"][name]["file"]).stem
    path = folder / "models" / dataset / "super-spatial-temporal-64-5-float64.pt"
    model = torch.jit.load(str(path), map_location=device).eval()
    lower, upper = [meta["fields"][name][key] for key in ("minimum", "maximum")]
    partial = folder / f"{name}.partial.npy"
    values = np.lib.format.open_memmap(partial, mode="w+", dtype="float64", shape=dims[::-1])
    flat = values.reshape(-1)
    with torch.inference_mode():
        for start in range(0, len(flat), 8192):
            stop = min(start+8192, len(flat))
            k, j, i = np.unravel_index(np.arange(start, stop), values.shape)
            xyz = np.column_stack((x[i], y[j], t[k]))
            # TorchScript takes normalized (time, y, x); VTK stores x fastest.
            normalized = (2*(xyz-lo)/(hi-lo)-1)[:, ::-1].copy()
            out = model(torch.as_tensor(normalized, device=device)).flatten().cpu().numpy()
            flat[start:stop] = lower + (out+1)*(upper-lower)/2
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite model predictions for {name}")
    value_range = [float(values.min()), float(values.max())]
    values.flush()
    del flat, values
    partial.replace(folder / f"{name}.npy")
    count = int(np.prod(dims))
    write_json(folder / f"sample_{name}.json", dict(field=name, grid_points=count,
               dimensions_xyz=list(dims), origin=lo.tolist(), upper=hi.tolist(),
               model_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), value_range=value_range))
    print(f"{name}: evaluated {count:,} model coordinates", flush=True)


def check_surface_on_disk(raw, path):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    saved = reader.GetOutput()
    if (saved.GetNumberOfCells(), saved.GetNumberOfPoints()) != (raw.GetNumberOfCells(), raw.GetNumberOfPoints()):
        raise ValueError("Surface size changed on serialization")
    def compare(a, b):
        if not np.array_equal(vtk_to_numpy(a), vtk_to_numpy(b)):
            raise ValueError(f"VTK array changed on serialization: {a.GetName()}")
    compare(raw.GetPoints().GetData(), saved.GetPoints().GetData())
    for get_cells in ("GetVerts", "GetLines", "GetPolys", "GetStrips"):
        a, b = getattr(raw, get_cells)(), getattr(saved, get_cells)()
        compare(a.GetOffsetsArray(), b.GetOffsetsArray())
        compare(a.GetConnectivityArray(), b.GetConnectivityArray())
    for a, b in ((raw.GetPointData(), saved.GetPointData()), (raw.GetCellData(), saved.GetCellData())):
        if a.GetNumberOfArrays() != b.GetNumberOfArrays():
            raise ValueError("Attribute count changed on serialization")
        for index in range(a.GetNumberOfArrays()):
            array = a.GetArray(index)
            compare(array, b.GetArray(array.GetName()))


def extract(folder, dims, meta, time_max):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
    from ttk_geometry import triangle_filter

    lo, hi = bounds(meta, time_max)
    spacing = (hi-lo)/(np.asarray(dims)-1)
    grid = vtk.vtkImageData()
    grid.SetDimensions(*dims)
    grid.SetOrigin(*lo)
    grid.SetSpacing(*spacing)
    for name in NAMES:
        values = np.load(folder / f"{name}.npy", mmap_mode="r")
        metadata = json.loads((folder / f"sample_{name}.json").read_text())
        if (values.shape != dims[::-1] or metadata["dimensions_xyz"] != list(dims)
                or metadata["origin"] != lo.tolist() or metadata["upper"] != hi.tolist()):
            raise ValueError(f"Mismatched grid for {name}")
        array = numpy_to_vtk(values.ravel(), deep=True)
        array.SetName(name)
        grid.GetPointData().AddArray(array)
    grid.GetPointData().SetActiveScalars(NAMES[0])
    volume_path = folder / "sampled_coordnet.partial.vti"
    vw = vtk.vtkXMLImageDataWriter()
    vw.SetFileName(str(volume_path))
    vw.SetInputData(grid)
    if not vw.Write():
        raise OSError(volume_path)
    volume_path.replace(folder / "sampled_coordnet.vti")

    filt = triangle_filter(grid, NAMES, QUERY, 3, True)
    filt.Update()
    raw = filt.GetOutput()
    n = raw.GetNumberOfPolys()
    if not n or raw.GetNumberOfCells() != n or not np.all(np.diff(vtk_to_numpy(raw.GetPolys().GetOffsetsArray())) == 3):
        raise ValueError("Expected a nonempty TTK triangle surface")
    path = folder / "fiber_surface.partial.vtp"
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputConnection(filt.GetOutputPort())
    if not writer.Write():
        raise OSError(path)
    check_surface_on_disk(raw, path)
    path.replace(folder / "fiber_surface.vtp")
    sides = np.bincount(vtk_to_numpy(raw.GetCellData().GetArray("EdgeIds")), minlength=3)
    report = dict(source_trajectory=meta["trajectory"], frame_start=meta["frame_start"],
                  frame_stop_exclusive=meta["frame_stop_exclusive"], training_grid=meta["dimensions_xyz"],
                  source_time=[meta["origin"][2], meta["upper"][2]], rescaled_time=[lo[2], hi[2]],
                  training_time_spacing=float((hi[2]-lo[2])/(meta["dimensions_xyz"][2]-1)),
                  dimensions_xyz=list(dims), grid_points=grid.GetNumberOfPoints(), grid_bounds=list(grid.GetBounds()),
                  grid_spacing=spacing.tolist(), triangle_AB=QUERY.tolist(),
                  triangles=n, points=raw.GetNumberOfPoints(), mesh_bounds=list(raw.GetBounds()),
                  triangles_by_edge_id=sides.tolist(), vtk=vtk.vtkVersion.GetVTKVersion(),
                  round_trip="Verified points, cell connectivity and offsets, point data and cell data against ttkFiberSurface.")
    write_json(folder / "surface.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--size", type=int, default=200)
    p.add_argument("--time-size", type=int, help="Time samples (defaults to --size)")
    p.add_argument("--time-max", type=float, default=2.0)
    p.add_argument("--metadata", type=Path, default=DATA / "gray_scott_bubbles.json")
    stages = p.add_subparsers(dest="stage", required=True)
    s = stages.add_parser("sample")
    s.add_argument("name", choices=NAMES)
    s.add_argument("--device", default="cpu")
    stages.add_parser("extract")
    args = p.parse_args()
    dims = (args.size, args.size, args.time_size if args.time_size is not None else args.size)
    if min(dims) < 2:
        p.error("Each axis needs at least two samples")
    if not np.isfinite(args.time_max) or args.time_max <= 0:
        p.error("--time-max must be finite and positive")
    settings = json.loads(args.metadata.read_text())
    if args.stage == "sample":
        sample_component(args.folder, args.name, dims, args.device, settings, args.time_max)
    else:
        extract(args.folder, dims, settings, args.time_max)
