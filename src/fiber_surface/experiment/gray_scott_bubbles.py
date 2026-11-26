"""Sample the two Gray–Scott CoordNets and extract their bivariate TTK surface."""

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "CoordNet/Data"
OUTPUT = ROOT / "build/gray-scott-bubbles"
QUERY = np.array([[0.32, 0.35], [0.77, 0.13], [0.32, 0.56]])
NAMES = ("A", "B")


def source_fields(meta):
    shape = tuple(reversed(meta["dimensions_xyz"]))
    return [np.memmap(DATA / meta["fields"][k]["file"], mode="r", dtype="<f4", shape=shape) for k in NAMES]


def membership(fields):
    a, b = fields
    return np.logical_and.reduce([(q[0]-p[0])*(b-p[1]) >= (q[1]-p[1])*(a-p[0])
                                  for p, q in zip(QUERY, np.roll(QUERY, -1, axis=0))])


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sample(meta, folder, device):
    import torch

    source = source_fields(meta)
    nt, ny, nx = source[0].shape
    path = folder / "coordnet_concentrations.npy"
    pred = np.lib.format.open_memmap(path.with_suffix(".partial.npy"), mode="w+", dtype="float64",
                                     shape=(2, nt, ny, nx))
    metrics = {}
    # Same ordering and independent per-field normalization as the training loader.
    with torch.inference_mode():
        for j, name in enumerate(NAMES):
            model_path = ROOT / f"CoordNet/Exp/gray_scott_bubbles_{name}/super-spatial-temporal-64-5-float64.pt"
            model = torch.jit.load(str(model_path), map_location=device).eval()
            values = pred[j].reshape(-1)
            low = meta["fields"][name]["minimum"]
            high = meta["fields"][name]["maximum"]
            for start in range(0, len(values), 8192):
                stop = min(start+8192, len(values))
                ids = np.arange(start, stop)
                ijk = np.column_stack(np.unravel_index(ids, (nt, ny, nx)))
                coords = 2.0 * ijk / np.array([nt-1, ny-1, nx-1]) - 1.0
                output = model(torch.as_tensor(coords, device=device)).cpu().numpy().ravel()
                values[start:stop] = low + (output+1.0)*(high-low)/2.0
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite model output for {name}")
            error = pred[j] - source[j]
            rmse = float(np.sqrt(np.mean(error**2)))
            metrics[name] = dict(rmse=rmse, range_normalized_rmse=rmse/(high-low),
                                 psnr_db=float(20*np.log10((high-low)/rmse)),
                                 max_abs_error=float(np.abs(error).max()))
    original, fitted = membership(source), membership(pred)
    intersection = (original & fitted).sum(axis=(1, 2))
    union = (original | fitted).sum(axis=(1, 2))
    metrics.update(triangle_AB=QUERY.tolist(), grid_points=original.size,
                   source_selected=int(original.sum()), model_selected=int(fitted.sum()),
                   intersection=int(intersection.sum()), union=int(union.sum()),
                   membership_iou=float(intersection.sum()/union.sum()),
                   per_frame_iou=np.divide(intersection, union, out=np.ones(nt), where=union>0).tolist())
    pred.flush()
    del pred
    path.with_suffix(".partial.npy").replace(path)
    write_json(folder / "fit.json", metrics)
    print(f"Native-grid membership IoU {metrics['membership_iou']:.6f}", flush=True)


def time_indices(fields, count):
    # Split the interval with the largest RMS temporal interpolation error.
    # Later, slowly evolving frames can then share wider intervals.
    nt = fields[0].shape[0]
    if count == nt:
        return np.arange(nt)
    samples = np.stack([a[:, ::4, ::4] / np.ptp(a) for a in fields], axis=-1)
    frames = [0, nt-1]
    while len(frames) < count:
        errors = np.full(nt, -1.0)
        for left, right in zip(frames, frames[1:]):
            ids = np.arange(left+1, right)
            if not len(ids):
                continue
            s = ((ids-left)/(right-left))[:, None, None, None]
            residual = (1-s)*samples[left] + s*samples[right] - samples[ids]
            errors[ids] = np.mean(residual**2, axis=(1, 2, 3))
        frames.append(int(errors.argmax()))
        frames.sort()
    return np.array(frames)


def extract(meta, fields, spatial_count, temporal_count, folder, source_name):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk
    from ttk_geometry import extract_triangle
    from evaluate_curved_prism import save_polydata

    nt, ny, nx = fields[0].shape
    xi = np.rint(np.linspace(0, nx-1, spatial_count)).astype(int)
    yi = np.rint(np.linspace(0, ny-1, spatial_count)).astype(int)
    ti = time_indices(fields, temporal_count)
    indices = np.ix_(ti, yi, xi)
    grid = vtk.vtkImageData()
    grid.SetDimensions(len(xi), len(yi), len(ti))
    for name, values in zip(NAMES, fields):
        data = numpy_to_vtk(np.asarray(values[indices], dtype=np.float64).ravel(), deep=True)
        data.SetName(name)
        grid.GetPointData().AddArray(data)
    xyz, faces, labels, edges, corners, collapsed = extract_triangle(grid, NAMES, QUERY, 3, True)
    # The sampled axes may have unequal intervals. Mapping the index-space TTK
    # geometry through its owning grid cell is affine on every tetrahedron.
    physical_axes = [np.linspace(lo, hi, n)[ids] for lo, hi, n, ids in
                     zip(meta["origin"], meta["upper"], [nx, ny, nt], [xi, yi, ti])]
    physical = np.column_stack([np.interp(xyz[:, j], np.arange(len(axis)), axis)
                                for j, axis in enumerate(physical_axes)])
    folder.mkdir(parents=True, exist_ok=True)
    save_polydata(folder / "fiber_surface.vtp", physical, faces, labels)
    if len(edges):
        save_polydata(folder / "corner_fibers.vtp", physical, edges, corners, lines=True)
    if source_name == "CoordNet" and (len(xi), len(yi), len(ti)) == (nx, ny, nt):
        grid.SetOrigin(*meta["origin"])
        grid.SetSpacing(*[(axis[-1]-axis[0])/(len(axis)-1) for axis in physical_axes])
        grid.GetPointData().SetActiveScalars("A")
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetInputData(grid)
        writer.SetFileName(str(folder / "sampled_coordnet.vti"))
        if not writer.Write():
            raise OSError("Failed to save the sampled CoordNet concentrations")
    counts = np.bincount(labels, minlength=4)[1:].tolist()
    report = dict(input=source_name, vtk=vtk.vtkVersion.GetVTKVersion(), triangle_AB=QUERY.tolist(),
                  domain_axes=["x", "y", "simulation_time"], x_indices=xi.tolist(), y_indices=yi.tolist(),
                  time_indices_in_window=ti.tolist(), simulation_times=physical_axes[2].tolist(),
                  triangles=len(faces), triangles_per_side=counts, vertices=len(np.unique(faces)),
                  corner_segments=np.bincount(corners, minlength=3).tolist(), collapsed_triangles=collapsed)
    write_json(folder / "surface.json", report)
    print(f"{source_name} {spatial_count} x {spatial_count} x {temporal_count}: {len(faces)} triangles, sides {counts}",
          flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest="stage", required=True)
    s = commands.add_parser("sample")
    s.add_argument("--device", default="cpu")
    e = commands.add_parser("extract")
    e.add_argument("--size", type=int, required=True)
    e.add_argument("--times", type=int, required=True)
    e.add_argument("--output", type=Path, required=True)
    e.add_argument("--raw", action="store_true", help="extract a source-data comparison")
    args = p.parse_args()
    meta = json.loads((DATA / "gray_scott_bubbles.json").read_text())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "sample":
        sample(meta, OUTPUT, args.device)
    else:
        if not 2 <= args.size <= 128 or not 2 <= args.times <= 128:
            p.error("Both sampled axis counts must be between 2 and 128")
        fields = source_fields(meta) if args.raw else np.load(OUTPUT / "coordnet_concentrations.npy", mmap_mode="r")
        extract(meta, fields, args.size, args.times, args.output, "original" if args.raw else "CoordNet")


if __name__ == "__main__":
    main()
