"""Sample the fitted wake, its TTK piecewise-linear approximation and any extracted caps or volume."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch
import vtk

from experiment.evaluate_curved_prism import read_mesh
from wake_reference import WAKE, LO, HI, WakeField


FACE_NAMES = ("xmin", "xmax", "ymin", "ymax", "tmin", "tmax")


def wilson(k, n):
    z = 1.959963984540054
    p = k/n
    center = p + z*z/(2*n)
    radius = z*np.sqrt((p*(1-p)+z*z/(4*n))/n)
    return [max(0., (center-radius)/(1+z*z/n)), min(1., (center+radius)/(1+z*z/n))]


def agreement(target, hit):
    n = len(target)
    both = int((target & hit).sum())
    missed = int((target & ~hit).sum())
    excess = int((~target & hit).sum())
    union = both+missed+excess
    return dict(both=both, missed=missed, excess=excess, neither=n-union,
                iou=both/union if union else None, missed_fraction=missed/n, excess_fraction=excess/n,
                missed_wilson95=wilson(missed, n), excess_wilson95=wilson(excess, n))


def linear_values(grid, xyz):
    # TTK's implicit grid uses the body diagonal B=(1,0,0) to G=(0,1,1).
    # Sort coordinates measured from B to find the containing one of six tets.
    dims = np.array(grid.shape[:3][::-1])
    pos = (xyz-LO)/(HI-LO)*(dims-1)
    origin = np.clip(np.floor(pos).astype(int), 0, dims-2)
    local = pos-origin
    local[:, 0] = 1-local[:, 0]
    order = np.argsort(-local, axis=1)
    sorted_f = np.take_along_axis(local, order, axis=1)
    weights = np.column_stack((1-sorted_f[:, 0], -np.diff(sorted_f, axis=1), sorted_f[:, 2]))
    rows = np.arange(len(xyz))
    index = origin.copy()
    index[:, 0] += 1
    result = weights[:, :1]*grid[index[:, 2], index[:, 1], index[:, 0]]
    for i in range(3):
        axis = order[:, i]
        index[rows, axis] += np.where(axis == 0, -1, 1)
        result += weights[:, i+1:i+2]*grid[index[:, 2], index[:, 1], index[:, 0]]
    return result


def sample_points(index, n=100000):
    seed = np.random.SeedSequence(2026, spawn_key=(index,))
    points = np.random.default_rng(seed).uniform(LO, HI, size=(n, 3))
    if index:
        face = index-1
        points[:, face//2] = (HI if face%2 else LO)[face//2]
    return points


def build_reference_samples(folder):
    field = WakeField()
    grid = np.load(folder / "grid_values.npy", mmap_mode="r")
    # Check the tetrahedral interpolation against TTK's exported corners.
    from vtk.util.numpy_support import vtk_to_numpy
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(folder / "grid0025/fibers.vtp"))
    reader.Update()
    curves = reader.GetOutput()
    points = vtk_to_numpy(curves.GetPoints().GetData()).astype(float)
    lines = vtk_to_numpy(curves.GetLines().GetConnectivityArray()).reshape(-1, 2)
    q = vtk_to_numpy(curves.GetCellData().GetScalars())
    difference = linear_values(grid, points)[lines]-field.query[q][:, None, :]
    err = float(np.linalg.norm(difference, axis=2).max())
    if err > 1e-4:
        raise ValueError(f"TTK corner and interpolated field disagree by {err}")
    report = dict(seed=2026, samples_per_region=100000, max_pl_residual_at_ttk_corners=err, regions={})
    labels = {}
    for i, name in enumerate(("domain", *FACE_NAMES)):
        xyz = sample_points(i)
        target = field.selected(field(xyz))
        labels[name] = target
        n = len(target)
        count = int(target.sum())
        measure = float(np.prod(HI-LO) if not i else np.prod(np.delete(HI-LO, (i-1)//2)))
        region = dict(in_target=count, target_fraction=count/n, target_measure=measure*count/n,
                      target_measure_wilson95=[measure*v for v in wilson(count, n)])
        for tag, g in (("pl_grid005", grid[::2, ::2, ::2]), ("pl_grid0025", grid)):
            region[tag] = agreement(target, field.selected(linear_values(g, xyz)))
        report["regions"][name] = region
    np.savez_compressed(folder / "sample_labels.npz", **labels)
    (folder / "membership.json").write_text(json.dumps(report, indent=2)+"\n")


def cap_membership(xyz, triangles, axis):
    in_plane = np.arange(3) != axis
    xy = xyz[:, in_plane]
    hit = np.zeros(len(xy), dtype=bool)
    for tri in triangles[:, :, in_plane]:
        use = np.flatnonzero(((xy >= tri.min(axis=0)) & (xy <= tri.max(axis=0))).all(axis=1))
        a, b = tri[1:]-tri[0]
        det = a[0]*b[1]-a[1]*b[0]
        if not len(use) or det == 0:
            continue
        d = xy[use]-tri[0]
        u = (d[:, 0]*b[1]-d[:, 1]*b[0])/det
        v = (a[0]*d[:, 1]-a[1]*d[:, 0])/det
        hit[use] |= (u >= 0) & (v >= 0) & (u+v <= 1)
    return hit


def assess_run(folder, reference, dest):
    labels = np.load(reference / "sample_labels.npz")
    raw = json.loads((folder / "volume_report.json").read_text())
    result = dict(complete=raw["complete"], volume_accuracy=None, caps={})
    if raw["complete"]:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(str(folder / "volume_volume.vtk"))
        reader.Update()
        mesh = reader.GetOutput()
        if not mesh.GetNumberOfCells():
            raise ValueError("Empty completed volume")
        locator = vtk.vtkStaticCellLocator()
        locator.SetDataSet(mesh)
        locator.BuildLocator()
        hit = np.array([locator.FindCell(p) >= 0 for p in sample_points(0)])
        result["volume_accuracy"] = agreement(labels["domain"], hit)
    surface = "volume_tet_boundary.vtk" if raw["complete"] else "volume_boundary.vtk"
    pts, faces, patch = read_mesh(folder / surface)
    for i, name in enumerate(FACE_NAMES):
        cap = pts[faces[patch == 4+i]]
        hit = cap_membership(sample_points(i+1), cap, i//2)
        result["caps"][name] = dict(triangles=len(cap), **agreement(labels[name], hit))
    (dest / "membership.json").write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=WAKE/"ttk")
    parser.add_argument("--run", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    torch.set_num_threads(2)
    if args.run is None:
        build_reference_samples(args.reference)
    else:
        assess_run(args.run, args.reference, args.output or args.run)
