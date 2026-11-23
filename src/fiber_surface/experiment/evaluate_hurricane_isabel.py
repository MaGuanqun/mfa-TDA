"""Check the seed-40 Isabel mesh against the time-30 grid and plot its boundary."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter
from vtkmodules.vtkIOLegacy import vtkPolyDataReader

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_isabel_grid import grid_selection, measure_grid


def read_boundary(folder, report):
    reader = vtkPolyDataReader()
    reader.SetFileName(str(folder / "volume_boundary.vtk"))
    reader.Update()
    mesh = reader.GetOutput()
    if mesh.GetNumberOfPolys() != report["boundary_triangles"]:
        raise ValueError("Boundary face count differs from the extraction report")
    return mesh


def plot_boundary(folder, report, mesh=None):
    if mesh is None:
        mesh = read_boundary(folder, report)
    points = vtk_to_numpy(mesh.GetPoints().GetData())
    faces = vtk_to_numpy(mesh.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    patches = vtk_to_numpy(mesh.GetCellData().GetArray("patch"))
    colors = ["#418bb7", "#d67742", "#568954"]
    corners = np.asarray(report["triangle"])

    fig = plt.figure(figsize=(12.0, 4.0))
    panels = fig.add_gridspec(1, 3, width_ratios=[.86, 1, 1], left=.07, right=.97,
                              bottom=.16, top=.89, wspace=.20)
    ax = fig.add_subplot(panels[0])
    ax.fill(*corners.T, color="#f3f3ed")
    for i, (p, q) in enumerate(zip(corners, np.roll(corners, -1, axis=0))):
        ax.plot(*np.vstack([p, q]).T, color=colors[i], linewidth=3)
        right, top = p >= corners.mean(axis=0)
        ax.annotate("ABC"[i], p, xytext=(5 if right else -5, 3 if top else -3),
                    ha="left" if right else "right", va="bottom" if top else "top",
                    textcoords="offset points", weight="bold")
    ax.set(xlabel="Pressure", ylabel="Temperature", title="Query in model range")
    ax.locator_params(axis="x", nbins=5)
    ax.locator_params(axis="y", nbins=6)
    ax.margins(.18)

    face_colors = [colors[p - 1] if 1 <= p <= 3 else "#a8aeb1" for p in patches]
    # Render the same opaque patches from opposite elevations.
    z_tick = round(points[:, 2].mean())
    for panel, elev, title in [(1, 35, "Boundary from above"), (2, -35, "Boundary from below")]:
        ax = fig.add_subplot(panels[panel], projection="3d", proj_type="ortho")
        ax.add_collection3d(Poly3DCollection(points[faces], facecolors=face_colors,
                                             edgecolors="#36444c", linewidths=.15))
        ax.auto_scale_xyz(*points.T)
        ax.set_box_aspect(np.ptp(points, axis=0))
        ax.locator_params(axis="x", nbins=4)
        ax.locator_params(axis="y", nbins=4)
        if np.ptp(points[:, 2]) >= 20:
            ax.locator_params(axis="z", nbins=3)
            ax.zaxis.set_major_formatter("z={x:g}")
        else:
            ax.set_zticks([z_tick], labels=[f"z={z_tick}"])
        ax.tick_params(labelsize=8, pad=1)
        ax.view_init(elev=elev, azim=-65)
        ax.set(xlabel="x", ylabel="y", title=title)
    if not report["boundary_valid"]:
        fig.suptitle("MFA boundary · validation failed")
    for ext in ("png", "pdf"):
        fig.savefig(folder / f"isabel_triangle.{ext}", dpi=200, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


def plot_footprint(folder, stats, footprint, fields, triangle, boundary=None, valid=False):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.0), layout="constrained",
                             sharex=True, sharey=True)
    ax = axes[0]
    image = ax.imshow(np.ma.masked_equal(footprint, 0), origin="lower", cmap="viridis",
                      interpolation="nearest")
    fig.colorbar(image, ax=ax, orientation="horizontal", shrink=.84, pad=.12,
                 label="Selected vertical grid levels")
    ax.set_title("Footprint across all levels")
    # Use the level with the largest selection for a reproducible horizontal section.
    z = int(np.argmax(stats["selected_by_z"]))
    p, t = (field[z] for field in fields)
    selected = grid_selection(p, t, triangle)
    segments = []
    if boundary is not None:
        plane = vtkPlane()
        plane.SetOrigin(0, 0, z)
        plane.SetNormal(0, 0, 1)
        cut = vtkCutter()
        cut.SetCutFunction(plane)
        cut.SetInputData(boundary)
        cut.Update()
        lines = cut.GetOutput()
        if lines.GetNumberOfLines():
            xy = vtk_to_numpy(lines.GetPoints().GetData())[:, :2]
            cells = lines.GetLines()
            offsets = vtk_to_numpy(cells.GetOffsetsArray())
            ids = vtk_to_numpy(cells.GetConnectivityArray())
            segments = [xy[ids[a:b]] for a, b in zip(offsets, offsets[1:])]

    for section, values, name in zip(axes[1:], (p, t), ("Pressure", "Temperature")):
        background = section.imshow(values, origin="lower", cmap="Greys_r", interpolation="nearest")
        fig.colorbar(background, ax=section, orientation="horizontal", shrink=.84, pad=.12,
                     label=name)
        if selected.any() and not selected.all():
            section.contour(selected, levels=[.5], colors=["#d34b13"], linewidths=1.4)
            section.plot([], [], color="#d34b13", label="Source-grid selection")
        if segments:
            section.add_collection(LineCollection(segments, colors="#0085a6", linewidths=.8,
                                   label="MFA mesh" if valid else "MFA mesh (unvalidated)"))
        section.set_title(f"{name} at z = {z}")
    if axes[1].get_legend_handles_labels()[0]:
        axes[1].legend(loc="lower left", fontsize=8, framealpha=.95)
    for axis in axes:
        axis.set(xlabel="x (grid units)", aspect="equal")
    ax.set_ylabel("y (grid units)")
    for ext in ("png", "pdf"):
        fig.savefig(folder / f"isabel_footprint.{ext}", dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("source_grid", type=Path)
    args = parser.parse_args()
    report = json.loads((args.folder / "volume_report.json").read_text())
    if (report["discovery_grid_shape"] != [202, 202, 40] or report["starting_points"] != 1632160
            or report["step"] != report["fiber_step"]):
        raise ValueError("Expected the full Isabel seed-40 lattice with equal tracing and marching steps")
    print(f"Corner fiber components at A, B, C: {report['corner_fiber_components']}")
    fields = [np.memmap(args.source_grid / name, dtype="<f4", mode="r", shape=(100, 500, 500))
              for name in ("Pf30_3d.bin", "TCf30_3d.bin")]
    agreement, footprint = measure_grid(args.folder, report, *fields)
    (args.folder / "source_grid.json").write_text(json.dumps(agreement, indent=2) + "\n")
    boundary = read_boundary(args.folder, report) if report["boundary_triangles"] else None
    plot_footprint(args.folder, agreement, footprint, fields, report["triangle"], boundary,
                   report["boundary_valid"])
    if boundary is not None:
        plot_boundary(args.folder, report, boundary)
    if not report["complete"]:
        print(f"Incomplete extraction at {report['failed_stage']}. Source footprint saved.")
        return 2
    print(f"Validated boundary: {report['boundary_triangles']} triangles, "
          f"{report['coverage_misses']} sampled coverage misses.")
    print(f"Validated fill: {report['tetrahedra']} tetrahedra, {report['volume']:.6f} cubic grid units.")
    print(f"Source-grid agreement: {agreement['both']} in both, {agreement['source_only']} source only, "
          f"{agreement['mesh_only']} mesh only. IoU {agreement['iou']:.6f}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
