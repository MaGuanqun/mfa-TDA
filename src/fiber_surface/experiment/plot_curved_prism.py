"""Render the curved-prism boundary, volume cutaway and four-panel teaser."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection, PolyCollection
from matplotlib.path import Path as MplPath
from matplotlib.text import Text
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from evaluate_curved_prism import CORNERS, DEFAULT_DIR, read_mesh


PATCH_COLORS = ["#74b7e5", "#e68aa5", "#bdcf70"]  # AB sky blue, BC rose, CA yellow-green
FIBER_COLORS = ["#12605b", "#9b420c", "#58429a"]
INK = "#263944"


def project(xyz):
    # Orthographic camera. x, y and z have the same spatial scale.
    az, el = np.deg2rad([-65, 24])
    right = [-np.sin(az), np.cos(az), 0]
    up = [-np.cos(az)*np.sin(el), -np.sin(az)*np.sin(el), np.cos(el)]
    return np.asarray(xyz) @ np.array([right, up, np.cross(right, up)]).T


def spatial_axes(ax, extent, arrow_length):
    ax.set_aspect("equal")
    ax.set_axis_off()
    bbox = project(extent)
    xmin, ymin = bbox[:, :2].min(axis=0)
    xmax, ymax = bbox[:, :2].max(axis=0)
    # Reserve more room on the left only when the axis indicator is visible.
    left_pad, right_pad = (.42, .16) if arrow_length else (.29, .29)
    ax.set_xlim(xmin-left_pad, xmax+right_pad)
    ax.set_ylim(ymin-.19, ymax+.24)
    if arrow_length == 0:
        return
    origin = np.array([xmin-.30, ymin+.10])
    for label, vec in zip("xyz", project(np.eye(3))):
        end = origin + arrow_length*vec[:2]
        ax.annotate("", xy=end, xytext=origin, arrowprops=dict(
            arrowstyle="-|>", lw=.75, color=INK, mutation_scale=6, shrinkA=0, shrinkB=0))
        pos = end + .065*vec[:2]
        ax.text(*pos, label, ha="center", va="center", fontsize=8, color=INK)


def draw_faces(ax, polygons, colors, edgecolors, linewidths):
    polygons = [project(p) for p in polygons]
    order = np.argsort([p[:, 2].mean() for p in polygons], kind="stable")
    paths = []
    for i in order:
        xy = polygons[i][:, :2]
        if len(xy) == 2:
            # A round join at the lower vertex meets the preceding segment;
            # the open path ends flush at the upper vertex, including the tip.
            paths.append(MplPath(xy[[1, 0, 1]]))
        else:
            paths.append(MplPath(np.vstack((xy, xy[0])), closed=True))
    ax.add_collection(PathCollection(paths,
                      facecolors=[colors[i] for i in order],
                      edgecolors=[edgecolors[i] for i in order],
                      linewidths=[linewidths[i] for i in order], joinstyle="round", capstyle="butt"))


def read_fibers(folder):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(folder / "corner_fibers.vtp"))
    reader.Update()
    mesh = reader.GetOutput()
    coords = vtk_to_numpy(mesh.GetPoints().GetData())
    curves = []
    for i in range(mesh.GetNumberOfCells()):
        ids = mesh.GetCell(i).GetPointIds()
        curves.append(coords[[ids.GetId(j) for j in range(ids.GetNumberOfIds())]])
    return curves


def save(fig, folder, name):
    for ext in ("pdf", "png"):
        fig.savefig(folder / f"{name}.{ext}", dpi=245, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)


def panel_caption(ax, text, fontsize):
    box = ax.get_position(original=True)
    ax.figure.text(box.x0 + box.width/2, .07, text, ha="center", va="baseline", fontsize=fontsize)


def draw_fiber_surface(ax, points, faces, labels, fibers, *, arrow_length=.17):
    spatial_axes(ax, np.concatenate(fibers), arrow_length)
    sides = labels <= 3
    polys = list(points[faces[sides]])
    colors = [matplotlib.colors.to_rgba(PATCH_COLORS[v-1], .70) for v in labels[sides]]
    edgecolors = [(0.14, .24, .31, .19)]*len(polys)
    linewidths = [.27]*len(polys)
    # Sort individual fiber segments with the mesh so nearer faces tint the fibers behind them.
    for curve, color in zip(fibers, FIBER_COLORS):
        segments = np.stack((curve[:-1], curve[1:]), axis=1)
        polys.extend(segments)
        colors.extend(["none"]*len(segments))
        edgecolors.extend([color]*len(segments))
        linewidths.extend([1.8]*len(segments))
    draw_faces(ax, polys, colors, edgecolors, linewidths)
    for i, curve in enumerate(fibers):
        xy = project(curve[-1])[:2]
        ax.annotate("ABC"[i], xy, xytext=((-4, 4), (0, 5), (0, 5))[i],
                    textcoords="offset points", color=FIBER_COLORS[i],
                    ha="center", weight="bold", fontsize=10,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=.4))


def draw_range_triangle(key):
    key.set_aspect("equal")
    key.set_xlim(-.2, 1.37)
    key.set_ylim(-.3, 1.38)
    key.set_axis_off()
    key.fill(*CORNERS.T, color="#f4f4f1")
    for i, (a, b) in enumerate(zip(CORNERS, np.roll(CORNERS, -1, axis=0))):
        key.plot(*np.array([a, b]).T, color=PATCH_COLORS[i], lw=3)
    for i, (x, y) in enumerate(CORNERS):
        key.scatter(x, y, s=28, color=FIBER_COLORS[i], zorder=3)
        key.annotate(f"{'ABC'[i]}  ({x:.0f}, {y:.0f})", (x, y),
                     xytext=((-4, -16), (-8, -16), (6, 4))[i], textcoords="offset points",
                     color=FIBER_COLORS[i], fontsize=9, weight="medium")


def cross_section(folder):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(folder / "volume_volume.vtk"))
    plane = vtk.vtkPlane()
    plane.SetNormal(0, 0, 1)
    plane.SetOrigin(0, 0, 0)
    cut = vtk.vtkCutter()
    cut.SetInputConnection(reader.GetOutputPort())
    cut.SetCutFunction(plane)
    cut.GenerateTrianglesOff()  # preserve quadrilateral sections of single cells
    cut.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    cut.Update()
    mesh = cut.GetOutput()
    points = vtk_to_numpy(mesh.GetPoints().GetData())
    cells = mesh.GetPolys()
    offsets = vtk_to_numpy(cells.GetOffsetsArray())
    ids = vtk_to_numpy(cells.GetConnectivityArray())
    return [points[ids[a:b]] for a, b in zip(offsets, offsets[1:])]


def lower_half(triangle):
    result = []
    for a, b in zip(triangle, np.roll(triangle, -1, axis=0)):
        if a[2] <= 0:
            result.append(a)
        if a[2]*b[2] < 0:
            result.append(a + (b-a)*(-a[2]/(b[2]-a[2])))
    return np.asarray(result)


def draw_cutaway(ax, points, faces, labels, fibers, sections, *, arrow_length=.17):
    spatial_axes(ax, np.concatenate(fibers), arrow_length)
    # Outline the removed half with the measured fiber coordinates.
    for curve, color in zip(fibers, FIBER_COLORS):
        at_cut = [np.interp(0., curve[:, 2], coord) for coord in curve.T]
        xy = project(np.vstack((at_cut, curve[curve[:, 2] > 0])))
        ax.plot(*xy[:, :2].T, color=color, alpha=.55, ls=(0, (2.5, 2)), lw=.8, zorder=0)
    upper = project([fibers[i][-1] for i in (0, 1, 2, 0)])
    for i, color in enumerate(PATCH_COLORS):
        ax.plot(*upper[i:i+2, :2].T, color=color, ls=(0, (2.5, 2)), lw=.8, zorder=0)
    polygons, colors = [], []
    # The cutaway is opaque. Cull outward normals pointing away from the
    # camera so hidden triangles of the bottom cap cannot overpaint the wall.
    a, b, c = points[faces].transpose(1, 0, 2)
    front = project(np.cross(b-a, c-a))[:, 2] > 0
    for face, label in zip(faces[front], labels[front]):
        poly = lower_half(points[face])
        if len(poly) >= 3:
            polygons.append(poly)
            colors.append("#709aa8" if label == 8 else PATCH_COLORS[label-1])
    n = len(polygons)
    draw_faces(ax, polygons + sections, colors + ["#f6ce91"]*len(sections),
               ["#708c99"]*n + ["#97642d"]*len(sections), [.24]*n + [.32]*len(sections))
    cut_center = project([*CORNERS.mean(axis=0), 0])[:2]
    ax.annotate("$z=0$", xy=cut_center, xytext=(-36, 1), textcoords="offset points",
                ha="right", va="center", color=INK,
                arrowprops=dict(arrowstyle="->", lw=.8, color=INK))


def draw_section(top, sections, *, show_exact=True):
    top.add_collection(PolyCollection([p[:, :2] for p in sections], facecolors="#f6ce91",
                                      edgecolors="#97642d", linewidths=.4, joinstyle="round"))
    if show_exact:
        outline = CORNERS[[0, 1, 2, 0]]
        top.plot(*outline.T, ls=(0, (3, 2)), lw=1.1, color=INK, label="Exact section", clip_on=False)
        top.legend(loc="upper right", frameon=False, fontsize=8, handlelength=1.4,
                   handletextpad=.5, borderaxespad=.2)
    top.set_aspect("equal")
    top.set_xlim(-.05, 1.06)
    top.set_ylim(-.05, 1.06)
    top.set_xticks([0, .5, 1])
    top.set_yticks([0, .5, 1])
    top.set_xlabel("$x$", labelpad=0)
    top.set_ylabel("$y$", rotation=0, labelpad=5)
    for spine in ("top", "right"):
        top.spines[spine].set_visible(False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, nargs="?", default=DEFAULT_DIR)
    folder = parser.parse_args().folder
    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "axes.edgecolor": "#8b989d", "xtick.color": INK, "ytick.color": INK,
                         "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    points, faces, labels = read_mesh(folder / "volume_tet_boundary.vtk")
    fibers = read_fibers(folder)
    sections = cross_section(folder)

    fig, (surface, key) = plt.subplots(1, 2, figsize=(6.3, 3.65), gridspec_kw={"width_ratios": [1.4, 1]})
    fig.subplots_adjust(left=.015, right=.97, top=.98, bottom=.10, wspace=.025)
    draw_fiber_surface(surface, points, faces, labels, fibers)
    draw_range_triangle(key)
    panel_caption(key, "Range triangle  $(u,v)$", fontsize=11)
    save(fig, folder, "boundary_fibers")

    fig, (volume, section) = plt.subplots(1, 2, figsize=(6.55, 3.7), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(left=.01, right=.98, top=.96, bottom=.20, wspace=.16)
    draw_cutaway(volume, points, faces, labels, fibers, sections)
    draw_section(section, sections)
    panel_caption(section, "Tetrahedral section at $z=0$", fontsize=10)
    save(fig, folder, "volume_cutaway")

    fig, axes = plt.subplots(1, 4, figsize=(6.15, 2.25))
    fig.subplots_adjust(left=.015, right=.985, bottom=.10, top=.916, wspace=.015)
    key, surface, volume, section = axes
    # Fixed aspect ratios set the plot widths from these heights.
    for ax, scale, anchor_y in ((key, 2/3, .22), (section, 1/2, .5)):
        box = ax.get_position(original=True)
        ax.set_position(box.shrunk(1, scale).anchored((.5, anchor_y), box))
    draw_range_triangle(key)
    draw_fiber_surface(surface, points, faces, labels, fibers, arrow_length=.23)
    draw_cutaway(volume, points, faces, labels, fibers, sections, arrow_length=0)
    draw_section(section, sections, show_exact=False)
    headings = ("(a) Range triangle $(u,v)$", "(b) Fiber surface", "(c) Volume cutaway", "(d) Section at $z=0$")
    for ax, heading in zip(axes, headings):
        panel_caption(ax, heading, fontsize=10)
    for label in fig.findobj(Text):
        label.set_fontsize(label.get_fontsize() * 2/3)
    save(fig, folder, "curved-prism-teaser")


if __name__ == "__main__":
    main()
