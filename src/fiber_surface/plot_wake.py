"""Plot the directly extracted wake boundary alongside its PL reference and show fitted-field slices."""

import argparse
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
import numpy as np
import torch

from evaluate_wake_geometry import read_reference, original_curves
from experiment.evaluate_curved_prism import read_mesh
from wake_reference import HI, LO, WAKE, WakeField


PATCH = ["#44a8a5", "#eba05b", "#a896cb"]
FIBER = ["#065955", "#99430a", "#58388d"]


def project(points):
    x, y, t = np.moveaxis(points, -1, 0)
    return np.stack((x+.48*y, t+.65*y), axis=-1)


def panel(ax, sides, curves, title):
    for index in product(range(2), repeat=3):
        start = np.where(index, HI, LO)
        for axis in range(3):
            if index[axis]:
                continue
            end = start.copy()
            end[axis] = HI[axis]
            ax.plot(*project(np.array([start, end])).T, color="#bec8cb", lw=.65, zorder=0)
    faces = np.concatenate(sides)
    order = np.argsort(faces[:, :, 1].mean(axis=1))
    colors = np.vstack([np.broadcast_to(to_rgb(c), (len(f), 3)) for c, f in zip(PATCH, sides)])
    ax.add_collection(PolyCollection(project(faces[order]), facecolors=colors[order], linewidths=0,
                                     rasterized=True, alpha=.85))
    for lines, color in zip(curves, FIBER):
        ax.add_collection(LineCollection(project(lines), colors=color, linewidths=.46, rasterized=True))
    ax.set_title(title, loc="left", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-.88, 7.88)
    ax.set_ylim(13.08, 15.42)
    ax.axis("off")
    ax.annotate("$x$", project(np.array([3.5,-.5,13.5])), xytext=(0,-8), textcoords="offset points")
    ax.annotate("$t$", project(np.array([-.5,-.5,14.25])), xytext=(-9,0), textcoords="offset points")
    ax.annotate("$y$", project(np.array([7.5,0,13.5])), xytext=(4,-4), textcoords="offset points")


def save(fig, path):
    for fmt in ("pdf", "png"):
        fig.savefig(path.with_suffix("."+fmt), dpi=210, bbox_inches="tight", pad_inches=.1)
    plt.close(fig)


def plot_geometry(run, trace, dest, reference, h, complete):
    src = "volume_tet_boundary.vtk" if complete else "volume_boundary.vtk"
    p, faces, labels = read_mesh(run/src)
    curves = original_curves(trace)
    parts = [np.concatenate([np.stack((xyz[:-1], xyz[1:]), axis=1) for xyz, _, _ in group]) for group in curves]
    fig, axes = plt.subplots(2, 1, figsize=(8.6,5.2), gridspec_kw={"hspace":.21})
    status = "closed" if complete else "partial"
    panel(axes[0], [p[faces[labels == j]] for j in (1,2,3)], parts,
          f"Continuous INR  ·  {status} boundary,  $h={h:g}$")
    panel(axes[1], read_reference(reference,False), read_reference(reference,True),
          "TTK  ·  piecewise-linear INR on the full 321 × 41 × 61 grid")
    legend = [Line2D([],[], color=c, lw=2, label=label) for c,label in zip(FIBER,("A  (0.60, 0)", "B  (0.90, 0)", "C  (0.75, 0.35)"))]
    fig.legend(handles=legend, title="Corner fibers in velocity coordinates", loc="lower center",
               bbox_to_anchor=(.5,-.025), ncol=3, frameon=False, fontsize=8, title_fontsize=8)
    save(fig, dest/"wake_boundary")


def plot_sections(dest):
    field = WakeField()
    x, y = np.meshgrid(np.linspace(LO[0],HI[0],801), np.linspace(LO[1],HI[1],101))
    fig, axes = plt.subplots(3,1, figsize=(8.5,4.1), sharex=True, layout="constrained")
    for ax, t in zip(axes,(13.8,14.1,14.4)):
        xyz = np.column_stack((x.ravel(), y.ravel(), np.full(x.size,t)))
        uv = field(xyz).reshape(*x.shape,2)
        tests = []
        for a,b in zip(field.query, np.roll(field.query,-1,axis=0)):
            d = b-a
            tests.append((uv-a) @ np.array([-d[1],d[0]]) / np.linalg.norm(d))
        margin = np.min(tests,axis=0)
        ax.contourf(x,y,margin,levels=[0,1],colors=["#2c918e"])
        ax.set_aspect("equal")
        ax.set_yticks([-.5,0,.5])
        ax.set_ylabel("$y$")
        ax.set_title(f"$t={t}$",loc="left",fontsize=9)
    axes[-1].set_xlabel("$x$")
    fig.suptitle("Upward-flow region evaluated directly on the fitted INR",fontsize=11)
    save(fig,dest/"wake_target_sections")


if __name__=="__main__":
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--trace", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    torch.set_num_threads(1)
    dest = args.output or args.run
    report = json.loads((args.run/"volume_report.json").read_text())
    plot_geometry(args.run, args.trace or args.run/"trace", dest, WAKE/"ttk/grid0025", report["step"], report["complete"])
    plot_sections(dest)
