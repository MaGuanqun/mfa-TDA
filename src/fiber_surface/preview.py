#!/usr/bin/env python3
"""Render an exported OBJ and its optional endpoint-fiber VTK (NumPy/Matplotlib)."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--title", default="Continuous-model fiber surface")
    args = parser.parse_args()
    vertices, triangles = [], []
    for line in args.mesh.read_text().splitlines():
        fields = line.split()
        if fields and fields[0] == "v":
            vertices.append([float(x) for x in fields[1:4]])
        elif fields and fields[0] == "f":
            triangles.append([int(x.split("/")[0])-1 for x in fields[1:4]])
    vertices, triangles = np.asarray(vertices), np.asarray(triangles, dtype=int)
    if not len(triangles):
        raise ValueError("The mesh has no triangles")
    fibers = []
    fiber_path = args.mesh.with_name(args.mesh.stem + "_fibers.vtk")
    if fiber_path.exists():
        tokens = fiber_path.read_text().split()
        start = tokens.index("LINES")
        count = int(tokens[start+1])
        for i in range(count):
            offset = start+3+3*i
            fibers.append([int(tokens[offset+1]), int(tokens[offset+2])])
    fig = plt.figure(figsize=(11, 5), constrained_layout=True)
    for i, (elevation, azimuth) in enumerate([(25, -55), (65, -55)], 1):
        ax = fig.add_subplot(1, 2, i, projection="3d")
        ax.add_collection3d(Poly3DCollection(vertices[triangles], facecolor="#79bdd1",
                                           edgecolor="#315c68", linewidth=.23, alpha=1))
        for edge in fibers:
            ax.plot(*vertices[edge].T, color="#e4572e", linewidth=1.8)
        low, high = vertices.min(axis=0), vertices.max(axis=0)
        extent = np.maximum(high-low, .03*max(high-low))
        ax.set(xlim=(low[0]-.03*extent[0], high[0]+.03*extent[0]),
               ylim=(low[1]-.03*extent[1], high[1]+.03*extent[1]),
               zlim=(low[2]-.03*extent[2], high[2]+.03*extent[2]),
               xlabel="x", ylabel="y", zlabel="z")
        ax.set_box_aspect(extent)
        ax.view_init(elev=elevation, azim=azimuth)
        ax.tick_params(labelsize=8)
    fig.suptitle(f"{args.title}\n{len(vertices):,} vertices · {len(triangles):,} triangles · orange: endpoint fibers")
    output = args.output or args.mesh.with_suffix(".png")
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
