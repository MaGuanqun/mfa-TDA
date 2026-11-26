"""Extract the paper's native-merge TTK baseline from sampled CoordNet grid nodes."""
import argparse
import json
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from join_isabel_refinements import write_ply
from sample_coordnet_grid import read_grid
from ttk_geometry import ttk_triangle_mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid", type=Path)
    parser.add_argument("output", type=Path, help="labeled binary boundary PLY")
    parser.add_argument("--threads", type=int, default=3)
    args = parser.parse_args()
    grid, fields = read_grid(args.grid)
    image = vtk.vtkImageData()
    image.SetDimensions(*grid["dimensions"])
    image.SetOrigin(*grid["box_min"])
    image.SetSpacing(*grid["spacing"])
    names = ("first", "second")
    for name, values in zip(names, fields):
        array = numpy_to_vtk(values.ravel())
        array.SetName(name)
        image.GetPointData().AddArray(array)
    query = np.array(grid["triangle"])
    merge_distance = 1e-6
    ttk = ttk_triangle_mesh(image, names, query, args.threads, range_octree=False, merge_distance=merge_distance)
    n = ttk.GetNumberOfPolys()
    if n:
        xyz = vtk_to_numpy(ttk.GetPoints().GetData())
        faces = vtk_to_numpy(ttk.GetPolys().GetConnectivityArray()).reshape(n, 3)
        labels = vtk_to_numpy(ttk.GetCellData().GetArray("EdgeIds")) + 1
        if len(labels) != n or not np.isin(labels, [1, 2, 3]).all():
            raise ValueError("Expected TTK cells labeled AB, BC and CA by EdgeIds")
    else:
        xyz, faces, labels = np.empty((0, 3)), np.empty((0, 3), dtype=int), np.empty(0, dtype=int)
    reference = args.output
    write_ply(reference, xyz, faces, labels, (0, 0, 0))
    meta = dict(ttk_native_point_merge=merge_distance, ttk_threads=args.threads,
                ttk_vtk_version=vtk.vtkVersion.GetVTKVersion(), ttk_triangles=n)
    reference.with_name(reference.stem + "_extraction.json").write_text(json.dumps(meta, indent=2) + "\n")


if __name__ == "__main__":
    main()
