"""Reduce the native TTK surface while retaining the three range-edge labels."""

import argparse
import ctypes
import json
from pathlib import Path
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from evaluate_curved_prism import save_polydata
from gray_scott_bubbles import DATA, write_json
from ttk_geometry import directed_distances


def read_surface(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    # meshoptimizer expects consistent winding across neighboring triangles.
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(reader.GetOutputPort())
    normals.ConsistencyOn()
    normals.SplittingOff()
    normals.NonManifoldTraversalOff()
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.Update()
    poly = normals.GetOutput()
    faces = vtk_to_numpy(poly.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    labels = vtk_to_numpy(poly.GetCellData().GetArray("patch"))
    if len(faces) != reader.GetOutput().GetNumberOfPolys() or not np.array_equal(np.unique(labels), [1, 2, 3]):
        raise ValueError("Expected all three TTK side patches")
    return vtk_to_numpy(poly.GetPoints().GetData()), faces, labels


def simplify(xyz, faces, labels, lo, hi, budget):
    # Duplicate only the vertex records shared by different patches. The
    # geometric positions coincide so meshoptimizer can collapse both sides of
    # a seam together while keeping their index buffers distinct.
    keys, inverse = np.unique(faces*3 + (labels[:, None]-1), return_inverse=True)
    vertex_ids = keys // 3
    patch_ids = keys % 3 + 1
    indices = np.ascontiguousarray(inverse.ravel(), dtype=np.uint32)
    positions = np.ascontiguousarray((xyz[vertex_ids]-lo)/(hi-lo), dtype=np.float32)
    destination = np.empty_like(indices)
    library = ctypes.CDLL(str(Path(sys.prefix) / "lib/libmeshoptimizer.so"))
    reduce_mesh = library.meshopt_simplifyWithAttributes
    ui = np.ctypeslib.ndpointer(np.uint32, ndim=1, flags="C_CONTIGUOUS")
    pts = np.ctypeslib.ndpointer(np.float32, ndim=2, flags="C_CONTIGUOUS")
    size = ctypes.c_size_t
    locks_type = np.ctypeslib.ndpointer(np.uint8, ndim=1, flags="C_CONTIGUOUS")
    reduce_mesh.argtypes = [ui, ui, size, pts, size, size, ctypes.c_void_p, size,
                            ctypes.c_void_p, size, locks_type, size, ctypes.c_float,
                            ctypes.c_uint, ctypes.POINTER(ctypes.c_float)]
    reduce_mesh.restype = size
    err = ctypes.c_float()
    locked = np.zeros(len(xyz), dtype=bool)
    for attempt in range(12):
        n = reduce_mesh(destination, indices, len(indices), positions, len(positions), positions.strides[0],
                        None, 0, None, 0, locked[vertex_ids].astype(np.uint8),
                        3*budget, 1.0, 0, ctypes.byref(err))
        wedges = destination[:n].reshape(-1, 3)
        patch = patch_ids[wedges]
        if np.any(patch != patch[:, :1]):
            raise ValueError("Reduction crossed a range-edge label")
        reduced = vertex_ids[wedges]
        triangles, multiplicity = np.unique(np.sort(reduced, axis=1), axis=0, return_counts=True)
        edge_keys = np.concatenate([np.minimum(reduced[:, i], reduced[:, j])*len(xyz)
                                    + np.maximum(reduced[:, i], reduced[:, j])
                                    for i, j in ((0, 1), (1, 2), (2, 0))])
        keys, count = np.unique(edge_keys, return_counts=True)
        bad = keys[count > 2]
        defects = np.unique(np.r_[triangles[multiplicity > 1].ravel(), bad//len(xyz), bad%len(xyz)])
        if not len(defects):
            return reduced, patch[:, 0], err.value, int(locked.sum())
        # Thin components can fold into coincident faces during aggressive
        # reduction. Retain their original one-ring and simplify again.
        affected = np.zeros(len(xyz), dtype=bool)
        affected[defects] = True
        locked[faces[affected[faces].any(axis=1)]] = True
    raise RuntimeError("Could not reduce this surface without duplicate faces or nonmanifold edges")


def inspect_edges(xyz, faces, labels, lo, hi):
    width = int(faces.max())+1
    e = np.concatenate([np.sort(faces[:, pair], axis=1) for pair in ((0, 1), (1, 2), (2, 0))])
    keys = e[:, 0]*width+e[:, 1]
    order = np.argsort(keys)
    unique, first, count = np.unique(keys[order], return_index=True, return_counts=True)
    start = first[count == 2]
    sides = np.tile(labels, 3)[order]
    a, b = sides[start], sides[start+1]
    seam = unique[count == 2][a != b]
    # For 1-based side IDs, the common corner of (1,2), (2,3), (1,3) is 1,2,0.
    corners = (2*(a[a != b]+b[a != b])+1) % 3
    edges = np.column_stack((seam//width, seam%width))
    exposed = unique[count == 1]
    endpoints = xyz[np.column_stack((exposed//width, exposed%width))]
    tol = 1e-6*(hi-lo)
    on_domain = ((np.abs(endpoints-lo) <= tol).all(axis=1) |
                 (np.abs(endpoints-hi) <= tol).all(axis=1)).any(axis=1)
    return edges, corners, dict(open_edges=len(exposed), open_on_domain=int(on_domain.sum()),
                                open_interior=int((~on_domain).sum()),
                                nonmanifold_edges=int(np.count_nonzero(count > 2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--budget", type=int, default=80000)
    args = parser.parse_args()
    meta = json.loads((DATA / "gray_scott_bubbles.json").read_text())
    lo, hi = np.array(meta["origin"]), np.array(meta["upper"])
    xyz, faces, labels = read_surface(args.native / "fiber_surface.vtp")
    reduced, patch, error, locked = simplify(xyz, faces, labels, lo, hi, args.budget)
    _, _, native_incidence = inspect_edges(xyz, faces, labels, lo, hi)
    edges, corners, incidence = inspect_edges(xyz, reduced, patch, lo, hi)
    print(f"{len(faces):,} -> {len(reduced):,} triangles", flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    save_polydata(args.output / "fiber_surface.vtp", xyz, reduced, patch)
    if len(edges):
        save_polydata(args.output / "corner_fibers.vtp", xyz, edges, corners, lines=True)
    report = dict(native_triangles=len(faces), target=args.budget, triangles=len(reduced),
                  triangles_per_side=np.bincount(patch, minlength=4)[1:].tolist(),
                  corner_segments=np.bincount(corners, minlength=3).tolist(), incidence=incidence,
                  native_incidence=native_incidence,
                  locked_native_vertices=locked, meshoptimizer_relative_error=error,
                  geometry_distance_units="each domain axis scaled to [0,1]",
                  sampled_distances_per_side={})
    for i in (1, 2, 3):
        a = xyz[faces[labels == i]]
        b = xyz[reduced[patch == i]]
        report["sampled_distances_per_side"][str(i)] = dict(
            native_to_reduced=directed_distances(a, b, 42+i, lo, hi),
            reduced_to_native=directed_distances(b, a, 84+i, lo, hi))
    write_json(args.output / "surface.json", report)


if __name__ == "__main__":
    main()
