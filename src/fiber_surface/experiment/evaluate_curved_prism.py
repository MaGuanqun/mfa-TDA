"""Evaluate the manuscript's curved prism from a completed fiber_volume_explicit run.

Run run_curved_prism.sh in this directory to reproduce the full experiment.
To evaluate existing results, pass the directory containing the extractor's
VTK files and report as the positional argument.
"""

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import vtk
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR / "critical_point_tracking"))
from evaluate_analytic_geometry import SegmentGraph


DEFAULT_DIR = SRC_DIR.parent / "build/fiber-volume-results/curved_prism_h010_seed10"
CORNERS = np.array([[0., 0.], [1., 0.], [0., 1.]])
LO = np.array([-2., -2., -1.])
HI = -LO
PROJECTION_EPS = 1e-12  # dimensionless segment parameter s


def field(points):
    return points[:, :2] + points[:, 2:]**2 / 2


def read_mesh(path, volume=False):
    reader = vtk.vtkUnstructuredGridReader() if volume else vtk.vtkPolyDataReader()
    reader.SetFileName(str(path))
    reader.ReadAllScalarsOn()
    reader.Update()
    mesh = reader.GetOutput()
    if not mesh.GetNumberOfCells():
        raise ValueError(f"Empty mesh in {path}")
    points = vtk_to_numpy(mesh.GetPoints().GetData())
    cells = mesh.GetCells() if volume else mesh.GetPolys()
    n = 4 if volume else 3
    if not np.all(np.diff(vtk_to_numpy(cells.GetOffsetsArray())) == n):
        raise ValueError(f"Unexpected cell type in {path}")
    if volume and not np.all(vtk_to_numpy(mesh.GetCellTypesArray()) == vtk.VTK_TETRA):
        raise ValueError(f"Expected tetrahedra in {path}")
    ids = vtk_to_numpy(cells.GetConnectivityArray()).reshape(-1, n)
    labels = None if volume else vtk_to_numpy(mesh.GetCellData().GetArray("patch"))
    return points, ids, labels


def incidence(cells, face_size):
    owners = defaultdict(list)
    for i, cell in enumerate(cells):
        for face in combinations(cell, face_size):
            owners[tuple(sorted(face))].append(i)
    return owners


def components(cells, face_size):
    # Components meet across full edges (surfaces) or full faces (volume).
    links = [(ids[0], other) for ids in incidence(cells, face_size).values() for other in ids[1:]]
    r, c = np.asarray(links, dtype=int).reshape(-1, 2).T
    graph = coo_matrix((np.ones(len(r)), (r, c)), shape=(len(cells), len(cells)))
    return int(connected_components(graph, directed=False, return_labels=False))


def save_polydata(path, points, cells, labels, lines=False, scalar_name=None):
    # Compact by vertex ID to preserve the shared patch seams.
    used = np.unique(np.concatenate(cells))
    xyz = vtk.vtkPoints()
    xyz.SetData(numpy_to_vtk(points[used], deep=True))
    offsets = np.r_[0, np.cumsum([len(cell) for cell in cells])].astype(np.int64)
    ids = np.searchsorted(used, np.concatenate(cells)).astype(np.int64)
    connectivity = vtk.vtkCellArray()
    connectivity.SetData(numpy_to_vtkIdTypeArray(offsets, deep=True),
                         numpy_to_vtkIdTypeArray(ids, deep=True))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(xyz)
    (mesh.SetLines if lines else mesh.SetPolys)(connectivity)
    scalars = numpy_to_vtk(np.asarray(labels, dtype=np.int32), deep=True)
    scalars.SetName(scalar_name or ("corner_id" if lines else "patch"))
    mesh.GetCellData().SetScalars(scalars)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(mesh)
    if not writer.Write():
        raise OSError(f"Could not write {path}")


def fiber_vertex_hausdorff(points, corner):
    recovered = (points - LO) / (HI - LO)
    z = np.linspace(-1, 1, 16385)
    reference = (np.column_stack((corner[0] - z*z/2, corner[1] - z*z/2, z)) - LO) / (HI - LO)
    reference_segments = SegmentGraph(reference[:-1], reference[1:])
    recovered_segments = SegmentGraph(recovered[:-1], recovered[1:])
    recovered_to_reference = float(reference_segments.nearest_many(recovered)[0].max())
    reference_to_recovered = float(recovered_segments.nearest_many(reference)[0].max())
    # The interpolation error is bounded by ||r''(z)|| * dz^2 / 8 in the unit cube.
    dz = z[1] - z[0]
    interpolation_bound = (math.sqrt(2) / 4) * dz**2 / 8
    return dict(distance=max(recovered_to_reference, reference_to_recovered),
                reference_vertices=len(reference), reference_interpolation_bound=interpolation_bound,
                recovered_to_reference=recovered_to_reference,
                reference_to_recovered=reference_to_recovered)


def seam_paths(points, faces, labels):
    owners = incidence(faces, 2)
    curves = []
    for i, corner in enumerate(CORNERS):
        pair = sorted((1 + (i-1) % 3, i+1))
        # Independently identify vertices shared by both patches, then check
        # every consecutive pair on the z-monotone fiber against face incidence.
        ids = np.intersect1d(faces[labels == pair[0]], faces[labels == pair[1]])
        ids = ids[np.argsort(points[ids, 2])]
        if len(ids) < 2 or not np.array_equal(points[ids[[0, -1]], 2], [-1, 1]):
            raise ValueError(f"Incomplete analytical fiber at {corner}")
        expected = {tuple(sorted(edge)) for edge in zip(ids[:-1], ids[1:])}
        actual = {edge for edge, incident in owners.items() if sorted(labels[incident]) == pair}
        if expected != actual:
            raise ValueError(f"Seam mismatch at {corner}: {len(expected ^ actual)} edges")
        curves.append(ids)
    return curves


def evaluate_fibers(folder, final_points, final_faces, final_labels):
    # The pre-fill boundary contains the original traced fibers. Check their
    # vertices and the seam incidence again on the actual tetrahedral exterior.
    points, faces, labels = read_mesh(folder / "volume_boundary.vtk")
    curves = seam_paths(points, faces, labels)
    filled = seam_paths(final_points, final_faces, final_labels)
    results = []
    with (folder / "fiber_vertices.csv").open("w", newline="") as stream:
        out = csv.writer(stream)
        out.writerow(["corner", "sample", "x", "y", "z", "F1", "F2", "range_residual"])
        for i, ids in enumerate(curves):
            xyz = points[ids]
            if not set(map(tuple, xyz)).issubset(map(tuple, final_points[filled[i]])):
                raise ValueError(f"Traced corner {i} vertices missing from filled mesh")
            uv = field(xyz)
            residuals = np.linalg.norm(uv - CORNERS[i], axis=1)
            out.writerows(("ABC"[i], j, *p, *v, r) for j, (p, v, r) in enumerate(zip(xyz, uv, residuals)))
            results.append(dict(corner="ABC"[i], q=CORNERS[i].tolist(), vertices=len(ids), segments=len(ids)-1,
                                components=1, z_span=xyz[[0, -1], 2].tolist(),
                                incident_patches=sorted((1 + (i-1) % 3, i+1)),
                                verified_two_patch_edges=len(filled[i])-1, max_vertex_residual=float(residuals.max()),
                                unit_cube_vertex_hausdorff=fiber_vertex_hausdorff(xyz, CORNERS[i])))
    save_polydata(folder / "corner_fibers.vtp", points, curves, range(3), lines=True)
    return results


def evaluate_patches(folder, points, faces, labels):
    results = []
    with (folder / "patch_samples.csv").open("w", newline="") as stream:
        out = csv.writer(stream)
        out.writerow(["patch", "sample", "x", "y", "z", "line_distance", "projection_s"])
        for i, corner in enumerate(CORNERS):
            patch = faces[labels == i+1]
            unique_edges = np.array(list(incidence(patch, 2)))
            direction = CORNERS[(i+1) % 3] - corner
            samples = {"vertices": points[np.unique(patch)], "edge_midpoints": points[unique_edges].mean(axis=1),
                       "centroids": points[patch].mean(axis=1)}
            stats = {}
            for kind, xyz in samples.items():
                delta = field(xyz) - corner
                s = (delta @ direction) / (direction @ direction)
                error = np.abs(delta @ [-direction[1], direction[0]]) / np.linalg.norm(direction)
                outside = np.maximum(0., np.maximum(-s, s-1))
                stats[kind] = dict(count=len(s), max_line_distance=float(error.max()),
                                   rms_line_distance=float(np.sqrt(np.mean(error**2))),
                                   projection_min=float(s.min()), projection_max=float(s.max()),
                                   outside_strict=int(np.count_nonzero(outside)),
                                   outside_beyond_roundoff=int(np.count_nonzero(outside > PROJECTION_EPS)),
                                   max_projection_overshoot_range_units=float(outside.max()*np.linalg.norm(direction)))
                out.writerows((i+1, kind, *p, d, t) for p, d, t in zip(xyz, error, s))
            results.append(dict(edge=("AB", "BC", "CA")[i], patch=i+1, triangles=len(patch),
                                components=components(patch, 2), samples=stats))
    sides = labels <= 3
    save_polydata(folder / "side_patches.vtp", points, faces[sides], labels[sides])
    return results


def cap_samples(points, faces, labels):
    # Classify each face directly in its two in-plane dimensions. A 3D volume
    # locator can give ambiguous membership for points exactly on a cap.
    n = 100000
    def wilson(k):
        z = 1.959963984540054
        f = k/n
        center = f + z*z/(2*n)
        radius = z * math.sqrt((f*(1-f) + z*z/(4*n))/n)
        return [max(0., (center-radius)/(1+z*z/n)), min(1., (center+radius)/(1+z*z/n))]

    runs = []
    for side, seed in enumerate(np.random.SeedSequence(2026).spawn(6)):
        axis = side // 2
        axes = [d for d in range(3) if d != axis]
        xyz = np.random.default_rng(seed).uniform(LO, HI, size=(n, 3))
        xyz[:, axis] = (HI if side % 2 else LO)[axis]
        uv = field(xyz)
        target = (uv.min(axis=1) >= 0) & (uv.sum(axis=1) <= 1)
        hit = np.zeros(n, dtype=bool)
        xy = xyz[:, axes]
        caps = faces[labels == 4+side]
        for tri in points[caps][:, :, axes]:
            subset = np.flatnonzero(((xy >= tri.min(axis=0)) & (xy <= tri.max(axis=0))).all(axis=1))
            if len(subset):
                a, b = (tri[1:] - tri[0])
                delta = xy[subset] - tri[0]
                det = a[0]*b[1] - a[1]*b[0]
                u = (delta[:, 0]*b[1] - delta[:, 1]*b[0]) / det
                v = (a[0]*delta[:, 1] - a[1]*delta[:, 0]) / det
                hit[subset] |= (u >= 0) & (v >= 0) & (u+v <= 1)
        missed, excess = int((target & ~hit).sum()), int((hit & ~target).sum())
        runs.append(dict(face=("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")[side],
                         triangles=len(caps), samples=n, both=int((target & hit).sum()),
                         neither=int((~target & ~hit).sum()), missed=missed, excess=excess,
                         missed_fraction=missed/n, excess_fraction=excess/n,
                         missed_wilson95=wilson(missed), excess_wilson95=wilson(excess)))
    return runs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path, nargs="?", default=DEFAULT_DIR)
    folder = parser.parse_args().folder
    report = json.loads((folder / "volume_report.json").read_text())
    if (not report["complete"] or report["triangle"] != CORNERS.tolist()
            or report["box_min"] != LO.tolist() or report["box_max"] != HI.tolist()):
        raise ValueError("Expected a completed full-domain curved prism with the manuscript triangle")
    vertices, cells, _ = read_mesh(folder / "volume_volume.vtk", volume=True)
    points, faces, labels = read_mesh(folder / "volume_tet_boundary.vtk")
    if not np.array_equal(points, vertices) or set(labels) != {1, 2, 3, 8, 9}:
        raise ValueError("Boundary and volume do not match the curved prism")
    tetra_faces = incidence(cells, 3)
    exterior = {face for face, ids in tetra_faces.items() if len(ids) == 1}
    if (len(faces) != len(exterior) or exterior != set(map(tuple, np.sort(faces, axis=1)))
            or any(len(ids) > 2 for ids in tetra_faces.values())):
        raise ValueError("Saved boundary does not match the exterior tetrahedral faces")
    fibers = evaluate_fibers(folder, points, faces, labels)
    patches = evaluate_patches(folder, points, faces, labels)
    tets = vertices[cells]
    volumes = np.linalg.det(tets[:, 1:] - tets[:, :1]) / 6
    volume = math.fsum(abs(volumes))
    if (not (np.all(volumes > 0) or np.all(volumes < 0))
            or not math.isclose(volume, report["volume"], rel_tol=1e-12)):
        raise ValueError("Cell orientation or independently summed volume failed")
    result = dict(triangle=CORNERS.tolist(), domain=[LO.tolist(), HI.tolist()],
                  tracing_step=report["fiber_step"], marching_step=report["step"],
                  spatial_distance="Euclidean after axis-wise normalization of the full domain to [0,1]^3",
                  projection_epsilon=PROJECTION_EPS, fibers=fibers, patches=patches,
                  side_surface_components=components(faces[labels <= 3], 2),
                  closed_boundary_components=components(faces, 2), closed_boundary_triangles=len(faces),
                  tetrahedra=len(cells), volume_components=components(cells, 3),
                  exact_volume=1., tetrahedral_volume=volume, relative_volume_error=abs(volume-1.),
                  face_sampling_seed=2026, cap_area_checks=cap_samples(points, faces, labels))
    (folder / "metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"{len(fibers)} fibers, {len(patches)} patches, {len(cells)} cells, volume {volume:.12f}")


if __name__ == "__main__":
    main()
