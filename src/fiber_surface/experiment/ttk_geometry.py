"""TTK side patches, endpoint fibers and sampled geometry distances for range triangles."""

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy


def extract_triangle(domain, names, query, threads, range_octree):
    from topologytoolkit import ttkFiberSurface

    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.column_stack((query, np.zeros(3))), deep=True))
    boundary = vtk.vtkUnstructuredGrid()
    boundary.SetPoints(pts)
    for i in range(3):
        boundary.InsertNextCell(vtk.VTK_LINE, 2, [i, (i+1) % 3])
    for i, name in enumerate(names):
        values = numpy_to_vtk(query[:, i].copy(), deep=True)
        values.SetName(name)
        boundary.GetPointData().AddArray(values)
    filt = ttkFiberSurface()
    filt.SetInputDataObject(0, domain)
    filt.SetInputDataObject(1, boundary)
    for i in range(4):
        filt.SetInputArrayToProcess(i, i//2, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, names[i%2])
    filt.SetUseAllCores(False)
    filt.SetThreadNumber(threads)
    filt.SetRangeOctree(range_octree)
    filt.SetRangeCoordinates(True)
    filt.SetEdgeIds(True)
    filt.SetPointMerge(False)
    filt.Update()
    raw = filt.GetOutput()
    if not raw.GetNumberOfPolys():
        raise ValueError("Empty ttkFiberSurface reference")
    ids = vtk_to_numpy(raw.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    patch = vtk_to_numpy(raw.GetCellData().GetArray("EdgeIds"))+1
    uv = np.column_stack([vtk_to_numpy(raw.GetPointData().GetArray(k)) for k in names])
    # Weld exactly coincident output coordinates. TTK emits duplicate vertices
    # per input tetrahedron and per polygon edge. Drop collapsed primitives.
    xyz, inverse = np.unique(vtk_to_numpy(raw.GetPoints().GetData()), axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    faces = inverse[ids]
    good = (np.diff(np.sort(faces, axis=1), axis=1) != 0).all(axis=1)
    edges, labels = [], []
    for i, q in enumerate(query):
        side = ids[patch == i+1]
        on = np.linalg.norm(uv-q, axis=1) < 1e-9
        ends = np.vstack([side[:, pair] for pair in ([0, 1], [1, 2], [2, 0])])
        ends = inverse[ends[on[ends].all(axis=1)]]
        ends = np.unique(np.sort(ends, axis=1), axis=0)
        ends = ends[ends[:, 0] != ends[:, 1]]
        edges.extend(ends)
        labels.extend([i]*len(ends))
    return (xyz, faces[good], patch[good], np.asarray(edges, dtype=int).reshape(-1, 2),
            np.asarray(labels, dtype=int), int((~good).sum()))


def original_curves(trace):
    result = []
    for i in range(3):
        records = (trace / f"segment_{i}/fibers.state").read_text().splitlines()
        n = int(records[1].split()[1])
        xyz = np.array([row.split()[:3] for row in records[2:n+2]], dtype=float)
        start = next(k for k, row in enumerate(records) if row.startswith("fibers "))
        count = int(records[start].split()[1])
        corner = []
        # Every corner is copied to both incident patches. Use its owning
        # patch's first endpoint once and keep the original trace connectivity.
        for record in records[start+1:start+1+count]:
            endpoint, closed, complete, length, *ids = map(int, record.split())
            if endpoint == 0:
                if length != len(ids):
                    raise ValueError("Invalid fiber checkpoint")
                curve = xyz[ids+[ids[0]]] if closed else xyz[ids]
                corner.append((curve, bool(closed), bool(complete)))
        result.append(corner)
    return result


def sizes(cells):
    v = cells[:, 1:] - cells[:, :1]
    if cells.shape[1] == 2:
        return np.linalg.norm(v[:, 0], axis=1)
    return np.linalg.norm(np.cross(v[:, 0], v[:, 1]), axis=1)/2


def directed_distances(src, dst, seed, lo, hi):
    # Uniform by source length or area in unit-cube coordinates. Distance is
    # to the continuous target segments/triangles, including their interiors.
    src = (src-lo)/(hi-lo)
    dst = (dst-lo)/(hi-lo)
    rng = np.random.default_rng(seed)
    weights = sizes(src)
    n = 10000
    primitives = src[rng.choice(len(src), n, p=weights/weights.sum())]
    u = rng.random((n, 1))
    if src.shape[1] == 2:
        samples = primitives[:, 0]*(1-u) + primitives[:, 1]*u
    else:
        u = np.sqrt(u)
        v = rng.random((n, 1))
        samples = primitives[:, 0]*(1-u) + primitives[:, 1]*u*(1-v) + primitives[:, 2]*u*v
    poly = vtk.vtkPolyData()
    p = vtk.vtkPoints()
    p.SetData(numpy_to_vtk(dst.reshape(-1, 3), deep=True))
    poly.SetPoints(p)
    cells = vtk.vtkCellArray()
    cells.SetData(dst.shape[1], numpy_to_vtkIdTypeArray(np.arange(len(dst)*dst.shape[1], dtype=np.int64), deep=True))
    (poly.SetLines if dst.shape[1] == 2 else poly.SetPolys)(cells)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()
    pos, cell, sub, d = [0.]*3, vtk.reference(0), vtk.reference(0), vtk.reference(0.)
    distances = np.empty(n)
    for i, x in enumerate(samples):
        locator.FindClosestPoint(x, pos, cell, sub, d)
        distances[i] = np.sqrt(float(d))
    return dict(median=float(np.median(distances)), p95=float(np.quantile(distances, .95)),
                sampled_max=float(distances.max()), samples=n, seed=seed)


def compare(cells, reference, seed, lo, hi):
    return dict(physical_measure=float(sizes(cells).sum()), ttk_physical_measure=float(sizes(reference).sum()),
                extracted_to_ttk=directed_distances(cells, reference, seed, lo, hi),
                ttk_to_extracted=directed_distances(reference, cells, seed+1, lo, hi))
