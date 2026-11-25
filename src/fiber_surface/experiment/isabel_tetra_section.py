"""Slice the tetrahedral VTK exported by the native fiber-volume pipeline."""
from itertools import combinations
import json

import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonDataModel import VTK_TETRA
from vtkmodules.vtkIOLegacy import vtkPolyDataReader, vtkUnstructuredGridReader


def mesh_section(mesh, x, volume_file):
    if not volume_file.name.endswith("_volume.vtk"):
        raise ValueError("Expected the fiber-volume pipeline's PREFIX_volume.vtk")
    prefix = volume_file.name[:-len("_volume.vtk")]
    report = json.loads(volume_file.with_name(prefix+"_report.json").read_text())
    if not report["complete"] or not report["boundary_valid"] or report["failed_stage"]:
        raise ValueError(f"Native volume did not pass validation: {report['failed_stage']}")

    # Confirm that both sides of the illustration come from the same repaired
    # surface. The native fill may add Steiner vertices to its tetrahedral mesh.
    reader = vtkPolyDataReader()
    reader.SetFileName(str(volume_file.with_name(prefix+"_boundary.vtk")))
    reader.ReadAllScalarsOn()
    reader.Update()
    boundary = reader.GetOutput()
    points, faces, labels = mesh
    if (boundary.GetNumberOfPolys() != len(faces) or boundary.GetNumberOfPoints() != len(points)
            or not np.array_equal(vtk_to_numpy(boundary.GetPoints().GetData()), points)
            or not np.array_equal(vtk_to_numpy(boundary.GetPolys().GetConnectivityArray()), faces.ravel())
            or not np.array_equal(vtk_to_numpy(boundary.GetCellData().GetArray("patch")), labels)):
        raise ValueError("Native volume boundary differs from the surface used in the figure")

    reader = vtkUnstructuredGridReader()
    reader.SetFileName(str(volume_file))
    reader.Update()
    grid = reader.GetOutput()
    n = grid.GetNumberOfCells()
    if (not n or n != report["tetrahedra"] or grid.GetNumberOfPoints() != report["tet_vertices"]
            or not np.all(vtk_to_numpy(grid.GetCellTypesArray()) == VTK_TETRA)):
        raise ValueError("The native VTK and its volume report disagree")
    vertices = vtk_to_numpy(grid.GetPoints().GetData())
    cells = vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(n, 4)

    # Preserve quadrilaterals when the plane crosses four edges of a tetrahedron.
    high = vertices[:, 0] > x
    polygons = []
    for t in vertices[cells[high[cells].any(axis=1) & ~high[cells].all(axis=1)]]:
        section = []
        for i, j in combinations(range(4), 2):
            a, b = t[i], t[j]
            if (a[0] > x) != (b[0] > x):
                section.append((a+(x-a[0])/(b[0]-a[0])*(b-a))[1:])
        poly = np.unique(section, axis=0)
        if len(poly) >= 3:
            dy, dz = (poly-poly.mean(axis=0)).T
            polygons.append(poly[np.argsort(np.arctan2(dz, dy))])
    return polygons, n
