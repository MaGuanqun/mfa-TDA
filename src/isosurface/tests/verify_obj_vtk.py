"""Optional independent OBJ verification; requires the vtk Python package.

Usage: python verify_obj_vtk.py path/to/sheet.obj
Shared-vertex triangle pairs are excluded from this independent intersection
check. The C++ tests cover the shared-vertex and coplanar predicates separately.
"""
import argparse
import json

from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator, vtkTriangle
from vtkmodules.vtkFiltersCore import (
    vtkFeatureEdges,
    vtkMassProperties,
    vtkPolyDataConnectivityFilter,
)
from vtkmodules.vtkIOGeometry import vtkOBJReader

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("obj")
args = parser.parse_args()
reader = vtkOBJReader()
reader.SetFileName(args.obj)
reader.Update()
mesh = reader.GetOutput()
if mesh.GetNumberOfCells() == 0:
    raise SystemExit("No faces found")

boundary = vtkFeatureEdges()
boundary.SetInputData(mesh)
boundary.BoundaryEdgesOn()
boundary.NonManifoldEdgesOn()
boundary.FeatureEdgesOff()
boundary.ManifoldEdgesOff()
boundary.Update()
connectivity = vtkPolyDataConnectivityFilter()
connectivity.SetInputData(mesh)
connectivity.SetExtractionModeToAllRegions()
connectivity.Update()
mass = vtkMassProperties()
mass.SetInputData(mesh)
mass.Update()
locator = vtkStaticCellLocator()
locator.SetDataSet(mesh)
locator.BuildLocator()
nearby = vtkIdList()
intersections = checked = 0
for i in range(mesh.GetNumberOfCells()):
    cell = mesh.GetCell(i)
    if cell.GetNumberOfPoints() != 3:
        raise SystemExit("Expected a triangular mesh")
    ids = [cell.GetPointId(k) for k in range(3)]
    points = [mesh.GetPoint(k) for k in ids]
    locator.FindCellsWithinBounds(cell.GetBounds(), nearby)
    for k in range(nearby.GetNumberOfIds()):
        j = nearby.GetId(k)
        if j <= i:
            continue
        other = mesh.GetCell(j)
        if other.GetNumberOfPoints() != 3:
            raise SystemExit("Expected a triangular mesh")
        other_ids = [other.GetPointId(k) for k in range(3)]
        if set(ids).intersection(other_ids):
            continue
        other_points = [mesh.GetPoint(k) for k in other_ids]
        checked += 1
        intersections += vtkTriangle.TrianglesIntersect(*points, *other_points)

report = {
    "vertices": mesh.GetNumberOfPoints(),
    "triangles": mesh.GetNumberOfCells(),
    "boundary_or_nonmanifold_edges": boundary.GetOutput().GetNumberOfCells(),
    "components": connectivity.GetNumberOfExtractedRegions(),
    "nonadjacent_triangle_intersections": intersections,
    "candidate_pairs_checked": checked,
    "volume": mass.GetVolume(),
    "area": mass.GetSurfaceArea(),
}
print(json.dumps(report, indent=2))
if report["boundary_or_nonmanifold_edges"] or report["components"] != 1 or intersections:
    raise SystemExit(1)
