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
    vtkAppendPolyData,
    vtkFeatureEdges,
    vtkMassProperties,
    vtkPolyDataConnectivityFilter,
)
from vtkmodules.vtkIOGeometry import vtkOBJReader

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("obj", nargs="+")
parser.add_argument("--allow-open", action="store_true")
parser.add_argument("--components", type=int, default=1)
args = parser.parse_args()
append = vtkAppendPolyData()
readers = []
for filename in args.obj:
    reader = vtkOBJReader()
    reader.SetFileName(filename)
    reader.Update()
    if reader.GetOutput().GetNumberOfCells() == 0:
        raise SystemExit(f"No faces found in {filename}")
    append.AddInputData(reader.GetOutput())
    readers.append(reader)
append.Update()
mesh = append.GetOutput()
if mesh.GetNumberOfCells() == 0:
    raise SystemExit("No faces found")

boundary = vtkFeatureEdges()
boundary.SetInputData(mesh)
boundary.BoundaryEdgesOn()
boundary.NonManifoldEdgesOff()
boundary.FeatureEdgesOff()
boundary.ManifoldEdgesOff()
boundary.Update()
nonmanifold = vtkFeatureEdges()
nonmanifold.SetInputData(mesh)
nonmanifold.BoundaryEdgesOff()
nonmanifold.NonManifoldEdgesOn()
nonmanifold.FeatureEdgesOff()
nonmanifold.ManifoldEdgesOff()
nonmanifold.Update()
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
    "boundary_edges": boundary.GetOutput().GetNumberOfCells(),
    "nonmanifold_edges": nonmanifold.GetOutput().GetNumberOfCells(),
    "components": connectivity.GetNumberOfExtractedRegions(),
    "nonadjacent_triangle_intersections": intersections,
    "candidate_pairs_checked": checked,
    "volume": mass.GetVolume() if not args.allow_open else None,
    "area": mass.GetSurfaceArea(),
}
print(json.dumps(report, indent=2))
if (report["boundary_edges"] and not args.allow_open) or report["nonmanifold_edges"] or report["components"] != args.components or intersections:
    raise SystemExit(1)
