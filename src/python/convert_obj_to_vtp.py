"""Convert an OBJ file containing 3D vertices and edges (lines) into a .vtp file.

The OBJ format is expected to contain:
  - vertex lines:  ``v x y z``
  - edge lines:    ``l i j [k ...]``  (1-based vertex indices, possibly a polyline)

OBJ indices are 1-based; they are converted to VTK's 0-based indexing.
"""

import argparse

import vtk


def read_obj(file_path):
    """Parse vertices and edges from an OBJ file.

    Returns a list of (x, y, z) vertices and a list of edges, where each edge
    is a tuple of 0-based vertex indices.
    """
    vertices = []
    edges = []
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split()
            if not parts:
                continue
            tag = parts[0]
            if tag == "v":  # Vertex (x, y, z)
                vertices.append(tuple(map(float, parts[1:4])))
            elif tag == "l":  # Edge / polyline (1-based indices)
                indices = [int(idx) - 1 for idx in parts[1:]]
                edges.append(tuple(indices))
    return vertices, edges


def build_polydata(vertices, edges):
    """Create a vtkPolyData from vertices and edges."""
    points = vtk.vtkPoints()
    for vertex in vertices:
        points.InsertNextPoint(vertex[0], vertex[1], vertex[2])

    lines = vtk.vtkCellArray()
    for edge in edges:
        if len(edge) < 2:
            continue
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(edge))
        for local_id, point_id in enumerate(edge):
            polyline.GetPointIds().SetId(local_id, point_id)
        lines.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    return polydata


def write_vtp(polydata, output_file):
    """Write a vtkPolyData to a .vtp file."""
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(polydata)
    writer.Write()


def main(obj_file, output_file):
    vertices, edges = read_obj(obj_file)

    polydata = build_polydata(vertices, edges)
    write_vtp(polydata, output_file)
    print(f"File saved: {output_file} ({len(vertices)} vertices, {len(edges)} edges)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert an OBJ file with vertices and edges to a .vtp file."
    )
    parser.add_argument(
        "-i", "--input_obj_name", type=str, default="file_name.obj", help="input obj file"
    )
    parser.add_argument(
        "-o", "--output_name", type=str, default="output.vtp", help="output vtp name"
    )
    args = parser.parse_args()

    main(args.input_obj_name, args.output_name)
