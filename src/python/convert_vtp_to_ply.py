"""Convert a .vtp file containing vertices and edges (lines) into a .ply file.

The input VTP is expected to store geometry as points plus line cells (each line
cell may be a single segment or a polyline). Polylines are split into individual
edges so the output PLY uses the standard ``edge`` element with two endpoints.

The output is an ASCII PLY with:
  - a ``vertex`` element holding the x/y/z coordinates
  - an ``edge`` element holding pairs of 0-based vertex indices (vertex1, vertex2)
"""

import argparse

import vtk


def read_vtp(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def extract_vertices_and_edges(polydata):
    """Return (vertices, edges) from a vtkPolyData.

    vertices is a list of (x, y, z) tuples.
    edges is a list of (i, j) tuples of 0-based vertex indices. Polylines are
    expanded into consecutive segments.
    """
    points = polydata.GetPoints()
    num_points = points.GetNumberOfPoints() if points is not None else 0
    vertices = [points.GetPoint(i) for i in range(num_points)]

    edges = []
    lines = polydata.GetLines()
    if lines is not None:
        id_list = vtk.vtkIdList()
        lines.InitTraversal()
        while lines.GetNextCell(id_list):
            nids = id_list.GetNumberOfIds()
            for k in range(nids - 1):
                edges.append((id_list.GetId(k), id_list.GetId(k + 1)))

    return vertices, edges


def write_ply(vertices, edges, output_file):
    """Write vertices and edges to an ASCII PLY file."""
    with open(output_file, "w") as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write(f"element vertex {len(vertices)}\n")
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write(f"element edge {len(edges)}\n")
        f.write("property int vertex1\n")
        f.write("property int vertex2\n")
        f.write("end_header\n")

        for x, y, z in vertices:
            f.write(f"{x} {y} {z}\n")
        for i, j in edges:
            f.write(f"{i} {j}\n")


def main(input_file, output_file):
    polydata = read_vtp(input_file)
    vertices, edges = extract_vertices_and_edges(polydata)
    write_ply(vertices, edges, output_file)
    print(f"File saved: {output_file} ({len(vertices)} vertices, {len(edges)} edges)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a .vtp file with vertices and edges to a .ply file."
    )
    parser.add_argument(
        "-i", "--input_vtp_name", type=str, default="input.vtp", help="input vtp file"
    )
    parser.add_argument(
        "-o", "--output_name", type=str, default="output.ply", help="output ply name"
    )
    args = parser.parse_args()

    main(args.input_vtp_name, args.output_name)
