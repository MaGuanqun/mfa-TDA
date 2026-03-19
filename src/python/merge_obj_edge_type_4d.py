from paraview.simple import *
import argparse

def read_obj(file_path):
    vertices = []
    edges = []
    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split()
            if not parts or parts[0].startswith("#"):
                continue

            if parts[0] == "v":  # Vertex (strict 4D)
                # OBJ 4D vertex should look like: v x y z w
                if len(parts) != 5:
                    raise ValueError(
                        f"Expected 4D OBJ vertex in form 'v x y z w', got: {line.strip()}"
                    )
                vertices.append(tuple(map(float, parts[1:])))  # (x, y, z, w)

            elif parts[0] == "l":  # Edge (strict 2-point line)
                # OBJ edge should look like: l i j  (1-based indices)
                if len(parts) != 3:
                    raise ValueError(
                        f"Expected edge in form 'l i j' (two indices), got: {line.strip()}"
                    )
                edges.append((int(parts[1]) - 1, int(parts[2]) - 1))  # Convert to 0-based indices
    return vertices, edges

def read_edge_values(file_path):
    with open(file_path, "r") as file:
        return [int(line.strip()) for line in file]
    
def create_programmable_source(vertices_xyz, vertices_w, edges, edge_values):
    """Create a Programmable Source in ParaView."""
    programmableSource = ProgrammableSource()

    # Set the dataset type
    programmableSource.OutputDataSetType = "vtkPolyData"

    # Generate the script for the Programmable Source
    vertices_str = str(vertices_xyz)
    vertices_w_str = str(vertices_w)
    edges_str = str(edges)
    edge_values_str = str(edge_values)

    programmableSource.Script = f"""
# Define vertices, edges, and edge values
vertices = {vertices_str}
vertices_w = {vertices_w_str}
edges = {edges_str}
edge_values = {edge_values_str}

# Import necessary modules from ParaView
from paraview.vtk import vtkPoints, vtkCellArray, vtkFloatArray, vtkLine

# Create points
points = vtkPoints()
for vertex in vertices:
    points.InsertNextPoint(vertex)

# Create lines
lines = vtkCellArray()
scalars = vtkFloatArray()
scalars.SetName("EdgeValues")
for i, edge in enumerate(edges):
    line = vtkLine()
    line.GetPointIds().SetId(0, edge[0])
    line.GetPointIds().SetId(1, edge[1])
    lines.InsertNextCell(line)
    scalars.InsertNextValue(edge_values[i])

# Attach the 4th coordinate as a separate point-data array.
# Note: VTK geometry itself is still 3D (x,y,z); we store w separately.
w_arr = vtkFloatArray()
w_arr.SetName("t")
w_arr.SetNumberOfComponents(1)
for w in vertices_w:
    w_arr.InsertNextValue(w)
output.GetPointData().AddArray(w_arr)

# Create the output polydata
output.SetPoints(points)
output.SetLines(lines)
output.GetCellData().SetScalars(scalars)
"""
    programmableSource.UpdatePipeline()
    return programmableSource


def main(obj_file, edge_values_file, output_file):
    # Read the input files
    vertices, edges = read_obj(obj_file)
    edge_values = read_edge_values(edge_values_file)

    assert len(edges) == len(edge_values), "Mismatch between edge count and edge values!"

    if not vertices:
        raise ValueError("OBJ file contains no vertices.")

    vertex_dim = len(vertices[0])
    if vertex_dim != 4 or any(len(v) != 4 for v in vertices):
        raise ValueError("This script expects ONLY strict 4D OBJ vertices: 'v x y z w'.")

    # VTK point coordinates are 3D, so we store (x, y, z) as geometry and (x, y, z, w)
    # as a point-data array for full 4D retention.
    vertices_xyz = [(v[0], v[1], v[2]) for v in vertices]
    vertices_w = [v[3] for v in vertices]

    # Create the dataset in ParaView
    programmableSource = create_programmable_source(
        vertices_xyz=vertices_xyz,
        vertices_w=vertices_w,
        edges=edges,
        edge_values=edge_values,
    )

    # Save the output as a .vtp file
    SaveData(output_file, proxy=programmableSource)
    print(f"File saved: {output_file}")


parser = argparse.ArgumentParser(description='merge obj edge types.')


parser.add_argument('-i', '--input_obj_name', type=str, default='file_name.obj', help='obj file')
parser.add_argument('-j', '--input_edge_type', type=str, default='file_name.txt', help='edge value file')
parser.add_argument('-o','--output_name', type=str, default='output.vtp', help='output csv name')


args = parser.parse_args()

obj_file_name=args.input_obj_name
type_file_name=args.input_edge_type
output_file=args.output_name

main(obj_file_name, type_file_name, output_file)


