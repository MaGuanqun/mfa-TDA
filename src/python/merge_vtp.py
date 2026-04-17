#!/usr/bin/env python

import argparse
import vtk


def read_vtp(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def merge_vtp(poly_a, poly_b):
    """
    Merge two vtkPolyData objects into one.

    - Keeps all points from both inputs (appended).
    - Keeps all line cells from both inputs (with ID offset for second input).
    - Preserves point and cell data arrays by taking the union of array names.
      Missing arrays in either input are filled with zeros.
    """
    out = vtk.vtkPolyData()

    points_a = poly_a.GetPoints()
    points_b = poly_b.GetPoints()
    npa = points_a.GetNumberOfPoints() if points_a is not None else 0
    npb = points_b.GetNumberOfPoints() if points_b is not None else 0

    merged_points = vtk.vtkPoints()

    for i in range(npa):
        merged_points.InsertNextPoint(points_a.GetPoint(i))
    for i in range(npb):
        merged_points.InsertNextPoint(points_b.GetPoint(i))

    out.SetPoints(merged_points)

    def collect_array_specs(data_a, data_b):
        specs = {}
        for data in (data_a, data_b):
            num = data.GetNumberOfArrays()
            for ai in range(num):
                arr = data.GetArray(ai)
                if arr is None:
                    continue
                name = arr.GetName()
                if not name:
                    continue
                ncomp = arr.GetNumberOfComponents()
                if name not in specs:
                    specs[name] = ncomp
                else:
                    specs[name] = max(specs[name], ncomp)
        return specs

    def build_array_map(data):
        arr_map = {}
        num = data.GetNumberOfArrays()
        for ai in range(num):
            arr = data.GetArray(ai)
            if arr is None:
                continue
            name = arr.GetName()
            if not name:
                continue
            arr_map[name] = arr
        return arr_map

    def pad_tuple(values, size):
        if len(values) >= size:
            return values[:size]
        return values + [0.0] * (size - len(values))

    # Merge point data arrays.
    pd_a = poly_a.GetPointData()
    pd_b = poly_b.GetPointData()
    point_specs = collect_array_specs(pd_a, pd_b)
    point_map_a = build_array_map(pd_a)
    point_map_b = build_array_map(pd_b)
    out_pd = out.GetPointData()

    for name, ncomp in point_specs.items():
        out_arr = vtk.vtkDoubleArray()
        out_arr.SetName(name)
        out_arr.SetNumberOfComponents(ncomp)

        arr_a = point_map_a.get(name)
        for i in range(npa):
            if arr_a is not None:
                vals = [arr_a.GetComponent(i, c) for c in range(arr_a.GetNumberOfComponents())]
            else:
                vals = []
            out_arr.InsertNextTuple(pad_tuple(vals, ncomp))

        arr_b = point_map_b.get(name)
        for i in range(npb):
            if arr_b is not None:
                vals = [arr_b.GetComponent(i, c) for c in range(arr_b.GetNumberOfComponents())]
            else:
                vals = []
            out_arr.InsertNextTuple(pad_tuple(vals, ncomp))

        out_pd.AddArray(out_arr)

    # Merge line cells and line cell data.
    merged_lines = vtk.vtkCellArray()
    in_lines_a = poly_a.GetLines()
    in_lines_b = poly_b.GetLines()

    cell_data_a = poly_a.GetCellData()
    cell_data_b = poly_b.GetCellData()
    cell_specs = collect_array_specs(cell_data_a, cell_data_b)
    cell_map_a = build_array_map(cell_data_a)
    cell_map_b = build_array_map(cell_data_b)
    out_cd = out.GetCellData()

    out_cell_arrays = {}
    for name, ncomp in cell_specs.items():
        arr = vtk.vtkDoubleArray()
        arr.SetName(name)
        arr.SetNumberOfComponents(ncomp)
        out_cell_arrays[name] = arr

    def append_lines_and_data(lines, pid_offset, arr_map):
        if lines is None:
            return
        id_list = vtk.vtkIdList()
        lines.InitTraversal()
        local_cid = 0
        while lines.GetNextCell(id_list):
            nids = id_list.GetNumberOfIds()
            if nids <= 0:
                local_cid += 1
                continue

            merged_lines.InsertNextCell(nids)
            for j in range(nids):
                merged_lines.InsertCellPoint(id_list.GetId(j) + pid_offset)

            for name, out_arr in out_cell_arrays.items():
                in_arr = arr_map.get(name)
                if in_arr is not None:
                    vals = [
                        in_arr.GetComponent(local_cid, c)
                        for c in range(in_arr.GetNumberOfComponents())
                    ]
                else:
                    vals = []
                out_arr.InsertNextTuple(pad_tuple(vals, out_arr.GetNumberOfComponents()))

            local_cid += 1

    append_lines_and_data(in_lines_a, 0, cell_map_a)
    append_lines_and_data(in_lines_b, npa, cell_map_b)

    out.SetLines(merged_lines)

    for arr in out_cell_arrays.values():
        out_cd.AddArray(arr)

    return out


def write_vtp(polydata, path):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(polydata)
    writer.Write()


def main():
    parser = argparse.ArgumentParser(
        description="Merge two input VTP files into one output VTP."
    )
    parser.add_argument("--input_a", "-a", required=True, help="First input .vtp file")
    parser.add_argument("--input_b", "-b", required=True, help="Second input .vtp file")
    parser.add_argument("--output", "-o", required=True, help="Merged output .vtp file")
    args = parser.parse_args()

    poly_a = read_vtp(args.input_a)
    poly_b = read_vtp(args.input_b)
    merged = merge_vtp(poly_a, poly_b)
    write_vtp(merged, args.output)

    print(f"Merged VTP written to: {args.output}")


if __name__ == "__main__":
    main()
