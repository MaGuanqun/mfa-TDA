#!/usr/bin/env python

import argparse

# Use VTK from ParaView (works in pvpython)
import vtk
import csv

def filter_polydata_in_circle(polydata, cx, cy, radius):
    """Remove points with (x,y) inside circle and edges incident to them."""

    r2 = radius * radius

    points = polydata.GetPoints()
    n_points = points.GetNumberOfPoints()

    # Map old point id -> new point id (-1 means removed)
    point_id_map = [-1] * n_points
    new_points = vtk.vtkPoints()

    for pid in range(n_points):
        x, y, z = points.GetPoint(pid)
        dx = x - cx
        dy = y - cy
        if dx * dx + dy * dy > r2:
            # Keep this point
            new_pid = new_points.InsertNextPoint(x, y, z)
            point_id_map[pid] = new_pid
        # else: point is inside (or on) circle, remove it

    # Prepare new cell arrays for lines
    new_lines = vtk.vtkCellArray()

    cell_data = polydata.GetCellData()
    num_cell_arrays = cell_data.GetNumberOfArrays()
    original_arrays = []
    new_cell_arrays = []

    # Copy metadata for cell data arrays
    for ai in range(num_cell_arrays):
        arr = cell_data.GetArray(ai)
        original_arrays.append(arr)
        new_arr = arr.NewInstance()
        new_arr.SetName(arr.GetName())
        new_arr.SetNumberOfComponents(arr.GetNumberOfComponents())
        new_cell_arrays.append(new_arr)

    num_cells = polydata.GetNumberOfCells()
    id_list = vtk.vtkIdList()

    for cid in range(num_cells):
        # Only consider line cells
        if polydata.GetCellType(cid) != vtk.VTK_LINE:
            continue

        polydata.GetCellPoints(cid, id_list)
        if id_list.GetNumberOfIds() != 2:
            continue

        p0 = id_list.GetId(0)
        p1 = id_list.GetId(1)

        new_p0 = point_id_map[p0]
        new_p1 = point_id_map[p1]

        # Skip edges that touch removed points
        if new_p0 == -1 or new_p1 == -1:
            continue

        # Add new line cell
        new_lines.InsertNextCell(2)
        new_lines.InsertCellPoint(new_p0)
        new_lines.InsertCellPoint(new_p1)

        # Copy cell data for this cell
        for ai, arr in enumerate(original_arrays):
            tup = arr.GetTuple(cid)
            new_cell_arrays[ai].InsertNextTuple(tup)

    # Build new polydata
    new_polydata = vtk.vtkPolyData()
    new_polydata.SetPoints(new_points)
    new_polydata.SetLines(new_lines)

    new_cell_data = new_polydata.GetCellData()
    for new_arr in new_cell_arrays:
        new_cell_data.AddArray(new_arr)

    # Optionally set the first array as active scalars
    if new_cell_arrays:
        new_cell_data.SetScalars(new_cell_arrays[0])

    return new_polydata



def filter_our_csv_in_circle(input_csv, output_csv, cx, cy, radius, for_our_csv=False):
    """
    Read a CSV file, interpret columns:
      PositionX -> x
      PositionY -> y
      PositionZ -> z
    Remove rows where (x, y) is inside the same circle and write to output_csv.
    """
    r2 = radius * radius

    with open(input_csv, "r", newline="") as fin:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames

        if fieldnames is None:
            raise ValueError("CSV file has no header / fieldnames.")

        required_cols = ["x0", "x1", "x2"]
        for col in required_cols:
            if col not in fieldnames:
                raise ValueError(
                    f"Required column '{col}' not found in CSV header: {fieldnames}"
                )

        kept_rows = 0
        removed_rows = 0

        with open(output_csv, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                try:
                    x = float(row["x0"])
                    y = float(row["x1"])
                    # z is read but not used in the circle test
                    _ = float(row["x2"])
                except ValueError:
                    # If conversion fails, keep the row unchanged (or you can choose to skip)
                    writer.writerow(row)
                    kept_rows += 1
                    continue

                dx = x - cx
                dy = y - cy

                # Keep if outside circle
                if dx * dx + dy * dy > r2:
                    writer.writerow(row)
                    kept_rows += 1
                else:
                    removed_rows += 1

    print(f"[CSV] Kept rows: {kept_rows}, removed rows: {removed_rows}")
    # print(f"[CSV] Filtered CSV written to: {output_csv}")
    
    
def filter_csv_in_circle(input_csv, output_csv, cx, cy, radius):
    """
    Read a CSV file, interpret columns:
      PositionX -> x
      PositionY -> y
      PositionZ -> z
    Remove rows where (x, y) is inside the same circle and write to output_csv.
    """
    r2 = radius * radius

    with open(input_csv, "r", newline="") as fin:
        reader = csv.DictReader(fin)
        fieldnames = reader.fieldnames

        if fieldnames is None:
            raise ValueError("CSV file has no header / fieldnames.")

        required_cols = ["PositionX", "PositionY", "PositionZ"]
        for col in required_cols:
            if col not in fieldnames:
                raise ValueError(
                    f"Required column '{col}' not found in CSV header: {fieldnames}"
                )

        kept_rows = 0
        removed_rows = 0

        with open(output_csv, "w", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                try:
                    x = float(row["PositionX"])
                    y = float(row["PositionY"])
                    # z is read but not used in the circle test
                    _ = float(row["PositionZ"])
                except ValueError:
                    # If conversion fails, keep the row unchanged (or you can choose to skip)
                    writer.writerow(row)
                    kept_rows += 1
                    continue

                dx = x - cx
                dy = y - cy

                # Keep if outside circle
                if dx * dx + dy * dy > r2:
                    writer.writerow(row)
                    kept_rows += 1
                else:
                    removed_rows += 1

    print(f"[CSV] Kept rows: {kept_rows}, removed rows: {removed_rows}")
    print(f"[CSV] Filtered CSV written to: {output_csv}")
    
    
def main(input_vtp, output_vtp, input_csv, output_csv, cx, cy, radius, for_our_csv=False):
    
    if input_vtp!="input.vtp":
        # Read input .vtp
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(input_vtp)
        reader.Update()
        polydata = reader.GetOutput()

        # Filter
        filtered = filter_polydata_in_circle(polydata, cx, cy, radius)

        # Write output .vtp
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_vtp)
        writer.SetInputData(filtered)
        writer.Write()

        print(f"Filtered VTP written to: {output_vtp}")
        
    if input_csv!="input.csv":
        if for_our_csv:
            filter_our_csv_in_circle(input_csv, output_csv, cx, cy, radius)
        else:
            filter_csv_in_circle(input_csv, output_csv, cx, cy, radius)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove points (and their edges) whose (x,y) lie inside a circle."
    )
    parser.add_argument(
        "-i", "--input_vtp", type=str, default="input.vtp",
        help="Input .vtp file"
    )
    parser.add_argument(
        "-o", "--output_vtp", type=str,
        default="output.vtp",
        help="Output .vtp file"
    )
    parser.add_argument(
        "-c",
        "--input_csv",
        type=str,
        default="input.csv",
        help="Input CSV file (must contain PositionX, PositionY, PositionZ columns)",
    )
    parser.add_argument(
        "-d",
        "--output_csv",
        default="output.csv",
        type=str,
        help="Output filtered CSV file",
    )
    parser.add_argument(
        "--for_our_csv", type=int, default=0,
        help="Whether the input CSV is for our CSV format"
    )
    parser.add_argument(
        "--data", type=str, required=True,
        help="dataset name"
    )

    args = parser.parse_args()

    if args.data == "vortex_street_3d":
        center = [0,0]
        radius =  0.11
    elif args.data == "boussinesq_3d":
        center = [0,-0.15]
        radius = 0.09
    else:
        raise ValueError(f"Unknown data set: {args.data}")
    
    main(
        input_vtp=args.input_vtp,
        output_vtp=args.output_vtp,
        input_csv=args.input_csv,
        output_csv=args.output_csv,
        cx=center[0],
        cy=center[1],
        radius=radius,
        for_our_csv=args.for_our_csv,
    )
