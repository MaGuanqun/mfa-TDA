#!/usr/bin/env python

import argparse

# Use VTK from ParaView (works in pvpython)
import vtk
import csv

def filter_polydata_in_circle(polydata, cx, cy, radius, return_removed=False):
    """Remove points on/inside the circle and lines incident to those points.

    Preserve point, cell and field attributes, including reference colors. Keep
    singleton detections as vertex cells; clipping a line may create new ones.
    The optional removed output contains inside points and fully inside lines.
    Crossing lines are discarded from both outputs, as in the original filter.
    """
    n_points = polydata.GetNumberOfPoints()
    r2 = radius * radius
    inside = []
    for pid in range(n_points):
        x, y, _ = polydata.GetPoint(pid)
        inside.append((x - cx)**2 + (y - cy)**2 <= r2)

    # Preserve the native VTK cell order: vertices first, followed by lines.
    vertices, lines = [], []
    incident_cell = [-1] * n_points
    for cid in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(cid)
        if cell.GetCellType() not in (vtk.VTK_VERTEX, vtk.VTK_LINE):
            raise ValueError("Obstacle filtering expects vertex cells and two-point lines")
        ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        (vertices if cell.GetCellType() == vtk.VTK_VERTEX else lines).append((cid, ids))
        for pid in ids:
            if incident_cell[pid] == -1:
                incident_cell[pid] = cid

    def select(remove_inside):
        selected = [pid for pid in range(n_points) if inside[pid] == remove_inside]
        mapping = {pid: index for index, pid in enumerate(selected)}
        result = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        if polydata.GetPoints() is not None:
            points.SetDataType(polydata.GetPoints().GetDataType())
        result.SetPoints(points)
        result.GetPointData().CopyAllocate(polydata.GetPointData(), len(selected))
        for new_id, old_id in enumerate(selected):
            points.InsertNextPoint(polydata.GetPoint(old_id))
            result.GetPointData().CopyData(polydata.GetPointData(), old_id, new_id)

        selected_lines = [(cid, ids) for cid, ids in lines if all(pid in mapping for pid in ids)]
        selected_vertices = [(cid, ids) for cid, ids in vertices if ids[0] in mapping]
        used = {pid for _, ids in selected_lines + selected_vertices for pid in ids}
        # An outside endpoint whose only incident line was clipped still needs
        # a cell, so that Surface representation displays every surviving point.
        orphan_vertices = [(incident_cell[pid], [pid]) for pid in selected if pid not in used]
        cells = selected_vertices + orphan_vertices + selected_lines
        result.GetCellData().CopyAllocate(polydata.GetCellData(), len(cells))
        for index in range(result.GetCellData().GetNumberOfArrays()):
            array = result.GetCellData().GetAbstractArray(index)
            array.SetNumberOfTuples(len(cells))
            # Orphan input points without an incident cell have no cell tuple.
            for value_index in range(array.GetNumberOfValues()):
                array.SetVariantValue(value_index, vtk.vtkVariant(0))
        for new_id, (old_id, _) in enumerate(cells):
            if old_id >= 0:
                result.GetCellData().CopyData(polydata.GetCellData(), old_id, new_id)
        # A newly isolated detection has no matched segment distance.
        distance = result.GetCellData().GetArray('MatchDistance')
        if distance is not None:
            for index in range(len(selected_vertices), len(selected_vertices) + len(orphan_vertices)):
                distance.SetTuple1(index, 0.)
        for collection, setter in [(selected_vertices + orphan_vertices, result.SetVerts),
                                   (selected_lines, result.SetLines)]:
            output_cells = vtk.vtkCellArray()
            for _, ids in collection:
                output_cells.InsertNextCell(len(ids))
                for pid in ids:
                    output_cells.InsertCellPoint(mapping[pid])
            setter(output_cells)
        result.GetFieldData().ShallowCopy(polydata.GetFieldData())
        return result

    kept = select(False)
    if return_removed:
        return kept, select(True)
    return kept



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
    
    
def main(
    input_vtp,
    output_vtp,
    input_csv,
    output_csv,
    cx,
    cy,
    radius,
    for_our_csv=False,
    removed_vtp="",
):
    
    if input_vtp!="input.vtp":
        # Read input .vtp
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(input_vtp)
        reader.Update()
        polydata = reader.GetOutput()

        # Filter
        if removed_vtp:
            filtered, removed_polydata = filter_polydata_in_circle(
                polydata, cx, cy, radius, return_removed=True
            )
        else:
            filtered = filter_polydata_in_circle(polydata, cx, cy, radius)

        # Write output .vtp
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(output_vtp)
        writer.SetInputData(filtered)
        writer.Write()

        print(f"Filtered VTP written to: {output_vtp}")

        if removed_vtp:
            removed_writer = vtk.vtkXMLPolyDataWriter()
            removed_writer.SetFileName(removed_vtp)
            removed_writer.SetInputData(removed_polydata)
            removed_writer.Write()
            print(f"Removed geometry VTP written to: {removed_vtp}")
        
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
        "--removed_vtp",
        type=str,
        default="",
        help="Optional output .vtp file for removed line geometry",
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
        radius =  0.12
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
        removed_vtp=args.removed_vtp,
    )
