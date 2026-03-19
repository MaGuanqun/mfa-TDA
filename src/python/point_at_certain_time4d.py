#!/usr/bin/env pvpython
"""
Extract points from 4D polyline data at a target time value.

Input:
  - Polyline geometry in 3D point coordinates (x, y, z)
  - 4th dimension stored in point-data array "t" (created by merge_obj_edge_type_4d.py)
  - Optional integer edge type stored in cell-data array (default auto-detect: EdgeValues)

Output point-data arrays:
  - origCellId (int)
  - segIndex   (int)
  - lineParam  (double)  # interpolation parameter along segment [0,1]
  - t          (double)  # target 4D time value
  - hitType    (int)     # 0=interior interpolation,1=exact vertex,2=co-time endpoint pair
  - edge_type  (int)     # integer code if unique, otherwise -1
  - edge_types (string)  # semicolon-separated list when multiple types are aggregated
"""

import argparse
import sys
import vtk

def find_type_array_name(cell_data, requested_name=None):
    """Return (name, vtkArray) if found, else (None, None)."""
    if requested_name:
        arr = cell_data.GetArray(requested_name)
        if arr:
            return requested_name, arr
    candidates = ["EdgeValues"]
    for name in candidates:
        arr = cell_data.GetArray(name)
        if arr:
            return name, arr
    return None, None

def variant_to_int_or_none(arr, idx):
    """
    Try to extract an integer from a vtk array at index idx.
    Returns an int on success, or None on failure.
    """
    try:
        # Prefer GetVariantValue if available (strings or variants)
        if hasattr(arr, "GetVariantValue"):
            v = arr.GetVariantValue(idx)
            # v may be a vtkVariant or similar; try to get integer
            try:
                return int(v)
            except Exception:
                try:
                    s = v.ToString()
                    return int(s)
                except Exception:
                    return None
        else:
            # numeric array types: GetValue
            v = arr.GetValue(idx)
            return int(v)
    except Exception:
        # fallback attempts
        try:
            v = arr.GetValue(idx)
            return int(v)
        except Exception:
            return None

def find_time_array_name(point_data, requested_name=None):
    """Return (name, vtkArray) for time array in point data, else (None, None)."""
    if requested_name:
        arr = point_data.GetArray(requested_name)
        if arr:
            return requested_name, arr
    for name in ["t", "time", "w"]:
        arr = point_data.GetArray(name)
        if arr:
            return name, arr
    return None, None


def compute_time_intersections_with_int_types(
    vtp_path, t_target, time_array_name="t", type_array_name=None, tol=1e-9, dedup_tol=1e-8
):
    # Read input
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    mesh = reader.GetOutput()

    n_cells = mesh.GetNumberOfCells()
    pts = mesh.GetPoints()
    point_data = mesh.GetPointData()
    cell_data = mesh.GetCellData()

    # Find time array in point data
    found_time_name, time_arr = find_time_array_name(point_data, time_array_name)
    if time_arr is None:
        raise ValueError(
            "Could not find point-data time array. Tried requested name '{}' plus candidates [t, time, w].".format(
                time_array_name
            )
        )
    time_array_name = found_time_name

    # Find type array
    if type_array_name is None:
        found_name, type_arr = find_type_array_name(cell_data, None)
        type_array_name = found_name
    else:
        type_arr = cell_data.GetArray(type_array_name)
        if type_arr is None:
            found_name, type_arr = find_type_array_name(cell_data, None)
            type_array_name = found_name

    # Collect raw intersections
    raw_points = []   # (x,y,z)
    raw_meta = []     # {'cell':cid,'seg':i,'lineParam':s,'hit':hit,'type_int':int_or_None}

    for cid in range(n_cells):
        cell = mesh.GetCell(cid)
        npts_cell = cell.GetNumberOfPoints()
        if npts_cell < 2:
            continue

        # get integer type if available
        type_int = None
        if type_arr is not None:
            type_int = variant_to_int_or_none(type_arr, cid)

        for i in range(npts_cell - 1):
            id0 = cell.GetPointId(i)
            id1 = cell.GetPointId(i+1)
            p0 = pts.GetPoint(id0)
            p1 = pts.GetPoint(id1)
            x0,y0,z0 = float(p0[0]), float(p0[1]), float(p0[2])
            x1,y1,z1 = float(p1[0]), float(p1[1]), float(p1[2])
            t0 = float(time_arr.GetComponent(id0, 0))
            t1 = float(time_arr.GetComponent(id1, 0))

            # vertex hit p0
            if abs(t0 - t_target) <= tol:
                raw_points.append((x0,y0,z0))
                raw_meta.append({'cell':cid,'seg':i,'lineParam':0.0,'hit':1,'type_int':type_int})

            # same-time segment
            if abs(t1 - t0) <= tol:
                if abs(t0 - t_target) <= tol:
                    raw_points.append((x0,y0,z0))
                    raw_meta.append({'cell':cid,'seg':i,'lineParam':0.0,'hit':2,'type_int':type_int})
                    raw_points.append((x1,y1,z1))
                    raw_meta.append({'cell':cid,'seg':i,'lineParam':1.0,'hit':2,'type_int':type_int})
                continue

            # crossing test in time dimension
            if (t_target - t0) * (t_target - t1) <= 0.0:
                s = (t_target - t0) / (t1 - t0)
                if s < -1e-12 or s > 1.0 + 1e-12:
                    continue
                s = max(0.0, min(1.0, s))
                x = x0 + s * (x1 - x0)
                y = y0 + s * (y1 - y0)
                z = z0 + s * (z1 - z0)
                if abs(s - 0.0) <= 1e-12 or abs(s - 1.0) <= 1e-12:
                    hit = 1
                else:
                    hit = 0
                raw_points.append((x,y,z))
                raw_meta.append({'cell':cid,'seg':i,'lineParam':s,'hit':hit,'type_int':type_int})

    # If no intersections
    if len(raw_points) == 0:
        out_poly = vtk.vtkPolyData()
        out_poly.SetPoints(vtk.vtkPoints())
        return out_poly

    # Deduplicate by rounded coords and aggregate integer types
    coord_tol = dedup_tol
    def key_for_point(x,y,z):
        return ( round(x/coord_tol)*coord_tol, round(y/coord_tol)*coord_tol, round(z/coord_tol)*coord_tol )

    dedup_map = {}
    uniq_points = []
    uniq_meta = []  # {'cell':first_cell, 'seg':first_seg, 'lineParam':first_s, 'hit':first_hit, 'types_set': set() }

    for pt,meta in zip(raw_points, raw_meta):
        k = key_for_point(*pt)
        if k not in dedup_map:
            dedup_map[k] = len(uniq_points)
            uniq_points.append(pt)
            agg = {'cell': meta['cell'], 'seg': meta['seg'], 'lineParam': meta['lineParam'], 'hit': meta['hit'], 'types_set': set()}
            if meta.get('type_int') is not None:
                agg['types_set'].add(int(meta['type_int']))
            uniq_meta.append(agg)
        else:
            idx = dedup_map[k]
            if meta.get('type_int') is not None:
                uniq_meta[idx]['types_set'].add(int(meta['type_int']))

    # Build output polydata
    out_pts = vtk.vtkPoints()
    out_verts = vtk.vtkCellArray()

    arr_origCellId = vtk.vtkIntArray(); arr_origCellId.SetName("origCellId")
    arr_segIndex = vtk.vtkIntArray(); arr_segIndex.SetName("segIndex")
    arr_lineParam = vtk.vtkDoubleArray(); arr_lineParam.SetName("lineParam")
    arr_t = vtk.vtkDoubleArray(); arr_t.SetName("t")
    arr_hit = vtk.vtkIntArray(); arr_hit.SetName("hitType")

    arr_edge_type = vtk.vtkIntArray(); arr_edge_type.SetName("edge_type")  # single int or -1 if multiple
    str_edge_types = vtk.vtkStringArray(); str_edge_types.SetName("edge_types")  # semicolon list for multiples or single

    for (x,y,z), agg in zip(uniq_points, uniq_meta):
        pid = out_pts.InsertNextPoint(float(x), float(y), float(z))
        out_verts.InsertNextCell(1)
        out_verts.InsertCellPoint(pid)

        arr_origCellId.InsertNextValue(int(agg['cell']))
        arr_segIndex.InsertNextValue(int(agg['seg']))
        arr_lineParam.InsertNextValue(float(agg['lineParam']))
        arr_t.InsertNextValue(float(t_target))
        arr_hit.InsertNextValue(int(agg['hit']))

        types_set = agg['types_set']
        if len(types_set) == 0:
            edge_type_val = -1
            edge_types_str = ""
        elif len(types_set) == 1:
            edge_type_val = next(iter(types_set))
            edge_types_str = str(edge_type_val)
        else:
            # multiple types -> mark -1 and provide semicolon list
            edge_type_val = -1
            edge_types_str = ";".join(str(tt) for tt in sorted(types_set))

        arr_edge_type.InsertNextValue(int(edge_type_val))
        str_edge_types.InsertNextValue(edge_types_str)

    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_poly.SetVerts(out_verts)
    out_poly.GetPointData().AddArray(arr_origCellId)
    out_poly.GetPointData().AddArray(arr_segIndex)
    out_poly.GetPointData().AddArray(arr_lineParam)
    out_poly.GetPointData().AddArray(arr_t)
    out_poly.GetPointData().AddArray(arr_hit)
    out_poly.GetPointData().AddArray(arr_edge_type)
    out_poly.GetPointData().AddArray(str_edge_types)
    out_poly.GetPointData().SetActiveScalars("t")

    return out_poly

def write_vtp(polydata, out_path):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(out_path)
    writer.SetDataModeToBinary()
    try:
        writer.SetInputData(polydata)
    except AttributeError:
        writer.SetInput(polydata)
    writer.Write()

def parse_args(argv):
    p = argparse.ArgumentParser(description="Find intersections of 4D polylines with a target point-data time value.")
    p.add_argument("--input", "-i", required=True, help="Input .vtp file")
    p.add_argument("--output", "-o", required=True, help="Output .vtp file")
    p.add_argument("--t", required=True, type=float, help="Target time value (4th dimension)")
    p.add_argument("--time-array", default="t", help="Name of point-data array storing time (default: t)")
    p.add_argument("--type-array", default=None, help="Name of cell-data array containing integer edge types (optional; auto-detect if omitted)")
    p.add_argument("--tol", type=float, default=1e-9, help="Time equality tolerance")
    p.add_argument("--dedup-tol", type=float, default=1e-8, help="coordinate dedupe tolerance")
    return p.parse_args(argv)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    in_vtp = args.input
    out_vtp = args.output
    t_target = args.t
    time_array_name = args.time_array
    type_array_name = args.type_array
    tol = args.tol
    dedup_tol = args.dedup_tol

    print("Input:", in_vtp)
    print("Output:", out_vtp)
    print("Target time t =", t_target)
    print("Point-data time array =", time_array_name)
    if type_array_name:
        print("Using cell-data type array:", type_array_name)
    else:
        print("Auto-detecting integer cell-data type array (if present).")

    out_poly = compute_time_intersections_with_int_types(
        in_vtp,
        t_target,
        time_array_name=time_array_name,
        type_array_name=type_array_name,
        tol=tol,
        dedup_tol=dedup_tol,
    )
    npts = out_poly.GetNumberOfPoints() if out_poly else 0
    print("Found {} point(s) at requested t.".format(npts))
    write_vtp(out_poly, out_vtp)
    print("Wrote:", out_vtp)
