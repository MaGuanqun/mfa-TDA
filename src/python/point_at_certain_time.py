#!/usr/bin/env pvpython
"""
pv_find_plane_intersections_with_int_types.py

Find intersections of polyline edges (cells) with plane z = z_target and write a .vtp of points.
Treats cell 'type' as integer codes (0,1,2,...). Writes:

 - origCellId (int)
 - segIndex   (int)
 - t          (double)
 - hitType    (int)  # 0=interior interpolation,1=exact vertex,2=co-planar endpoint
 - edge_type  (int)  # integer code if unique, otherwise -1 (ambiguous/multiple)
 - edge_types (string) # semicolon-separated list of types (present when multiple types aggregated)
 - ColorId    (int)  # integer color id if unique, otherwise -1 (ambiguous/multiple)

Usage:
  pvpython pv_find_plane_intersections_with_int_types.py --input in.vtp --output out.vtp --z 123.45 --type-array edge_type
"""

import sys
import argparse
import vtk
import math

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

def find_colorid_array_name(cell_data):
    """Return (name, vtkArray) for ColorId-like arrays if found, else (None, None)."""
    candidates = ["ColorId", "colorId", "color_id"]
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

def compute_plane_intersections_with_int_types(vtp_path, z_target, type_array_name=None, tol=1e-9, dedup_tol=1e-8):
    # Read input
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(vtp_path)
    reader.Update()
    mesh = reader.GetOutput()

    n_cells = mesh.GetNumberOfCells()
    pts = mesh.GetPoints()
    cell_data = mesh.GetCellData()

    # Find type array
    if type_array_name is None:
        found_name, type_arr = find_type_array_name(cell_data, None)
        type_array_name = found_name
    else:
        type_arr = cell_data.GetArray(type_array_name)
        if type_arr is None:
            found_name, type_arr = find_type_array_name(cell_data, None)
            type_array_name = found_name

    # Find ColorId array (optional)
    colorid_array_name, colorid_arr = find_colorid_array_name(cell_data)

    # Collect raw intersections
    raw_points = []   # (x,y,z)
    raw_meta = []     # dicts: {'cell':cid,'seg':i,'t':t,'hit':hit,'type_int':int_or_None,'colorid_int':int_or_None}

    for cid in range(n_cells):
        cell = mesh.GetCell(cid)
        npts_cell = cell.GetNumberOfPoints()
        if npts_cell < 2:
            continue

        # get integer type if available
        type_int = None
        if type_arr is not None:
            type_int = variant_to_int_or_none(type_arr, cid)
        colorid_int = None
        if colorid_arr is not None:
            colorid_int = variant_to_int_or_none(colorid_arr, cid)

        for i in range(npts_cell - 1):
            id0 = cell.GetPointId(i)
            id1 = cell.GetPointId(i+1)
            p0 = pts.GetPoint(id0)
            p1 = pts.GetPoint(id1)
            x0,y0,z0 = float(p0[0]), float(p0[1]), float(p0[2])
            x1,y1,z1 = float(p1[0]), float(p1[1]), float(p1[2])

            # vertex hit p0
            if abs(z0 - z_target) <= tol:
                raw_points.append((x0,y0,z0))
                raw_meta.append({'cell':cid,'seg':i,'t':0.0,'hit':1,'type_int':type_int,'colorid_int':colorid_int})

            # co-planar horizontal segment
            if abs(z1 - z0) <= tol:
                if abs(z0 - z_target) <= tol:
                    # endpoints
                    raw_points.append((x0,y0,z0))
                    raw_meta.append({'cell':cid,'seg':i,'t':0.0,'hit':2,'type_int':type_int,'colorid_int':colorid_int})
                    raw_points.append((x1,y1,z1))
                    raw_meta.append({'cell':cid,'seg':i,'t':1.0,'hit':2,'type_int':type_int,'colorid_int':colorid_int})
                continue

            # crossing test
            if (z_target - z0) * (z_target - z1) <= 0.0:
                t = (z_target - z0) / (z1 - z0)
                if t < -1e-12 or t > 1.0 + 1e-12:
                    continue
                t = max(0.0, min(1.0, t))
                x = x0 + t * (x1 - x0)
                y = y0 + t * (y1 - y0)
                z = z_target
                if abs(t - 0.0) <= 1e-12 or abs(t - 1.0) <= 1e-12:
                    hit = 1
                else:
                    hit = 0
                raw_points.append((x,y,z))
                raw_meta.append({'cell':cid,'seg':i,'t':t,'hit':hit,'type_int':type_int,'colorid_int':colorid_int})

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
    uniq_meta = []  # {'cell':first_cell, 'seg':first_seg, 't':first_t, 'hit':first_hit, 'types_set': set(), 'colorids_set': set() }

    for pt,meta in zip(raw_points, raw_meta):
        k = key_for_point(*pt)
        if k not in dedup_map:
            dedup_map[k] = len(uniq_points)
            uniq_points.append(pt)
            agg = {'cell': meta['cell'], 'seg': meta['seg'], 't': meta['t'], 'hit': meta['hit'], 'types_set': set(), 'colorids_set': set()}
            if meta.get('type_int') is not None:
                agg['types_set'].add(int(meta['type_int']))
            if meta.get('colorid_int') is not None:
                agg['colorids_set'].add(int(meta['colorid_int']))
            uniq_meta.append(agg)
        else:
            idx = dedup_map[k]
            if meta.get('type_int') is not None:
                uniq_meta[idx]['types_set'].add(int(meta['type_int']))
            if meta.get('colorid_int') is not None:
                uniq_meta[idx]['colorids_set'].add(int(meta['colorid_int']))

    # Build output polydata
    out_pts = vtk.vtkPoints()
    out_verts = vtk.vtkCellArray()

    arr_origCellId = vtk.vtkIntArray(); arr_origCellId.SetName("origCellId")
    arr_segIndex = vtk.vtkIntArray(); arr_segIndex.SetName("segIndex")
    arr_t = vtk.vtkDoubleArray(); arr_t.SetName("t")
    arr_hit = vtk.vtkIntArray(); arr_hit.SetName("hitType")
    arr_colorid = vtk.vtkIntArray(); arr_colorid.SetName("ColorId")

    arr_edge_type = vtk.vtkIntArray(); arr_edge_type.SetName("edge_type")  # single int or -1 if multiple
    str_edge_types = vtk.vtkStringArray(); str_edge_types.SetName("edge_types")  # semicolon list for multiples or single

    for (x,y,z), agg in zip(uniq_points, uniq_meta):
        pid = out_pts.InsertNextPoint(float(x), float(y), float(z))
        out_verts.InsertNextCell(1)
        out_verts.InsertCellPoint(pid)

        arr_origCellId.InsertNextValue(int(agg['cell']))
        arr_segIndex.InsertNextValue(int(agg['seg']))
        arr_t.InsertNextValue(float(agg['t']))
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

        colorids_set = agg['colorids_set']
        if len(colorids_set) == 1:
            colorid_val = next(iter(colorids_set))
        else:
            # no ColorId found or ambiguous due to dedup across differently colored edges
            colorid_val = -1
        arr_colorid.InsertNextValue(int(colorid_val))

    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_poly.SetVerts(out_verts)
    out_poly.GetPointData().AddArray(arr_origCellId)
    out_poly.GetPointData().AddArray(arr_segIndex)
    out_poly.GetPointData().AddArray(arr_t)
    out_poly.GetPointData().AddArray(arr_hit)
    out_poly.GetPointData().AddArray(arr_colorid)
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
    p = argparse.ArgumentParser(description="Find intersections of polylines with plane z and store integer edge types.")
    p.add_argument("--input", "-i", required=True, help="Input .vtp file")
    p.add_argument("--output", "-o", required=True, help="Output .vtp file")
    p.add_argument("--z", required=True, type=float, help="Target z coordinate")
    p.add_argument("--type-array", "-t", default=None, help="Name of cell-data array containing integer edge types (optional; auto-detect if omitted)")
    p.add_argument("--tol", type=float, default=1e-9, help="z equality tolerance")
    p.add_argument("--dedup-tol", type=float, default=1e-8, help="coordinate dedupe tolerance")
    return p.parse_args(argv)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    in_vtp = args.input
    out_vtp = args.output
    z_target = args.z
    type_array_name = args.type_array
    tol = args.tol
    dedup_tol = args.dedup_tol

    print("Input:", in_vtp)
    print("Output:", out_vtp)
    print("Plane z =", z_target)
    if type_array_name:
        print("Using cell-data type array:", type_array_name)
    else:
        print("Auto-detecting integer cell-data type array (if present).")

    out_poly = compute_plane_intersections_with_int_types(in_vtp, z_target, type_array_name=type_array_name, tol=tol, dedup_tol=dedup_tol)
    npts = out_poly.GetNumberOfPoints() if out_poly else 0
    print("Found {} intersection point(s).".format(npts))
    write_vtp(out_poly, out_vtp)
    print("Wrote:", out_vtp)
