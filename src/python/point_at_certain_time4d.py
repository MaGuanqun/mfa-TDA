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
  - Optional scalars copied from input point data (default: ColorId, RegionId) when present:
    one int if unique after dedupe, otherwise -1. Cell-data ColorId is used only if point-data
    ColorId is absent but --copy-point-arrays includes ColorId.
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


def find_colorid_array_name(cell_data):
    """Return (name, vtkArray) for ColorId-like arrays on cell data, else (None, None)."""
    for name in ("ColorId", "colorId", "color_id"):
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


def _parse_copy_point_arrays(s):
    if not s or not str(s).strip():
        return []
    return [x.strip() for x in str(s).split(",") if x.strip()]


def _tracked_scalar_names(copy_names, point_data, cell_data):
    """
    Which extra scalar array names to write: point-data arrays that exist, plus
    ColorId from cell data only if ColorId is requested and missing on points.
    """
    tracked = []
    for n in copy_names:
        if point_data.GetArray(n) is not None:
            tracked.append(n)
        elif n == "ColorId":
            _cn, arr = find_colorid_array_name(cell_data)
            if arr is not None:
                tracked.append("ColorId")
    seen = set()
    out = []
    for n in tracked:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def point_scalar_int_or_none(point_data, name, pid):
    arr = point_data.GetArray(name)
    if arr is None:
        return None
    try:
        return int(round(float(arr.GetComponent(int(pid), 0))))
    except Exception:
        return None


def compute_time_intersections_with_int_types(
    vtp_path,
    t_target,
    time_array_name="t",
    type_array_name=None,
    tol=1e-9,
    dedup_tol=1e-8,
    copy_point_arrays="ColorId,RegionId",
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

    copy_names = _parse_copy_point_arrays(copy_point_arrays)
    tracked_extra = _tracked_scalar_names(copy_names, point_data, cell_data)
    use_cell_color_only = (
        "ColorId" in tracked_extra
        and point_data.GetArray("ColorId") is None
    )
    _, cell_color_arr = (
        find_colorid_array_name(cell_data) if use_cell_color_only else (None, None)
    )

    # Find type array
    if type_array_name is None:
        found_name, type_arr = find_type_array_name(cell_data, None)
        type_array_name = found_name
    else:
        type_arr = cell_data.GetArray(type_array_name)
        if type_arr is None:
            found_name, type_arr = find_type_array_name(cell_data, None)
            type_array_name = found_name

    def _append_hit(x, y, z, cid, seg_i, s_lp, hit, type_int, pids_for_extras):
        raw_points.append((x, y, z))
        cell_color_int = None
        if use_cell_color_only and cell_color_arr is not None:
            cell_color_int = variant_to_int_or_none(cell_color_arr, cid)
        raw_meta.append(
            {
                "cell": cid,
                "seg": seg_i,
                "lineParam": s_lp,
                "hit": hit,
                "type_int": type_int,
                "pids": pids_for_extras,
                "cell_color_int": cell_color_int,
            }
        )

    # Collect raw intersections
    raw_points = []   # (x,y,z)
    raw_meta = []     # dicts per hit

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
                _append_hit(x0, y0, z0, cid, i, 0.0, 1, type_int, [id0])

            # same-time segment
            if abs(t1 - t0) <= tol:
                if abs(t0 - t_target) <= tol:
                    _append_hit(x0, y0, z0, cid, i, 0.0, 2, type_int, [id0])
                    _append_hit(x1, y1, z1, cid, i, 1.0, 2, type_int, [id1])
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
                    pids_ex = [id0] if abs(s - 0.0) <= 1e-12 else [id1]
                else:
                    hit = 0
                    pids_ex = [id0, id1]
                _append_hit(x, y, z, cid, i, s, hit, type_int, pids_ex)

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
    uniq_meta = []  # agg dict incl. types_set, extras_sets

    def _feed_extras_into(agg, meta):
        ex = agg.setdefault("extras_sets", {})
        for name in tracked_extra:
            if name == "ColorId" and use_cell_color_only:
                if meta.get("cell_color_int") is not None:
                    ex.setdefault("ColorId", set()).add(int(meta["cell_color_int"]))
                continue
            for pid in meta.get("pids") or []:
                v = point_scalar_int_or_none(point_data, name, pid)
                if v is not None:
                    ex.setdefault(name, set()).add(int(v))

    for pt, meta in zip(raw_points, raw_meta):
        k = key_for_point(*pt)
        if k not in dedup_map:
            dedup_map[k] = len(uniq_points)
            uniq_points.append(pt)
            agg = {
                "cell": meta["cell"],
                "seg": meta["seg"],
                "lineParam": meta["lineParam"],
                "hit": meta["hit"],
                "types_set": set(),
                "extras_sets": {},
            }
            if meta.get("type_int") is not None:
                agg["types_set"].add(int(meta["type_int"]))
            _feed_extras_into(agg, meta)
            uniq_meta.append(agg)
        else:
            idx = dedup_map[k]
            if meta.get("type_int") is not None:
                uniq_meta[idx]["types_set"].add(int(meta["type_int"]))
            _feed_extras_into(uniq_meta[idx], meta)

    # Build output polydata
    out_pts = vtk.vtkPoints()
    out_verts = vtk.vtkCellArray()

    arr_origCellId = vtk.vtkIntArray()
    arr_origCellId.SetName("origCellId")
    arr_segIndex = vtk.vtkIntArray()
    arr_segIndex.SetName("segIndex")
    arr_lineParam = vtk.vtkDoubleArray()
    arr_lineParam.SetName("lineParam")
    arr_t = vtk.vtkDoubleArray()
    arr_t.SetName("t")
    arr_hit = vtk.vtkIntArray()
    arr_hit.SetName("hitType")

    arr_edge_type = vtk.vtkIntArray()
    arr_edge_type.SetName("edge_type")
    str_edge_types = vtk.vtkStringArray()
    str_edge_types.SetName("edge_types")

    extra_write = []
    for name in tracked_extra:
        arr_i = vtk.vtkIntArray()
        arr_i.SetName(name)
        arr_s = vtk.vtkStringArray()
        arr_s.SetName(name + "s")
        extra_write.append((name, arr_i, arr_s))

    for (x, y, z), agg in zip(uniq_points, uniq_meta):
        pid = out_pts.InsertNextPoint(float(x), float(y), float(z))
        out_verts.InsertNextCell(1)
        out_verts.InsertCellPoint(pid)

        arr_origCellId.InsertNextValue(int(agg["cell"]))
        arr_segIndex.InsertNextValue(int(agg["seg"]))
        arr_lineParam.InsertNextValue(float(agg["lineParam"]))
        arr_t.InsertNextValue(float(t_target))
        arr_hit.InsertNextValue(int(agg["hit"]))

        types_set = agg["types_set"]
        if len(types_set) == 0:
            edge_type_val = -1
            edge_types_str = ""
        elif len(types_set) == 1:
            edge_type_val = next(iter(types_set))
            edge_types_str = str(edge_type_val)
        else:
            edge_type_val = -1
            edge_types_str = ";".join(str(tt) for tt in sorted(types_set))

        arr_edge_type.InsertNextValue(int(edge_type_val))
        str_edge_types.InsertNextValue(edge_types_str)

        ex_sets = agg.get("extras_sets") or {}
        for name, arr_i, arr_s in extra_write:
            sset = ex_sets.get(name, set())
            if len(sset) == 0:
                arr_i.InsertNextValue(-1)
                arr_s.InsertNextValue("")
            elif len(sset) == 1:
                v = int(next(iter(sset)))
                arr_i.InsertNextValue(v)
                arr_s.InsertNextValue(str(v))
            else:
                arr_i.InsertNextValue(-1)
                arr_s.InsertNextValue(";".join(str(tt) for tt in sorted(sset)))

    out_poly = vtk.vtkPolyData()
    out_poly.SetPoints(out_pts)
    out_poly.SetVerts(out_verts)
    pd_out = out_poly.GetPointData()
    pd_out.AddArray(arr_origCellId)
    pd_out.AddArray(arr_segIndex)
    pd_out.AddArray(arr_lineParam)
    pd_out.AddArray(arr_t)
    pd_out.AddArray(arr_hit)
    pd_out.AddArray(arr_edge_type)
    pd_out.AddArray(str_edge_types)
    for _name, arr_i, arr_s in extra_write:
        pd_out.AddArray(arr_i)
        pd_out.AddArray(arr_s)
    pd_out.SetActiveScalars("t")

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
    p.add_argument(
        "--copy-point-arrays",
        type=str,
        default="ColorId,RegionId",
        help=(
            "Comma-separated point-data scalars to copy (first component, int). "
            "Missing arrays are skipped. ColorId falls back to cell data if absent on points. "
            "Use empty string to disable."
        ),
    )
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
    copy_point_arrays = args.copy_point_arrays

    print("Input:", in_vtp)
    print("Output:", out_vtp)
    print("Target time t =", t_target)
    print("Point-data time array =", time_array_name)
    if type_array_name:
        print("Using cell-data type array:", type_array_name)
    else:
        print("Auto-detecting integer cell-data type array (if present).")
    print("Copy point arrays:", copy_point_arrays or "(none)")

    out_poly = compute_time_intersections_with_int_types(
        in_vtp,
        t_target,
        time_array_name=time_array_name,
        type_array_name=type_array_name,
        tol=tol,
        dedup_tol=dedup_tol,
        copy_point_arrays=copy_point_arrays,
    )
    npts = out_poly.GetNumberOfPoints() if out_poly else 0
    print("Found {} point(s) at requested t.".format(npts))
    write_vtp(out_poly, out_vtp)
    print("Wrote:", out_vtp)
