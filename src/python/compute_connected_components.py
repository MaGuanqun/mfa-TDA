#!/usr/bin/env pvpython

import argparse
import os
from typing import Dict, List, Optional, Tuple

import vtk

# try:
from paraview.simple import *
from paraview import servermanager
# except Exception:  # pragma: no cover
#     Connectivity = None  # type: ignore
#     FindSource = None  # type: ignore
#     GetActiveViewOrCreate = None  # type: ignore
#     XMLPolyDataReader = None  # type: ignore
#     paraview = None  # type: ignore
#     servermanager = None  # type: ignore


def plugin_log(is_server: int):
    """
    Load TTK plugin on local machine (0) or remote server (1).
    Adjust paths to your ParaView installations.
    """
    if is_server == 0:
        LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals()) 
    else:
        LoadPlugin("/home/u1435513-gma/apps/ParaView-5.11.2-osmesa-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",remote=False,ns=globals())

def _get_edge_type_array(polydata: vtk.vtkPolyData, edge_type_array: str):
    cell_data = polydata.GetCellData()
    arr = cell_data.GetArray(edge_type_array)
    if arr is not None:
        return arr

    scalars = cell_data.GetScalars()
    if scalars is not None:
        return scalars

    raise ValueError(
        f"Could not find cell edge type array '{edge_type_array}' "
        "and no active cell scalars are present."
    )


def _read_vtp_via_paraview(
    input_vtp: Optional[str],
    source_name: Optional[str],
    region_id_array: str,
) -> vtk.vtkPolyData:
    """
    Use ParaView's Connectivity filter to label regions and return a vtkPolyData
    that contains point-data `region_id_array` (default: RegionId).
    """
    if Connectivity is None or servermanager is None:
        raise RuntimeError(
            "ParaView is not available. Run this script with pvpython."
        )

    # Match ParaView trace defaults.

    if source_name is not None:
        if FindSource is None:
            raise RuntimeError("ParaView FindSource is unavailable (unexpected).")
        src = FindSource(source_name)
        if src is None:
            raise ValueError(f"Could not FindSource('{source_name}').")
    else:
        if input_vtp is None:
            raise ValueError("Provide --input_vtp when --source_name is not set.")
        if not os.path.isfile(input_vtp):
            raise FileNotFoundError(f"Input VTP not found: {input_vtp}")
        if XMLPolyDataReader is None:
            raise RuntimeError("ParaView XMLPolyDataReader is unavailable (unexpected).")
        src = XMLPolyDataReader(
            registrationName=os.path.basename(input_vtp),
            FileName=[os.path.abspath(input_vtp)],
        )
        try:
            src.TimeArray = "None"
        except Exception:
            pass

    # Create a view like the trace does (not strictly required for Fetch, but keeps parity).
    try:
        if GetActiveViewOrCreate is not None:
            GetActiveViewOrCreate("RenderView")
    except Exception:
        pass

    conn = Connectivity(Input=src)
    # Ensure RegionId is generated and all regions are labeled.
    try:
        conn.ColorRegions = 1
    except Exception:
        pass
    try:
        conn.RegionIdAssignmentMode = "Cell"
    except Exception:
        pass

    conn.UpdatePipeline()
    poly = servermanager.Fetch(conn)
    if poly is None:
        raise RuntimeError("Failed to fetch ParaView Connectivity output.")
    if poly.GetPointData().GetArray(region_id_array) is None:
        # ParaView's Connectivity typically uses "RegionId"
        raise ValueError(
            f"Connectivity output is missing point-data array '{region_id_array}'. "
            "Try --region_id_array RegionId."
        )
    return poly


def _build_line_edges(
    poly: vtk.vtkPolyData, edge_arr
) -> List[Tuple[int, int, int]]:
    """Return list of edges as (p0, p1, edge_type_int) in increasing cell-id order."""
    poly.BuildCells()
    edges: List[Tuple[int, int, int]] = []
    pts = vtk.vtkIdList()
    n_cells = poly.GetNumberOfCells()
    for cid in range(n_cells):
        if poly.GetCellType(cid) != vtk.VTK_LINE:
            continue
        poly.GetCellPoints(cid, pts)
        if pts.GetNumberOfIds() != 2:
            continue
        p0 = int(pts.GetId(0))
        p1 = int(pts.GetId(1))
        edge_type_val = edge_arr.GetTuple1(cid)
        edge_type_int = int(round(float(edge_type_val)))
        edges.append((p0, p1, edge_type_int))
    return edges


def _connected_components_by_region_id(poly: vtk.vtkPolyData, region_id_array: str):
    arr = poly.GetPointData().GetArray(region_id_array)
    if arr is None:
        raise ValueError(f"Missing point data array '{region_id_array}'.")
    comps: Dict[int, List[int]] = {}
    n_points = poly.GetNumberOfPoints()
    for pid in range(n_points):
        rid = int(round(float(arr.GetTuple1(pid))))
        comps.setdefault(rid, []).append(pid)
    return comps


def _write_component_obj_and_types(
    poly: vtk.vtkPolyData,
    comp_point_ids: List[int],
    edges: List[Tuple[int, int, int]],
    output_obj: str,
    output_edge_types_txt: str,
    w_arr,
):
    comp_point_set = set(comp_point_ids)
    comp_point_ids_sorted = sorted(comp_point_ids)
    old_to_new = {old: (i + 1) for i, old in enumerate(comp_point_ids_sorted)}

    # Keep edge order as in the original cell iteration order (via `edges` list).
    comp_edges_new: List[Tuple[int, int, int]] = []
    for p0, p1, edge_type_int in edges:
        if p0 in comp_point_set and p1 in comp_point_set:
            new_i = old_to_new[p0]
            new_j = old_to_new[p1]
            if new_i == new_j:
                continue
            comp_edges_new.append((new_i, new_j, edge_type_int))

    os.makedirs(os.path.dirname(os.path.abspath(output_obj)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(output_edge_types_txt)), exist_ok=True)

    with open(output_obj, "w", encoding="utf-8") as fobj:
        fobj.write(f"# Connected component: n_points={len(comp_point_ids_sorted)}\n")
        for old_pid in comp_point_ids_sorted:
            x, y, z = poly.GetPoint(old_pid)
            if w_arr is not None:
                w = w_arr.GetTuple1(old_pid)
                fobj.write(f"v {x} {y} {z} {w}\n")
            else:
                fobj.write(f"v {x} {y} {z}\n")
        for new_i, new_j, _edge_type_int in comp_edges_new:
            fobj.write(f"l {new_i} {new_j}\n")

    with open(output_edge_types_txt, "w", encoding="utf-8") as ft:
        for _new_i, _new_j, edge_type_int in comp_edges_new:
            ft.write(f"{edge_type_int}\n")

    return len(comp_point_ids_sorted), len(comp_edges_new)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Read a .vtp graph (VTK_LINE cells) with per-edge types and optional 4th coordinate "
            "as point-data array, compute all connected components, then export the largest K "
            "components as OBJ + edge type txt. Uses ParaView's Connectivity filter to label "
            "regions (RegionId)."
        )
    )
    parser.add_argument(
        "-i",
        "--input_vtp",
        type=str,
        required=False,
        help="Input .vtp (e.g. created by merge_obj_edge_type_4d.py).",
    )
    parser.add_argument(
        "--source_name",
        type=str,
        default=None,
        help=(
            "If running inside ParaView/pvpython with an existing pipeline source, "
            "use FindSource(source_name) instead of opening --input_vtp."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        type=str,
        required=True,
        help=(
            "Output prefix. Will write: <prefix>_cc1.obj, <prefix>_cc1_edge_types.txt, ... "
            "(cc1 is the largest component)."
        ),
    )
    parser.add_argument(
        "--num_k",
        type=int,
        default=3,
        help="How many largest components to export (default: 3).",
    )
    parser.add_argument(
        "--start_k",
        type=int,
        default=50,
        help="The starting index of the components to export (default: 50).",
    )
    parser.add_argument(
        "--edge_type_array",
        type=str,
        default="EdgeValues",
        help="Cell data array name for edge types (default: EdgeValues).",
    )
    parser.add_argument(
        "--w_array",
        type=str,
        default="t",
        help="Point data array name for 4th coordinate (default: t).",
    )
    parser.add_argument(
        "--region_id_array",
        type=str,
        default="RegionId",
        help="Point-data array name for connectivity region id (default: RegionId).",
    )
    parser.add_argument(
        "--require_w",
        action="store_true",
        help="Fail if the point-data w_array is missing (useful when you expect 4D).",
    )
    parser.add_argument(
        "--server",
        type=int,
        default=0,
        help="0 for local machine, 1 for remote server.",
    )

    args = parser.parse_args()


    plugin_log(is_server=args.server)

    poly = _read_vtp_via_paraview(
        input_vtp=args.input_vtp,
        source_name=args.source_name,
        region_id_array=args.region_id_array,
    )
    n_points = poly.GetNumberOfPoints()
    if n_points == 0:
        raise ValueError("Input polydata has 0 points.")

    w_arr = poly.GetPointData().GetArray(args.w_array)
    if args.require_w and w_arr is None:
        raise ValueError(
            f"Missing point data array '{args.w_array}' but --require_w was set."
        )

    edge_arr = _get_edge_type_array(poly, edge_type_array=args.edge_type_array)
    edges = _build_line_edges(poly, edge_arr=edge_arr)
    if not edges:
        raise ValueError("No VTK_LINE edges found in input polydata.")

    comps = _connected_components_by_region_id(poly, region_id_array=args.region_id_array)

    # Sort components by number of vertices (descending).
    comp_items = [(root, pids) for root, pids in comps.items()]
    comp_items.sort(key=lambda x: len(x[1]), reverse=True)

    # k = max(args.start_k, min(int(args.num_k)+args.start_k, len(comp_items)))
    # if k == 0:
    #     raise ValueError("No components found to export (top_k resolved to 0).")

    for idx in range(args.start_k, args.num_k+args.start_k):
        _root, pids = comp_items[idx]
        cc_name = f"{args.output_prefix}_cc{idx+1}"
        obj_out = cc_name + ".obj"
        edge_types_out = cc_name + "_edge_types.txt"
        n_v, n_e = _write_component_obj_and_types(
            poly=poly,
            comp_point_ids=pids,
            edges=edges,
            output_obj=obj_out,
            output_edge_types_txt=edge_types_out,
            w_arr=w_arr,
        )
        print(
            f"Wrote component {idx+1}/{args.num_k}: n_vertices={n_v} n_edges={n_e} "
            f"obj={obj_out} edge_types={edge_types_out}"
        )


if __name__ == "__main__":
    main()