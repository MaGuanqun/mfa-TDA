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


def _read_vtp_raw(input_vtp: str) -> vtk.vtkPolyData:
    """Load .vtp without ParaView Connectivity — keeps input point-data RegionId (or any id array)."""
    if not os.path.isfile(input_vtp):
        raise FileNotFoundError(f"Input VTP not found: {input_vtp}")
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(os.path.abspath(input_vtp))
    reader.Update()
    poly = reader.GetOutput()
    if poly is None:
        raise RuntimeError(f"VTK failed to read: {input_vtp}")
    return poly


def _build_line_edges(
    poly: vtk.vtkPolyData, edge_arr
) -> List[Tuple[int, int, int, int]]:
    """Return list of edges as (cell_id, p0, p1, edge_type_int) in increasing cell-id order."""
    poly.BuildCells()
    edges: List[Tuple[int, int, int, int]] = []
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
        edges.append((cid, p0, p1, edge_type_int))
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
    edges: List[Tuple[int, int, int, int]],
    output_obj: str,
    output_edge_types_txt: str,
    w_arr,
):
    comp_point_set = set(comp_point_ids)
    comp_point_ids_sorted = sorted(comp_point_ids)
    old_to_new = {old: (i + 1) for i, old in enumerate(comp_point_ids_sorted)}

    # Keep edge order as in the original cell iteration order (via `edges` list).
    comp_edges_new: List[Tuple[int, int, int]] = []
    for _cid, p0, p1, edge_type_int in edges:
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


def _write_component_vtp(
    poly: vtk.vtkPolyData,
    comp_point_ids: List[int],
    edges: List[Tuple[int, int, int, int]],
    output_vtp: str,
    trajectory_order_idx_0based: int,
):
    comp_point_set = set(comp_point_ids)
    comp_point_ids_sorted = sorted(comp_point_ids)
    old_to_new = {old: i for i, old in enumerate(comp_point_ids_sorted)}

    out_poly = vtk.vtkPolyData()
    out_points = vtk.vtkPoints()
    out_points.SetNumberOfPoints(len(comp_point_ids_sorted))
    for old_pid, new_pid in old_to_new.items():
        out_points.SetPoint(new_pid, poly.GetPoint(old_pid))
    out_poly.SetPoints(out_points)

    # Copy all point-data arrays ("vortex values" and any other point arrays).
    src_point_data = poly.GetPointData()
    out_point_data = out_poly.GetPointData()
    for arr_idx in range(src_point_data.GetNumberOfArrays()):
        src_arr = src_point_data.GetArray(arr_idx)
        if src_arr is None:
            continue
        dst_arr = src_arr.NewInstance()
        dst_arr.DeepCopy(src_arr)
        dst_arr.SetNumberOfTuples(len(comp_point_ids_sorted))
        for old_pid, new_pid in old_to_new.items():
            dst_arr.SetTuple(new_pid, old_pid, src_arr)
        out_point_data.AddArray(dst_arr)
    if src_point_data.GetScalars() is not None:
        active_scalars_name = src_point_data.GetScalars().GetName()
        if active_scalars_name:
            out_point_data.SetActiveScalars(active_scalars_name)

    # 0-based index in global size-sorted order (0 = largest); same on every point.
    order_arr = vtk.vtkIntArray()
    order_arr.SetName("TrajectoryOrderIndex")
    order_arr.SetNumberOfTuples(len(comp_point_ids_sorted))
    for new_pid in range(len(comp_point_ids_sorted)):
        order_arr.SetValue(new_pid, int(trajectory_order_idx_0based))
    out_point_data.AddArray(order_arr)

    # Build line cells belonging to this connected component, while preserving
    # all cell-data arrays ("edge values" and any other cell arrays).
    out_lines = vtk.vtkCellArray()
    kept_cell_ids: List[int] = []
    line = vtk.vtkLine()
    for cid, p0, p1, _edge_type in edges:
        if p0 in comp_point_set and p1 in comp_point_set:
            line.GetPointIds().SetId(0, old_to_new[p0])
            line.GetPointIds().SetId(1, old_to_new[p1])
            out_lines.InsertNextCell(line)
            kept_cell_ids.append(cid)
    out_poly.SetLines(out_lines)

    src_cell_data = poly.GetCellData()
    out_cell_data = out_poly.GetCellData()
    n_out_cells = len(kept_cell_ids)
    for arr_idx in range(src_cell_data.GetNumberOfArrays()):
        src_arr = src_cell_data.GetArray(arr_idx)
        if src_arr is None:
            continue
        dst_arr = src_arr.NewInstance()
        dst_arr.DeepCopy(src_arr)
        dst_arr.SetNumberOfTuples(n_out_cells)
        for new_cid, old_cid in enumerate(kept_cell_ids):
            dst_arr.SetTuple(new_cid, old_cid, src_arr)
        out_cell_data.AddArray(dst_arr)
    if src_cell_data.GetScalars() is not None:
        active_cell_scalars_name = src_cell_data.GetScalars().GetName()
        if active_cell_scalars_name:
            out_cell_data.SetActiveScalars(active_cell_scalars_name)

    os.makedirs(os.path.dirname(os.path.abspath(output_vtp)), exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(output_vtp)
    writer.SetInputData(out_poly)
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write VTP file: {output_vtp}")

    return len(comp_point_ids_sorted), n_out_cells


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Read a .vtp graph (VTK_LINE cells) with per-edge types and optional 4th coordinate "
            "as point-data array, group points by RegionId, then export selected components as "
            "OBJ + edge type txt + VTP. By default, ParaView's Connectivity filter relabels regions "
            "(0..C-1). Use --no-connectivity to keep RegionId values already stored in the file "
            "(e.g. ids from branch splitting)."
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
            "Output prefix. Will write: <prefix>_cc1.obj, <prefix>_cc1_edge_types.txt, "
            "<prefix>_cc1.vtp, ... "
            "(cc1 is the largest component), plus <prefix>_trajectory_ranks.txt "
            "(TrajectoryOrderIndex per exported .vtp, 0=largest)."
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
        help=(
            "Starting index. In --selection_mode sorted, this is the rank start "
            "(0-based). In --selection_mode region_id, this is the first RegionId "
            "to export (default: 50)."
        ),
    )
    parser.add_argument(
        "--selection_mode",
        type=str,
        choices=["sorted", "region_id"],
        default="sorted",
        help=(
            "How to select components to export: "
            "'sorted' = by component size rank; "
            "'region_id' = exact RegionId values from [start_k, start_k + num_k)."
        ),
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
        "--no-connectivity",
        action="store_true",
        help=(
            "Do not run ParaView's Connectivity filter; read the .vtp with VTK only and "
            "split components by the RegionId (or --region_id_array) already on the file. "
            "Use this when your ids (e.g. 1081) come from branch splitting or another step — "
            "the default Connectivity pass renumbers regions to 0..C-1 and discards those labels."
        ),
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

    if args.no_connectivity and args.source_name:
        raise ValueError("--no-connectivity cannot be used with --source_name (use --input_vtp).")

    if args.no_connectivity:
        if not args.input_vtp:
            raise ValueError("--no-connectivity requires --input_vtp.")
    else:
        plugin_log(is_server=args.server)

    if args.no_connectivity:
        poly = _read_vtp_raw(args.input_vtp)
        if poly.GetPointData().GetArray(args.region_id_array) is None:
            raise ValueError(
                f"Input has no point-data array '{args.region_id_array}'. "
                "Add it to the .vtp or set --region_id_array."
            )
    else:
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

    # 0-based index in the global ordering by component size (largest = 0).
    comp_items_sorted = sorted(comps.items(), key=lambda x: len(x[1]), reverse=True)
    region_to_order_idx: Dict[int, int] = {
        rid: i for i, (rid, _) in enumerate(comp_items_sorted)
    }

    selected_components: List[Tuple[int, List[int], int, int]] = []
    if args.selection_mode == "sorted":
        # Sort components by number of vertices (descending).
        comp_items = comp_items_sorted
        end_idx = args.start_k + args.num_k
        if args.start_k < 0 or end_idx > len(comp_items):
            raise ValueError(
                f"Requested sorted indices [{args.start_k}, {end_idx}) but only "
                f"{len(comp_items)} components are available."
            )
        for idx in range(args.start_k, end_idx):
            region_id, pids = comp_items[idx]
            selected_components.append((idx + 1, pids, region_id, idx))
    else:
        # Exact RegionId selection in [start_k, start_k + num_k).
        end_region_id = args.start_k + args.num_k
        for rid in range(args.start_k, end_region_id):
            if rid not in comps:
                raise ValueError(
                    f"RegionId {rid} was requested but does not exist in input."
                )
            order_idx = region_to_order_idx[rid]
            selected_components.append((rid, comps[rid], rid, order_idx))

    ranks_manifest = os.path.abspath(f"{args.output_prefix}_trajectory_ranks.txt")
    ranks_dir = os.path.dirname(ranks_manifest)
    if ranks_dir:
        os.makedirs(ranks_dir, exist_ok=True)

    manifest_lines: List[str] = [
        "# vtp_filename\tTrajectoryOrderIndex(0=largest)\tRegionId\tselection_mode\n"
    ]
    for export_id, pids, region_id, order_idx in selected_components:
        cc_name = f"{args.output_prefix}_cc{export_id}"
        obj_out = cc_name + ".obj"
        edge_types_out = cc_name + "_edge_types.txt"
        vtp_out = cc_name + ".vtp"
        n_v, n_e = _write_component_obj_and_types(
            poly=poly,
            comp_point_ids=pids,
            edges=edges,
            output_obj=obj_out,
            output_edge_types_txt=edge_types_out,
            w_arr=w_arr,
        )
        n_v_vtp, n_e_vtp = _write_component_vtp(
            poly=poly,
            comp_point_ids=pids,
            edges=edges,
            output_vtp=vtp_out,
            trajectory_order_idx_0based=order_idx,
        )
        manifest_lines.append(
            f"{os.path.basename(vtp_out)}\t{order_idx}\t{region_id}\t"
            f"{args.selection_mode}\n"
        )
        print(
            f"Wrote component cc={export_id} (region_id={region_id}, "
            f"TrajectoryOrderIndex={order_idx}): "
            f"n_vertices={n_v} n_edges={n_e} "
            f"obj={obj_out} edge_types={edge_types_out} vtp={vtp_out} "
            f"(vtp_n_vertices={n_v_vtp}, vtp_n_edges={n_e_vtp})"
        )

    with open(ranks_manifest, "w", encoding="utf-8") as fr:
        fr.writelines(manifest_lines)
    print(f"Wrote trajectory rank manifest: {ranks_manifest}")


if __name__ == "__main__":
    main()