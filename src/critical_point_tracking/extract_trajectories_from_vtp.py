#!/usr/bin/env pvpython

import argparse
from collections import deque
import os
from typing import Optional

import vtk


class DSU:
    """Disjoint-set union (union-find) for fast connected components."""

    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def _get_edge_type_array(polydata: vtk.vtkPolyData, edge_type_array: str):
    cell_data = polydata.GetCellData()
    arr = cell_data.GetArray(edge_type_array)
    if arr is not None:
        return arr

    # Common fallback: "EdgeValues" is usually stored as active scalars.
    scalars = cell_data.GetScalars()
    if scalars is not None:
        return scalars

    raise ValueError(
        f"Could not find cell edge type array '{edge_type_array}' "
        "and no active cell scalars are present."
    )


def extract_component_and_write_obj(
    input_vtp: str,
    seed_point_id: int,
    seed_point_id_base: int,
    output_obj: str,
    output_edge_types_txt: str,
    edge_type_array: str,
    w_array: str,
    method: str,
    output_point_t_txt: Optional[str],
):
    if seed_point_id_base not in (0, 1):
        raise ValueError("--seed_point_id_base must be 0 or 1.")

    # Convert to VTK's 0-based point ids.
    seed_pid = seed_point_id - seed_point_id_base
    if seed_pid < 0:
        raise ValueError("Seed point id becomes negative after base conversion.")

    if not os.path.isfile(input_vtp):
        raise FileNotFoundError(f"Input VTP not found: {input_vtp}")

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(input_vtp)
    reader.Update()
    poly = reader.GetOutput()

    n_points = poly.GetNumberOfPoints()
    if seed_pid >= n_points:
        raise ValueError(f"Seed point id {seed_pid} out of range [0, {n_points - 1}].")

    # Detect w (4th coordinate) if present, so merge_obj_edge_type.py can embed 4D->3D.
    point_data = poly.GetPointData()
    w_arr = point_data.GetArray(w_array)
    has_w = w_arr is not None

    edge_arr = _get_edge_type_array(poly, edge_type_array=edge_type_array)

    # Speed-up: build internal connectivity once so GetPointCells(pid) is fast.
    # Without links, VTK may need to do extra work per query.
    poly.BuildCells()
    poly.BuildLinks()

    method = method.lower().strip()
    if method not in ("bfs", "dsu"):
        raise ValueError("--method must be one of: bfs, dsu")

    # Extract only the component reachable from seed_point_id.
    # `bfs` uses VTK's GetPointCells() to expand the graph locally, avoiding scanning all points.
    # `dsu` computes components globally (union-find) and then filters.
    if method == "bfs":
        visited_points: set[int] = {seed_pid}
        component_points_ordered: list[int] = [seed_pid]
        seen_edge_cell_ids: set[int] = set()
        component_edges_old: list[tuple[int, int, int]] = []  # (p0_old, p1_old, edge_type_int)

        q: deque[int] = deque([seed_pid])
        cell_ids = vtk.vtkIdList()
        pts = vtk.vtkIdList()

        while q:
            pid = q.popleft()

            # Incident cells for this point (includes line edges).
            poly.GetPointCells(pid, cell_ids)
            for k in range(cell_ids.GetNumberOfIds()):
                cid = cell_ids.GetId(k)
                if cid in seen_edge_cell_ids:
                    continue
                if poly.GetCellType(cid) != vtk.VTK_LINE:
                    continue

                poly.GetCellPoints(cid, pts)
                if pts.GetNumberOfIds() != 2:
                    continue

                p0_old = pts.GetId(0)
                p1_old = pts.GetId(1)

                seen_edge_cell_ids.add(cid)

                # Record edge type value for this exact VTK cell id.
                edge_type_val = edge_arr.GetTuple1(cid)
                edge_type_int = int(round(float(edge_type_val)))
                component_edges_old.append((p0_old, p1_old, edge_type_int))

                # Expand BFS frontier.
                if p0_old not in visited_points:
                    visited_points.add(p0_old)
                    component_points_ordered.append(p0_old)
                    q.append(p0_old)
                if p1_old not in visited_points:
                    visited_points.add(p1_old)
                    component_points_ordered.append(p1_old)
                    q.append(p1_old)

        if len(component_points_ordered) == 0:
            raise RuntimeError("Connected component is empty (unexpected).")

        old_pid_to_new = {
            old_pid: (new_i + 1) for new_i, old_pid in enumerate(component_points_ordered)
        }

        # Final edges written as (new_i, new_j, edge_type_int)
        component_edges = [
            (old_pid_to_new[p0_old], old_pid_to_new[p1_old], edge_type_int)
            for (p0_old, p1_old, edge_type_int) in component_edges_old
            if old_pid_to_new[p0_old] != old_pid_to_new[p1_old]
        ]

    else:
        # DSU fallback (global union-find). Useful if you need stable ordering (we sort vertices).
        dsu = DSU(n_points)
        edge_cells = []  # (p0, p1, cell_id) in the original cell iteration order

        id_list = vtk.vtkIdList()
        n_cells = poly.GetNumberOfCells()
        for cid in range(n_cells):
            if poly.GetCellType(cid) != vtk.VTK_LINE:
                continue
            poly.GetCellPoints(cid, id_list)
            if id_list.GetNumberOfIds() != 2:
                continue
            p0 = id_list.GetId(0)
            p1 = id_list.GetId(1)
            dsu.union(p0, p1)
            edge_cells.append((p0, p1, cid))

        seed_root = dsu.find(seed_pid)

        # Collect points in the connected component.
        component_points = [pid for pid in range(n_points) if dsu.find(pid) == seed_root]
        component_point_set = set(component_points)

        if not component_points:
            raise RuntimeError("Connected component is empty (unexpected).")

        component_points_sorted = sorted(component_points)
        old_pid_to_new = {
            old_pid: (new_i + 1) for new_i, old_pid in enumerate(component_points_sorted)
        }

        # Filter edges that lie fully within the component, and keep order aligned with OBJ lines.
        component_edges = []  # (new_i, new_j, edge_type_int)
        for p0, p1, cid in edge_cells:
            if p0 in component_point_set and p1 in component_point_set:
                new_i = old_pid_to_new[p0]
                new_j = old_pid_to_new[p1]
                if new_i == new_j:
                    continue
                edge_type_val = edge_arr.GetTuple1(cid)
                edge_type_int = int(round(float(edge_type_val)))
                component_edges.append((new_i, new_j, edge_type_int))

        component_points_ordered = component_points_sorted

    os.makedirs(os.path.dirname(os.path.abspath(output_obj)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(output_edge_types_txt)), exist_ok=True)

    # Write OBJ.
    with open(output_obj, "w", encoding="utf-8") as fobj:
        fobj.write(f"# Extracted connected component containing VTK point id {seed_pid}\n")
        if method == "dsu":
            component_points_ordered = component_points_ordered  # sorted
        for old_pid in component_points_ordered:
            x, y, z = poly.GetPoint(old_pid)
            if has_w:
                w = w_arr.GetTuple1(old_pid)
                fobj.write(f"v {x} {y} {z} {w}\n")
            else:
                fobj.write(f"v {x} {y} {z}\n")
        for new_i, new_j, _edge_type_int in component_edges:
            fobj.write(f"l {new_i} {new_j}\n")

    # Write edge types in the exact same order as OBJ edge lines.
    with open(output_edge_types_txt, "w", encoding="utf-8") as ft:
        for _new_i, _new_j, edge_type_int in component_edges:
            ft.write(f"{edge_type_int}\n")

    # Optionally write point->t mapping (old VTK point id -> t value).
    if output_point_t_txt is not None:
        if not has_w:
            raise ValueError(
                f"--output_point_t_txt requested but point data array '{w_array}' does not exist."
            )
        os.makedirs(os.path.dirname(os.path.abspath(output_point_t_txt)), exist_ok=True)
        with open(output_point_t_txt, "w", encoding="utf-8") as fpt:
            fpt.write("# old_vtk_point_id t\n")
            for old_pid in component_points_ordered:
                t_val = w_arr.GetTuple1(old_pid)
                fpt.write(f"{old_pid} {t_val}\n")

    print(
        "Component extracted.",
        f"seed_vtk_point_id={seed_pid}",
        f"n_component_points={len(component_points_ordered)}",
        f"n_component_edges={len(component_edges)}",
        f"output_obj={output_obj}",
        f"output_edge_types_txt={output_edge_types_txt}",
    )


def get_seed_point_id(data: str):
    if data == "vortex":
        return [177845,337745,186904, 155672]
    else:
        raise ValueError(f"Invalid data: {data}")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract connected component(s) from a .vtp edge graph (VTK_LINE cells). "
            "You can provide one seed (--seed_point_id) or multiple (--seed_point_ids), "
            "or use a preset seed list via --data."
        )
    )
    parser.add_argument(
        "-i",
        "--input_vtp",
        type=str,
        required=True,
        help="Input .vtp containing line edges (VTK_LINE cells) and cell edge-type values.",
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="Optional preset name for a list of seed point ids (e.g. 'vortex').",
    )
    parser.add_argument(
        "--seed_point_id",
        type=int,
        default=None,
        help="Single seed point id (use --seed_point_ids for multiple).",
    )
    parser.add_argument(
        "--seed_point_ids",
        type=int,
        nargs="+",
        default=None,
        help="Multiple seed point ids. Writes one OBJ/TXT per seed.",
    )
    parser.add_argument(
        "--seed_point_id_base",
        type=int,
        default=0,
        choices=[0, 1],
        help="0 if seed ids are VTK 0-based point ids; 1 if they are 1-based.",
    )
    parser.add_argument(
        "-o",
        "--output_obj",
        type=str,
        required=True,
        help=(
            "Output .obj file. If multiple seeds are used, filenames are generated as "
            "<root>_seed<id>.obj."
        ),
    )
    parser.add_argument(
        "-t",
        "--output_edge_types_txt",
        type=str,
        default=None,
        help=(
            "Output txt for edge types (one int per OBJ 'l' edge line). For multiple seeds, "
            "this is ignored and per-seed files are generated next to each OBJ."
        ),
    )
    parser.add_argument(
        "--edge_type_array",
        type=str,
        default="EdgeValues",
        help="Cell data array name storing per-edge values (default: EdgeValues).",
    )
    parser.add_argument(
        "--w_array",
        type=str,
        default="t",
        help="Point data array name storing the 4th coordinate (default: t).",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="bfs",
        choices=["bfs", "dsu"],
        help="Connected-component extraction method (default: bfs; faster for large graphs).",
    )
    parser.add_argument(
        "--output_point_t_txt",
        type=str,
        default=None,
        help=(
            "Optional txt file: old_vtk_point_id -> t value. For multiple seeds, per-seed "
            "files are generated as <root>_seed<id>_point_t.txt."
        ),
    )

    args = parser.parse_args()

    # Resolve seeds in priority order: explicit list -> preset --data -> single seed
    if args.seed_point_ids is not None:
        seed_ids = list(args.seed_point_ids)
    elif args.data is not None:
        seed_ids = list(get_seed_point_id(args.data))
    elif args.seed_point_id is not None:
        seed_ids = [args.seed_point_id]
    else:
        raise ValueError("Provide --seed_point_id, --seed_point_ids, or --data.")

    base_root, base_ext = os.path.splitext(args.output_obj)
    if not base_ext:
        base_ext = ".obj"

    def _obj_out_for_seed(sid: int) -> str:
        if len(seed_ids) == 1:
            return args.output_obj
        return f"{base_root}_seed{sid}{base_ext}"

    def _edge_types_out_for_obj(obj_path: str) -> str:
        if len(seed_ids) == 1 and args.output_edge_types_txt is not None:
            return args.output_edge_types_txt
        r, _e = os.path.splitext(obj_path)
        return r + "_edge_types.txt"

    def _point_t_out_for_obj(obj_path: str) -> Optional[str]:
        if args.output_point_t_txt is None:
            return None
        if len(seed_ids) == 1:
            return args.output_point_t_txt
        r, _e = os.path.splitext(obj_path)
        return r + "_point_t.txt"

    for sid in seed_ids:
        obj_out = _obj_out_for_seed(sid)
        extract_component_and_write_obj(
            input_vtp=args.input_vtp,
            seed_point_id=sid,
            seed_point_id_base=args.seed_point_id_base,
            output_obj=obj_out,
            output_edge_types_txt=_edge_types_out_for_obj(obj_out),
            edge_type_array=args.edge_type_array,
            w_array=args.w_array,
            method=args.method,
            output_point_t_txt=_point_t_out_for_obj(obj_out),
        )


if __name__ == "__main__":
    main()