#!/usr/bin/env python3
"""Extract degenerate points of a critical-point trajectory stored in a .vtp.

The trajectory is a graph living in space-time: every point is ``(x, y, t)``
with the third coordinate ``z`` interpreted as time, and the line/polyline cells
of the .vtp give the edges that connect consecutive tracked positions.

A point is considered *degenerate* when one of the following holds:

  1. It has degree exactly 2 and BOTH of its neighbors have a time (z) that is
     larger than its own, or BOTH have a time smaller than its own.  This is a
     turning point in time (a birth / death fold of the tracked feature).
  2. It has degree greater than 2 (a junction where trajectories merge / split).
  3. It has degree exactly 1, i.e. a trajectory that starts or ends inside the
     domain (a feature appearing / disappearing away from the time boundary).

Points that sit on the temporal domain boundary (z at the global min or max
time) are excluded, because their one-sided neighborhood is an artifact of the
finite time window rather than a real degenerate event -- in particular a
degree-1 endpoint on the time boundary is just the trajectory leaving the time
window, not a genuine start / end.  Boundary proximity is controlled per
dimension via ``--boundary-threshold TX TY TZ``: a point within the given
threshold of the bounding box min or max on ANY dimension is treated as on the
boundary (a negative threshold disables that dimension).

The result is written as a CSV with three columns (x, y, z) giving the location
of every degenerate point.
"""

import argparse
import csv
import os

import vtk


def read_vtp(path):
    """Read a .vtp file and return its vtkPolyData."""
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def build_adjacency(mesh):
    """Return (points, neighbors).

    ``points`` is a list of (x, y, z) tuples indexed by VTK point id.
    ``neighbors`` is a list of sets; ``neighbors[i]`` holds the point ids that
    share an edge with point ``i``.  Polyline cells are decomposed into their
    consecutive segments so the adjacency is per-edge.
    """
    pts = mesh.GetPoints()
    n_points = mesh.GetNumberOfPoints()
    points = [tuple(pts.GetPoint(i)) for i in range(n_points)]
    neighbors = [set() for _ in range(n_points)]

    # Edges are carried both by line cells (GetLines) and, defensively, by any
    # cell that is a (poly)line.  Iterating over all cells covers both cases.
    for cid in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(cid)
        n_cell_pts = cell.GetNumberOfPoints()
        if n_cell_pts < 2:
            continue
        ids = [cell.GetPointId(k) for k in range(n_cell_pts)]
        for a, b in zip(ids[:-1], ids[1:]):
            if a == b:
                continue
            neighbors[a].add(b)
            neighbors[b].add(a)

    return points, neighbors


def compute_bounds(points):
    """Return (mins, maxs) bounding-box tuples over all points."""
    mins = [min(p[d] for p in points) for d in range(3)]
    maxs = [max(p[d] for p in points) for d in range(3)]
    return mins, maxs


def is_on_boundary(point, mins, maxs, thresholds):
    """True if ``point`` is within the per-dimension threshold of the bounding box.

    ``thresholds`` is a 3-tuple ``(tx, ty, tz)``.  A point is treated as on the
    boundary when, for ANY dimension ``d``, it lies within ``thresholds[d]`` of
    that dimension's min or max.  A negative threshold disables the check on that
    dimension.
    """
    for d in range(3):
        t = thresholds[d]
        if t < 0:
            continue
        if point[d] - mins[d] <= t or maxs[d] - point[d] <= t:
            return True
    return False


def find_degenerate_points(points, neighbors, boundary_thresholds, time_tol):
    """Return a list of (x, y, z) degenerate points following the rules above."""
    mins, maxs = compute_bounds(points)
    results = []

    for i, nbrs in enumerate(neighbors):
        degree = len(nbrs)
        if degree < 1:
            continue  # isolated points carry no trajectory information

        if is_on_boundary(points[i], mins, maxs, boundary_thresholds):
            continue

        z = points[i][2]
        degenerate = False

        if degree == 1:
            # Rule 3: trajectory start / end inside the domain.
            degenerate = True
        elif degree > 2:
            # Rule 2: junction (merge / split).
            degenerate = True
        else:
            # Rule 1: degree exactly 2 -> turning point in time when both
            # neighbors lie strictly on the same side in time.
            nz = [points[j][2] for j in nbrs]
            both_larger = all(t - z > time_tol for t in nz)
            both_smaller = all(z - t > time_tol for t in nz)
            degenerate = both_larger or both_smaller

        if degenerate:
            results.append(points[i])

    return results


def deduplicate(points, tol):
    """Collapse points that coincide within ``tol`` (keeps first occurrence)."""
    seen = {}
    unique = []
    for p in points:
        key = tuple(round(c / tol) for c in p) if tol > 0 else p
        if key not in seen:
            seen[key] = True
            unique.append(p)
    return unique


def write_csv(points, path, header):
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as fout:
        writer = csv.writer(fout)
        if header:
            writer.writerow(["x", "y", "z"])
        for p in points:
            writer.writerow([repr(float(p[0])), repr(float(p[1])), repr(float(p[2]))])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract degenerate points (time turning points and junctions) "
                    "from a critical-point trajectory .vtp and save their locations to CSV."
    )
    parser.add_argument("-i", "--input", required=True, help="input .vtp trajectory file")
    parser.add_argument("-o", "--output", required=True, help="output .csv file (x, y, z)")
    parser.add_argument(
        "--boundary-threshold", type=float, nargs=3, metavar=("TX", "TY", "TZ"),
        default=[-1.0, -1.0, 1e-6],
        help="per-dimension distance thresholds (x y z) to the bounding box; a "
             "point within the threshold of the min or max on ANY dimension is "
             "treated as on the boundary and excluded. A negative threshold "
             "disables that dimension. Default: -1 -1 1e-6 (time boundary only)."
    )
    parser.add_argument(
        "--time-tol", type=float, default=0.0,
        help="time difference below which a neighbor is treated as level (default 0)."
    )
    parser.add_argument(
        "--dedup-tol", type=float, default=1e-9,
        help="coordinate tolerance for deduplicating output points (0 disables)."
    )
    parser.add_argument(
        "--no-header", action="store_true", help="do not write the x,y,z CSV header."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    mesh = read_vtp(args.input)
    points, neighbors = build_adjacency(mesh)
    print(f"read {len(points)} points and "
          f"{mesh.GetNumberOfCells()} cells from {args.input}")

    degenerate = find_degenerate_points(
        points, neighbors,
        boundary_thresholds=args.boundary_threshold,
        time_tol=args.time_tol,
    )
    degenerate = deduplicate(degenerate, args.dedup_tol)

    write_csv(degenerate, args.output, header=not args.no_header)
    print(f"wrote {len(degenerate)} degenerate point(s) to {args.output}")


if __name__ == "__main__":
    main()
