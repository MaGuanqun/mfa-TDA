#!/usr/bin/env python3
"""Compare saved analytic tracking graphs with exact critical-point curves.

Requires NumPy, SciPy and VTK. Run from the repository root, for example:
  python src/critical_point_tracking/evaluate_analytic_geometry.py

Defaults to the current ori_NAME_4.obj outputs and matching event files.
Use --graph-format vtp or --ratios to evaluate other saved outputs explicitly.

Distances are Euclidean in the original (x,y,t) coordinates. The reference
domain is [-2,2] x [-2,2] x [0,4], without clipping to the recovered time range.
The Hausdorff calculation includes line interiors in BOTH directions. Analytic
curves are approximated with a second-derivative interpolation error bound;
adaptive segment subdivision bounds the remaining maximization error. Bounds
are mathematical bounds in exact arithmetic, with ordinary float64 rounding
in the implementation (not interval arithmetic).
"""

import argparse
import csv
import hashlib
import heapq
import json
import math
from pathlib import Path
import struct
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy


FUNCTIONS = ('quartic_potential_2', 'rotating_quartic_multiwell')
RATIOS = (4,)


def ground_truth_events(name):
    if name == 'rotating_quartic_multiwell':
        return np.empty((0, 3))
    if name != 'quartic_potential_2':
        raise ValueError(name)
    t = math.pi / 2
    x = math.sqrt(t - 1)
    return np.array([[0, 0, 1], [-x, 0, t], [0, 0, t], [x, 0, t]])


def read_events(path):
    """Read utility::writeMatrixVector's int32 / size_t64 / Eigen float64 file."""
    raw = Path(path).read_bytes()
    if len(raw) < 20 or struct.unpack_from('<i', raw)[0] != 1:
        raise ValueError(f'{path}: expected one matrix')
    rows, cols = struct.unpack_from('<QQ', raw, 4)
    if cols != 3 or len(raw) != 20 + rows * cols * 8:
        raise ValueError(f'{path}: invalid event matrix dimensions or length')
    points = np.frombuffer(raw, dtype='<f8', offset=20).reshape((rows, cols), order='F').copy()
    if not np.isfinite(points).all():
        raise ValueError(f'{path}: nonfinite event coordinates')
    return points


def point_segment_distances(point, starts, ends):
    delta = ends - starts
    length2 = np.einsum('...i,...i->...', delta, delta)
    fraction = np.divide(np.sum((point - starts) * delta, axis=-1), length2,
                         out=np.zeros_like(length2), where=length2 > 0)
    nearest = starts + np.clip(fraction, 0, 1)[..., None] * delta
    return np.linalg.norm(point - nearest, axis=-1)


class SegmentGraph:
    def __init__(self, starts, ends):
        self.starts = np.ascontiguousarray(starts, dtype=np.float64).reshape((-1, 3))
        self.ends = np.ascontiguousarray(ends, dtype=np.float64).reshape((-1, 3))
        if not len(self.starts) or self.starts.shape != self.ends.shape:
            raise ValueError('Expected a nonempty set of 3D line segments')
        if not np.isfinite(self.starts).all() or not np.isfinite(self.ends).all():
            raise ValueError('Nonfinite graph coordinates')
        points = np.stack((self.starts, self.ends), axis=1).reshape((-1, 3))
        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_to_vtk(points, deep=True))
        cells = vtk.vtkCellArray()
        cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, len(points) + 1, 2), deep=True),
                      numpy_to_vtkIdTypeArray(np.arange(len(points)), deep=True))
        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(vtk_points)
        self.poly.SetLines(cells)
        self.locator = vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(self.poly)
        self.locator.BuildLocator()
        self.closest = [0., 0., 0.]
        self.cell = vtk.reference(0)
        self.sub = vtk.reference(0)
        self.distance2 = vtk.reference(0.)
        self.queries = 0

    def nearest(self, point):
        self.locator.FindClosestPoint(point, self.closest, self.cell, self.sub, self.distance2)
        self.queries += 1
        return math.sqrt(max(0., float(self.distance2))), int(self.cell)

    def nearest_many(self, points):
        distances = np.empty(len(points))
        cells = np.empty(len(points), dtype=np.int64)
        for i, point in enumerate(points):
            distances[i], cells[i] = self.nearest(point)
        return distances, cells


def read_graph(path):
    if Path(path).suffix == '.obj':
        vertices, edges = [], []
        for line in Path(path).read_text().splitlines():
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            if parts[0] == 'v':
                vertices.append(list(map(float, parts[1:4])))
            elif parts[0] == 'l':
                ids = [int(p.split('/')[0]) for p in parts[1:]]
                ids = [i-1 if i > 0 else len(vertices)+i for i in ids]
                edges.extend(zip(ids[:-1], ids[1:]))
            elif parts[0] not in ('o', 'g', 's', 'usemtl', 'mtllib'):
                raise ValueError(f'{path}: unsupported OBJ record {parts[0]}')
        points, edges = np.asarray(vertices, dtype=float), np.asarray(edges, dtype=np.int64)
        if not len(edges) or edges.min() < 0 or edges.max() >= len(points):
            raise ValueError(f'{path}: invalid or empty graph connectivity')
        used = np.unique(edges)
        graph = SegmentGraph(points[edges[:, 0]], points[edges[:, 1]])
        return graph, dict(points=len(points), segments=len(edges), unused_points=len(points)-len(used),
                           coordinate_type='float64 parsed from OBJ',
                           bounds_min=points[used].min(axis=0).tolist(),
                           bounds_max=points[used].max(axis=0).tolist())
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    poly = reader.GetOutput()
    if poly.GetNumberOfPolys() or poly.GetNumberOfStrips() or not poly.GetNumberOfPoints():
        raise ValueError(f'{path}: expected points and lines only')
    points = vtk_to_numpy(poly.GetPoints().GetData()).astype(np.float64)
    starts, ends, used = [], [], set()
    for i in range(poly.GetNumberOfCells()):
        ids = poly.GetCell(i).GetPointIds()
        indices = [ids.GetId(j) for j in range(ids.GetNumberOfIds())]
        used.update(indices)
        if len(indices) == 1:
            starts.append(points[indices[0]])
            ends.append(points[indices[0]])
        else:
            for a, b in zip(indices[:-1], indices[1:]):
                starts.append(points[a])
                ends.append(points[b])
    graph = SegmentGraph(starts, ends)
    info = dict(points=len(points), segments=len(starts), unused_points=len(points)-len(used),
                coordinate_type=poly.GetPoints().GetData().GetDataTypeAsString(),
                bounds_min=points[list(used)].min(axis=0).tolist(),
                bounds_max=points[list(used)].max(axis=0).tolist())
    return graph, info


def analytic_reference(name, error):
    """Linear interpolation has error <= M * du**2 / 8 for ||curve''|| <= M.

    Quartic birth curves use square-root time parameters, removing the infinite
    time derivative at birth. For y != 0 let s=sqrt(t-pi/2), b=4-pi/2 and
    g(s)=sinc(s*s)=integral_0^1 cos(s*s*u) du. Then y=+/-s*sqrt(g),
    g>=sin(b)/b=m>0, |g'|<=s, |g''|<=1+4*s*s/3. These inequalities bound y''
    throughout [0,sqrt(b)]. If x != 0, |x''|<=1/sqrt(pi/2-1); t''=2.
    """
    branches, descriptions = [], []

    def append(label, lo, hi, second_derivative_bound, evaluate):
        count = max(1, math.ceil((hi-lo) * math.sqrt(second_derivative_bound/(8*error))))
        u = np.linspace(lo, hi, count+1)
        branch = evaluate(u)
        branches.append(branch)
        descriptions.append(dict(label=label, segments=count,
                                 second_derivative_bound=second_derivative_bound,
                                 interpolation_error_bound=second_derivative_bound*((hi-lo)/count)**2/8))

    if name == 'rotating_quartic_multiwell':
        for u in (-1, 0, 1):
            for v in (-1, 0, 1):
                append(f'u={u},v={v}', 0, 4, math.hypot(u, v),
                       lambda t, u=u, v=v: np.column_stack((u*np.cos(t)-v*np.sin(t),
                                                            u*np.sin(t)+v*np.cos(t), t)))
    elif name == 'quartic_potential_2':
        append('x=0,y=0', 0, 4, 0, lambda t: np.column_stack((t*0, t*0, t)))
        for sx in (-1, 1):
            append(f'x_sign={sx},y=0', 0, math.sqrt(3), 2,
                   lambda s, sx=sx: np.column_stack((sx*s, s*0, 1+s*s)))
        a, b = math.pi/2-1, 4-math.pi/2
        smax, m = math.sqrt(b), math.sin(b)/b
        y2 = smax/math.sqrt(m) + smax*(1+4*b/3)/(2*math.sqrt(m)) + smax**3/(4*m**1.5)
        for sx in (-1, 0, 1):
            for sy in (-1, 1):
                bound = math.sqrt(y2*y2 + 4 + (1/a if sx else 0))
                append(f'x_sign={sx},y_sign={sy}', 0, smax, bound,
                       lambda s, sx=sx, sy=sy: np.column_stack((sx*np.sqrt(a+s*s),
                           sy*s*np.sqrt(np.sinc(s*s/np.pi)), math.pi/2+s*s)))
    else:
        raise ValueError(name)
    graph = SegmentGraph(np.concatenate([p[:-1] for p in branches]),
                         np.concatenate([p[1:] for p in branches]))
    return graph, descriptions


def interval_upper(a, b, da, db, ca, cb, target):
    # Distance to a fixed segment is convex along [a,b]. Thus the larger of
    # its two endpoint distances bounds the full interval's distance to the
    # target union. Take the best such bound from both nearest target segments.
    cells = np.array([ca, cb])
    fixed_segment = np.maximum(point_segment_distances(a, target.starts[cells], target.ends[cells]),
                               point_segment_distances(b, target.starts[cells], target.ends[cells])).min()
    lipschitz = (da + db + np.linalg.norm(b-a))/2
    return max(da, db, min(float(fixed_segment), float(lipschitz)))


def directed_hausdorff(source, target, tolerance, global_lower=0.):
    """Bound max_{p in source segments} distance(p, target segments).

    An optional lower bound from the other direction avoids refining intervals
    that cannot determine the symmetric distance. In that case this direction's
    own interval can be wider than tolerance; the symmetric interval is not.
    """
    endpoints, inverse = np.unique(np.concatenate((source.starts, source.ends)), axis=0,
                                   return_inverse=True)
    distances, cells = target.nearest_many(endpoints)
    best = int(np.argmax(distances))
    lower, witness = float(distances[best]), endpoints[best].copy()
    n = len(source.starts)
    left, right = inverse[:n], inverse[n:]
    # Vectorize the initial upper bounds; only unresolved intervals enter heap.
    ca, cb = cells[left], cells[right]
    bound_a = np.maximum(point_segment_distances(source.starts, target.starts[ca], target.ends[ca]),
                         point_segment_distances(source.ends, target.starts[ca], target.ends[ca]))
    bound_b = np.maximum(point_segment_distances(source.starts, target.starts[cb], target.ends[cb]),
                         point_segment_distances(source.ends, target.starts[cb], target.ends[cb]))
    upper = np.minimum(np.minimum(bound_a, bound_b),
                       (distances[left]+distances[right]+np.linalg.norm(source.ends-source.starts, axis=1))/2)
    upper = np.maximum(upper, np.maximum(distances[left], distances[right]))
    queue, serial = [], 0
    for i in np.flatnonzero(upper > max(lower, global_lower) + tolerance):
        heapq.heappush(queue, (-float(upper[i]), serial, source.starts[i], source.ends[i],
                              distances[left[i]], distances[right[i]], ca[i], cb[i]))
        serial += 1
    iterations = 0
    while queue and -queue[0][0] > max(lower, global_lower) + tolerance:
        _, _, a, b, da, db, ca, cb = heapq.heappop(queue)
        mid = (a+b)/2
        dm, cm = target.nearest(mid)
        if dm > lower:
            lower, witness = dm, mid.copy()
        for p, q, dp, dq, cp, cq in ((a, mid, da, dm, ca, cm), (mid, b, dm, db, cm, cb)):
            ub = interval_upper(p, q, dp, dq, cp, cq, target)
            if ub > max(lower, global_lower) + tolerance:
                heapq.heappush(queue, (-ub, serial, p, q, dp, dq, cp, cq))
                serial += 1
        iterations += 1
        if iterations > 5_000_000:
            raise RuntimeError('Hausdorff refinement limit exceeded')
    # All omitted intervals lie below the global cutoff plus tolerance.
    return dict(lower=lower, upper=max(lower, global_lower)+tolerance,
                witness=witness.tolist(), subdivisions=iterations, global_lower_cutoff=global_lower)


def match_events(detected, reference):
    """One-to-one minimum-total Euclidean assignment; retain unmatched counts."""
    n, m = len(detected), len(reference)
    summary = dict(detected_events=n, reference_events=m, matched_events=min(n, m),
                   unmatched_detected_events=max(0, n-m), unmatched_reference_events=max(0, m-n),
                   event_max_distance=None, event_mean_distance=None)
    if not n or not m:
        return summary, []
    distances = np.linalg.norm(detected[:, None, :] - reference[None, :, :], axis=2)
    rows, cols = linear_sum_assignment(distances)
    errors = distances[rows, cols]
    summary.update(event_max_distance=float(errors.max()), event_mean_distance=float(errors.mean()))
    details = []
    for row, col, error in zip(rows, cols, errors):
        details.append(dict(detected_index=int(row), reference_index=int(col),
                            detected_x=detected[row, 0], detected_y=detected[row, 1], detected_t=detected[row, 2],
                            reference_x=reference[col, 0], reference_y=reference[col, 1], reference_t=reference[col, 2],
                            distance=float(error), spatial_error=float(np.linalg.norm(detected[row, :2]-reference[col, :2])),
                            time_error=float(abs(detected[row, 2]-reference[col, 2]))))
    return summary, details


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path):
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-root', type=Path, default=Path('build/src'))
    parser.add_argument('--output-dir', type=Path, default=Path('build/analytic_geometry_validation/ratio_4'))
    parser.add_argument('--graph-format', choices=('vtp', 'obj'), default='obj',
                        help='vtp: NAME_original_RATIO.vtp; obj: ori_NAME_RATIO.obj')
    parser.add_argument('--functions', nargs='+', choices=FUNCTIONS, default=list(FUNCTIONS))
    parser.add_argument('--ratios', nargs='+', type=int, default=list(RATIOS))
    parser.add_argument('--tolerance', type=float, default=1e-8,
                        help='Maximum width of final Hausdorff interval in original coordinate units')
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error('--tolerance must be positive and finite')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    script_info = fingerprint(Path(__file__))
    rows, event_rows, details, references = [], [], [], {}
    for name in args.functions:
        reference, branches = analytic_reference(name, args.tolerance/4)
        truth = ground_truth_events(name)
        references[name] = dict(branches=branches, events=truth.tolist())
        for ratio in args.ratios:
            folder = args.input_root/name
            path = folder/(f'{name}_original_{ratio}.vtp' if args.graph_format == 'vtp'
                           else f'ori_{name}_{ratio}.obj')
            recovered, info = read_graph(path)
            backward = directed_hausdorff(reference, recovered, args.tolerance/2)
            forward = directed_hausdorff(recovered, reference, args.tolerance/2, backward['lower'])
            lower = max(0., max(forward['lower'], backward['lower'])-args.tolerance/4)
            upper = max(forward['upper'], backward['upper'])+args.tolerance/4
            event_path = folder/f'{name}_degenerate_{ratio}.dat'
            if len(truth) and not event_path.exists():
                raise FileNotFoundError(event_path)
            detected = read_events(event_path) if event_path.exists() else np.empty((0, 3))
            event_summary, matches = match_events(detected, truth)
            # Missing files for event-free functions are not evidence of detector output.
            event_status = 'file_evaluated' if event_path.exists() else 'not_applicable_no_reference_events_no_event_file'
            graph_event_distances = recovered.nearest_many(detected)[0] if len(detected) else []
            graph_vertices = np.concatenate((recovered.starts, recovered.ends))
            vertex_distances = np.linalg.norm(detected[:, None, :]-graph_vertices[None, :, :], axis=2).min(axis=1) if len(detected) else []
            vertex_discrepancy = max(vertex_distances) if len(detected) else None
            if vertex_discrepancy is not None and vertex_discrepancy > 1e-6:
                event_status = 'file_evaluated_graph_vertices_differ'
            row = dict(function=name, step_ratio=ratio, nominal_step=0.4/ratio, graph_file=str(path),
                       hausdorff_estimate=(lower+upper)/2, hausdorff_lower=lower, hausdorff_upper=upper,
                       **event_summary, event_file_status=event_status,
                       event_file_max_distance_to_graph_vertex=vertex_discrepancy,
                       recovered_time_min=info['bounds_min'][2], recovered_time_max=info['bounds_max'][2])
            rows.append(row)
            event_rows.extend(dict(function=name, step_ratio=ratio, **match) for match in matches)
            details.append(dict(**row, graph=fingerprint(path), graph_info=info,
                                event_file=fingerprint(event_path) if event_path.exists() else None,
                                detected_event_distance_to_graph=list(map(float, graph_event_distances)),
                                recovered_to_reference_polyline=forward, reference_polyline_to_recovered=backward))
            write_csv(args.output_dir/'summary.csv', rows)
            write_csv(args.output_dir/'event_matches.csv', event_rows)
            print(f'{name} ratio={ratio}: H in [{lower:.10g}, {upper:.10g}], '
                  f'event max={event_summary["event_max_distance"]}, mean={event_summary["event_mean_distance"]}', flush=True)
    metadata = dict(command=sys.argv, metric='Euclidean (x,y,t), original coordinates, no rescaling',
                    time_domain=[0, 4], tolerance=args.tolerance, reference_interpolation_bound=args.tolerance/4,
                    directed_polyline_maximization_tolerance=args.tolerance/2,
                    arithmetic='float64; geometric error bounds do not include floating-point rounding',
                    reference_definitions=references, results=details,
                    versions=dict(numpy=np.__version__, vtk=vtk.vtkVersion.GetVTKVersion()),
                    script=script_info)
    (args.output_dir/'details.json').write_text(json.dumps(metadata, indent=2)+'\n')
    report = [
        'Analytic tracking geometry evaluation',
        f'Inputs: saved {args.graph_format} graphs; original Euclidean (x,y,t) coordinates.',
        'Entire exact reference on t=[0,4]; recovered graphs are not extended or clipped.',
        'Both line interiors and both Hausdorff directions are evaluated.',
        f'Hausdorff is reported as a bounded numerical estimate; interval width <= {args.tolerance:g}.',
        'Bounds include analytic curve interpolation and polyline maximization errors, but not float64 rounding.',
        'Input coordinate export error is retained as part of the saved graph error.',
        'Degenerate points: one-to-one minimum-total-distance assignment to the four exact quartic events.',
        'Mean/max use all matched events; unmatched counts are recorded separately.',
        'Event files and graph vertices are cross-checked; mismatches > 1e-6 are flagged in summary.csv.',
        'Rotating multiwell has no exact degenerate points, so event distances are N/A.',
        'An absent rotating event file does not establish that the event detector found no false positives.',
        'Ratios label saved runs; nominal step=0.4/ratio follows the explicit tracker configuration.',
        'The separate no_newton_correction runs are not included.',
        '', 'function,ratio,H_estimate,H_lower,H_upper,event_max,event_mean',
    ]
    for row in rows:
        report.append(','.join(str(row[key]) for key in ('function', 'step_ratio', 'hausdorff_estimate',
                                                        'hausdorff_lower', 'hausdorff_upper',
                                                        'event_max_distance', 'event_mean_distance')))
    (args.output_dir/'report.txt').write_text('\n'.join(report)+'\n')


if __name__ == '__main__':
    main()
