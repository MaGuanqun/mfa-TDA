#!/usr/bin/env python3
"""Symmetric whole-graph Hausdorff distance for quartic_potential_3d, ratio 4.

Use all four original coordinates (x,y,z,t), with no normalization or time
trimming. The analytic reference covers t=[0,4]. Reuse the existing adaptive
Hausdorff bounds, but replace the three-coordinate VTK locator with an exact
point-to-segment search in four dimensions. Analytic interpolation and
maximization bounds exclude ordinary float64 rounding.
"""

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

from evaluate_analytic_geometry import (directed_hausdorff, fingerprint,
                                       point_segment_distances, write_csv)


class SegmentGraphND:
    """Nearest points on segments, using midpoint balls to prune candidates.

    If a candidate point-to-segment distance is d, any closer segment has its
    midpoint within d + half_length of the query. Within each length group,
    searching a ball of radius d + max_half_length includes every possible
    closer segment. Length groups prevent a long straight segment from making
    candidate balls unnecessarily large for the finely sampled curved pieces.
    Projection onto each candidate segment supplies the actual distance.
    """

    def __init__(self, starts, ends):
        self.starts = np.ascontiguousarray(starts, dtype=np.float64)
        self.ends = np.ascontiguousarray(ends, dtype=np.float64)
        if (self.starts.ndim != 2 or not len(self.starts)
                or self.starts.shape != self.ends.shape):
            raise ValueError('Expected matching, nonempty arrays of segment endpoints')
        if not np.isfinite(self.starts).all() or not np.isfinite(self.ends).all():
            raise ValueError('Nonfinite segment coordinates')
        self.dimension = self.starts.shape[1]
        midpoints = (self.starts+self.ends)/2
        self.tree = cKDTree(midpoints)
        half_lengths = np.linalg.norm(self.ends-self.starts, axis=1)/2
        # frexp groups positive lengths within a factor of two. Put point
        # segments in a separate group with zero search-radius enlargement.
        keys = np.frexp(half_lengths)[1]
        keys[half_lengths == 0] = np.iinfo(keys.dtype).min
        self.groups = []
        for key in np.unique(keys):
            ids = np.flatnonzero(keys == key)
            self.groups.append((cKDTree(midpoints[ids]), ids, float(half_lengths[ids].max())))
        self.scale = max(1., float(np.abs(self.starts).max()), float(np.abs(self.ends).max()))
        self.queries = 0

    def nearest_many(self, points):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != self.dimension or not np.isfinite(points).all():
            raise ValueError('Invalid query coordinates')
        distances = np.empty(len(points))
        cells = np.empty(len(points), dtype=np.int64)
        for lo in range(0, len(points), 2048):
            batch = points[lo:lo+2048]
            _, initial = self.tree.query(batch)
            upper = point_segment_distances(batch, self.starts[initial], self.ends[initial])
            # Enlarge the search radius slightly to avoid pruning at a rounded boundary.
            pad = 64*np.finfo(float).eps*(self.scale+np.abs(batch).max(axis=1)+upper)
            best_cells = initial.copy()
            for tree, group_ids, max_half_length in self.groups:
                radii = np.nextafter(upper+max_half_length+pad, np.inf)
                midpoint_distances, _ = tree.query(batch)
                active = np.flatnonzero(midpoint_distances <= radii)
                if not len(active):
                    continue
                candidates = tree.query_ball_point(batch[active], radii[active], return_sorted=False)
                counts = np.fromiter(map(len, candidates), dtype=np.int64, count=len(active))
                if np.any(counts == 0):
                    raise RuntimeError('Midpoint search unexpectedly excluded a nearby candidate')
                ids = group_ids[np.concatenate(candidates).astype(np.int64, copy=False)]
                owners = np.repeat(np.arange(len(active)), counts)
                values = point_segment_distances(batch[active[owners]], self.starts[ids], self.ends[ids])
                offsets = np.r_[0, np.cumsum(counts)[:-1]]
                minima = np.minimum.reduceat(values, offsets)
                winners = np.minimum.reduceat(
                    np.where(values == minima[owners], np.arange(len(values)), len(values)), offsets)
                improved = minima < upper[active]
                upper[active[improved]] = minima[improved]
                best_cells[active[improved]] = ids[winners[improved]]
            distances[lo:lo+len(batch)] = upper
            cells[lo:lo+len(batch)] = best_cells
        self.queries += len(points)
        return distances, cells

    def nearest(self, point):
        distances, cells = self.nearest_many(np.asarray(point, dtype=float)[None, :])
        return float(distances[0]), int(cells[0])


def read_obj_4d(path):
    """The tracker writes four-coordinate v records; the fourth value is time."""
    vertices, edges = [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith('#'):
            continue
        if parts[0] == 'v':
            if len(parts) != 5:
                raise ValueError(f'{path}: expected exactly four vertex coordinates')
            vertices.append(list(map(float, parts[1:])))
        elif parts[0] == 'l':
            ids = [int(p.split('/')[0]) for p in parts[1:]]
            ids = [i-1 if i > 0 else len(vertices)+i for i in ids]
            edges.extend(zip(ids[:-1], ids[1:]))
        elif parts[0] not in ('o', 'g', 's', 'usemtl', 'mtllib'):
            raise ValueError(f'{path}: unsupported OBJ record {parts[0]}')
    points, edges = np.asarray(vertices, dtype=float), np.asarray(edges, dtype=np.int64)
    if not len(edges) or edges.min() < 0 or edges.max() >= len(points):
        raise ValueError('Invalid or empty graph connectivity')
    used = np.unique(edges)
    graph = SegmentGraphND(points[edges[:, 0]], points[edges[:, 1]])
    return graph, dict(points=len(points), used_points=len(used), segments=len(edges),
                       unused_points=len(points)-len(used),
                       bounds_min=points[used].min(axis=0).tolist(),
                       bounds_max=points[used].max(axis=0).tolist())


def analytic_reference_4d(error):
    """All spatial roots of f on [-2,2]^3 x [0,4], with bounded interpolation.

    x=0 or +/-sqrt(t-1); independently y,z=0 or +/-sqrt(-cos(t)).
    Nonzero x starts at t=1, nonzero y/z at t=pi/2. These give 27 analytic
    trajectory pieces whose union is one connected graph.

    For births at pi/2, set t=pi/2+s^2 and y=sy*s*sqrt(g(s)), likewise z,
    where g(s)=sinc(s^2)=integral_0^1 cos(s^2*v)dv. On [0,sqrt(b)], b=4-pi/2,
    g>=m=sin(b)/b>0, |g'|<=s, |g''|<=1+4*s^2/3. Differentiating s*sqrt(g)
    yields the bound y2 below. Also |x''|<=1/sqrt(pi/2-1) and t''=2.
    Thus the norm bound includes separate contributions from BOTH y and z.
    """
    if not math.isfinite(error) or error <= 0:
        raise ValueError('Reference error must be positive and finite')
    branches, descriptions = [], []

    def append(label, lo, hi, bound, evaluate):
        count = max(1, math.ceil((hi-lo)*math.sqrt(bound/(8*error))))
        parameter = np.linspace(lo, hi, count+1)
        points = evaluate(parameter)
        branches.append(points)
        descriptions.append(dict(label=label, segments=count, second_derivative_bound=bound,
                                 interpolation_error_bound=bound*((hi-lo)/count)**2/8,
                                 t_min=float(points[0, -1]), t_max=float(points[-1, -1])))

    append('x=0,y=0,z=0', 0, 4, 0,
           lambda t: np.column_stack((t*0, t*0, t*0, t)))
    for sx in (-1, 1):
        append(f'x_sign={sx},y=0,z=0', 0, math.sqrt(3), 2,
               lambda s, sx=sx: np.column_stack((sx*s, s*0, s*0, 1+s*s)))
    a, b = math.pi/2-1, 4-math.pi/2
    smax, m = math.sqrt(b), math.sin(b)/b
    y2 = smax/math.sqrt(m) + smax*(1+4*b/3)/(2*math.sqrt(m)) + smax**3/(4*m**1.5)
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for sz in (-1, 0, 1):
                if sy == sz == 0:
                    continue
                bound = math.sqrt((sy*sy+sz*sz)*y2*y2 + 4 + (1/a if sx else 0))

                def evaluate(s, sx=sx, sy=sy, sz=sz):
                    yz = s*np.sqrt(np.sinc(s*s/np.pi))
                    return np.column_stack((sx*np.sqrt(a+s*s), sy*yz, sz*yz, math.pi/2+s*s))

                append(f'x_sign={sx},y_sign={sy},z_sign={sz}', 0, smax, bound, evaluate)
    graph = SegmentGraphND(np.concatenate([p[:-1] for p in branches]),
                           np.concatenate([p[1:] for p in branches]))
    return graph, descriptions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('build/src/quartic_potential_3d/ori_quartic_potential_3d_4.obj'))
    parser.add_argument('--output-dir', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/quartic_potential_3d/saved_metrics'))
    parser.add_argument('--tolerance', type=float, default=1e-8)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error('--tolerance must be positive and finite')
    recovered, info = read_obj_4d(args.input)
    print(f'Recovered: {info}', flush=True)
    reference, descriptions = analytic_reference_4d(args.tolerance/4)
    error = max(d['interpolation_error_bound'] for d in descriptions)
    print(f'Reference: {len(descriptions)} analytic pieces, {len(reference.starts)} segments, '
          f'interpolation bound {error:.12g}', flush=True)
    backward = directed_hausdorff(reference, recovered, args.tolerance/2)
    print(f'Reference to recovered: {backward}', flush=True)
    forward = directed_hausdorff(recovered, reference, args.tolerance/2, backward['lower'])
    print(f'Recovered to reference: {forward}', flush=True)
    lower = max(0., max(forward['lower'], backward['lower'])-error)
    upper = max(forward['upper'], backward['upper'])+error
    row = dict(function='quartic_potential_3d', ratio=4, comparison='whole_graph',
               coordinates='x,y,z,t', hausdorff_estimate=(lower+upper)/2,
               hausdorff_lower=lower, hausdorff_upper=upper, tolerance=args.tolerance,
               reference_interpolation_bound=error, reference_segments=len(reference.starts),
               recovered_vertices=info['used_points'], recovered_segments=info['segments'],
               graph_file=str(args.input.resolve()))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir/'summary.csv', [row])
    metadata = dict(command=sys.argv, input=fingerprint(args.input),
                    script=fingerprint(Path(__file__)),
                    geometry_evaluator=fingerprint(Path(__file__).with_name('evaluate_analytic_geometry.py')),
                    metric='Symmetric Hausdorff of complete graphs in original Euclidean (x,y,z,t)',
                    reference='x in {0,+/-sqrt(t-1)}, y,z independently in {0,+/-sqrt(-cos(t))}; '
                              'include nonzero roots only when real, t in [0,4]',
                    time_domain=[0, 4], coordinate_normalization='none',
                    arithmetic='float64; geometric bounds exclude roundoff',
                    locator='Exact segment projection with exhaustive midpoint-ball candidate bounds',
                    recovered=info, reference_pieces=descriptions,
                    reference_to_recovered=backward, recovered_to_reference=forward, summary=row)
    (args.output_dir/'details.json').write_text(json.dumps(metadata, indent=2)+'\n')
    report = [f'quartic_potential_3d, ratio 4: symmetric whole-graph Hausdorff = {(lower+upper)/2:.12g}',
              f'Bounds: [{lower:.16g}, {upper:.16g}] (excluding float64 rounding).',
              'Original Euclidean (x,y,z,t); all four coordinates retained, no normalization.',
              'Full analytic reference on t=[0,4]; segment interiors included in both directions.',
              f'Input graph: {args.input.resolve()}',
              f'Recovered used vertices: {info["used_points"]}; segments: {info["segments"]}.',
              f'Reference analytic pieces: {len(descriptions)}; segments: {len(reference.starts)}.',
              f'Reference interpolation bound: {error:.16g}.']
    (args.output_dir/'report.txt').write_text('\n'.join(report)+'\n')
    print('\n'.join(report), flush=True)


if __name__ == '__main__':
    main()
