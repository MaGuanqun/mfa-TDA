#!/usr/bin/env python3
"""Hausdorff errors for corresponding rotating-quartic connected components.

Components use OBJ vertex IDs and line connectivity, with no geometric welding.
Each component is retained in full, without splitting at vertices or by type.
The exact component label (u,v) is constant in the rotating coordinate system.
Compare each matched component over the FULL reference interval [0,4]; missing
end portions contribute to the error. Means give all nine components equal weight.
The default input is ratio 4. This utility does not change the quartic evaluation.
"""

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from evaluate_analytic_geometry import (SegmentGraph, analytic_reference,
    directed_hausdorff, fingerprint, write_csv)


ROOTS = np.array([(u, v) for u in (-1, 0, 1) for v in (-1, 0, 1)])


def read_obj_connectivity(path):
    points, edges = [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith('#'):
            continue
        if parts[0] == 'v':
            points.append(list(map(float, parts[1:4])))
        elif parts[0] == 'l':
            ids = [int(p.split('/')[0]) for p in parts[1:]]
            ids = [i-1 if i > 0 else len(points)+i for i in ids]
            edges.extend(zip(ids[:-1], ids[1:]))
        elif parts[0] not in ('o', 'g', 's', 'usemtl', 'mtllib'):
            raise ValueError(f'{path}: unsupported OBJ record {parts[0]}')
    points, edges = np.asarray(points, dtype=float), np.asarray(edges, dtype=np.int64)
    if not len(edges) or edges.min() < 0 or edges.max() >= len(points):
        raise ValueError('Invalid or empty graph connectivity')
    if points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('Expected finite (x,y,t) coordinates')
    return points, edges


def component_indices(points, edges):
    used = np.unique(edges)
    adjacency = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                           shape=(len(points), len(points))).tocsr()[used][:, used]
    count, labels = connected_components(adjacency, directed=False)
    vertex_labels = np.full(len(points), -1, dtype=int)
    vertex_labels[used] = labels
    return [(used[labels == i], np.flatnonzero(vertex_labels[edges[:, 0]] == i))
            for i in range(count)]


def match_components(points, edges):
    components = component_indices(points, edges)
    if len(components) != len(ROOTS):
        raise ValueError(f'Expected nine connected components, found {len(components)}; '
                         'cannot form a complete one-to-one correspondence')
    coordinates, costs = [], []
    for vertices, _ in components:
        x, y, t = points[vertices].T
        uv = np.column_stack((x*np.cos(t)+y*np.sin(t), -x*np.sin(t)+y*np.cos(t)))
        coordinates.append(uv)
        costs.append(np.mean(np.sum((uv[:, None, :]-ROOTS[None, :, :])**2, axis=2), axis=0))
    rows, columns = linear_sum_assignment(costs)
    matches = []
    for row, col in zip(rows, columns):
        uv = coordinates[row]
        nearest = np.argmin(np.sum((uv[:, None, :]-ROOTS[None, :, :])**2, axis=2), axis=1)
        if not np.all(nearest == col):
            raise ValueError(f'Component {row} does not consistently match one analytic component')
        vertices, edge_indices = components[row]
        matches.append(dict(component_id=int(row), reference_id=int(col),
                            u=int(ROOTS[col, 0]), v=int(ROOTS[col, 1]),
                            vertices=vertices, edge_indices=edge_indices,
                            invariant_max_error=float(np.linalg.norm(uv-ROOTS[col], axis=1).max())))
    return sorted(matches, key=lambda item: item['reference_id'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('build/src/rotating_quartic_multiwell/ori_rotating_quartic_multiwell_4.obj'))
    parser.add_argument('--output-dir', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/rotating_components'))
    parser.add_argument('--tolerance', type=float, default=1e-8)
    args = parser.parse_args()
    if not math.isfinite(args.tolerance) or args.tolerance <= 0:
        parser.error('--tolerance must be positive and finite')
    points, edges = read_obj_connectivity(args.input)
    matches = match_components(points, edges)
    full_reference, descriptions = analytic_reference('rotating_quartic_multiwell', args.tolerance/4)
    offsets = np.r_[0, np.cumsum([d['segments'] for d in descriptions])]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, details = [], []
    for match in matches:
        reference_id = match['reference_id']
        lo, hi = offsets[reference_id:reference_id+2]
        reference = SegmentGraph(full_reference.starts[lo:hi], full_reference.ends[lo:hi])
        component_edges = edges[match['edge_indices']]
        recovered = SegmentGraph(points[component_edges[:, 0]], points[component_edges[:, 1]])
        backward = directed_hausdorff(reference, recovered, args.tolerance/2)
        forward = directed_hausdorff(recovered, reference, args.tolerance/2, backward['lower'])
        error = descriptions[reference_id]['interpolation_error_bound']
        lower = max(0., max(forward['lower'], backward['lower'])-error)
        upper = max(forward['upper'], backward['upper'])+error
        u, v = match['u'], match['v']
        kind = 'minimum' if u and v else ('saddle' if u or v else 'maximum (stationary)')
        ts = points[match['vertices'], 2]
        row = dict(u=u, v=v, critical_type=kind, component_id=match['component_id'],
                   vertices=len(match['vertices']), segments=len(component_edges),
                   recovered_t_min=float(ts.min()), recovered_t_max=float(ts.max()),
                   reference_t_min=0., reference_t_max=4.,
                   hausdorff_estimate=(lower+upper)/2, hausdorff_lower=lower, hausdorff_upper=upper,
                   invariant_max_error=match['invariant_max_error'])
        rows.append(row)
        details.append(dict(**row, vertex_ids=match['vertices'].tolist(),
                            source_edge_indices=match['edge_indices'].tolist(),
                            recovered_to_reference=forward, reference_to_recovered=backward))
        print(f'(u,v)=({u},{v}): t=[{ts.min():.9g},{ts.max():.9g}], '
              f'H=[{lower:.10g},{upper:.10g}]', flush=True)
    aggregates = []
    for group, included in [('all_nine_connected_components', rows)]:
        aggregates.append(dict(group=group, connected_components=len(included),
            mean_component_hausdorff=float(np.mean([r['hausdorff_estimate'] for r in included])),
            mean_lower=float(np.mean([r['hausdorff_lower'] for r in included])),
            mean_upper=float(np.mean([r['hausdorff_upper'] for r in included])),
            max_component_hausdorff=max(r['hausdorff_estimate'] for r in included),
            max_lower=max(r['hausdorff_lower'] for r in included),
            max_upper=max(r['hausdorff_upper'] for r in included)))
    write_csv(args.output_dir/'component_errors.csv', rows)
    write_csv(args.output_dir/'summary.csv', aggregates)
    metadata = dict(command=sys.argv, input=fingerprint(args.input), script=fingerprint(Path(__file__)),
                    geometry_evaluator=fingerprint(Path(__file__).with_name('evaluate_analytic_geometry.py')),
                    metric='Symmetric Hausdorff in original Euclidean (x,y,t) coordinates',
                    reference='(x,y,t)=(u*cos(t)-v*sin(t),u*sin(t)+v*cos(t),t), t in [0,4]',
                    matching='OBJ connected components; bijective assignment by constant rotated (u,v); '
                             'every vertex must agree with its component label',
                    component_count=len(matches), exact_component_count=9,
                    tolerance=args.tolerance, mean_weighting='Equal weight per connected component',
                    arithmetic='float64; geometric bounds exclude roundoff',
                    aggregates=aggregates, components=details)
    (args.output_dir/'details.json').write_text(json.dumps(metadata, indent=2)+'\n')
    report = [__doc__.strip(), '', f'Input: {args.input}',
              'Original Euclidean (x,y,t), exact reference on [0,4], no trimming of time intervals.',
              f'Numerical Hausdorff bound width <= {args.tolerance:g}, excluding float64 rounding.',
              'Nine recovered connected components match nine exact connected components.',
              'The stationary center is included as one connected component.',
              'Mean component Hausdorff is the average of the nine component Hausdorff distances.',
              '', 'u,v,type,recovered_t_min,recovered_t_max,H_estimate']
    for r in rows:
        report.append(','.join(str(r[k]) for k in ('u','v','critical_type','recovered_t_min',
                                                   'recovered_t_max','hausdorff_estimate')))
    report.extend(['', json.dumps(aggregates, indent=2)])
    (args.output_dir/'report.txt').write_text('\n'.join(report)+'\n')
    print(json.dumps(aggregates, indent=2))


if __name__ == '__main__':
    main()
