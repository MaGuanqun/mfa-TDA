#!/usr/bin/env python3
"""Directed vertex-to-exact-curve errors, not symmetric graph Hausdorff.

Use the graph files in ratio_4/comparison_summary.csv. Rotating components keep
their analytic (u,v) correspondence. Quartic vertices compare with the entire
exact critical-point graph. Every used OBJ vertex receives equal weight.
Continuous nearest points are obtained by scalar stationary-point searches,
restricted to a time window bounded by an existing candidate distance. This is
a float64 numerical calculation, not a certified interval computation.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

from evaluate_analytic_geometry import fingerprint, write_csv
from evaluate_rotating_components import read_obj_connectivity, match_components


def quartic_curves():
    # Each item contains an evaluation, derivative, and time-to-parameter map.
    curves = [('x=0,y=0', 0., 4., lambda s: np.array([0., 0., s]),
               lambda s: np.array([0., 0., 1.]), lambda t: t)]
    for sx in (-1, 1):
        curves.append((f'x={sx},y=0', 1., 4.,
            lambda s, sx=sx: np.array([sx*s, 0., 1+s*s]),
            lambda s, sx=sx: np.array([sx, 0., 2*s]),
            lambda t: math.sqrt(max(0., t-1))))
    a, t0 = math.pi/2-1, math.pi/2
    for sx in (-1, 0, 1):
        for sy in (-1, 1):
            def evaluate(s, sx=sx, sy=sy):
                g = np.sinc(s*s/np.pi)
                return np.array([sx*math.sqrt(a+s*s), sy*s*math.sqrt(g), t0+s*s])
            def derivative(s, sx=sx, sy=sy):
                g = np.sinc(s*s/np.pi)
                return np.array([sx*s/math.sqrt(a+s*s), sy*math.cos(s*s)/math.sqrt(g), 2*s])
            curves.append((f'x={sx},y={sy}', t0, 4., evaluate, derivative,
                           lambda t: math.sqrt(max(0., t-t0))))
    return curves


def rotating_curve(u, v):
    def evaluate(t):
        return np.array([u*math.cos(t)-v*math.sin(t), u*math.sin(t)+v*math.cos(t), t])
    def derivative(t):
        q = evaluate(t)
        return np.array([-q[1], q[0], 1.])
    return [(f'u={u},v={v}', 0., 4., evaluate, derivative, lambda t: t)]


def nearest(point, curves, brackets):
    candidates = []
    # An explicit point on a reference curve gives an upper bound U. Any closer
    # point must have |t_reference-t_vertex| <= U in Euclidean spacetime distance.
    for label, tmin, tmax, curve, _, parameter in curves:
        q = curve(parameter(float(np.clip(point[2], tmin, tmax))))
        candidates.append((float(np.linalg.norm(point-q)), q, label))
    upper = min(c[0] for c in candidates)
    if upper == 0:
        return min(candidates, key=lambda item: item[0])
    pad = 32*np.finfo(float).eps*max(1., abs(point[2]))
    for label, tmin, tmax, curve, derivative, parameter in curves:
        lo_t, hi_t = max(tmin, point[2]-upper-pad), min(tmax, point[2]+upper+pad)
        if hi_t < lo_t:
            continue
        lo, hi = parameter(lo_t), parameter(hi_t)
        def stationary(s):
            return float(np.dot(curve(s)-point, derivative(s)))
        grid = np.linspace(lo, hi, brackets+1)
        values = [stationary(s) for s in grid]
        roots = [lo, hi]
        for i, (left, right) in enumerate(zip(grid[:-1], grid[1:])):
            if values[i] == 0:
                roots.append(left)
            if np.signbit(values[i]) != np.signbit(values[i+1]):
                roots.append(brentq(stationary, left, right,
                                    xtol=np.nextafter(0., 1.), rtol=4*np.finfo(float).eps))
        for s in roots:
            q = curve(s)
            candidates.append((float(np.linalg.norm(point-q)), q, label))
    return min(candidates, key=lambda item: item[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/comparison_summary.csv'))
    parser.add_argument('--output-dir', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/vertex_errors'))
    parser.add_argument('--quartic-input', type=Path)
    parser.add_argument('--brackets', type=int, default=32)
    args = parser.parse_args()
    if args.brackets < 2:
        parser.error('--brackets must be at least 2')
    inputs = list(csv.DictReader(args.selection.open()))
    summary, records, sources = [], [], []
    for selected in inputs:
        name = selected['function']
        path = args.quartic_input if name == 'quartic_potential_2' and args.quartic_input else Path(selected['graph_file'])
        points, edges = read_obj_connectivity(path)
        used = np.unique(edges)
        if name == 'quartic_potential_2':
            groups = [(-1, used, quartic_curves())]
        else:
            groups = [(m['component_id'], m['vertices'], rotating_curve(m['u'], m['v']))
                      for m in match_components(points, edges)]
        errors = []
        for component, ids, curves in groups:
            for vertex in ids:
                point = points[vertex]
                distance, q, label = nearest(point, curves, args.brackets)
                errors.append(distance)
                records.append(dict(function=name, component_id=component, vertex_id=int(vertex),
                    x=point[0], y=point[1], t=point[2], reference_label=label,
                    reference_x=q[0], reference_y=q[1], reference_t=q[2], distance=distance))
        row = dict(function=name, ratio=4, vertices=len(used), graph_file=str(path),
                   max_vertex_to_reference_distance=max(errors),
                   mean_vertex_to_reference_distance=math.fsum(errors)/len(errors))
        summary.append(row)
        sources.append(fingerprint(path))
        print(json.dumps(row), flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir/'summary.csv', summary)
    write_csv(args.output_dir/'vertex_errors.csv', records)
    metadata = dict(command=sys.argv, metric='Directed vertex-to-continuous-reference Euclidean (x,y,t) distance',
        exclusions='Does not evaluate line interiors or the reverse reference-to-output direction',
        weighting='Equal weight for each used OBJ vertex; rotating references are component-matched',
        numerical_method='Scalar stationary-point search in time windows restricted by an upper distance bound; float64',
        brackets=args.brackets, inputs=sources, script=fingerprint(Path(__file__)), results=summary)
    (args.output_dir/'details.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    main()
