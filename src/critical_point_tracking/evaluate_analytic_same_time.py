#!/usr/bin/env python3
"""Spatial vertex-to-analytic-reference errors at exactly the recovered time.

Evaluate the selected ratio-4 graphs for quartic_potential_2,
rotating_quartic_multiwell and quartic_potential_3d. Each used OBJ vertex has
equal weight, including degenerate vertices. Quartic vertices use the closest
spatial root at the same t; rotating vertices retain the previously established
connected-component correspondence. No curve sampling, time adjustment, or
segment-interior comparison is performed. This is not graph Hausdorff distance.
"""

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import sys

import numpy as np

from evaluate_analytic_geometry import fingerprint, write_csv
from evaluate_rotating_components import match_components


SPATIAL_DIMENSIONS = {'quartic_potential_2': 2, 'rotating_quartic_multiwell': 2,
                      'quartic_potential_3d': 3}


def read_vertices(path, spatial_dimensions):
    """Read the tracker's OBJ convention: spatial coordinates followed by t."""
    points, edges = [], []
    for line in Path(path).read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith('#'):
            continue
        if parts[0] == 'v':
            if len(parts) != spatial_dimensions+2:
                raise ValueError(f'{path}: wrong number of vertex coordinates')
            points.append(list(map(float, parts[1:])))
        elif parts[0] == 'l':
            ids = [int(p.split('/')[0]) for p in parts[1:]]
            if 0 in ids:
                raise ValueError(f'{path}: OBJ vertex indices cannot be zero')
            ids = [i-1 if i > 0 else len(points)+i for i in ids]
            edges.extend(zip(ids[:-1], ids[1:]))
        elif parts[0] not in ('o', 'g', 's', 'usemtl', 'mtllib'):
            raise ValueError(f'{path}: unsupported OBJ record {parts[0]}')
    points, edges = np.asarray(points, dtype=float), np.asarray(edges, dtype=np.int64)
    if not len(edges) or edges.min() < 0 or edges.max() >= len(points):
        raise ValueError(f'{path}: invalid or empty graph connectivity')
    if not np.isfinite(points).all():
        raise ValueError(f'{path}: nonfinite vertex coordinates')
    return points, edges


def same_time_roots(name, t, component=None):
    """Enumerate analytic spatial roots; the final coordinate t is held fixed."""
    if name not in SPATIAL_DIMENSIONS:
        raise ValueError(f'Unsupported analytic function: {name}')
    if not math.isfinite(t) or not 0 <= t <= 4:
        raise ValueError(f'Vertex time {t} is outside the analytic domain [0,4]')
    if name == 'rotating_quartic_multiwell':
        if component is None or tuple(component) not in itertools.product((-1, 0, 1), repeat=2):
            raise ValueError('A matched rotating component (u,v) is required')
        u, v = component
        return np.array([[u*math.cos(t)-v*math.sin(t), u*math.sin(t)+v*math.cos(t)]])
    x_roots = (0.,)
    if t > 1:
        x = math.sqrt(t-1)
        x_roots = (-x, 0., x)
    yz_roots = (0.,)
    squared_yz = -math.cos(t)
    if squared_yz > 0:
        yz = math.sqrt(squared_yz)
        yz_roots = (-yz, 0., yz)
    axes = [x_roots]+[yz_roots]*(SPATIAL_DIMENSIONS[name]-1)
    return np.array(list(itertools.product(*axes)), dtype=float)


def nearest_same_time(point, name, component=None):
    point = np.asarray(point, dtype=float)
    if point.shape != (SPATIAL_DIMENSIONS[name]+1,) or not np.isfinite(point).all():
        raise ValueError('Invalid vertex coordinates')
    roots = same_time_roots(name, float(point[-1]), component)
    distances = np.linalg.norm(roots-point[:-1], axis=1)
    best = int(np.argmin(distances))
    return float(distances[best]), np.r_[roots[best], point[-1]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/comparison_summary.csv'))
    parser.add_argument('--output-dir', type=Path,
                        default=Path('build/analytic_geometry_validation/ratio_4/same_time_vertex_errors'))
    args = parser.parse_args()
    with args.selection.open(newline='') as stream:
        inputs = list(csv.DictReader(stream))
    if (len(inputs) != 3 or {r['function'] for r in inputs} != set(SPATIAL_DIMENSIONS)
            or any(int(r['ratio']) != 4 for r in inputs)):
        parser.error('--selection must contain exactly the three analytic functions at ratio 4')
    summary, records, sources, component_summary = [], [], [], []
    for selected in inputs:
        name, path = selected['function'], Path(selected['graph_file'])
        dimension = SPATIAL_DIMENSIONS[name]
        points, edges = read_vertices(path, dimension)
        used = np.unique(edges)
        if name == 'rotating_quartic_multiwell':
            groups = [(m['component_id'], m['vertices'], (m['u'], m['v']))
                      for m in match_components(points, edges)]
        else:
            groups = [(-1, used, None)]
        errors = []
        for component_id, ids, label in groups:
            group_errors = []
            for vertex in ids:
                point = points[vertex]
                distance, reference = nearest_same_time(point, name, label)
                errors.append(distance)
                group_errors.append(distance)
                records.append(dict(function=name, ratio=4, component_id=component_id,
                    vertex_id=int(vertex), x=point[0], y=point[1], z=point[2] if dimension == 3 else '',
                    t=point[-1], reference_x=reference[0], reference_y=reference[1],
                    reference_z=reference[2] if dimension == 3 else '', reference_t=reference[-1],
                    spatial_distance=distance, time_difference=float(reference[-1]-point[-1])))
            if label is not None:
                component_summary.append(dict(function=name, component_id=component_id,
                    u=label[0], v=label[1], vertices=len(ids),
                    mean_spatial_distance=math.fsum(group_errors)/len(group_errors),
                    max_spatial_distance=max(group_errors)))
        row = dict(function=name, ratio=4, spatial_dimensions=dimension, vertices=len(used),
                   mean_spatial_distance=math.fsum(errors)/len(errors), max_spatial_distance=max(errors),
                   graph_file=str(path.resolve()))
        summary.append(row)
        sources.append(fingerprint(path))
        print(json.dumps(row), flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir/'summary.csv', summary)
    write_csv(args.output_dir/'vertex_errors.csv', records)
    write_csv(args.output_dir/'rotating_component_errors.csv', component_summary)
    metadata = dict(command=sys.argv,
        metric='Spatial Euclidean distance from each recovered vertex to the closest analytic root at exactly the same t',
        time_constraint='reference_t == vertex_t; no clamping, temporal shift, or time resampling',
        matching='Quartic: all real spatial roots at t. Rotating: matched connected component (u,v).',
        weighting='Equal weight per used OBJ vertex, including degenerate vertices; each vertex ID counted once',
        exclusions='Unused OBJ vertices, segment interiors, and reverse reference-to-output coverage are not evaluated',
        numerical_method='Direct analytic root evaluation and spatial minimum in float64; no interpolated reference',
        reference_formulas=dict(quartic='x=0 or +/-sqrt(t-1); y (and z in 3D)=0 or +/-sqrt(-cos(t)); retain real roots only',
                                rotating='(x,y)=(u*cos(t)-v*sin(t),u*sin(t)+v*cos(t)) for the matched component'),
        inputs=sources, selection=fingerprint(args.selection), script=fingerprint(Path(__file__)),
        component_matcher=fingerprint(Path(__file__).with_name('evaluate_rotating_components.py')),
        max_absolute_time_difference=max(abs(r['time_difference']) for r in records), results=summary)
    (args.output_dir/'details.json').write_text(json.dumps(metadata, indent=2)+'\n')
    lines = ['Ratio 4: spatial vertex-to-analytic-reference distance at exactly the same time.',
             'Original spatial coordinates (x,y) or (x,y,z), without normalization.',
             'Mean weights every used OBJ vertex equally, including degenerate vertices.',
             'Rotating comparisons remain restricted to the corresponding connected component.',
             'This evaluates vertex spatial accuracy, not graph Hausdorff distance.', '',
             'function,vertices,mean_spatial_distance,max_spatial_distance']
    for row in summary:
        lines.append(f'{row["function"]},{row["vertices"]},{row["mean_spatial_distance"]:.12g},{row["max_spatial_distance"]:.12g}')
    (args.output_dir/'report.txt').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
