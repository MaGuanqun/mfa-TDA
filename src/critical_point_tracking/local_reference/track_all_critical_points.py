"""Track every saved TTK scalar-field critical point with no persistence filtering.

This is a direct, type-preserving assignment tracker, not TTKTrackingFromFields.
Adjacent frames are matched by minimum total squared XY displacement, allowing
birth/death. Unmatched detections remain explicit vertex cells. No point is
discarded because of persistence, matching cost, or trajectory length.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import random
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'python'))
from seperate_diffferent_branches import compute_connectivity_regions, vtk_read_polydata, vtk_write_polydata
from seperate_branch_assign_color2 import assign_color_ids_from_reference


def values(poly, name):
    array = poly.GetPointData().GetArray(name)
    if array is None or array.GetNumberOfTuples() != poly.GetNumberOfPoints():
        raise ValueError(f'Missing or malformed point array: {name}')
    return vtk_to_numpy(array).reshape(-1)


def coordinates(poly):
    return (vtk_to_numpy(poly.GetPoints().GetData()) if poly.GetNumberOfPoints()
            else np.empty((0, 3), dtype=np.float64))


def add_array(data, name, array):
    result = numpy_to_vtk(np.asarray(array), deep=True)
    result.SetName(name)
    data.AddArray(result)


def select_points(poly, indices):
    """Copy original point coordinates and attributes, without any old cells."""
    result = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(coordinates(poly)[indices], deep=True))
    result.SetPoints(points)
    result.GetPointData().CopyAllocate(poly.GetPointData(), len(indices))
    for new_id, old_id in enumerate(indices):
        result.GetPointData().CopyData(poly.GetPointData(), int(old_id), new_id)
    result.GetFieldData().ShallowCopy(poly.GetFieldData())
    return result


def assign_frame(left, right, max_displacement):
    """Allow a match only if cheaper than a birth+death (each costs 0.5)."""
    n, m = len(left), len(right)
    if not n or not m:
        return []
    distance2 = cdist(left[:, :2], right[:, :2], metric='sqeuclidean') / max_displacement**2
    costs = np.full((n+m, n+m), np.inf)
    costs[:n, :m] = np.where(distance2 <= 1., distance2, np.inf)
    costs[np.arange(n), m+np.arange(n)] = .5
    costs[n+np.arange(m), np.arange(m)] = .5
    costs[n:, m:] = 0.
    rows, columns = linear_sum_assignment(costs)
    return [(int(a), int(b)) for a, b in zip(rows, columns) if a < n and b < m]


def track_points(source, max_displacement):
    if not np.isfinite(max_displacement) or max_displacement <= 0:
        raise ValueError('Maximum displacement must be finite and positive')
    xyz = coordinates(source)
    times = values(source, 'TimeIndex')
    kinds = values(source, 'CriticalType')
    vertices = values(source, 'VertexId')
    if (not np.isfinite(xyz).all() or not np.isfinite(times).all()
            or np.any(times < 0) or np.any(times != np.floor(times))):
        raise ValueError('Invalid point coordinates or time indices')
    if len(set(zip(times, vertices))) != len(times):
        raise ValueError('Duplicate detection at the same time/grid vertex')
    groups = defaultdict(list)
    for point, (time, kind) in enumerate(zip(times, kinds)):
        groups[(int(time), int(kind))].append(point)
    for group in groups.values():
        group.sort(key=lambda point: int(vertices[point]))
    edges = []
    for (time, kind), left in sorted(groups.items()):
        right = groups.get((time+1, kind), [])
        for a, b in assign_frame(xyz[left], xyz[right], max_displacement):
            edges.append((left[a], right[b]))

    incoming = np.zeros(len(xyz), dtype=np.int32)
    outgoing = np.zeros(len(xyz), dtype=np.int32)
    for a, b in edges:
        outgoing[a] += 1
        incoming[b] += 1
        if times[b]-times[a] != 1 or kinds[a] != kinds[b]:
            raise RuntimeError('Invalid cross-time/type match')
    if np.any(incoming > 1) or np.any(outgoing > 1):
        raise RuntimeError('Assignments must be one-to-one between adjacent frames')
    singletons = np.flatnonzero((incoming+outgoing) == 0)
    result = select_points(source, np.arange(len(xyz)))
    verts, lines = vtk.vtkCellArray(), vtk.vtkCellArray()
    for point in singletons:
        verts.InsertNextCell(1)
        verts.InsertCellPoint(int(point))
    for a, b in edges:
        lines.InsertNextCell(2)
        lines.InsertCellPoint(a)
        lines.InsertCellPoint(b)
    result.SetVerts(verts)
    result.SetLines(lines)
    distances = np.asarray([np.linalg.norm(xyz[b, :2]-xyz[a, :2]) for a, b in edges])
    add_array(result.GetCellData(), 'MatchDistance', np.r_[np.zeros(len(singletons)), distances])
    add_array(result.GetPointData(), 'HasPrevious', incoming)
    add_array(result.GetPointData(), 'HasNext', outgoing)
    add_array(result.GetPointData(), 'IsSingleton', ((incoming+outgoing) == 0).astype(np.int32))
    result = compute_connectivity_regions(result, 'RegionId')
    labels = values(result, 'RegionId')
    lengths = Counter(labels)
    add_array(result.GetPointData(), 'TrackLength', np.asarray([lengths[label] for label in labels], dtype=np.int32))
    verify_coverage(source, result)
    stats = {'points': len(xyz), 'segments': len(edges), 'tracks': len(lengths),
             'linked_points': len(xyz)-len(singletons), 'singleton_points': len(singletons),
             'points_by_type': {str(int(k)): int(v) for k, v in sorted(Counter(kinds).items())}}
    return result, stats


def verify_coverage(source, result):
    if not np.array_equal(coordinates(source), coordinates(result)):
        raise RuntimeError('Tracking changed input point positions or ordering')
    for index in range(source.GetPointData().GetNumberOfArrays()):
        array = source.GetPointData().GetArray(index)
        if not np.array_equal(vtk_to_numpy(array), values(result, array.GetName())):
            raise RuntimeError(f'Tracking changed input point data: {array.GetName()}')
    used = set()
    for index in range(result.GetNumberOfCells()):
        cell = result.GetCell(index)
        used.update(cell.GetPointId(j) for j in range(cell.GetNumberOfPoints()))
    if used != set(range(source.GetNumberOfPoints())):
        raise RuntimeError('Some detections are absent from output cells')


def point_keys(poly, info=None):
    if info is None:
        return set(zip(values(poly, 'TimeIndex').astype(int), values(poly, 'VertexId').astype(int),
                       values(poly, 'CriticalType').astype(int)))
    # Legacy TTK's TimeStep marks the segment start on both endpoints. Use Z.
    ijk = np.rint((coordinates(poly)-info['sample_min'])/info['spacing']).astype(int)
    return set(zip(ijk[:, 2], ijk[:, 1]*info['dimensions'][0]+ijk[:, 0],
                   values(poly, 'CriticalType').astype(int)))


def write_point_csv(poly, path):
    names = [poly.GetPointData().GetArrayName(i) for i in range(poly.GetPointData().GetNumberOfArrays())]
    arrays = [values(poly, name) for name in names]
    with Path(path).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['PositionX', 'PositionY', 'PositionZ']+names)
        for index, xyz in enumerate(coordinates(poly)):
            writer.writerow([*xyz, *(array[index] for array in arrays)])


def run(sample, reference_vtp, max_displacement, spatial_step, time_step, color_threshold):
    sample = Path(sample).resolve()
    info = json.loads(Path(str(sample)+'.json').read_text())
    input_path = Path(str(sample)+'_critical_points_all.vtp')
    if not input_path.is_file():
        raise FileNotFoundError(f'{input_path}: run --critical-points-only first')
    source = vtk_read_polydata(str(input_path))
    reference = vtk_read_polydata(str(reference_vtp))
    counts = json.loads(Path(str(sample)+'_critical_points.json').read_text())
    nt = info['dimensions'][2]
    if source.GetNumberOfPoints() != counts['raw_points'] or counts['slices'] != nt:
        raise ValueError('Critical-point reference/metadata counts disagree')
    times = values(source, 'TimeIndex').astype(int)
    if np.any(times < 0) or np.any(times >= nt):
        raise ValueError('Critical-point time indices outside sample grid')
    np.testing.assert_array_equal(coordinates(source)[:, 2], np.linspace(info['sample_min'][2], info['sample_max'][2], nt)[times])
    prefix = str(sample)+'_tracking_all_critical_points'
    raw, raw_stats = track_points(source, max_displacement)
    vtk_write_polydata(raw, prefix+'.vtp', 'binary')

    # Apply the requested six-face exclusion before matching the interior set.
    # Boundary detections must not compete with interior detections for matches.
    lo, hi = np.asarray(info['sample_min']), np.asarray(info['sample_max'])
    tolerance = np.minimum(8*np.finfo(np.float32).eps*np.maximum(1., np.maximum(abs(lo), abs(hi))),
                           np.asarray(info['spacing'])/4)
    xyz = coordinates(source)
    boundary = np.any((abs(xyz-lo) <= tolerance) | (abs(xyz-hi) <= tolerance), axis=1)
    interior_source = select_points(source, np.flatnonzero(~boundary))
    if interior_source.GetNumberOfPoints() != counts['interior_points']:
        raise ValueError('Boundary exclusion disagrees with the per-slice interior reference')
    interior, interior_stats = track_points(interior_source, max_displacement)
    vtk_write_polydata(interior, prefix+'_interior.vtp', 'binary')
    matches, _ = assign_color_ids_from_reference(
        interior, reference, 'RegionId', 'RegionId', 'ColorId', 'Time', '3d', 'mean',
        48, 32, color_threshold, random.Random(0),
        distance_scales=(spatial_step, spatial_step, time_step),
        matches_csv=prefix+'_color_matches.csv')
    for data in (interior.GetPointData(), interior.GetCellData()):
        labels = vtk_to_numpy(data.GetArray('RegionId'))
        add_array(data, 'ReferenceColorMatched',
                  np.asarray([matches[int(label)][2] for label in labels], dtype=np.int32))
    verify_coverage(interior_source, interior)
    vtk_write_polydata(interior, prefix+'_regionid_colorid.vtp', 'binary')
    write_point_csv(interior, prefix+'_points.csv')
    interior_stats['reference_color_matches'] = sum(match[2] for match in matches.values())
    interior_stats['fallback_colors'] = len(matches)-interior_stats['reference_color_matches']

    ref_keys = point_keys(source)
    legacy_path = Path(str(sample)+'_tracking.vtp')
    legacy_keys = point_keys(vtk_read_polydata(str(legacy_path)), info) if legacy_path.exists() else set()
    legacy_metadata = Path(str(sample)+'_tracking.json')
    legacy_settings = (json.loads(legacy_metadata.read_text()).get('ttk_settings', {})
                       if legacy_metadata.exists() else {})
    missing = ref_keys-legacy_keys
    interior_keys = point_keys(interior_source)
    old_by_time = Counter(key[0] for key in ref_keys & legacy_keys)
    raw_by_time = Counter(times)
    inside_by_time = Counter(values(interior_source, 'TimeIndex').astype(int))
    new_by_time = Counter(values(raw, 'TimeIndex').astype(int))
    new_inside_by_time = Counter(values(interior, 'TimeIndex').astype(int))
    with Path(prefix+'_coverage.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['TimeIndex', 'ReferencePoints', 'LegacyTTKPoints', 'NewPoints',
                         'InteriorReferencePoints', 'InteriorNewPoints'])
        for i in range(nt):
            if raw_by_time[i] != new_by_time[i] or inside_by_time[i] != new_inside_by_time[i]:
                raise RuntimeError(f'Incomplete point coverage at time slice {i}')
            writer.writerow([i, raw_by_time[i], old_by_time[i], new_by_time[i], inside_by_time[i], new_inside_by_time[i]])
    summary = {
        'sample': str(sample), 'input': str(input_path),
        'method': 'Direct per-type adjacent-frame assignment on TTK scalar-field critical points (SciPy linear_sum_assignment)',
        'persistence_filtering': False, 'minimum_track_length_filtering': False,
        'maximum_xy_displacement': max_displacement,
        'assignment_cost': 'squared XY displacement / maximum_xy_displacement^2; birth/death cost 0.5 each',
        'unmatched_points': 'retained as vertex cells', 'slices': nt,
        'interior_boundary_policy': 'Remove all six sample-box faces before matching interior detections',
        'boundary_tolerance': tolerance.tolist(),
        'raw': raw_stats, 'interior': interior_stats,
        'legacy_ttk': {'settings': legacy_settings,
                       'covered_points': len(ref_keys & legacy_keys), 'missing_points': len(missing),
                       'extra_points': len(legacy_keys-ref_keys),
                       'missing_by_type': {str(int(k)): int(v) for k, v in sorted(Counter(key[2] for key in missing).items())},
                       'covered_interior_points': len(interior_keys & legacy_keys)},
        'color_reference': str(Path(reference_vtp).resolve()),
        'color_distance_scales': [spatial_step, spatial_step, time_step],
        'color_match_threshold': color_threshold,
        'colored_vtp': prefix+'_regionid_colorid.vtp', 'coverage_csv': prefix+'_coverage.csv',
    }
    Path(prefix+'.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--reference-vtp', required=True)
    parser.add_argument('--max-displacement', type=float, required=True)
    parser.add_argument('--spatial-step-size', type=float, required=True)
    parser.add_argument('--t-step-size', type=float, required=True)
    parser.add_argument('--match-distance-threshold', type=float, default=1.)
    args = parser.parse_args()
    run(args.sample, args.reference_vtp, args.max_displacement, args.spatial_step_size,
        args.t_step_size, args.match_distance_threshold)


if __name__ == '__main__':
    main()
