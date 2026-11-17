"""Recheck reference colors on surviving connected paths after obstacle removal."""
import argparse
import csv
import json
from pathlib import Path
import random

import numpy as np
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from seperate_diffferent_branches import compute_connectivity_regions, vtk_read_polydata, vtk_write_polydata
from seperate_branch_assign_color2 import assign_color_ids_from_reference


def add_array(data, name, values):
    array = numpy_to_vtk(np.asarray(values), deep=True)
    array.SetName(name)
    data.AddArray(array)


def recolor(source, reference, scales, threshold, matches_csv, color_overrides=None):
    color_overrides = color_overrides or {}
    # RegionId remains the original trajectory label. Clipping may split it
    # into several disconnected paths that need independent color matches.
    result = compute_connectivity_regions(source, 'ColorMatchComponentId')
    before = {}
    for association, data in [('point', source.GetPointData()), ('cell', source.GetCellData())]:
        before[association] = {
            name: vtk_to_numpy(data.GetArray(name)).copy()
            for name in ('ColorId', 'ReferenceColorMatched', 'RegionId')}
    assign_color_ids_from_reference(
        result, reference, 'ColorMatchComponentId', 'RegionId', 'ColorId', 'Time',
        '3d', 'mean', 48, 32, threshold, random.Random(0),
        distance_scales=tuple(scales), matches_csv=str(matches_csv))
    with Path(matches_csv).open() as stream:
        rows = list(csv.DictReader(stream))
    components = vtk_to_numpy(result.GetPointData().GetArray('ColorMatchComponentId'))
    for row in rows:
        component = int(row.pop('RegionId'))
        ids = np.flatnonzero(components == component)
        for values in before['point'].values():
            if len(np.unique(values[ids])) != 1:
                raise ValueError('A connected path has inconsistent original labels/colors')
        row['ComponentId'] = component
        row['OriginalRegionId'] = int(before['point']['RegionId'][ids[0]])
        row['PreviousColorId'] = int(before['point']['ColorId'][ids[0]])
        row['PreviousMatchedReference'] = int(before['point']['ReferenceColorMatched'][ids[0]])
        # Keep existing fallback colors stable when the geometry still fails
        # the original threshold, instead of reshuffling random colors.
        if not int(row['MatchedReference']):
            row['ColorId'] = row['PreviousColorId']
        row['AutomaticColorId'] = int(row['ColorId'])
        row['ManualColorOverride'] = int(row['OriginalRegionId'] in color_overrides)
        if row['ManualColorOverride']:
            row['ColorId'] = color_overrides[row['OriginalRegionId']]
        row['ColorChanged'] = int(int(row['ColorId']) != row['PreviousColorId'])

    missing = set(color_overrides) - {row['OriginalRegionId'] for row in rows}
    if missing:
        raise ValueError(f'Color override RegionId not present: {sorted(missing)}')

    fields = ['ComponentId', 'OriginalRegionId', 'ReferenceRegionId', 'ReferenceColorId',
              'ColorId', 'Distance', 'MatchedReference', 'Points', 'PreviousColorId',
              'PreviousMatchedReference', 'ColorChanged', 'AutomaticColorId', 'ManualColorOverride']
    with Path(matches_csv).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for association, data in [('point', result.GetPointData()), ('cell', result.GetCellData())]:
        ids = vtk_to_numpy(data.GetArray('ColorMatchComponentId'))
        color = np.zeros(len(rows), dtype=np.int32)
        automatic_color = np.zeros(len(rows), dtype=np.int32)
        manual = np.zeros(len(rows), dtype=np.int32)
        matched = np.zeros(len(rows), dtype=np.int32)
        ref_ids = np.full(len(rows), -1, dtype=np.int32)
        distances = np.zeros(len(rows), dtype=np.float64)
        for row in rows:
            component = row['ComponentId']
            color[component] = int(row['ColorId'])
            automatic_color[component] = row['AutomaticColorId']
            manual[component] = row['ManualColorOverride']
            matched[component] = int(row['MatchedReference'])
            distances[component] = float(row['Distance'])
            if matched[component]:
                ref_ids[component] = int(row['ReferenceRegionId'])
        add_array(data, 'ColorIdBeforeObstacleRematch', before[association]['ColorId'])
        add_array(data, 'ReferenceColorMatchedBeforeObstacleRematch', before[association]['ReferenceColorMatched'])
        add_array(data, 'ColorId', color[ids])
        add_array(data, 'AutomaticColorId', automatic_color[ids])
        add_array(data, 'ManualColorOverride', manual[ids])
        add_array(data, 'ReferenceColorMatched', matched[ids])
        add_array(data, 'ReferenceRegionId', ref_ids[ids])
        add_array(data, 'ReferenceMatchDistance', distances[ids])
    result.GetPointData().SetActiveScalars('ColorId')
    return result, rows


def write_points(poly, path):
    data = poly.GetPointData()
    names = [data.GetArrayName(i) for i in range(data.GetNumberOfArrays())]
    arrays = [vtk_to_numpy(data.GetArray(name)) for name in names]
    with Path(path).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['PositionX', 'PositionY', 'PositionZ'] + names)
        for index in range(poly.GetNumberOfPoints()):
            writer.writerow([*poly.GetPoint(index), *(array[index] for array in arrays)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-vtp', 'output-vtp', 'output-csv', 'tracking-json', 'matches-csv', 'report-json'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--reference-vtp', help='Override the color reference stored in tracking metadata')
    parser.add_argument('--color-override', nargs=2, type=int, action='append', default=[],
                        metavar=('REGION_ID', 'COLOR_ID'), help='Override an original RegionId color after automatic matching')
    args = parser.parse_args()
    settings = json.loads(Path(args.tracking_json).read_text())
    reference_path = str(Path(args.reference_vtp or settings['color_reference']).resolve())
    if not Path(reference_path).is_file():
        raise FileNotFoundError(reference_path)
    source = vtk_read_polydata(args.input_vtp)
    reference = vtk_read_polydata(reference_path)
    result, rows = recolor(source, reference, settings['color_distance_scales'],
                          settings['color_match_threshold'], args.matches_csv, dict(args.color_override))
    # This operation changes color metadata only, never coordinates or cells.
    np.testing.assert_array_equal(vtk_to_numpy(source.GetPoints().GetData()),
                                  vtk_to_numpy(result.GetPoints().GetData()))
    for before_cells, after_cells in [(source.GetVerts(), result.GetVerts()),
                                      (source.GetLines(), result.GetLines())]:
        np.testing.assert_array_equal(vtk_to_numpy(before_cells.GetConnectivityArray()),
                                      vtk_to_numpy(after_cells.GetConnectivityArray()))
        np.testing.assert_array_equal(vtk_to_numpy(before_cells.GetOffsetsArray()),
                                      vtk_to_numpy(after_cells.GetOffsetsArray()))
    vtk_write_polydata(result, args.output_vtp, 'binary')
    write_points(result, args.output_csv)
    report = dict(
        input=args.input_vtp, output=args.output_vtp, reference=reference_path,
        color_distance_scales=settings['color_distance_scales'],
        color_match_threshold=settings['color_match_threshold'],
        method='Original mean nearest-point matcher, 48 samples and 32 candidates, on surviving connected components',
        unmatched_policy='Preserve previous fallback color unless manually overridden; ReferenceRegionId=-1',
        manual_color_overrides=dict(args.color_override),
        manual_override_components=sum(row['ManualColorOverride'] for row in rows),
        manual_override_points=sum(int(row['Points']) for row in rows if row['ManualColorOverride']),
        components=len(rows), matched_components=sum(int(row['MatchedReference']) for row in rows),
        changed_components=sum(row['ColorChanged'] for row in rows),
        changed_points=sum(int(row['Points']) for row in rows if row['ColorChanged']),
        unmatched_components=sum(not int(row['MatchedReference']) for row in rows),
        fallback_components=sum(not int(row['MatchedReference']) and not row['ManualColorOverride'] for row in rows),
        geometry_unchanged=True, matches_csv=args.matches_csv)
    Path(args.report_json).write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
