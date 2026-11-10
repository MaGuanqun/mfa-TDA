"""Validate reference colors and export one VTP per SFFF RegionId (pvpython)."""
import argparse
import csv
import json
from pathlib import Path
import re
import struct

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_poly(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    poly = reader.GetOutput()
    if not poly.GetNumberOfPoints() or not poly.GetNumberOfLines():
        raise ValueError(f'Empty trajectory: {path}')
    return poly


def export(folder, reference):
    poly = read_poly(folder / 'tracking_colorid.vtp')
    ref = read_poly(reference)
    with (folder / 'color_matches.csv').open() as stream:
        matches = list(csv.DictReader(stream))
    ref_regions = vtk_to_numpy(ref.GetPointData().GetArray('RegionId'))
    ref_colors = vtk_to_numpy(ref.GetPointData().GetArray('ColorId'))
    for data in (poly.GetPointData(), poly.GetCellData()):
        regions = vtk_to_numpy(data.GetArray('RegionId'))
        colors = vtk_to_numpy(data.GetArray('ColorId'))
        if set(regions) != {int(row['RegionId']) for row in matches}:
            raise ValueError('Color report does not match trajectory branches')
        for row in matches:
            region, color = int(row['RegionId']), int(row['ColorId'])
            expected = ref_colors[ref_regions == int(row['ReferenceRegionId'])]
            if (not int(row['MatchedReference']) or not len(expected)
                    or not np.all(expected == color) or not np.all(colors[regions == region] == color)):
                raise ValueError(f'Branch {region} does not carry its reference color')

    branch_dir = folder / 'branches'
    branch_dir.mkdir(exist_ok=True)
    rows = []
    for match in matches:
        region, color = int(match['RegionId']), int(match['ColorId'])
        threshold = vtk.vtkThreshold()
        threshold.SetInputData(poly)
        threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, 'RegionId')
        threshold.SetLowerThreshold(region)
        threshold.SetUpperThreshold(region)
        threshold.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
        geometry = vtk.vtkGeometryFilter()
        geometry.SetInputConnection(threshold.GetOutputPort())
        geometry.Update()
        branch = geometry.GetOutput()
        if not branch.GetNumberOfLines():
            raise ValueError(f'Empty branch {region}')
        filename = f'branch_{region:04d}_color_{color}.vtp'
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(str(branch_dir / filename))
        writer.SetInputData(branch)
        if writer.Write() != 1:
            raise OSError(f'Could not write {filename}')
        rows.append(dict(match, Points=branch.GetNumberOfPoints(),
                         Segments=branch.GetNumberOfLines(), File=f'branches/{filename}'))
    if sum(row['Segments'] for row in rows) != poly.GetNumberOfLines():
        raise ValueError('Branch export changed the total segment count')
    # Remove only this exporter's stale branch outputs after a successful rerun.
    current = {Path(row['File']).name for row in rows}
    for path in branch_dir.glob('branch_*_color_*.vtp'):
        if path.name not in current:
            path.unlink()
    with (folder / 'branches.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with (folder / 'gradient.vff').open('rb') as stream:
        if stream.read(4) != b'VFF1' or struct.unpack('<III', stream.read(12)) != (0, 3, 3):
            raise ValueError('Expected a float64 3D full-gradient VFF1 file')
        dimensions = struct.unpack('<3I', stream.read(12))
        lower = struct.unpack('<3d', stream.read(24))
        upper = struct.unpack('<3d', stream.read(24))
    points = vtk_to_numpy(poly.GetPoints().GetData()).astype(float)
    if not np.isfinite(points).all() or np.any(points < np.array(lower)-1e-6) or np.any(points > np.array(upper)+1e-6):
        raise ValueError('Trajectory has invalid or out-of-domain coordinates')
    x, y, t = points.T
    residual = np.hypot(x**3 + (1-t)*x, y**3 + np.cos(t)*y)
    with (folder / 'seeds.csv').open() as stream:
        seed_count = sum(1 for _ in csv.DictReader(stream))
    tracking_log = (folder / 'tracking.log').read_text()
    trace_stats = re.search(
        r'stream lines built: (\d+), closed loops: (\d+), seeds skipped \(already covered\): (\d+)',
        tracking_log,
    )
    if trace_stats is None:
        raise ValueError('Missing tracking statistics in tracking.log')
    report = dict(
        function='quartic_potential_2', dimensions=dimensions,
        domain_min=lower, domain_max=upper,
        spacing=((np.array(upper)-lower)/(np.array(dimensions)-1)).tolist(),
        gradient_components=['df/dx', 'df/dy', 'df/dt'], seeds=seed_count,
        points=poly.GetNumberOfPoints(), segments=poly.GetNumberOfLines(),
        branches=len(rows), matched_reference_regions=sorted({int(row['ReferenceRegionId']) for row in rows}),
        reference_regions=sorted(int(r) for r in set(ref_regions)),
        streamlines=int(trace_stats[1]), closed_loops=int(trace_stats[2]), skipped_seeds=int(trace_stats[3]),
        color_matching='Closest reference branch by mean sampled point distance; no distance cutoff',
        maximum_branch_mean_match_distance=max(float(row['Distance']) for row in rows),
        analytical_spatial_gradient_norm=dict(mean=float(residual.mean()), max=float(residual.max())),
        reference_vtp=str(reference.resolve()),
        colored_vtp=str((folder / 'tracking_colorid.vtp').resolve()),
    )
    (folder / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Validated colors and exported {len(rows)} branches to {branch_dir}')
    print(f'Analytical spatial-gradient residual: mean={residual.mean():.6g}, max={residual.max():.6g}')
    print(f'Summary: {folder / "summary.json"}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment-dir', type=Path, required=True)
    parser.add_argument('--reference-vtp', type=Path, required=True)
    args = parser.parse_args()
    export(args.experiment_dir, args.reference_vtp)
