"""Extract per-slice TTK critical-point references from saved local-field-v1 samples.

Run with pvpython. Reuses mfa_bin_critical_point's preconditioning,
tetrahedralization and ScalarFieldCriticalPoints pipeline in bounded chunks.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mfa_bin_critical_point import critical_points_from_vtkImageData_chunk
from tracking_script import load_ttk


FIELDS = ['PositionX', 'PositionY', 'PositionZ', 'TimeIndex', 'Time',
          'CriticalType', 'IsOnBoundary', 'OnTimeBoundary', 'VertexId', 'Scalar']


def open_samples(prefix):
    info = json.loads(Path(str(prefix) + '.json').read_text())
    if (info.get('format') != 'local-field-v1' or info.get('dtype') != '<f8'
            or info.get('layout') != 't,y,x'):
        raise ValueError('Unsupported sample metadata format')
    dims = info['dimensions']
    if len(dims) != 3 or any(type(n) is not int or n < 2 for n in dims):
        raise ValueError('Each sample dimension must have at least two points')
    lo, hi = np.asarray(info['sample_min']), np.asarray(info['sample_max'])
    if lo.shape != (3,) or hi.shape != (3,) or not np.isfinite([lo, hi]).all() or np.any(hi <= lo):
        raise ValueError('Invalid sample bounds')
    spacing = (hi-lo)/(np.asarray(dims)-1)
    if not np.allclose(info['spacing'], spacing, rtol=1e-12, atol=1e-14):
        raise ValueError('Sample spacing does not match bounds/dimensions')
    nx, ny, nt = dims
    binary = Path(str(prefix) + '.bin')
    if binary.stat().st_size != nx*ny*nt*8:
        raise ValueError('Binary sample size does not match metadata')
    data = np.memmap(binary, dtype='<f8', mode='r', shape=(nt, ny, nx))
    return info, data, np.linspace(lo[2], hi[2], nt)


def make_chunk(info, data, times, first, last):
    nx, ny, _ = info['dimensions']
    image = vtk.vtkImageData()
    image.SetDimensions(nx, ny, 1)
    image.SetOrigin(*info['sample_min'])
    image.SetSpacing(*info['spacing'][:2], 1.)
    for index in range(first, last):
        values = data[index].reshape(-1)
        if not np.isfinite(values).all():
            raise ValueError(f'Nonfinite sample values at slice {index}')
        name = f'step_{index:09d}'
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        image.GetPointData().AddArray(array)
        time_array = vtk.vtkDoubleArray()
        time_array.SetName(name)
        time_array.InsertNextValue(float(times[index]))
        image.GetFieldData().AddArray(time_array)
    return image


def write_reference(prefix, rows):
    """CSV coordinates first for matching_ratio.py; VTP vertex cells for ParaView."""
    with Path(str(prefix) + '.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(FIELDS)
        writer.writerows(rows)
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.asarray([row[:3] for row in rows], dtype=np.float64).reshape(-1, 3), deep=True))
    poly.SetPoints(points)
    vertices = vtk.vtkCellArray()
    for index in range(len(rows)):
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(index)
    poly.SetVerts(vertices)
    for column, name in enumerate(FIELDS[3:], 3):
        dtype = np.float64 if name in ('Time', 'Scalar') else np.int64
        array = numpy_to_vtk(np.asarray([row[column] for row in rows], dtype=dtype), deep=True)
        array.SetName(name)
        poly.GetPointData().AddArray(array)
    poly.GetPointData().SetActiveScalars('CriticalType')
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(prefix) + '.vtp')
    writer.SetInputData(poly)
    if writer.Write() != 1:
        raise OSError(f'Could not write {prefix}.vtp')


def extract_reference(prefix, threads=3, chunk_size=8, plugin=None):
    if threads < 1 or chunk_size < 1:
        raise ValueError('Threads and chunk size must be positive')
    prefix = Path(prefix).resolve()
    info, data, times = open_samples(prefix)
    nx, ny, nt = info['dimensions']
    load_ttk(plugin)
    all_rows, slice_rows, seen = [], [], set()
    started = time.monotonic()

    def collect(name, mesh):
        index = int(name.removeprefix('step_'))
        if index in seen or not 0 <= index < nt:
            raise RuntimeError(f'Unexpected or repeated time slice: {name}')
        seen.add(index)
        pd = mesh.GetPointData()
        count = mesh.GetNumberOfPoints()
        names = ('CriticalType', 'IsOnBoundary', 'ttkVertexScalarField', name)
        arrays = []
        for key in names:
            array = pd.GetArray(key)
            if array is None or array.GetNumberOfTuples() != count:
                raise RuntimeError(f'Missing/malformed TTK array {key} at slice {index}')
            arrays.append(vtk_to_numpy(array).reshape(-1))
        types, boundary, vertex_ids, scalars = arrays
        vertex_ids = vertex_ids.astype(np.int64)
        if np.any(vertex_ids < 0) or np.any(vertex_ids >= nx*ny):
            raise RuntimeError(f'Invalid source vertex ids at slice {index}')
        if not np.array_equal(scalars, data[index].reshape(-1)[vertex_ids]):
            raise RuntimeError(f'TTK scalar values differ from sampled slice {index}')
        coords = (vtk_to_numpy(mesh.GetPoints().GetData()).astype(np.float64, copy=True)
                  if count else np.empty((0, 3)))
        expected_xy = np.column_stack((info['sample_min'][0] + (vertex_ids % nx)*info['spacing'][0],
                                       info['sample_min'][1] + (vertex_ids // nx)*info['spacing'][1]))
        tolerance = 8*np.finfo(np.float32).eps*np.maximum(1., abs(expected_xy))
        if np.any(abs(coords[:, :2]-expected_xy) > tolerance):
            raise RuntimeError(f'TTK positions do not match the sampled XY grid at slice {index}')
        coords[:, 2] = times[index]
        time_boundary = int(index in (0, nt-1))
        reference_types = types[boundary != 1]
        for point, kind, edge, vertex, scalar in zip(coords, types, boundary, vertex_ids, scalars):
            all_rows.append([*map(float, point), index, float(times[index]), int(kind),
                             int(edge), time_boundary, int(vertex), float(scalar)])
        reference_count = int((boundary != 1).sum())
        slice_rows.append([index, float(times[index]), count, int((boundary == 1).sum()),
                           reference_count, 0 if time_boundary else reference_count,
                           int((reference_types == 0).sum()),
                           int(np.isin(reference_types, [1, 2]).sum()),
                           int((reference_types == 3).sum()),
                           int((~np.isin(reference_types, [0, 1, 2, 3])).sum())])

    for first in range(0, nt, chunk_size):
        last = min(first+chunk_size, nt)
        image = make_chunk(info, data, times, first, last)
        critical_points_from_vtkImageData_chunk(image, None, threads=threads,
                                               on_slice=collect, quiet=True)
        del image
        print(f'Critical points: slices {last}/{nt}; {len(all_rows)} raw points; '
              f'{time.monotonic()-started:.1f}s elapsed', flush=True)
    if seen != set(range(nt)):
        raise RuntimeError('Critical point extraction did not cover every time slice')
    reference_rows = [row for row in all_rows if row[6] != 1]
    interior_rows = [row for row in reference_rows if not row[7]]
    output = Path(str(prefix) + '_critical_points')
    write_reference(output, reference_rows)
    write_reference(str(output) + '_all', all_rows)
    write_reference(str(output) + '_interior', interior_rows)
    with Path(str(output) + '_slices.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['TimeIndex', 'Time', 'RawPoints', 'SpatialBoundaryPoints',
                         'ReferencePoints', 'InteriorPoints', 'Minima', 'Saddles', 'Maxima', 'OtherTypes'])
        writer.writerows(slice_rows)
    summary = {
        'sample': str(prefix), 'slices': nt, 'dimensions': info['dimensions'],
        'time_range': [float(times[0]), float(times[-1])],
        'raw_points': len(all_rows), 'reference_points': len(reference_rows),
        'interior_points': len(interior_rows),
        'slices_with_reference_points': sum(row[4] > 0 for row in slice_rows),
        'threads': threads, 'chunk_size': chunk_size,
        'pipeline': 'TTKArrayPreconditioning -> Tetrahedralize -> TTKScalarFieldCriticalPoints',
        'reference_boundary_policy': 'Remove XY boundary points; retain all time slices',
        'interior_boundary_policy': 'Remove XY boundary points and first/last time slices',
        'reference_csv': str(output) + '.csv', 'reference_vtp': str(output) + '.vtp',
        'elapsed_seconds': time.monotonic()-started,
    }
    Path(str(output) + '.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True, help='Prefix of existing samples.bin/.json')
    parser.add_argument('--threads', type=int, default=3)
    parser.add_argument('--chunk-size', type=int, default=8)
    parser.add_argument('--ttk-plugin')
    args = parser.parse_args()
    extract_reference(args.sample, threads=args.threads, chunk_size=args.chunk_size, plugin=args.ttk_plugin)


if __name__ == '__main__':
    main()
