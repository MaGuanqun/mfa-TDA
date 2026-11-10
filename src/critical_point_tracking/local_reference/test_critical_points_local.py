"""Run with pvpython: analytic critical points and compatibility with the old extractor."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from critical_points_local import extract_reference, make_chunk, open_samples


class CriticalPointReferenceTests(unittest.TestCase):
    def test_analytic_points_and_original_extractor(self):
        with tempfile.TemporaryDirectory() as temp:
            prefix = Path(temp) / 'samples'
            lo, hi = np.array([-1., -.5, 13.5]), np.array([1., .5, 14.1])
            dims = np.array([41, 21, 5])
            x, y = np.linspace(lo[0], hi[0], dims[0]), np.linspace(lo[1], hi[1], dims[1])
            values = []
            for i in range(5):
                dx = x[None, :] - (-.4+.2*i)
                value = dx**2 + 2*y[:, None]**2
                if i in (1, 4):
                    value = -value
                elif i == 2:
                    value = dx**2 - 2*y[:, None]**2
                values.append(value)
            np.asarray(values, dtype='<f8').tofile(str(prefix) + '.bin')
            metadata = dict(format='local-field-v1', dtype='<f8', layout='t,y,x',
                            dimensions=dims.tolist(), sample_min=lo.tolist(), sample_max=hi.tolist(),
                            spacing=((hi-lo)/(dims-1)).tolist())
            Path(str(prefix) + '.json').write_text(json.dumps(metadata))
            result = extract_reference(prefix, threads=1, chunk_size=2)
            self.assertEqual(result['slices'], 5)
            self.assertEqual(result['reference_points'], 5)
            self.assertEqual(result['interior_points'], 3)
            with Path(str(prefix) + '_critical_points.csv').open() as stream:
                rows = list(csv.DictReader(stream))
            times = np.linspace(lo[2], hi[2], 5)
            for i, row in enumerate(rows):
                self.assertEqual(int(row['TimeIndex']), i)
                self.assertEqual(float(row['Time']), times[i])
                self.assertEqual(float(row['PositionZ']), times[i])
                self.assertEqual(int(row['IsOnBoundary']), 0)
                self.assertEqual(int(row['OnTimeBoundary']), int(i in (0, 4)))
                self.assertIn(int(row['CriticalType']), (1, 2) if i == 2 else (3,) if i in (1, 4) else (0,))
                np.testing.assert_allclose([float(row['PositionX']), float(row['PositionY'])],
                                           [-.4+.2*i, 0.], rtol=0, atol=1e-6)

            # Invoke the existing all-slices script without changing its code.
            info, data, times = open_samples(prefix)
            image = make_chunk(info, data, times, 0, 5)
            vti = str(prefix) + '.vti'
            writer = vtk.vtkXMLImageDataWriter()
            writer.SetFileName(vti)
            writer.SetInputData(image)
            self.assertEqual(writer.Write(), 1)
            old_csv = Path(temp) / 'original.csv'
            script = Path(__file__).resolve().parents[1] / 'extract_all_critical_points.py'
            pvpython = Path(sys.executable).resolve().parent.parent / 'bin' / 'pvpython'
            run = subprocess.run([str(pvpython), '--no-mpi', '--force-offscreen-rendering',
                                  str(script), '-i', vti, '-o', str(old_csv)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            with old_csv.open() as stream:
                old_rows = list(csv.DictReader(stream))
            expected = [[float(row[name]) for name in ('PositionX', 'PositionY', 'PositionZ')] for row in old_rows]
            actual = [[float(row[name]) for name in ('PositionX', 'PositionY', 'PositionZ')] for row in rows]
            np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-6)

            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(str(prefix) + '_critical_points_interior.vtp')
            reader.Update()
            interior = reader.GetOutput()
            self.assertEqual(interior.GetNumberOfVerts(), 3)
            np.testing.assert_array_equal(vtk_to_numpy(interior.GetPointData().GetArray('TimeIndex')), [1, 2, 3])
            with Path(str(prefix) + '_critical_points_slices.csv').open() as stream:
                slices = list(csv.DictReader(stream))
            self.assertEqual([int(row['TimeIndex']) for row in slices], list(range(5)))
            self.assertEqual([int(row['ReferencePoints']) for row in slices], [1]*5)


if __name__ == '__main__':
    unittest.main()
