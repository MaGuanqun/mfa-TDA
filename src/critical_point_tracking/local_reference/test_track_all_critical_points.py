"""Run with pvpython: coverage, saddles, births/deaths, gaps and assignments."""
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from track_all_critical_points import assign_frame, track_points, values


def detections(rows):
    # Each row: time index, x, y, critical type; scalar values may be zero.
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    coords = np.asarray([[r[1], r[2], 13.5+r[0]*.01] for r in rows], dtype=float).reshape(-1, 3)
    points.SetData(numpy_to_vtk(coords, deep=True))
    poly.SetPoints(points)
    for name, array in [('TimeIndex', np.asarray([r[0] for r in rows], dtype=np.int64)),
                        ('CriticalType', np.asarray([r[3] for r in rows], dtype=np.int32)),
                        ('VertexId', np.arange(len(rows), dtype=np.int64)),
                        ('Scalar', np.zeros(len(rows)))]:
        data = numpy_to_vtk(array, deep=True)
        data.SetName(name)
        poly.GetPointData().AddArray(data)
    return poly


class TrackAllCriticalPointsTests(unittest.TestCase):
    def test_types_saddles_and_coverage(self):
        source = detections([(t, x+.001*t, 0., kind)
                             for t in range(4) for x, kind in [(0., 0), (.1, 1), (.2, 3)]])
        result, stats = track_points(source, .00625)
        self.assertEqual(stats['points'], 12)
        self.assertEqual(stats['segments'], 9)
        self.assertEqual(stats['tracks'], 3)
        self.assertEqual(stats['singleton_points'], 0)
        self.assertEqual(stats['points_by_type'], {'0': 4, '1': 4, '3': 4})
        np.testing.assert_array_equal(vtk_to_numpy(result.GetPoints().GetData()),
                                      vtk_to_numpy(source.GetPoints().GetData()))
        for kind in [0, 1, 3]:
            labels = values(result, 'RegionId')[values(result, 'CriticalType') == kind]
            self.assertEqual(len(set(labels)), 1)

    def test_unmatched_and_short_lived_points_are_retained(self):
        source = detections([(0, 0., 0., 0), (1, .1, 0., 0),  # too far to connect
                             (1, .2, 0., 1), (2, .201, 0., 1),  # short saddle track
                             (2, .3, 0., 3)])  # single-frame maximum
        result, stats = track_points(source, .00625)
        self.assertEqual(stats['points'], 5)
        self.assertEqual(stats['segments'], 1)
        self.assertEqual(stats['singleton_points'], 3)
        self.assertEqual(result.GetNumberOfVerts(), 3)
        np.testing.assert_array_equal(values(result, 'TrackLength'), [1, 1, 2, 2, 1])

    def test_no_type_changes_or_gap_bridging(self):
        source = detections([(0, 0., 0., 0), (1, 0., 0., 3), (2, 0., 0., 0)])
        _, stats = track_points(source, .00625)
        self.assertEqual(stats['singleton_points'], 3)
        self.assertEqual(stats['segments'], 0)

    def test_global_assignment_instead_of_greedy_nearest(self):
        left = np.array([[0., 0.], [.08, 0.]])
        right = np.array([[.06, 0.], [.14, 0.]])
        self.assertEqual(assign_frame(left, right, .1), [(0, 0), (1, 1)])

    def test_empty_and_invalid_input(self):
        result, stats = track_points(detections([]), .00625)
        self.assertEqual(result.GetNumberOfPoints(), 0)
        self.assertEqual(stats['tracks'], 0)
        with self.assertRaises(ValueError):
            track_points(detections([]), 0)


if __name__ == '__main__':
    unittest.main()
