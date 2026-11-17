"""Independent geometric checks for the analytic comparison utility."""

import math
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from evaluate_analytic_geometry import (SegmentGraph, analytic_reference,
    directed_hausdorff, ground_truth_events, match_events, read_events, read_graph)
from evaluate_rotating_components import component_indices, match_components


class GeometryTests(unittest.TestCase):
    def assert_bounded(self, actual, expected, tolerance):
        self.assertLessEqual(actual['lower'], expected + 1e-12)
        self.assertGreaterEqual(actual['upper'], expected - 1e-12)
        self.assertLessEqual(actual['upper']-actual['lower'], tolerance * 1.000001)

    def test_translated_lines(self):
        a = SegmentGraph([[0, 0, 0]], [[1, 0, 0]])
        b = SegmentGraph([[0, 0, 0.25]], [[1, 0, 0.25]])
        for source, target in ((a, b), (b, a)):
            self.assert_bounded(directed_hausdorff(source, target, 1e-8), .25, 1e-8)

    def test_line_interior_and_asymmetry(self):
        flat = SegmentGraph([[-1, 0, 0]], [[1, 0, 0]])
        tent = SegmentGraph([[-1, 0, 0], [0, 1, 0]], [[0, 1, 0], [1, 0, 0]])
        self.assert_bounded(directed_hausdorff(flat, tent, 1e-8), math.sqrt(.5), 1e-8)
        self.assert_bounded(directed_hausdorff(tent, flat, 1e-8), 1., 1e-8)

    def test_missing_endpoint_is_not_clipped(self):
        short = SegmentGraph([[0, 0, 0]], [[0, 0, 3.9]])
        full = SegmentGraph([[0, 0, 0]], [[0, 0, 4]])
        self.assert_bounded(directed_hausdorff(short, full, 1e-8), 0., 1e-8)
        self.assert_bounded(directed_hausdorff(full, short, 1e-8), .1, 1e-8)

    def test_symmetric_cutoff_retains_bounds(self):
        flat = SegmentGraph([[-1, 0, 0]], [[1, 0, 0]])
        tent = SegmentGraph([[-1, 0, 0], [0, 1, 0]], [[0, 1, 0], [1, 0, 0]])
        backward = directed_hausdorff(tent, flat, 1e-8)
        forward = directed_hausdorff(flat, tent, 1e-8, backward['lower'])
        lower = max(forward['lower'], backward['lower'])
        upper = max(forward['upper'], backward['upper'])
        self.assert_bounded(dict(lower=lower, upper=upper), 1., 1e-8)

    def test_obj_polyline_preserves_segments(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'graph.obj'
            path.write_text('# example\nv -1 0 0\nv 0 1 0\nv 1 0 0\nl 1 2 3\n')
            graph, info = read_graph(path)
            self.assertEqual(info['segments'], 2)
            self.assertEqual(info['unused_points'], 0)
            distance, _ = graph.nearest([0, 0, 0])
            self.assertAlmostEqual(distance, math.sqrt(.5))

    def test_analytic_reference_roots(self):
        for name in ('quartic_potential_2', 'rotating_quartic_multiwell'):
            graph, _ = analytic_reference(name, 1e-6)
            x, y, t = graph.starts.T
            if name == 'quartic_potential_2':
                fx, fy = x*(x*x+1-t), y*(y*y+np.cos(t))
            else:
                u, v = x*np.cos(t)+y*np.sin(t), -x*np.sin(t)+y*np.cos(t)
                fx, fy = u*(u*u-1), v*(v*v-1)
            self.assertLess(np.max(np.abs(fx)), 1e-13)
            self.assertLess(np.max(np.abs(fy)), 1e-13)

    def test_event_assignment_and_empty_reference(self):
        truth = ground_truth_events('quartic_potential_2')
        detected = truth[[3, 1, 0, 2]] + np.array([0, .003, .004])
        summary, details = match_events(detected, truth)
        self.assertAlmostEqual(summary['event_max_distance'], .005)
        self.assertAlmostEqual(summary['event_mean_distance'], .005)
        self.assertEqual([row['reference_index'] for row in details], [3, 1, 0, 2])
        summary, _ = match_events(detected[:3], truth)
        self.assertEqual(summary['unmatched_reference_events'], 1)
        summary, _ = match_events(detected, ground_truth_events('rotating_quartic_multiwell'))
        self.assertIsNone(summary['event_max_distance'])
        self.assertEqual(summary['unmatched_detected_events'], 4)

    def test_eigen_binary_layout(self):
        points = np.array([[1., 2., 3.], [4., 5., 6.]])
        raw = struct.pack('<iQQ', 1, 2, 3) + points.tobytes(order='F')
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'events.dat'
            path.write_bytes(raw)
            np.testing.assert_array_equal(read_events(path), points)

    def test_components_preserve_distinct_coincident_vertices(self):
        points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 0., 0.], [0., 1., 0.]])
        edges = np.array([[0, 1], [2, 3]])
        self.assertEqual(len(component_indices(points, edges)), 2)

    def test_rotating_correspondence_and_fragmentation(self):
        labels = [(u, v) for u in (1, -1, 0) for v in (0, 1, -1)]
        times = np.array([.1, .4, 1.3, 3.8])
        points, edges = [], []
        for u, v in labels:
            offset = len(points)
            points.extend(np.column_stack((u*np.cos(times)-v*np.sin(times),
                                            u*np.sin(times)+v*np.cos(times), times)))
            edges.extend((offset+i, offset+i+1) for i in range(3))
        points, edges = np.array(points), np.array(edges)
        matches = match_components(points, edges)
        self.assertEqual(len(matches), 9)
        for match in matches:
            self.assertEqual((match['u'], match['v']), labels[match['component_id']])
            self.assertLess(match['invariant_max_error'], 1e-14)
        with self.assertRaisesRegex(ValueError, 'found 10'):
            match_components(points, np.delete(edges, 1, axis=0))


if __name__ == '__main__':
    unittest.main()
