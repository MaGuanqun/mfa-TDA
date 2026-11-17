"""Independent checks for the four-coordinate analytic Hausdorff evaluation."""

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from evaluate_analytic_geometry import directed_hausdorff
from evaluate_quartic_3d_geometry import SegmentGraphND, analytic_reference_4d, read_obj_4d


class FourDimensionalGeometryTests(unittest.TestCase):
    def test_locator_matches_exhaustive_projection(self):
        rng = np.random.default_rng(723)
        starts = rng.normal(size=(150, 4))
        ends = starts+rng.normal(size=(150, 4))*rng.lognormal(size=(150, 1))
        ends[0] = starts[0]  # zero-length segments remain valid point targets
        points = rng.normal(size=(2500, 4))
        graph = SegmentGraphND(starts, ends)
        distances, cells = graph.nearest_many(points)
        delta = ends-starts
        length2 = np.sum(delta*delta, axis=1)
        fraction = np.sum((points[:, None, :]-starts)*delta, axis=2)/np.where(length2 > 0, length2, 1.)
        projected = starts+np.clip(fraction, 0, 1)[:, :, None]*delta
        brute = np.linalg.norm(points[:, None, :]-projected, axis=2)
        np.testing.assert_allclose(distances, brute.min(axis=1), atol=1e-14, rtol=1e-14)
        np.testing.assert_allclose(distances, brute[np.arange(len(points)), cells], atol=1e-14)

    def test_far_midpoint_of_long_segment_is_considered(self):
        graph = SegmentGraphND([[0, 2, 0, 0], [-100, 0, 0, 0]],
                               [[0, 3, 0, 0], [1, 0, 0, 0]])
        distance, cell = graph.nearest([0, 0, 0, 0])
        self.assertAlmostEqual(distance, 0.)
        self.assertEqual(cell, 1)

    def test_time_coordinate_and_interior_maximum(self):
        line = SegmentGraphND([[-1, 0, 0, 0]], [[1, 0, 0, 0]])
        tent = SegmentGraphND([[-1, 0, 0, 0], [0, 0, 0, 1]],
                              [[0, 0, 0, 1], [1, 0, 0, 0]])
        for source, target, expected in ((line, tent, math.sqrt(.5)), (tent, line, 1.)):
            bounds = directed_hausdorff(source, target, 1e-8)
            self.assertLessEqual(bounds['lower'], expected+1e-12)
            self.assertGreaterEqual(bounds['upper'], expected-1e-12)
            self.assertLessEqual(bounds['upper']-bounds['lower'], 1.00001e-8)

    def test_obj_preserves_fourth_coordinate_and_connectivity(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'graph.obj'
            path.write_text('v 0 0 0 0\nv 0 0 0 2\nv 0 0 0 4\nl 1 2 3\n')
            graph, info = read_obj_4d(path)
            self.assertEqual(info['segments'], 2)
            self.assertEqual(info['bounds_max'], [0, 0, 0, 4])
            self.assertAlmostEqual(graph.nearest([0, 0, 0, 5])[0], 1.)

    def test_complete_reference_and_midpoint_interpolation_bounds(self):
        graph, pieces = analytic_reference_4d(1e-5)
        self.assertEqual(len(pieces), 27)
        xyz_t = np.r_[graph.starts, graph.ends]
        x, y, z, t = xyz_t.T
        gradient = np.column_stack((x*(x*x+1-t), y*(y*y+np.cos(t)), z*(z*z+np.cos(t))))
        self.assertLess(np.max(np.abs(gradient)), 1e-13)
        terminal = np.unique(xyz_t[np.abs(t-4) < 1e-14], axis=0)
        self.assertEqual(len(terminal), 27)
        self.assertEqual(float(t.min()), 0.)
        offset = 0
        for piece in pieces:
            n = piece['segments']
            a, b = graph.starts[offset:offset+n], graph.ends[offset:offset+n]
            middle = (a+b)/2
            if offset == 0:
                exact = middle
            elif piece['t_min'] == 1.:
                parameter = (np.abs(a[:, 0])+np.abs(b[:, 0]))/2
                exact = middle.copy()
                exact[:, -1] = 1+parameter**2
            else:
                # Reconstruct the actual uniform sampling parameter independently
                # of rounded times near the degenerate starting endpoint.
                parameter = (np.arange(n)+.5)*math.sqrt(4-math.pi/2)/n
                times = math.pi/2+parameter**2
                signs = np.sign(b[-1, :3])
                yz = np.sqrt(np.maximum(0., -np.cos(times)))
                exact = np.column_stack((signs[0]*np.sqrt(times-1), signs[1]*yz, signs[2]*yz, times))
            self.assertLessEqual(float(np.linalg.norm(middle-exact, axis=1).max()),
                                 piece['interpolation_error_bound']+1e-11)
            offset += n


if __name__ == '__main__':
    unittest.main()
