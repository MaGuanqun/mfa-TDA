"""Checks that spatial projection respects fixed time, dimension, and component."""

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from evaluate_analytic_same_time import nearest_same_time, read_vertices, same_time_roots


class SameTimeTests(unittest.TestCase):
    def test_time_cannot_move_to_make_a_vertex_exact(self):
        # x=.2 is an exact root at t=1.04, but only x=0 is real at t=1.
        error, reference = nearest_same_time([.2, 0., 1.], 'quartic_potential_2')
        self.assertAlmostEqual(error, .2)
        np.testing.assert_array_equal(reference, [0., 0., 1.])

    def test_only_real_roots_exist_before_a_birth(self):
        error, reference = nearest_same_time([0., .2, 1.5], 'quartic_potential_2')
        self.assertAlmostEqual(error, .2)
        np.testing.assert_array_equal(reference, [0., 0., 1.5])
        error, reference = nearest_same_time([0., .8, math.pi], 'quartic_potential_2')
        self.assertAlmostEqual(error, .2)
        np.testing.assert_array_equal(reference, [0., 1., math.pi])

    def test_third_spatial_coordinate_is_included(self):
        error, reference = nearest_same_time([0., 0., .4, .5], 'quartic_potential_3d')
        self.assertAlmostEqual(error, .4)
        np.testing.assert_array_equal(reference, [0., 0., 0., .5])
        error, reference = nearest_same_time([0., .9, -.8, math.pi], 'quartic_potential_3d')
        self.assertAlmostEqual(error, math.sqrt(.05))
        np.testing.assert_array_equal(reference, [0., 1., -1., math.pi])

    def test_rotating_reference_cannot_switch_components(self):
        error, reference = nearest_same_time([1., 0., 0.], 'rotating_quartic_multiwell', (-1, 0))
        self.assertAlmostEqual(error, 2.)
        np.testing.assert_array_equal(reference, [-1., 0., 0.])
        with self.assertRaisesRegex(ValueError, 'matched rotating component'):
            nearest_same_time([1., 0., 0.], 'rotating_quartic_multiwell')

    def test_outside_time_domain_is_rejected_without_clamping(self):
        for t in (-1e-8, 4+1e-8):
            with self.assertRaisesRegex(ValueError, 'outside the analytic domain'):
                same_time_roots('quartic_potential_3d', t)

    def test_obj_retains_time_and_original_vertex_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'graph.obj'
            path.write_text('v 0 0 .4 .5\nv 0 0 .6 .7\nv 9 9 9 9\nl 1 2\n')
            points, edges = read_vertices(path, 3)
            np.testing.assert_array_equal(points[0], [0., 0., .4, .5])
            np.testing.assert_array_equal(np.unique(edges), [0, 1])
            with self.assertRaisesRegex(ValueError, 'wrong number of vertex coordinates'):
                read_vertices(path, 2)


if __name__ == '__main__':
    unittest.main()
