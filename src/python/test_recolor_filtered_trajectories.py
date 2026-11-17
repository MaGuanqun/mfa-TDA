"""Run with pvpython: split paths must be matched independently after clipping."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from recolor_filtered_trajectories import add_array, recolor


def polydata(xyz, lines, verts, regions, colors):
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.asarray(xyz, dtype=float), deep=True))
    poly.SetPoints(points)
    for cells, setter in [(verts, poly.SetVerts), (lines, poly.SetLines)]:
        array = vtk.vtkCellArray()
        for ids in cells:
            array.InsertNextCell(len(ids))
            for pid in ids:
                array.InsertCellPoint(pid)
        setter(array)
    for name, values in [('RegionId', regions), ('ColorId', colors),
                         ('ReferenceColorMatched', [0]*len(xyz))]:
        add_array(poly.GetPointData(), name, np.array(values, dtype=np.int32))
        add_array(poly.GetCellData(), name, np.array([values[ids[0]] for ids in verts+lines], dtype=np.int32))
    return poly


class RecolorTests(unittest.TestCase):
    def setUp(self):
        self.source = polydata([(0.1,0,0), (0.1,0,1), (10.1,0,0), (10.1,0,1), (100,0,0)],
                               [[0,1],[2,3]], [[4]], [7,7,7,7,8], [39,39,39,39,5])
        self.reference = polydata([(0,0,0), (0,0,1), (10,0,0), (10,0,1)],
                                  [[0,1],[2,3]], [], [0,0,1,1], [10,10,20,20])

    def test_split_components_and_stable_unmatched_color(self):
        with tempfile.TemporaryDirectory() as folder:
            result, rows = recolor(self.source, self.reference, [1,1,1], 1, Path(folder)/'matches.csv')
        pd = result.GetPointData()
        np.testing.assert_array_equal(vtk_to_numpy(pd.GetArray('ColorId')), [10,10,20,20,5])
        np.testing.assert_array_equal(vtk_to_numpy(pd.GetArray('RegionId')), [7,7,7,7,8])
        np.testing.assert_array_equal(vtk_to_numpy(pd.GetArray('ReferenceRegionId')), [0,0,1,1,-1])
        np.testing.assert_array_equal(vtk_to_numpy(pd.GetArray('ColorIdBeforeObstacleRematch')), [39,39,39,39,5])
        np.testing.assert_array_equal(vtk_to_numpy(result.GetCellData().GetArray('ColorId')), [5,10,20])
        np.testing.assert_array_equal(vtk_to_numpy(self.source.GetPointData().GetArray('ColorId')), [39,39,39,39,5])
        np.testing.assert_array_equal(vtk_to_numpy(result.GetLines().GetConnectivityArray()), [0,1,2,3])
        self.assertEqual(sum(row['ColorChanged'] for row in rows), 2)
        self.assertEqual(result.GetNumberOfVerts(), 1)

    def test_same_threshold_and_scales_are_respected(self):
        with tempfile.TemporaryDirectory() as folder:
            result, rows = recolor(self.source, self.reference, [.01,.01,.01], 1, Path(folder)/'matches.csv')
        np.testing.assert_array_equal(vtk_to_numpy(result.GetPointData().GetArray('ColorId')), [39,39,39,39,5])
        self.assertTrue(all(not int(row['MatchedReference']) for row in rows))


if __name__ == '__main__':
    unittest.main()
