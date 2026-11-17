"""Run with pvpython: clipping must retain detection attributes and singletons."""
import unittest

import vtk

from remove_obstacle import filter_polydata_in_circle


class ObstacleFilterTests(unittest.TestCase):
    def test_mixed_cells_attributes_and_new_singletons(self):
        source = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        points.SetDataTypeToDouble()
        for xyz in [(1, 0, 14), (0, 0, 14), (.5, 0, 14), (2, 0, 14),
                    (2, 1, 14), (.1, 0, 14), (0, 1, 14), (0, .1, 14),
                    (3, 0, 14), (0, -1, 14)]:
            points.InsertNextPoint(xyz)
        source.SetPoints(points)
        verts, lines = vtk.vtkCellArray(), vtk.vtkCellArray()
        for pid in [6, 7]:
            verts.InsertNextCell(1)
            verts.InsertCellPoint(pid)
        for a, b in [(0, 3), (3, 4), (4, 1), (0, 2), (1, 5), (9, 1)]:
            lines.InsertNextCell(2)
            lines.InsertCellPoint(a)
            lines.InsertCellPoint(b)
        source.SetVerts(verts)
        source.SetLines(lines)
        for data, count in [(source.GetPointData(), 10), (source.GetCellData(), 8)]:
            ids = vtk.vtkLongLongArray()
            ids.SetName('ColorId')
            names = vtk.vtkStringArray()
            names.SetName('Name')
            for i in range(count):
                ids.InsertNextValue(2**60 + i)
                names.InsertNextValue(f'value-{i}')
            data.SetScalars(ids)
            data.AddArray(names)
        field = vtk.vtkDoubleArray()
        field.SetName('Time')
        field.InsertNextValue(14.)
        source.GetFieldData().AddArray(field)

        kept, removed = filter_polydata_in_circle(source, 0, 0, .5, return_removed=True)
        expected = [0, 3, 4, 6, 8, 9]
        self.assertEqual(kept.GetNumberOfPoints(), 6)
        self.assertEqual(kept.GetNumberOfVerts(), 3)
        self.assertEqual(kept.GetNumberOfLines(), 2)
        self.assertEqual(removed.GetNumberOfPoints(), 4)
        self.assertEqual(removed.GetNumberOfLines(), 1)
        self.assertEqual(kept.GetPoints().GetDataType(), vtk.VTK_DOUBLE)
        for i, pid in enumerate(expected):
            self.assertEqual(kept.GetPoint(i), source.GetPoint(pid))
            self.assertEqual(kept.GetPointData().GetScalars().GetValue(i), 2**60 + pid)
            self.assertEqual(kept.GetPointData().GetAbstractArray('Name').GetValue(i), f'value-{pid}')
        # Existing vertex, original orphan, newly isolated endpoint, two lines.
        for new_id, old_id in [(0, 0), (2, 7), (3, 2), (4, 3)]:
            self.assertEqual(kept.GetCellData().GetScalars().GetValue(new_id), 2**60 + old_id)
            self.assertEqual(kept.GetCellData().GetAbstractArray('Name').GetValue(new_id), f'value-{old_id}')
        used = {kept.GetCell(i).GetPointId(j) for i in range(kept.GetNumberOfCells())
                for j in range(kept.GetCell(i).GetNumberOfPoints())}
        self.assertEqual(used, set(range(6)))
        self.assertEqual(kept.GetFieldData().GetArray('Time').GetValue(0), 14.)
        self.assertEqual(source.GetNumberOfPoints(), 10)

    def test_empty_input(self):
        kept, removed = filter_polydata_in_circle(vtk.vtkPolyData(), 0, 0, .12, return_removed=True)
        self.assertEqual(kept.GetNumberOfPoints(), 0)
        self.assertEqual(removed.GetNumberOfCells(), 0)


if __name__ == '__main__':
    unittest.main()
