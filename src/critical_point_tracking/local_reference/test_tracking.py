"""Run with pvpython: test file layout, physical time, and a known moving maximum."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from paraview import simple as pv, servermanager
from track_local_region import load_sample, annotate, save_tracking_outputs
from tracking_script import track_image, compute_tracking, load_ttk
from matching_ratio import load_csv_points, load_vtp_points


def original_pipeline(image, times):
    """Independent reproduction of the checked-in original tracking_script.py.

    Only its invalid output port 1 is corrected to mesh port 0 for this TTK.
    Do not call the new tracking helper here: this is a compatibility baseline.
    """
    load_ttk()
    proxies = []
    try:
        producer = pv.TrivialProducer()
        producer.GetClientSideObject().SetOutput(image)
        proxies.append(producer)
        precond = pv.TTKArrayPreconditioning(Input=producer)
        proxies.append(precond)
        precond.PointDataArrays = [image.GetPointData().GetArrayName(i)
                                  for i in range(image.GetPointData().GetNumberOfArrays())]
        precond.UpdatePipeline()
        tet = pv.Tetrahedralize(Input=precond)
        proxies.append(tet)
        tracker = pv.TTKTrackingFromFields(Input=tet)
        proxies.append(tracker)
        tracker.ForceZtranslation = 1
        tracker.ZTranslation = float(times[1]-times[0])
        tracker.Persistencethreshold = 0.
        mesh = vtk.vtkUnstructuredGrid()
        mesh.DeepCopy(servermanager.Fetch(tracker))
        return mesh
    finally:
        for proxy in reversed(proxies):
            pv.Delete(proxy)


def assert_same_mesh(test, actual, expected):
    test.assertEqual(actual.GetNumberOfPoints(),expected.GetNumberOfPoints())
    test.assertEqual(actual.GetNumberOfCells(),expected.GetNumberOfCells())
    if actual.GetNumberOfPoints():
        test.assertEqual(actual.GetPoints().GetDataType(),expected.GetPoints().GetDataType())
        np.testing.assert_array_equal(vtk_to_numpy(actual.GetPoints().GetData()),
                                      vtk_to_numpy(expected.GetPoints().GetData()))
    for c in range(actual.GetNumberOfCells()):
        test.assertEqual(actual.GetCellType(c),expected.GetCellType(c))
        a,b = actual.GetCell(c),expected.GetCell(c)
        test.assertEqual([a.GetPointId(i) for i in range(a.GetNumberOfPoints())],
                         [b.GetPointId(i) for i in range(b.GetNumberOfPoints())])
    for a,b in ((actual.GetPointData(),expected.GetPointData()),
                (actual.GetCellData(),expected.GetCellData()),
                (actual.GetFieldData(),expected.GetFieldData())):
        test.assertEqual(a.GetNumberOfArrays(),b.GetNumberOfArrays())
        for i in range(a.GetNumberOfArrays()):
            aa,bb = a.GetAbstractArray(i),b.GetAbstractArray(i)
            test.assertEqual((aa.GetName(),aa.GetDataType(),aa.GetNumberOfComponents()),
                             (bb.GetName(),bb.GetDataType(),bb.GetNumberOfComponents()))
            np.testing.assert_array_equal(vtk_to_numpy(aa),vtk_to_numpy(bb))
        if hasattr(a,"GetScalars"):
            test.assertEqual(a.GetScalars().GetName() if a.GetScalars() else None,
                             b.GetScalars().GetName() if b.GetScalars() else None)


class TrackingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.prefix = str(Path(self.temp.name)/"sample")
        self.lo = [-1.,-.5,13.5]
        self.hi = [1.,.5,14.1]
        self.dims = [41,21,5]
        self.times = np.linspace(self.lo[2],self.hi[2],self.dims[2])
        x = np.linspace(self.lo[0],self.hi[0],self.dims[0])
        y = np.linspace(self.lo[1],self.hi[1],self.dims[1])
        self.values = np.array([-((x[None,:]-(-.4+.2*i))**2+2*y[:,None]**2)
                                for i in range(self.dims[2])],dtype="<f8")
        self.values.tofile(self.prefix+".bin")
        self.info = dict(format="local-field-v1",dtype="<f8",layout="t,y,x",
                         dimensions=self.dims,sample_min=self.lo,sample_max=self.hi,
                         region_min=[-.6,-.2,13.6],region_max=[.6,.2,14.],
                         spacing=((np.array(self.hi)-self.lo)/(np.array(self.dims)-1)).tolist())
        with open(self.prefix+".json","w") as stream:
            json.dump(self.info,stream)

    def tearDown(self):
        self.temp.cleanup()

    def test_layout_and_corrupt_size(self):
        _, mapped, image, times = load_sample(self.prefix)
        np.testing.assert_allclose(times,self.times)
        self.assertEqual(image.GetDimensions(),(41,21,1))
        self.assertEqual(image.GetOrigin(),tuple(self.lo))
        for i in range(5):
            np.testing.assert_array_equal(vtk_to_numpy(image.GetPointData().GetArray(i)),self.values[i].ravel())
        del image,mapped
        with open(self.prefix+".bin","ab") as stream:
            stream.write(b"x")
        with self.assertRaisesRegex(ValueError,"size"):
            load_sample(self.prefix)

    def test_moving_maximum_and_vti_roundtrip(self):
        info, mapped, image, times = load_sample(self.prefix)
        mesh, _ = track_image(image,times,threads=1)
        self.assertGreater(mesh.GetNumberOfCells(),0)
        self.assertGreater(save_tracking_outputs(mesh,self.prefix,info,times),0)
        self.assertIsNone(mesh.GetPointData().GetArray("TimeIndex"))
        self.assertIsNone(mesh.GetPointData().GetArray("InsideRequestedRegion"))
        coords = vtk_to_numpy(mesh.GetPoints().GetData())
        time_index = np.rint((coords[:,2]-times[0])/(times[1]-times[0])).astype(int)
        np.testing.assert_allclose(coords[:,2],times[time_index],rtol=0,atol=2e-6)
        # At every time the exact maximum is a grid vertex; require it in the
        # tracked mesh. A transpose or lost physical time offset fails here.
        for i, time in enumerate(times):
            frame = coords[time_index == i]
            self.assertGreater(len(frame),0)
            error = np.linalg.norm(frame[:,:2]-[-.4+.2*i,0],axis=1)
            self.assertLess(error.min(),1e-6)
        maxima = np.linalg.norm(coords[:,:2]-np.column_stack((-.4+.2*time_index,
                                                            np.zeros(len(time_index)))),axis=1) < 1e-6
        links = set()
        for c in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(c)
            if cell.GetNumberOfPoints() == 2:
                a,b = cell.GetPointId(0),cell.GetPointId(1)
                if maxima[a] and maxima[b]:
                    links.add(tuple(sorted((int(time_index[a]),int(time_index[b])))))
        self.assertEqual(links,{(0,1),(1,2),(2,3),(3,4)})
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.prefix+"_tracking.vtp")
        reader.Update()
        assert_same_mesh(self,reader.GetOutput(),mesh)
        # Exercise the existing analysis readers without modifying them.
        np.testing.assert_array_equal(load_vtp_points(self.prefix+"_tracking.vtp"),coords)
        np.testing.assert_array_equal(load_csv_points(self.prefix+"_tracking_points.csv"),coords)
        reference = original_pipeline(image,times)
        assert_same_mesh(self,mesh,reference)
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(self.prefix+".vti")
        writer.SetInputData(image)
        self.assertEqual(writer.Write(),1)
        compute_tracking(self.prefix+".vti",self.prefix+"_roundtrip.vtu",threads=1)
        reader2 = vtk.vtkXMLUnstructuredGridReader()
        reader2.SetFileName(self.prefix+"_roundtrip.vtu")
        reader2.Update()
        other = reader2.GetOutput()
        assert_same_mesh(self,other,reference)

    def test_original_pipeline_with_zero_plane_origin(self):
        _, mapped, image, times = load_sample(self.prefix)
        image.SetOrigin(self.lo[0],self.lo[1],0.)
        mesh, _ = track_image(image,times,threads=1)
        assert_same_mesh(self,mesh,original_pipeline(image,times))

    def test_bad_time_order_and_empty_annotation(self):
        info, mapped, image, times = load_sample(self.prefix)
        with self.assertRaisesRegex(ValueError,"increasing"):
            track_image(image,times[::-1])
        mesh = vtk.vtkUnstructuredGrid()
        self.assertEqual(annotate(mesh,info,times),0)


if __name__ == "__main__":
    unittest.main()
