"""Track an XY image containing one scalar array per physical time, using pvpython.

Also imported by local_reference/track_local_region.py. The default pipeline
matches the original script: array preconditioning, explicit triangulation,
TTK's matching defaults, and Z translation by the physical time step. TTK's
native output geometry, connectivity and arrays are preserved.
"""
import argparse
import os
from pathlib import Path
import sys

import numpy as np
import vtk
from paraview import servermanager
from paraview import simple as pv


def load_ttk(plugin=None):
    if hasattr(pv, "TTKTrackingFromFields"):
        return
    plugin = plugin or os.environ.get("TTK_PLUGIN")
    if not plugin:
        root = Path(sys.executable).resolve().parent.parent
        candidates = sorted(root.glob("lib/paraview*/plugins/TopologyToolKit/TopologyToolKit.so"))
        if candidates:
            plugin = str(candidates[0])
    if plugin:
        pv.LoadPlugin(plugin, remote=False, ns=vars(pv))
    else:
        pv.LoadDistributedPlugin("TopologyToolKit", remote=False, ns=vars(pv))
    if not hasattr(pv, "TTKTrackingFromFields"):
        raise RuntimeError("TTKTrackingFromFields unavailable; supply --ttk-plugin /path/to/TopologyToolKit.so")


def track_image(image, times, plugin=None, threads=3, persistence=0.,
                destruction_cost=None, tetrahedralize=True, spatial_weight=None):
    if image.GetDimensions()[2] != 1 or len(times) < 3:
        raise ValueError("Tracking requires an XY image and at least three time-step arrays for legacy TTK compatibility")
    times = np.asarray(times, dtype=float)
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Time values must be finite and strictly increasing")
    delta_t = float(times[1] - times[0])
    if not np.allclose(np.diff(times), delta_t, rtol=1e-9, atol=1e-12):
        raise ValueError("The original TTK Z-translation convention requires uniformly spaced time values")
    if threads < 1 or not np.isfinite(persistence) or not 0 <= persistence <= 100:
        raise ValueError("Threads must be positive; persistence must be in [0,100]")
    load_ttk(plugin)
    source_image = vtk.vtkImageData()
    source_image.ShallowCopy(image)
    producer = pv.TrivialProducer()
    producer.GetClientSideObject().SetOutput(source_image)
    proxies, errors = [producer], []
    try:
        # Preserve the original pipeline, including its ordering/tie breaking.
        source = pv.TTKArrayPreconditioning(Input=producer)
        proxies.append(source)
        source.PointDataArrays = [image.GetPointData().GetArrayName(i)
                                  for i in range(image.GetPointData().GetNumberOfArrays())
                                  if not image.GetPointData().GetArrayName(i).endswith("_Order")]
        if "UseAllCores" in source.ListProperties():
            source.UseAllCores = 0
            source.ThreadNumber = threads
        source.UpdatePipeline()
        if tetrahedralize:
            source = pv.Tetrahedralize(Input=source)
            proxies.append(source)
        tracker = pv.TTKTrackingFromFields(Input=source)
        proxies.append(tracker)
        properties = tracker.ListProperties()
        tracker.ForceZtranslation = 1
        tracker.ZTranslation = delta_t
        tracker.Persistencethreshold = persistence
        # Explicitly use every available frame. A zero persistence threshold
        # still does not make the legacy persistence-pair backend emit saddles.
        for name, value in (("Firsttimestep", 0), ("Lasttimestep", -1), ("Timesampling", 1)):
            if name in properties:
                setattr(tracker, name, value)
        if spatial_weight is not None:
            if not np.isfinite(spatial_weight) or not 0 <= spatial_weight <= 1:
                raise ValueError("Spatial matching weight must be in [0,1]")
            tracker.Xweight = spatial_weight
            tracker.Yweight = spatial_weight
        if "UseAllCores" in properties:
            tracker.UseAllCores = 0
            tracker.ThreadNumber = threads
        if destruction_cost is not None:
            if "Relativedestructioncost" not in properties:
                raise ValueError("This TTK build has no Relativedestructioncost property; omit that option or use a newer TTK")
            if not np.isfinite(destruction_cost) or not 0 <= destruction_cost <= 1:
                raise ValueError("Relative destruction cost must be in [0,1]")
            tracker.Relativedestructioncost = destruction_cost
        tracker.GetClientSideObject().AddObserver("ErrorEvent", lambda *args: errors.append("TTK ErrorEvent"))
        tracker.UpdatePipeline()
        if errors:
            raise RuntimeError("TTK reported an error; inspect its diagnostic output")
        # The installed TTK exposes its tracking mesh on port 0, not port 1.
        fetched = servermanager.Fetch(tracker)
        if not isinstance(fetched, vtk.vtkUnstructuredGrid):
            raise RuntimeError("Expected a TTK tracking mesh on output port 0")
        output = vtk.vtkUnstructuredGrid()
        output.DeepCopy(fetched)
        settings = {name: str(getattr(tracker, name)) for name in properties
                    if name in ("Backend", "Assignmentmethod", "AssignmentMethod",
                                "Persistencethreshold", "Relativedestructioncost",
                                "Firsttimestep", "Lasttimestep", "Timesampling",
                                "Xweight", "Yweight", "Zweight", "Fweight",
                                "Extremumweight", "Saddleweight", "pparameter")}
        return output, settings
    finally:
        for proxy in reversed(proxies):
            pv.Delete(proxy)


def write_mesh(mesh, filename):
    extension = Path(filename).suffix.lower()
    if extension == ".vtu":
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputData(mesh)
    elif extension == ".vtp":
        surface = vtk.vtkDataSetSurfaceFilter()
        surface.SetInputData(mesh)
        surface.Update()
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetInputData(surface.GetOutput())
    else:
        raise ValueError("Tracking output must end in .vtu or .vtp")
    writer.SetFileName(str(filename))
    writer.SetHeaderTypeToUInt64()
    if not writer.Write():
        raise OSError("Could not write " + str(filename))


def compute_tracking(input, output, **options):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(input))
    reader.Update()
    raw = reader.GetOutput()
    names = sorted(raw.GetPointData().GetArrayName(i)
                   for i in range(raw.GetPointData().GetNumberOfArrays())
                   if not raw.GetPointData().GetArrayName(i).endswith("_Order"))
    image = vtk.vtkImageData()
    image.ShallowCopy(raw)
    times = []
    # Preserve input scalar/FieldData names, as in the original reader pipeline.
    for name in names:
        time = raw.GetFieldData().GetArray(name)
        if time is None or time.GetNumberOfTuples() != 1:
            raise ValueError("Missing physical time in FieldData for " + name)
        times.append(time.GetTuple1(0))
    mesh, settings = track_image(image, times, **options)
    write_mesh(mesh, output)
    print(f"Saved {mesh.GetNumberOfCells()} trajectory segments to {output}; settings: {settings}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_name", required=True)
    parser.add_argument("-o", "--output_name", required=True)
    parser.add_argument("--ttk-plugin")
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--persistence", type=float, default=0.)
    parser.add_argument("--destruction-cost", type=float)
    parser.add_argument("--spatial-weight", type=float)
    geometry = parser.add_mutually_exclusive_group()
    geometry.add_argument("--tetrahedralize", dest="tetrahedralize", action="store_true",
                          help="Use the original explicit triangulation (default)")
    geometry.add_argument("--implicit-grid", dest="tetrahedralize", action="store_false",
                          help="Opt into implicit triangulation; results can differ")
    parser.set_defaults(tetrahedralize=True)
    args = parser.parse_args()
    compute_tracking(args.input_name, args.output_name, plugin=args.ttk_plugin,
                     threads=args.threads, persistence=args.persistence,
                     destruction_cost=args.destruction_cost, tetrahedralize=args.tetrahedralize,
                     spatial_weight=args.spatial_weight)


if __name__ == "__main__":
    main()
