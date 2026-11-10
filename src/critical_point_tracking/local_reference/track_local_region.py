"""Read local-field-v1 samples, save a VTI, and run TTK in physical coordinates."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tracking_script import track_image, write_mesh


def load_sample(prefix):
    with open(str(prefix) + ".json") as stream:
        info = json.load(stream)
    if (info.get("format") != "local-field-v1" or info.get("dtype") != "<f8"
            or info.get("layout") != "t,y,x"):
        raise ValueError("Unsupported sample metadata format")
    dims = info["dimensions"]
    if len(dims) != 3 or any(type(n) is not int or n < 2 for n in dims):
        raise ValueError("Each sample dimension must contain at least two points")
    nx, ny, nt = dims
    if nt < 3:
        raise ValueError("At least three frames are needed for legacy TTK tracking")
    lo, hi = np.array(info["sample_min"]), np.array(info["sample_max"])
    if lo.shape != (3,) or hi.shape != (3,) or not np.isfinite([lo, hi]).all() or np.any(hi <= lo):
        raise ValueError("Invalid physical sample bounds")
    spacing = (hi-lo)/(np.array(dims)-1)
    if not np.allclose(info["spacing"], spacing, rtol=1e-12, atol=1e-14):
        raise ValueError("Sample spacing is inconsistent with dimensions/bounds")
    binary = Path(str(prefix) + ".bin")
    if binary.stat().st_size != nx*ny*nt*8:
        raise ValueError("Binary sample size does not match metadata")
    data = np.memmap(binary, dtype="<f8", mode="r", shape=(nt,ny,nx))
    image = vtk.vtkImageData()
    image.SetDimensions(nx,ny,1)
    # The original script adds i*dt to the input plane's z coordinate.
    image.SetOrigin(float(lo[0]),float(lo[1]),float(lo[2]))
    image.SetSpacing(float(spacing[0]),float(spacing[1]),1.)
    times = np.linspace(lo[2],hi[2],nt)
    for i, time in enumerate(times):
        flat = data[i].reshape(-1)
        if not np.isfinite(flat).all():
            raise ValueError(f"Nonfinite scalar samples at time index {i}")
        # Retain the memmap while VTK shares its buffers.
        arr = numpy_to_vtk(flat, deep=False)
        name = f"step_{i:09d}"
        arr.SetName(name)
        image.GetPointData().AddArray(arr)
        time_arr = vtk.vtkDoubleArray()
        time_arr.SetName(name)
        time_arr.InsertNextValue(float(time))
        image.GetFieldData().AddArray(time_arr)
    return info, data, image, times


def annotate(mesh, info, times):
    if mesh.GetNumberOfPoints():
        coords = vtk_to_numpy(mesh.GetPoints().GetData())
        lo, hi = np.array(info["region_min"]), np.array(info["region_max"])
        # TTK may represent spatial positions internally using float32.
        tol = 8*np.finfo(np.float32).eps*np.maximum(1.,np.abs(coords))
        inside = np.all((coords >= lo-tol) & (coords <= hi+tol),axis=1).astype(np.uint8)
        bounds_lo, bounds_hi = np.array(info["sample_min"]), np.array(info["sample_max"])
        boundary = np.any((np.abs(coords[:,:2]-bounds_lo[:2]) <= tol[:,:2]) |
                          (np.abs(coords[:,:2]-bounds_hi[:2]) <= tol[:,:2]),axis=1).astype(np.uint8)
        indices = np.rint((coords[:,2]-times[0])/(times[1]-times[0])).astype(np.int64)
        if np.any(indices < 0) or np.any(indices >= len(times)):
            raise ValueError("TTK trajectory time lies outside the sampled interval")
    else:
        inside = boundary = np.empty(0,dtype=np.uint8)
        indices = np.empty(0,dtype=np.int64)
    for name, values in (("InsideRequestedRegion",inside),("OnSampleSpatialBoundary",boundary),
                         ("TimeIndex",indices),("Time",np.asarray(times)[indices])):
        arr = numpy_to_vtk(values,deep=True)
        arr.SetName(name)
        mesh.GetPointData().AddArray(arr)
    return int(inside.sum())


def write_points(mesh, filename):
    pd = mesh.GetPointData()
    names = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())
             if pd.GetArray(i).GetNumberOfComponents() == 1]
    with open(filename,"w",newline="") as stream:
        writer = csv.writer(stream)
        # Existing matching_ratio.py reads columns 0,1,2 as coordinates.
        writer.writerow(["PositionX","PositionY","PositionZ","PointId"]+names)
        for i in range(mesh.GetNumberOfPoints()):
            writer.writerow([*mesh.GetPoint(i),i]+[pd.GetArray(name).GetTuple1(i) for name in names])


def save_tracking_outputs(mesh, prefix, info, times):
    """Write native trajectories first; diagnostics belong to a separate CSV."""
    write_mesh(mesh,prefix+"_tracking.vtu")
    write_mesh(mesh,prefix+"_tracking.vtp")
    diagnostics = vtk.vtkUnstructuredGrid()
    diagnostics.ShallowCopy(mesh)
    interior = annotate(diagnostics,info,times)
    write_points(diagnostics,prefix+"_tracking_points.csv")
    return interior


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample",required=True,help="Prefix of .bin/.json files")
    parser.add_argument("--ttk-plugin")
    parser.add_argument("--threads",type=int,default=3)
    parser.add_argument("--persistence",type=float,default=0.)
    parser.add_argument("--destruction-cost",type=float)
    parser.add_argument("--spatial-weight",type=float,help="Override TTK X/Y matching weights (default: original plugin defaults)")
    geometry = parser.add_mutually_exclusive_group()
    geometry.add_argument("--tetrahedralize",dest="tetrahedralize",action="store_true")
    geometry.add_argument("--implicit-grid",dest="tetrahedralize",action="store_false",
                          help="Opt into implicit triangulation; results can differ")
    parser.set_defaults(tetrahedralize=True)
    parser.add_argument("--skip-vti",action="store_true",help="Omit the optional VTI copy")
    args = parser.parse_args()
    info, mapped_data, image, times = load_sample(args.sample)
    if not args.skip_vti:
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(args.sample+".vti")
        writer.SetInputData(image)
        writer.SetHeaderTypeToUInt64()
        if not writer.Write():
            raise OSError("Could not write sampled VTI")
    mesh, settings = track_image(image,times,plugin=args.ttk_plugin,threads=args.threads,
        persistence=args.persistence,destruction_cost=args.destruction_cost,
        tetrahedralize=args.tetrahedralize,spatial_weight=args.spatial_weight)
    interior = save_tracking_outputs(mesh,args.sample,info,times)
    summary = {"sample": str(Path(args.sample).resolve()), "ttk_settings": settings,
               "tetrahedralize": args.tetrahedralize, "threads": args.threads,
               "points": mesh.GetNumberOfPoints(), "segments": mesh.GetNumberOfCells(),
               "points_inside_requested_region": interior,
               "time_range": [float(times[0]),float(times[-1])]}
    with open(args.sample+"_tracking.json","w") as stream:
        json.dump(summary,stream,indent=2)
    print(json.dumps(summary,indent=2))
    if not mesh.GetNumberOfCells():
        print("No trajectory segments were found in this region with these TTK settings.")


if __name__ == "__main__":
    main()
