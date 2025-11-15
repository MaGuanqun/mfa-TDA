from paraview.simple import *
from paraview.servermanager import Fetch

import argparse
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

from scipy.ndimage import median_filter, gaussian_filter

LoadPlugin(
    "/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/"
    "lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so",
    remote=False, ns=globals()
)

def read_vtk_dataset(filename):
    """Read .vtk (legacy) or .vti (XML)."""
    if filename.endswith(".vti"):
        reader = vtk.vtkXMLImageDataReader()
    else:
        reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def smooth_xy_slices_numpy(image_data, array_name,
                           median_radius_xy=1,
                           sigma_xy=2.0):
    """
    For a scalar array on a regular 3D grid, apply per-(x,y) slice filtering:
      1) 2D median filter (radius median_radius_xy) on each z slice
      2) 2D Gaussian filter with sigma_xy on each z slice

    image_data : vtkImageData (regular grid)
    array_name : name of the point-data array to smooth
    median_radius_xy : radius in x,y for median (kernel = 2*radius+1); 0 to disable
    sigma_xy   : Gaussian sigma in x,y (voxel units); 0 to disable
    """
    pd = image_data.GetPointData()
    vtk_arr = pd.GetArray(array_name)
    if vtk_arr is None:
        print(f"[WARN] Array '{array_name}' not found, skipping.")
        return

    nx, ny, nz = image_data.GetDimensions()  # (nx, ny, nz)
    npts = image_data.GetNumberOfPoints()

    flat = vtk_to_numpy(vtk_arr)  # shape: (npts,)
    if flat.size != npts:
        raise RuntimeError(
            f"Array '{array_name}' size {flat.size} != number of points {npts}"
        )

    # Reshape to (z, y, x); consistent with your earlier scripts
    vol = flat.reshape((nz, ny, nx))

    # ---- 1) Median per slice (no mixing along z) ----
    if median_radius_xy > 0:
        k = 2 * median_radius_xy + 1
        vol = median_filter(
            vol,
            size=(1, k, k),     # (z,y,x) -> 1 in z, k×k in y,x
            mode='nearest'
        )

    # ---- 2) Gaussian per slice (no mixing along z) ----
    if sigma_xy > 0.0:
        vol = gaussian_filter(
            vol,
            sigma=(0.0, sigma_xy, sigma_xy),   # (z,y,x)
            mode='nearest'
        )

    # Flatten back and put into VTK
    flat_smoothed = vol.reshape(-1).astype(flat.dtype)
    vtk_smoothed = numpy_to_vtk(flat_smoothed, deep=True)
    vtk_smoothed.SetName(array_name)

    # Replace the original array with the smoothed one
    pd.RemoveArray(array_name)
    pd.AddArray(vtk_smoothed)


# ----------------- main -----------------

parser = argparse.ArgumentParser(
    description='Median+Gaussian smoothing on each 2D xy slice of a regular 3D dataset'
)

parser.add_argument(
    '-i', '--input_vti', type=str, default='file_name.vti',
    help='input file to compute critical points tracking'
)
parser.add_argument(
    '-o', '--output_vtk', type=str, default='file_name.vtk',
    help='output file to compute critical points tracking'
)
parser.add_argument(
    '-b', '--out_bin', type=str, default='output_volume.bin',
    help='output binary file for 3D volume'
)
parser.add_argument(
    '-s', '--sigma', type=float, default=2.0,
    help='Gaussian sigma in xy (in voxel units) for per-slice smoothing; 0 disables Gaussian'
)
parser.add_argument(
    '--median_radius', type=int, default=1,
    help='median filter radius in xy (kernel = 2*radius+1); 0 disables median'
)

args = parser.parse_args()

# === 1. Read the regular-grid dataset as pure VTK ===
dataset = read_vtk_dataset(args.input_vti)
if not isinstance(dataset, vtk.vtkImageData):
    raise RuntimeError(
        f"Input file '{args.input_vti}' is not vtkImageData; "
        "per-slice xy smoothing assumes a regular grid."
    )

print(f"Loaded dataset from '{args.input_vti}'")
dims = dataset.GetDimensions()
print(f"Dimensions = {dims}, number of points = {dataset.GetNumberOfPoints()}")

# === 2. Smooth each point-data array per xy slice ===
pd = dataset.GetPointData()
n_arrays = pd.GetNumberOfArrays()
print(f"Number of point-data arrays = {n_arrays}")

sigma_xy = args.sigma
median_radius_xy = args.median_radius

array_names = []
for i in range(n_arrays):
    arr = pd.GetArray(i)
    if arr is None:
        continue
    name = arr.GetName() or f"array_{i}"
    array_names.append(name)

print("Arrays before filtering:")
for name in array_names:
    print(" ", name)

# Decide which arrays to smooth:
arrays_to_smooth = []
for name in array_names:
    # Skip "order" and validity masks, etc.
    if name.endswith("_Order"):
        print(f"  Skipping '{name}' (order array)")
        continue
    if name == 'vtkValidPointMask':
        print(f"  Skipping '{name}' (valid point mask)")
        continue
    arrays_to_smooth.append(name)

print("\nSmoothing these arrays (per xy slice):")
for name in arrays_to_smooth:
    print(" ", name)

for name in arrays_to_smooth:
    print(
        f"  Smoothing '{name}' with "
        f"median_radius_xy={median_radius_xy}, sigma_xy={sigma_xy} ..."
    )
    smooth_xy_slices_numpy(
        dataset, name,
        median_radius_xy=median_radius_xy,
        sigma_xy=sigma_xy
    )

print("\nFinished smoothing all selected arrays.")

# === 3. Write out the smoothed dataset ===
if args.output_vtk.endswith(".vti"):
    writer = vtk.vtkXMLImageDataWriter()
else:
    writer = vtk.vtkDataSetWriter()

writer.SetFileName(args.output_vtk)
writer.SetInputData(dataset)
writer.Write()
print(f"Wrote smoothed dataset to '{args.output_vtk}'")

# === 4. Export concatenated binary volume ===
dataset = read_vtk_dataset(args.output_vtk)
pd = dataset.GetPointData()
n_arrays = pd.GetNumberOfArrays()

if n_arrays == 0:
    raise RuntimeError("No point-data arrays found in the VTK file.")

all_data = []
num = 0
for i in range(n_arrays):
    arr = pd.GetArray(i)
    if arr is None:
        continue
    name = arr.GetName() or f"array_{i}"

    # Skip arrays like "0000_Order"
    if name.endswith("_Order"):
        print(f"  Skipping '{name}'")
        continue
    if name == 'vtkValidPointMask':
        continue

    flat = vtk_to_numpy(arr).astype(np.float32)

    print(f"  Saving array '{name}' length={flat.size}")
    all_data.append(flat)
    num += 1

all_data = np.concatenate(all_data).astype(np.float32)
print(num * flat.size, num)

# Save whole [z,y,x] volume as pure float32 binary
all_data.tofile(args.out_bin)
print(f"\nSaved 3D volume to '{args.out_bin}'")
print(f"Shape (total,) = {all_data.shape}, dtype = {all_data.dtype}")