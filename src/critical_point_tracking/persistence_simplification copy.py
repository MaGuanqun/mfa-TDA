
from paraview.simple import *
from paraview.servermanager import Fetch


import argparse
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals()) 


def simplify_all_arrays_to_vtk(tetra, array_names, threshold, output_vtk, image_data):
    """
    image_data : vtkImageData with all arrays already added (0000, 0001, …)
    array_names : list of array names to simplify, e.g. ['0000', '0001', ...]
    threshold : same PersistenceThreshold you used in GUI (0.1 etc)
    output_vtk : path to .vtk file (unstructured grid) to save
    """

    # Wrap the in-memory vtkImageData as a ParaView source

    # 2) Chain TTKTopologicalSimplificationByPersistence filters, one per array
    current = tetra
    for name in array_names:
        ttk = TTKTopologicalSimplificationByPersistence(
            registrationName=f"TTKTopologicalSimplificationByPersistence_{name}",
            Input=current
        )
        ttk.InputArray = ['POINTS', name]
        ttk.PersistenceThreshold = threshold        # 0.1 etc
        ttk.ThresholdIsAbsolute = 0                 # same as in GUI trace
        ttk.UpdatePipeline()
        current = ttk   # next iteration works on already-simplified data

    resampled = ResampleWithDataset(
        current,
        registrationName="ResampleWithDataset_SimplifiedToImage",
        DestinationMesh=image_data   # destination grid / sampling grid
    )
    resampled.UpdatePipeline()
    # 3) Save the final dataset (with all arrays simplified) to .vtk
    SaveData(
        output_vtk,
        proxy=resampled,
        PointDataArrays=array_names  # or add *_Order etc if you want them too
    )
    
def persistence_simplification(tetra, array_name, threshold=0.1):
    """
    Run TTKTopologicalSimplificationByPersistence on a POINTS scalar array
    in a vtkImageData and return a new vtkDataArray with simplified values.
    """

    # Wrap the vtkImageData in a ParaView source
    # producer = TrivialProducer()
    # producer.GetClientSideObject().SetOutput(image_data)

    
    ttk_simpl = TTKTopologicalSimplificationByPersistence(
        Input=tetra
    )

    # Same as in your GUI: use point data array by name
    ttk_simpl.InputArray = ['POINTS', array_name]
    ttk_simpl.PairType = 'Extremum-Saddle'
    # Relative persistence threshold (0 = relative, 1 = absolute)
    ttk_simpl.PersistenceThreshold = threshold
    ttk_simpl.ThresholdIsAbsolute = 0

    # Execute the pipeline
    ttk_simpl.UpdatePipeline()

    # Fetch result back as VTK dataset
    simplified_vtk = servermanager.Fetch(ttk_simpl)

    # Extract simplified scalar
    arr = simplified_vtk.GetPointData().GetArray(array_name)
    if arr is None:
        raise RuntimeError(f"Simplified array '{array_name}' not found in TTK output")

    # Deep copy so we can safely attach it to the original image_data
    out_array = arr.NewInstance()
    out_array.DeepCopy(arr)
    out_array.SetName(array_name)

    return out_array


def read_vtk_dataset(filename):
    """Read .vtk (legacy) or .vti (XML)."""
    if filename.endswith(".vti"):
        reader = vtk.vtkXMLImageDataReader()
    else:
        reader = vtk.vtkDataSetReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def infer_nx_ny(dataset, npts):
    """
    Infer nx, ny from vtkImageData if possible.
    Otherwise, try sqrt(npts) as fallback.
    """
    if isinstance(dataset, vtk.vtkImageData):
        nx, ny, nz = dataset.GetDimensions()
        # Your case: each array is an xy slice, so we expect npts == nx * ny
        if nx * ny == npts:
            return nx, ny

    # Fallback: assume square
    root = int(round(npts ** 0.5))
    if root * root == npts:
        return root, root

    raise RuntimeError(
        f"Cannot infer (nx, ny). npts = {npts} is not nx*ny from image "
        "and not a perfect square. Please modify script to pass nx, ny."
    )


parser = argparse.ArgumentParser(description='Apply TTK persistence simplification to every 2D slice')

parser.add_argument('-i', '--input_vti', type=str, default='file_name.vti', help='input file to compute critical points tracking')
parser.add_argument('-o', '--output_vtk', type=str, default='file_name.vtk', help='output file to compute critical points tracking')
parser.add_argument('-b', '--out_bin', type=str, default='output_volume.bin', help='output binary file for 3D volume')
parser.add_argument('-t', '--threshold', type=float, default=0.1, help='persistence threshold for simplification')
args = parser.parse_args()

 # === 1. Read VTI as a ParaView proxy (not pure VTK) ===
reader = XMLImageDataReader(FileName=[args.input_vti])
reader.UpdatePipeline()

# === 2. Tetrahedralize (proxy) ===
tetra = Tetrahedralize(Input=reader)
tetra.UpdatePipeline()

# === 3. Inspect arrays via Fetch to decide which ones to simplify ===
vtk_tetra = Fetch(tetra)
point_data = vtk_tetra.GetPointData()
num_arrays = point_data.GetNumberOfArrays()


# Cache original array names before we start modifying them
array_names = []
for i in range(num_arrays):
    name = point_data.GetArrayName(i)
    if not name:
        continue
    arr = point_data.GetArray(i)
    if arr is None:
        continue
    if arr.GetNumberOfComponents() != 1:
        continue  # only scalar arrays
    array_names.append(name)

print("Will simplify the following scalar point-data arrays:")
for name in array_names:
    print("  ", name)

# === 4. Run TTK for all arrays and write a single .vtk ===
simplify_all_arrays_to_vtk(
    tetra=tetra,
    array_names=array_names,
    threshold=args.threshold,
    output_vtk=args.output_vtk,
    image_data = reader
    
)


dataset = read_vtk_dataset(args.output_vtk)
pd = dataset.GetPointData()
n_arrays = pd.GetNumberOfArrays()

if n_arrays == 0:
    raise RuntimeError("No point-data arrays found in the VTK file.")

all_data = []

num = 0
for i in range(n_arrays):
    arr = pd.GetArray(i)
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
print(num*flat.size,num)


# Save whole [z,y,x] volume as pure float32 binary
all_data.tofile(args.out_bin)
print(f"\nSaved 3D volume to '{args.out_bin}'")
print(f"Shape (z, y, x) = {all_data.shape}, dtype = {all_data.dtype}")