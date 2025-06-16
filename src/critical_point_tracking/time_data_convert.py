import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import argparse




def convert_structured_grid_to_vti(input_vtk, output_vti):
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(input_vtk)
    reader.Update()
    grid = reader.GetOutput()

    dims = grid.GetDimensions()
    points = grid.GetPoints()
    pd = grid.GetPointData()

    array = pd.GetArray(0)
    array_name = array.GetName()
    # Estimate spacing assuming uniform grid
    def get_spacing(index_a, index_b):
        return np.linalg.norm(np.array(points.GetPoint(index_b)) - np.array(points.GetPoint(index_a)))

    dx = get_spacing(0, 1) if dims[0] > 1 else 1.0
    dy = get_spacing(0, dims[0]) if dims[1] > 1 else 1.0
    dz = get_spacing(0, dims[0]*dims[1]) if dims[2] > 1 else 1.0

    spacing = (dx, dy, dz)
    origin = points.GetPoint(0)

    # # Assume each array is one time slice: "000", "001", ...
    # array_names = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
    # array_names = sorted([name for name in array_names if name is not None])
    np_array = vtk_to_numpy(array).reshape((dims[2], dims[1], dims[0]))
    

    # slices = []
    # for name in array_names:
    #     array = pd.GetArray(name)
    #     np_array = vtk_to_numpy(array).reshape((dims[2], dims[1], dims[0]))  # (z, y, x)
    #     for z in range(dims[2]):
    #         slice_2d = np_array[z, :, :]  # (y, x)
    #         slices.append(slice_2d)

    # Stack into (x, y, t)
    # volume = np.stack(slices, axis=-1)  # shape: (ny, nx, nt)
    # volume = np.transpose(volume, (1, 0, 2))  # -> (nx, ny, nt)

    # Create vtkImageData
    image_data = vtk.vtkImageData()
    image_data.SetDimensions((dims[0], dims[1], 1))
    image_data.SetSpacing((dx, dy, 1.0))
    image_data.SetOrigin(origin)

    # For each array (which contains a full 3D volume):
    for z in range(dims[2]):
        slice_2d = np_array[z, :, :]  # shape (ny, nx)
        flat = slice_2d.T.ravel(order='F')  # transpose to (x, y) then flatten

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(f"{z:03d}")  # "000", "001", ...
        image_data.GetPointData().AddArray(vtk_arr)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_vti)
    writer.SetInputData(image_data)
    writer.SetCompressorTypeToZLib()                # Add this line
    writer.SetDataModeToAppended()                  # Optional, makes sure appended format
    writer.SetEncodeAppendedData(0)                 # 0 = raw (not base64), matches your sample
    writer.SetHeaderTypeToUInt64()
    writer.Write()
    
    
parser = argparse.ArgumentParser(description='TTK-critical points.')


parser.add_argument('-i', '--input_name', type=str, default='file_name.ply', help='file to compute critical points')
parser.add_argument('-o','--output_name', type=str, default='ttk.csv', help='output csv name')


args = parser.parse_args()


input_file=args.input_name
output_file=args.output_name

convert_structured_grid_to_vti(input_file, output_file)