import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import argparse




def convert_binary_file_to_vti(input_bin, output_vti, dims,min,max):
    data =np.fromfile(input_bin, dtype=np.float64)
    if data.size != dims[0] * dims[1] * dims[2]:
        print(data.size)
        raise ValueError("Data size does not match the provided dimensions.")


    np_array = data.reshape((dims[2], dims[1], dims[0]))
    

    distance=np.array(max)-np.array(min)
    
    dx= distance[0]/(dims[0]-1) if dims[0] > 1 else 1.0
    dy= distance[1]/(dims[1]-1) if dims[1] > 1 else 1.0
    z_value = distance[2]/(dims[2]-1) if dims[2] > 1 else 1.0
    # Create vtkImageData
    image_data = vtk.vtkImageData()
    image_data.SetDimensions((dims[0], dims[1], 1))
    image_data.SetSpacing((dx, dy, 1.0))
    image_data.SetOrigin(min)

    # For each array (which contains a full 3D volume):
    for z in range(dims[2]):
        slice_2d = np_array[z, :, :]  # shape (ny, nx)
        flat = slice_2d.T.ravel(order='F')  # transpose to (x, y) then flatten

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(f"{z:03d}")  # "000", "001", ...
        image_data.GetPointData().AddArray(vtk_arr)
        
        z_array = vtk.vtkDoubleArray()
        z_array.SetName(f"{z:03d}")
        z_array.SetNumberOfComponents(1)
        # z_array.InsertNextValue(3.0)
        z_array.InsertNextValue(min[2]+z_value*z)
        image_data.GetFieldData().AddArray(z_array)



    

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_vti)
    writer.SetInputData(image_data)
    writer.SetCompressorTypeToZLib()                # Add this line
    writer.SetDataModeToAppended()                  # Optional, makes sure appended format
    writer.SetEncodeAppendedData(0)                 # 0 = raw (not base64), matches your sample
    writer.SetHeaderTypeToUInt64()
    writer.Write()
    
    
parser = argparse.ArgumentParser(description='TTK-critical points.')


parser.add_argument('-i', '--input_name', type=str, default='file_name.dat', help='file to compute critical points')
parser.add_argument('-o','--output_name', type=str, default='ttk.csv', help='output csv name')
parser.add_argument('-f', '--function', type=str, default='rotating_gaussian', help='function')

args = parser.parse_args()


input_file=args.input_name
output_file=args.output_name
function=args.function

if function == 'rotating_gaussian':
    min = np.array([-2.0, -2.0, 0.0])
    max =np.array([2.0, 2.0, 4.0])
    dim = np.array([100,100,100])
elif function == 'quartic_potential' or function == 'quartic_potential_2':
    min = np.array([-2.0, -2.0, 0.0])
    max =np.array([2.0, 2.0, 4.0])
    dim = np.array([200,200,100])
elif function == 'vortex_street':
    dim = np.array([100, 80, 50])
    min = np.array([0.0, 0.0, 0.0])
    max = dim - 1
    
convert_binary_file_to_vti(input_file, output_file,dim,min,max)