
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import argparse




def convert_binary_file_to_vti(input_bin, output_vti, dims,min,max,args):
    
    dtype = np.float32 if args.float_type == 'float32' else np.float64
    data = np.fromfile(input_bin, dtype=dtype)
    if data.size != dims[0] * dims[1] * dims[2]:
        print(data.size)
        print(dims[0] * dims[1] * dims[2])
        raise ValueError("Data size does not match the provided dimensions.")


    # Original binary is already stored with x varying fastest (x, then y, then z).
    # Reshape directly to (nx, ny, nz) and pass to VTK.
    volume_xyz = data.reshape((int(dims[0]), int(dims[1]), int(dims[2])), order='F')

    distance = np.array(max) - np.array(min)

    dx = distance[0] / (dims[0] - 1) if dims[0] > 1 else 1.0
    dy = distance[1] / (dims[1] - 1) if dims[1] > 1 else 1.0
    dz = distance[2] / (dims[2] - 1) if dims[2] > 1 else 1.0

    # Create full 3D vtkImageData
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(int(dims[0]), int(dims[1]), int(dims[2]))
    image_data.SetSpacing(dx, dy, dz)
    image_data.SetOrigin(min)

    # VTK expects data ordered with x varying fastest, then y, then z.
    flat = volume_xyz.ravel(order='F')

    vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
    vtk_arr.SetName("scalars")
    image_data.GetPointData().SetScalars(vtk_arr)


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
parser.add_argument('--float_type', type=str, default='float32', help='float32 or float64')
parser.add_argument('--step_size', type=int, default=1, help='step size between frames')

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
    dim = np.array([100,100,100])
elif function == 'vortex_street':
    dim = np.array([100, 80, 50])
    min = np.array([0.0, 0.0, 0.0])
    max = np.array([99, 79, 49])
elif function =='vortex_street_3d':
    dim = np.array([80,10,15])
    min = np.array([-0.5, -0.5, 13.5])
    max = np.array([7.5, 0.5, 15.0])
    # min = np.array([-0.5, -0.5, 0.0])
    # max = np.array([0.5, 7.5, 15])
elif function =='boussinesq_3d':
    dim = np.array([10,30,15])
    min = np.array([-0.5, -0.5, 0.0])
    max = np.array([0.5, 2.5, 1.5])
elif function =='fluid':
    dim = np.array([10,10,10])
    min = np.array([0.0, 0.0, 0.0])
    max = np.array([1.0, 1.0, 1.0])
elif function =='cylinder':
    dim = np.array([40,10,10])
    min = np.array([1.5, 0.5, 0.0])
    max = np.array([5.5,1.5,1])
elif function =='cylinder2':
    dim = np.array([23,10,10])
    min = np.array([3.2, 0.5, 0.0])
    max = np.array([5.5,1.5,1])
elif function =='cylinder3':
    dim = np.array([20,10,10])
    min = np.array([3.5, 0.5, 0.0])
    max = np.array([5.5,1.5,5.0])
elif function =='vortex':
    dim = np.array([15,15,15])
    min = np.array([0.0, 0.0, 0.0])
    max = np.array([127.0, 127.0, 127.0])

size = args.step_size*dim
    
convert_binary_file_to_vti(input_file, output_file,size,min,max,args)