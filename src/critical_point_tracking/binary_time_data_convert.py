
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import argparse



def persistence_simplification(image_data, array_name, threshold=0.5):

    # create a new 'Tetrahedralize'
    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(image_data)
    # tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=producer)
    
    # tetrahedralize1.UpdatePipeline()
    
    # np_arr = vtk_to_numpy(image_data.GetPointData().GetArray(array_name))
    # print("1",np_arr,array_name)

    # create a new 'TTK TopologicalSimplificationByPersistence'
    tTKTopologicalSimplificationByPersistence1 = TTKTopologicalSimplificationByPersistence(registrationName='TTKTopologicalSimplificationByPersistence1', Input=producer)
    tTKTopologicalSimplificationByPersistence1.InputArray = ['POINTS', array_name]

    # Properties modified on tTKTopologicalSimplificationByPersistence1
    tTKTopologicalSimplificationByPersistence1.PersistenceThreshold = threshold
    tTKTopologicalSimplificationByPersistence1.ThresholdIsAbsolute = 0
    
    tTKTopologicalSimplificationByPersistence1.UpdatePipeline()
    # create a new 'Clean to Grid'
    
    simplified_vtk = Fetch(tTKTopologicalSimplificationByPersistence1)
    
    arr = simplified_vtk.GetPointData().GetArray(array_name)
    if arr is None:
        raise RuntimeError(
            f"Simplified array '{array_name}' not found in TTK output"
        )
        
    # np_arr = vtk_to_numpy(arr)
    # print("2",np_arr,array_name)

    # Deep copy into a standalone vtkDataArray so we can attach to image_data
    out_array = arr.NewInstance()
    out_array.DeepCopy(arr)
    out_array.SetName(array_name)
    
    return out_array


def convert_binary_file_to_vti(input_bin, output_vti, dims,min,max,args):
    
    dtype = np.float32 if args.float_type == 'float32' else np.float64
    data = np.fromfile(input_bin, dtype=dtype)
    if data.size != dims[0] * dims[1] * dims[2]:
        print(data.size)
        print(dims[0] * dims[1] * dims[2])
        raise ValueError("Data size does not match the provided dimensions.")


    np_array = data.reshape((dims[2], dims[1], dims[0]))
    

    distance=np.array(max)-np.array(min)
    
    dx= distance[0]/(dims[0]-1) if dims[0] > 1 else 1.0
    dy= distance[1]/(dims[1]-1) if dims[1] > 1 else 1.0
    z_value = distance[2]/(dims[2]-1) if dims[1] > 1 else 1.0
    # Create vtkImageData
    # image_data = vtk.vtkImageData()
    # image_data.SetDimensions((dims[0], dims[1], 1))
    # image_data.SetSpacing((dx, dy, 1.0))
    # image_data.SetOrigin(min)
    
    image_data = vtk.vtkImageData()
    image_data.SetDimensions((dims[0], dims[1], 1))
    image_data.SetSpacing((dx, dy, 1.0))
    image_data.SetOrigin(min)

    # For each array (which contains a full 3D volume):
    for z in range(dims[2]):
        slice_2d = np_array[z, :, :]  # shape (ny, nx)
        flat = slice_2d.T.ravel(order='F')  # transpose to (x, y) then flatten

        vtk_arr = numpy_to_vtk(flat, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_arr.SetName(f"{z:04d}")  # "000", "001", ...
        image_data.GetPointData().AddArray(vtk_arr)
        
        # simplified=persistence_simplification(image_data, f"{z:04d}", threshold=0.1)
        # image_data.GetPointData().RemoveArray(f"{z:04d}")
        # image_data.GetPointData().AddArray(simplified)
        
        # np_arr = vtk_to_numpy(image_data.GetPointData().GetArray(f"{z:04d}"))
        # print(np_arr)
        
        # print("=====",f"{z:04d}")
        
        z_array = vtk.vtkDoubleArray()
        z_array.SetName(f"{z:04d}")
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
    min = np.array([0.0, 0.0, 1350.0])
    max = np.array([639.0, 79.0, 1500.0])
    # min = np.array([-0.5, -0.5, 0.0])
    # max = np.array([0.5, 7.5, 15])
elif function =='boussinesq_3d':
    dim = np.array([10,30,20])
    min = np.array([0.0, 0.0, 0.0])
    max = np.array([149, 449, 199])
elif function =='fluid':
    dim = np.array([10,10,10])
    min = np.array([0.0, 0.0, 0.0])
    max = np.array([1.0, 1.0, 1.0])
elif function =='cylinder':
    dim = np.array([40,10,10])
    min = np.array([1.5, 0.5, 0.0])
    max = np.array([5.5,1.5,1])

size = args.step_size*dim
    
convert_binary_file_to_vti(input_file, output_file,size,min,max,args)