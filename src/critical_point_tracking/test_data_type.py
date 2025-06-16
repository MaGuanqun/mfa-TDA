import vtk

def get_vti_scalar_array_names(filename):
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(filename)
    reader.Update()
    image_data = reader.GetOutput()
    
    point_data = image_data.GetPointData()
    array_names = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]
    return array_names

# Compare array names in both files
sample_file = "/mnt/data/timeTracking.vti"
user_file = "/mnt/data/rotating_gaussian.mfa.vti"

sample_arrays = get_vti_scalar_array_names(sample_file)
user_arrays = get_vti_scalar_array_names(user_file)

sample_arrays, user_arrays