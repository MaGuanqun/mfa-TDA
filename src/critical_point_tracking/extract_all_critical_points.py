from paraview.simple import *
from paraview import servermanager as sm
import vtk.util.numpy_support as VN
import argparse
import numpy as np
import csv

LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals()) 


def compute_tracking(input, output):
    timeTrackingvti = XMLImageDataReader(FileName=[input])
    # timeTrackingvti.CellArrayStatus = []

    all_point_arrays = list(timeTrackingvti.PointData.keys())
    timeTrackingvti.PointArrayStatus = all_point_arrays
    timeTrackingvti.TimeArray = 'None'
    
    

    timeTrackingvti.UpdatePipeline()
    precond = TTKArrayPreconditioning(Input=timeTrackingvti)
    # precond.PointDataArrays = all_point_arrays
    precond.UpdatePipeline()
    tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=precond)
    print(all_point_arrays)
    
    critical_points_pos = np.empty((0, 3))
    critical_points_type = np.array([]) 
    critical_points_boundary = np.array([])

    
    
    for array_id in range(len(all_point_arrays)):#
        tTKScalarFieldCriticalPoints1 = TTKScalarFieldCriticalPoints(registrationName='TTKScalarFieldCriticalPoints', Input=tetrahedralize1)
        tTKScalarFieldCriticalPoints1.ScalarField = ['POINTS', all_point_arrays[array_id]]
        tTKScalarFieldCriticalPoints1.InputOffsetField = ['POINTS', all_point_arrays[array_id]]
        # tTKScalarFieldCriticalPoints1.Withvertexscalars = 0
        # Update the pipeline to compute the critical points
        tTKScalarFieldCriticalPoints1.UpdatePipeline()
        
        Points_fetched = sm.Fetch(tTKScalarFieldCriticalPoints1)

        criPointArray =  VN.vtk_to_numpy(Points_fetched.GetPoints().GetData())
        # Read the field data for the current array
        field_data_array = VN.vtk_to_numpy(Points_fetched.GetFieldData().GetArray(all_point_arrays[array_id]))

        
        criPointArray[:, 2] = field_data_array[0]
        

        critical_points_pos = np.concatenate((critical_points_pos, criPointArray), axis=0)
        
        criTypeArray = VN.vtk_to_numpy(Points_fetched.GetPointData().GetArray("CriticalType"))   
             
        critical_points_type = np.concatenate((critical_points_type, criTypeArray), axis=0)
        
        criBoundaryArray = VN.vtk_to_numpy(Points_fetched.GetPointData().GetArray("IsOnBoundary")) 
        critical_points_boundary = np.concatenate((critical_points_boundary, criBoundaryArray), axis=0)


    mask = critical_points_boundary != 1
    filtered_critical_points_pos = critical_points_pos[mask]
    filtered_critical_points_type = critical_points_type[mask]
    with open(output, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["CriticalType", "PositionX", "PositionY", "PositionZ"])
        for i in range(filtered_critical_points_pos.shape[0]):
            csvwriter.writerow([int(filtered_critical_points_type[i]), 
                                filtered_critical_points_pos[i, 0], 
                                filtered_critical_points_pos[i, 1], 
                                filtered_critical_points_pos[i, 2]])

    
   
    
    
parser = argparse.ArgumentParser(description='TTK-critical points.')


parser.add_argument('-i', '--input_name', type=str, default='file_name.vti', help='input file to compute critical points tracking')
parser.add_argument('-o', '--output_name', type=str, default='file_name.vtp', help='output file to compute critical points tracking')


args = parser.parse_args()


input_file=args.input_name
output_file=args.output_name

compute_tracking(input_file,output_file)