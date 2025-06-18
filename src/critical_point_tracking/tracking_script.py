from paraview.simple import *
from paraview import servermanager as sm
import argparse
import vtk.util.numpy_support as VN

LoadPlugin("/home/guanqunma/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/lib/paraview-5.11/plugins/TopologyToolKit/TopologyToolKit.so", remote=False, ns=globals()) 


def compute_tracking(input, output):
    timeTrackingvti = XMLImageDataReader(FileName=[input])
    # timeTrackingvti.CellArrayStatus = []

    all_point_arrays = list(timeTrackingvti.PointData.keys())
    timeTrackingvti.PointArrayStatus = all_point_arrays
    timeTrackingvti.TimeArray = 'None'
    # ["{:0>3}".format(i) for i in range(0, 20, 1)]
    timeTrackingvti.UpdatePipeline()
    precond = TTKArrayPreconditioning(Input=timeTrackingvti)
    precond.PointDataArrays = all_point_arrays
    precond.UpdatePipeline()
    tetrahedralize1 = Tetrahedralize(registrationName='Tetrahedralize1', Input=precond)
    
    Points_fetched = sm.Fetch(tetrahedralize1)
    value0 = VN.vtk_to_numpy(Points_fetched.GetFieldData().GetArray(all_point_arrays[0]))[0]
    value1 = VN.vtk_to_numpy(Points_fetched.GetFieldData().GetArray(all_point_arrays[1]))[0]
    
    tTKTrackingFromFields1 = TTKTrackingFromFields(Input=tetrahedralize1)
    tTKTrackingFromFields1.ForceZtranslation = 1
    tTKTrackingFromFields1.ZTranslation = value1-value0
    tTKTrackingFromFields1.Persistencethreshold = 0.0
    # tTKTrackingFromFields1.Fweight = 1
    print(f"Output data type: {tTKTrackingFromFields1.GetDataInformation().GetDataSetTypeAsString()}")
    
    SaveData(output, proxy=OutputPort(tTKTrackingFromFields1, 1))  # Table
    # SaveData(output, proxy=OutputPort(tTKTrackingFromFields1, 1), options={'WriteLines': True})
    
    
parser = argparse.ArgumentParser(description='TTK-critical points.')


parser.add_argument('-i', '--input_name', type=str, default='file_name.vti', help='input file to compute critical points tracking')
parser.add_argument('-o', '--output_name', type=str, default='file_name.vtp', help='output file to compute critical points tracking')


args = parser.parse_args()


input_file=args.input_name
output_file=args.output_name

compute_tracking(input_file,output_file)