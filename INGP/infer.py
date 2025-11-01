from __future__ import absolute_import, division, print_function
import argparse
import os
from utils import PSNR, create_path, make_coord_grid
from NGP import NGP_TCNN
import json
from dataset import Dataset
import torch
import numpy as np
import torch.nn.functional as F
from time import time
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
project_folder_path = os.path.dirname(os.path.abspath(__file__))
output_folder = os.path.join(project_folder_path, "Output")

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def sample_grid(model, grid, align_corners:bool = False,
                device:str="cuda", data_device:str="cuda", max_points:int = 100000):
    coord_grid = make_coord_grid(grid, 
        data_device, flatten=False,
        align_corners=align_corners)
    coord_grid_shape = list(coord_grid.shape)
    coord_grid = coord_grid.view(-1, coord_grid.shape[-1])
    vals = forward_maxpoints(model, coord_grid, 
                             max_points = max_points,
                             data_device=data_device,
                             device=device
                             )
    coord_grid_shape[-1] = -1
    vals = vals.reshape(coord_grid_shape)
    return vals

def forward_maxpoints(model, coords, out_dim=1, max_points=100000, 
                      data_device="cuda", device="cuda"):
    output_shape = list(coords.shape)
    output_shape[-1] = out_dim
    output = torch.empty(output_shape, 
        dtype=torch.float32, 
        device=data_device)
    
    for start in range(0, coords.shape[0], max_points):
        output[start:min(start+max_points, coords.shape[0])] = \
            model(coords[start:min(start+max_points, coords.shape[0])].to(device), scale_output=True).to(data_device)
    return output

def load_options(load_location):
    #print(load_location)
    if not os.path.exists(load_location):
        print("%s doesn't exist, load failed" % load_location)
        return
        
    if os.path.exists(os.path.join(load_location, "options.json")):
        with open(os.path.join(load_location, "options.json"), 'r') as fp:
            opt2 = json.load(fp)
    else:
        print("%s doesn't exist, load failed" % "options.json")
        return

    return opt2

def tensor_to_raw(tensor, path):
    v = tensor.squeeze().cpu()
    v = np.asarray(v,dtype='<f')
    v = v.flatten('F')
    v.tofile(path,format='<f')

def model_reconstruction(model, opt):
    
    # Load the reference data
    with torch.no_grad():
        result = sample_grid(model, opt['full_shape'], max_points=1000000,
                             align_corners=opt['align_corners'],
                             device=opt['device'],
                             data_device=opt['data_device'])
    result = result.to(opt['data_device'])
    result = result.permute(3, 0, 1, 2).unsqueeze(0)
    create_path(os.path.join(output_folder, "Reconstruction"))
    # ic(result.shape)
    # exit()
    tensor_to_raw(result, 
        os.path.join(output_folder, 
        "Reconstruction", opt['save_name']+".raw"))

def model_reconstruction_chunked(model, opt):
    
    chunk_size = 512
    full_shape = list(map(int, opt['data_dims'].split(',')))
    
    output = torch.empty(full_shape, 
        dtype=torch.float32, 
        device=opt['data_device']).unsqueeze(0).unsqueeze(0)
    
    with torch.no_grad():
        for z_ind in range(0, full_shape[0], chunk_size):
            z_ind_end = min(full_shape[0], z_ind+chunk_size)
            z_range = z_ind_end-z_ind
            for y_ind in range(0, full_shape[1], chunk_size):
                y_ind_end = min(full_shape[1], y_ind+chunk_size)
                y_range = y_ind_end-y_ind            
                for x_ind in range(0, full_shape[2], chunk_size):
                    x_ind_end = min(full_shape[2], x_ind+chunk_size)
                    x_range = x_ind_end-x_ind
                    
                    opt['extents'] = f"{z_ind},{z_ind_end},{y_ind},{y_ind_end},{x_ind},{x_ind_end}"
                    print(f"Extents: {z_ind},{z_ind_end},{y_ind},{y_ind_end},{x_ind},{x_ind_end}")
                                                                
                    grid = [z_range, y_range, x_range]
                    coord_grid = make_coord_grid(grid, 
                        opt['data_device'], flatten=True,
                        align_corners=True,
                        use_half=False)
                    
                    coord_grid += 1.0
                    coord_grid /= 2.0
                    
                    coord_grid[:,0] *= (x_range-1) / (full_shape[2]-1)
                    coord_grid[:,1] *= (y_range-1) / (full_shape[1]-1)
                    coord_grid[:,2] *= (z_range-1) / (full_shape[0]-1)
                    
                    coord_grid[:,0] += x_ind / (full_shape[2]-1)
                    coord_grid[:,1] += y_ind / (full_shape[1]-1)
                    coord_grid[:,2] += z_ind / (full_shape[0]-1)
                    
                    coord_grid *= 2.0
                    coord_grid -= 1.0
                    
                    out_tmp = forward_maxpoints(model, 
                                                coord_grid, max_points=2**20, 
                                                data_device=opt['data_device'],
                                                device=opt['device'])
                    out_tmp = out_tmp.permute(1,0)
                    out_tmp = out_tmp.view([out_tmp.shape[0]] + grid)
                    output[0,:,z_ind:z_ind_end,y_ind:y_ind_end,x_ind:x_ind_end] = out_tmp

                    print(f"Chunk {z_ind},{z_ind_end},{y_ind},{y_ind_end},{x_ind},{x_ind_end}")
        
    create_path(output_folder)
    tensor_to_raw(output, os.path.join(output_folder, opt['dataset_name']+".raw"))


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate a model on some tests')

    parser.add_argument('--load_from',default="./SavedModels/vortex_street",type=str,help="Model name to load")
    parser.add_argument('--device',default="cuda:0",type=str,
                        help="Device to load model to")
    parser.add_argument('--data_device',default="cuda:0",type=str,
                        help="Device to load data to")
    args = vars(parser.parse_args())
    
    # Load the model
    opt = load_options(args['load_from'])
    opt['device'] = args['device']
    opt['data_device'] = args['data_device']
   
    model = NGP_TCNN(opt)
    ckpt = torch.load(os.path.join(args['load_from'], 'model.ckpt.tar'), map_location = opt['device'], weights_only=False)  
    
    model.load_state_dict(ckpt['state_dict'])
    model = model.to(opt['device'])
    model.train(False)
    model.eval()
    # print(model)
    # exit()
    # Perform tests
    tic = time()
    model_reconstruction_chunked(model, opt),
    toc = time()
    print(f"Reconstruction time: {toc-tic}")
    
        
    
        



        

