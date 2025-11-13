from __future__ import absolute_import, division, print_function
import argparse
from dataset import Dataset
import datetime
from utils import str2bool, save_model
from NGP import NGP_TCNN
from fVSRN import fVSRN
import torch
import torch.optim as optim
from torch.nn import functional as F
import time
import os
# from option import *
import numpy as np
from torch.utils.data import DataLoader
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
project_folder_path = os.path.dirname(os.path.abspath(__file__))
output_folder = os.path.join(project_folder_path, "Output")

def second_derivative_xy_penalty(model, x, base_output=None, h=1e-2):
    """
    Approximate second derivatives f_xx and f_yy in the xy-plane using finite differences,
    and penalize their squared magnitude.

    - x: [B, 3] coordinates in [-1, 1]^3 (assumed order [x, y, z])
    - model: fVSRN model
    - base_output: optional precomputed f(x) of shape [B, 1] or [B]
    - h: finite-difference step size in normalized coordinate space
    """
    # x is just input coordinates; it does not require gradients.
    device = x.device
    dtype = x.dtype

    if base_output is None:
        base_output = model(x)

    # ensure h is a tensor on the right device / dtype
    h_tensor = torch.tensor(h, device=device, dtype=dtype)

    # Unit offsets in x and y (normalized domain)
    e_x = torch.zeros_like(x)
    e_y = torch.zeros_like(x)
    e_x[:, 0] = h_tensor      # perturb x
    e_y[:, 1] = h_tensor      # perturb y

    # Shifted coordinates, clamped to [-1, 1] to stay in domain
    x_plus  = (x + e_x).clamp(-1.0, 1.0)
    x_minus = (x - e_x).clamp(-1.0, 1.0)
    y_plus  = (x + e_y).clamp(-1.0, 1.0)
    y_minus = (x - e_y).clamp(-1.0, 1.0)

    # Batch all shifted coords into one forward pass for efficiency
    coords_other = torch.cat([x_plus, x_minus, y_plus, y_minus], dim=0)  # [4B, 3]
    vals_other = model(coords_other)  # [4B, 1] or [4B]

    B = x.shape[0]
    f_xplus, f_xminus, f_yplus, f_yminus = vals_other.split(B, dim=0)
    f = base_output

    # Finite-difference approximations
    # f_xx ≈ (f(x+h, y, z) - 2 f(x, y, z) + f(x-h, y, z)) / h^2
    # f_yy ≈ (f(x, y+h, z) - 2 f(x, y, z) + f(x, y-h, z)) / h^2
    denom = h_tensor * h_tensor
    f_xx = (f_xplus - 2.0 * f + f_xminus) / denom
    f_yy = (f_yplus - 2.0 * f + f_yminus) / denom

    # You can also add f_xy using cross terms if needed, but this is cheaper
    # and already encourages smoothness.
    penalty = (f_xx ** 2 + f_yy ** 2).mean()

    return penalty


def log_to_writer(iteration, losses, writer, opt, preconditioning=None):
    with torch.no_grad():   
        print_str = f"Iteration {iteration}/{opt['iterations']}, "
        for key in losses.keys():
            if(losses[key] is not None):    
                print_str = print_str + str(key) + f": {losses[key].mean().item() : 0.07f} " 
                writer.add_scalar(str(key), losses[key].mean().item(), iteration)
        print(print_str)
        if("cuda" in opt['device']):
            GBytes = (torch.cuda.max_memory_allocated(device=opt['device']) \
                / (1024**3))
            if preconditioning is None:
                writer.add_scalar('GPU memory (GB)', GBytes, iteration)
            elif "model" in preconditioning:
                writer.add_scalar('Preconditioning model GPU memory (GB)', GBytes, iteration)
            elif "grid" in preconditioning:
                writer.add_scalar('Preconditioning grid GPU memory (GB)', GBytes, iteration)

def train_step_vanilla(opt, iteration, batch, dataset, model, optimizer, scheduler,
                       early_stopping_data=None):
    opt['iteration_number'] = iteration
    optimizer.zero_grad()
       
    x, y = batch
    x = x.to(opt['device'])
    y = y.to(opt['device'])
    
    model_output = model(x)
    
    data_loss = F.mse_loss(model_output, y, reduction='none').mean()
    
    # Second-derivative regularization in xy-plane (only for 3D inputs)
    lam = float(opt.get('second_deriv_weight', 0.0))
    # reg_loss = torch.tensor(0.0, device=opt['device'])

    # if lam > 0.0 and int(opt.get('n_dims', 3)) == 3:
    #     reg_loss = second_derivative_xy_penalty(
    #         model,
    #         x,
    #         base_output=model_output,
    #         h=float(opt.get('second_deriv_h', 1e-2))
    #     )

    total_loss = data_loss #+ lam * reg_loss
    total_loss.backward()
    optimizer.step()
    scheduler.step()

    if lam > 0.0:
        print(
            f"Iteration {iteration} "
            f"data_loss: {data_loss.item():0.07f} "
            # f"reg_loss: {reg_loss.item():0.07f} "
            f"total_loss: {total_loss.item():0.07f}"
        )
    else:
        print(f"Iteration {iteration} loss: {data_loss.item():0.07f}")


    # loss = F.mse_loss(model_output, y, reduction='none')
    # loss.mean().backward()                   
    # optimizer.step()
    # scheduler.step()   
    # print(f"Iteration {iteration} loss: {loss.mean().item():0.07f}")

def train( model, dataset, opt):
    model = model.to(opt['device'])
    print(model)
    print("Training on %s" % (opt["device"]), os.path.join(opt['save_folder'], opt["dataset_name"]))
    
    dataloader = DataLoader(dataset, 
                            batch_size=None, 
                            num_workers=4 if ("cpu" in opt['data_device'] and "cuda" in opt['device']) else 0,
                            pin_memory=True if ("cpu" in opt['data_device'] and "cuda" in opt['device']) else False,
                            pin_memory_device=opt['device'] if ("cpu" in opt['data_device'] and "cuda" in opt['device']) else "")
    
    model.train(True)

    # choose the specific training iteration function based on the model
    train_step = train_step_vanilla
    optimizer = optim.Adam(model.parameters(), lr=opt["lr"], 
        betas=[opt['beta_1'], opt['beta_2']]) 
    scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer,
        [opt['iterations']*(2/5), opt['iterations']*(3/5), opt['iterations']*(4/5)],
        gamma=0.33)
    early_stopping_data = (False,
        torch.zeros([opt['iterations']], 
        dtype=torch.float32, device=opt['device'])
    )
    
    start_time = time.time()
    for (iteration, batch) in enumerate(dataloader):
        early_stopping_data = train_step(opt,
                iteration,
                batch,
                dataset,
                model,
                optimizer,
                scheduler,
                early_stopping_data=early_stopping_data)
    end_time = time.time()
    sec_passed = end_time-start_time
    mins = sec_passed / 60
    
    print(f"Model completed training after {int(mins)}m {sec_passed%60:0.02f}s")
    save_model(model, opt)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trains an implicit model on data.')
    # part of Hash Grid (NGP model) hyperparameters
    parser.add_argument('--model_type',default='fVSRN', type=str, choices=['NGP','fVSRN'],
        help='Type of model to use - NGP or fVSRN')
    parser.add_argument('--n_dims',default=3, type=int,
        help='Number of dimensions in the data')
    parser.add_argument('--n_outputs',default=1,type=int,
        help='Number of output channels for the data (ex. 1 for scalar field, 3 for vector field)')
    parser.add_argument('--n_grids',default=18,type=int,
        help='Number of grids')
    parser.add_argument('--hash_log2_size',default=15,type=int,
        help='Size of hash table')
    parser.add_argument('--hash_base_resolution',default=16,type=int,
        help='Minimum resolution of a single dimension')
    parser.add_argument('--hash_max_resolution',default=512,type=int,
        help='Maximum resolution of a single dimension') 
    
      
    # part of fVSRN model hyperparameters
    parser.add_argument('--feature_grid_shape',default='64,64,64',type=str,
        help='Feature grid shape for fVSRN model, increase to capture more details, but may not be smooth')
    parser.add_argument('--requires_padded_feats',default=False,type=bool,
        help='Whether to use padded features for fVSRN model, I recommend always set to False')
    parser.add_argument('--num_positional_encoding_terms',default=6,type=int,
        help='Number of positional encoding terms for fVSRN model, lower values for smooth reconstruction, but may miss high frequency details')
    
        # Second-derivative (smoothness) regularization in xy-plane
    parser.add_argument(
        '--second_deriv_weight',
        default=0.0,
        type=float,
        help='Weight for second-derivative smoothness penalty in the xy-plane (0 disables it)'
    )
    parser.add_argument(
        '--second_deriv_h',
        default=1e-2,
        type=float,
        help='Finite-difference step size in normalized coordinate space for second-derivative penalty'
    )
    
    # hyperparameters shared by both fVSRN and NGP models
    parser.add_argument('--n_features',default=2,type=int,
        help='Number of features in the feature grid') 
    parser.add_argument('--n_layers',default=2,type=int,
        help='Number of layers in the model') 
    parser.add_argument('--nodes_per_layer',default=64,type=int,
        help='Nodes per layer in the model')  
    
    # Training and Saving hyperparameters
    parser.add_argument('--data_path',default='./Data/vortex_street.bin',type=str,
        help='Data file path')
    parser.add_argument('--save_folder',default='./SavedModels',type=str,
        help='Data file path')
    parser.add_argument('--dataset_name',default='vortex_street',type=str,
        help='Save name for the model')  
    parser.add_argument('--data_dims',default='100,80,50',type=str,help='Data dimensions')
    
        
    
    parser.add_argument('--device',default='cuda:0',type=str,
        help='Which device to train on')
    parser.add_argument('--data_device',default='cuda:0',type=str,
        help='Which device to keep the data on')

    parser.add_argument('--iterations',default=50000, type=int,
        help='Number of iterations to train')
    parser.add_argument('--points_per_iteration',default=100000, type=int,
        help='Number of points to sample per training loop update')
    parser.add_argument('--lr',default=0.01, type=float,
        help='Learning rate for the adam optimizer')
    parser.add_argument('--beta_1',default=0.9, type=float,
        help='Beta1 for the adam optimizer')
    parser.add_argument('--beta_2',default=0.99, type=float,
        help='Beta2 for the adam optimizer')

    args = vars(parser.parse_args())
    
    os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"
    torch.manual_seed(42)
    torch.backends.cuda.matmul.allow_tf32 = True
    
    dataset = Dataset(args)
    args['data_min'] = dataset.min().item()
    args['data_max'] = dataset.max().item()
    
    #opt['data_min'] = max(dataset.min(), dataset.data.mean() - dataset.data.std()*3).item()
    #opt['data_max'] = min(dataset.max(), dataset.data.mean() + dataset.data.std()*3).item()
    #opt['data_min'] = dataset.data.mean().item()
    #opt['data_max'] = max(dataset.data.mean()-dataset.data.min(), dataset.data.max() -dataset.data.mean()).item()
    if args['model_type'] == 'NGP':
        model = NGP_TCNN(args)
    elif args['model_type'] == 'fVSRN':
        model = fVSRN(args)
    model = model.to(args['device'])
    

    now = datetime.datetime.now()
    start_time = time.time()
    
    train(model, dataset, args)
