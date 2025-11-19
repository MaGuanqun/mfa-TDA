import os
import torch
from utils import make_coord_grid, dat_to_tensor, bin_to_tensor
import torch.nn.functional as F
import time
import sys

class Dataset(torch.utils.data.Dataset):
    def __init__(self, opt):
        
        self.opt = opt
        self.min_ = None
        self.max_ = None
        self.mean_ = None
        self.full_coord_grid = None

        t1 = time.time()
        if opt['dataset_name'] in ['vortex_street', 'quartic_potential_2','boussinesq_3d','fluid','vortex_street_3d','cylinder','cylinder2','cylinder3']:
            d, full_shape, d_min, d_max = bin_to_tensor(opt['data_path'], opt)
        if opt['dataset_name'] in ['vortex']:
            d, full_shape, d_min, d_max = dat_to_tensor(opt['data_path'], opt)
        self.full_shape = full_shape
        self.opt['points_per_iteration'] = min(
            opt['points_per_iteration'],
            self.full_shape[0]*self.full_shape[1]*self.full_shape[2]
        )
        self.min_ = d_min
        self.max_ = d_max
        
        d = d.to(self.opt['data_device'])
        t2 = time.time()
        print(f"Data: {d.shape} from full extents {full_shape}. IO time loading data: {t2-t1 : 0.04f}")
        self.data = d
    
    def min(self):
        if self.min_ is not None:
            return self.min_
        else:
            self.min_ = self.data.min()
            return self.min_

    def max(self):
        if self.max_ is not None:
            return self.max_
        else:
            self.max_ = self.data.max()
            return self.max_

    def get_2D_slice(self):
        if(len(self.data.shape) == 4):
            return self.data[0].clone()
        else:
            return self.data[0,:,:,:,int(self.data.shape[4]/2)].clone()

    def sample_rect(self, starts, widths, samples):
        positions = []
        for i in range(len(starts)):
            positions.append(
                torch.arange(starts[i], starts[i] + widths[i], widths[i] / samples[i], 
                    dtype=torch.float32, device=self.opt['data_device'])
            )
            positions[i] -= 0.5
            positions[i] *= 2
        grid_to_sample = torch.stack(torch.meshgrid(*positions), dim=-1).unsqueeze(0)

        vals = F.grid_sample(self.data, 
                grid_to_sample, mode='bilinear', 
                align_corners=True)
        #print('dataset sample rect vals shape')
        print(vals.shape)
        return vals

    
    
    def total_points(self):
        t = 1
        for i in range(2, len(self.data.shape)):
            t *= self.data.shape[i]
        return t

    def get_full_coord_grid(self):
        if self.full_coord_grid is None:
            self.full_coord_grid = make_coord_grid(self.data.shape[2:], 
                    self.opt['data_device'], flatten=True, 
                    align_corners=True)
        return self.full_coord_grid
        
    def get_random_points(self, n_points):        
        
        x = torch.rand([1, 1, 1, n_points, self.opt['n_dims']], 
                device=self.opt['data_device']) * 2 - 1

        y = F.grid_sample(self.data,
            x, mode='bilinear', 
            align_corners=True)
        
        x = x.squeeze()
        y = y.squeeze()
        if(len(y.shape) == 1):
            y = y.unsqueeze(0)    
        
        y = y.permute(1,0)
        return x, y

    def __len__(self):
        return self.opt['iterations']
    
    def __getitem__(self, idx):
        return self.get_random_points(
            self.opt['points_per_iteration']
        )