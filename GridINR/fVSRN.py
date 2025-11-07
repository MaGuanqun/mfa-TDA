import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import List
import math

class ReLULayer(nn.Module):
    def __init__(self, in_features, out_features, 
                bias=True, dtype=torch.float32):
        super().__init__()
        
        self.in_features = in_features
        self.linear = nn.Linear(in_features, out_features, 
            bias=bias, dtype=dtype)
        
        self.init_weights()
    
    def init_weights(self):
        with torch.no_grad():
            nn.init.xavier_normal_(self.linear.weight)
            if(self.linear.bias is not None):
                nn.init.zeros_(self.linear.bias)

    def forward(self, input):
        return F.relu(self.linear(input))

class PositionalEncoding(nn.Module):
    def __init__(self, num_terms:int, n_dims:int):
        super(PositionalEncoding, self).__init__()  
        self.n_dims = n_dims      
        
        self.L = num_terms
        
        # Frequencies: [L], each term = 2^i * π
        freqs = torch.pow(2.0, torch.arange(0, num_terms, dtype=torch.float32)) * math.pi
        self.register_buffer("freqs", freqs, persistent=False)
        
        # L_terms = torch.arange(0, num_terms, 
        #     dtype=torch.float32).repeat_interleave(2*n_dims)
        # L_terms = torch.pow(2, L_terms) * torch.pi
        # self.L_t = L_terms
        # self.register_buffer("L_terms", L_terms, persistent=False)

    def forward(self, locations):
        loc = locations.unsqueeze(-1)                               # [..., n_dims, 1]
        # broadcast freqs: [1, 1, L] to match [..., n_dims, L]
        freqs = self.freqs.view(*([1] * (loc.dim() - 1)), -1)       # [..., 1, L]

        # phase: [..., n_dims, L]
        phase = loc * freqs                                         # out-of-place

        sin_part = torch.sin(phase)                                 # [..., n_dims, L]
        cos_part = torch.cos(phase)                                 # [..., n_dims, L]

        # Stack sin/cos along "channel" dim -> [..., 2*n_dims, L]
        enc = torch.cat([sin_part, cos_part], dim=-2)

        # Flatten last two dims -> [..., L * 2 * n_dims]
        enc = enc.reshape(*locations.shape[:-1], -1)
        return enc
    
        # repeats = len(locations.shape) * [1]
        # repeats[-1] = self.L*2
        # locations = locations.repeat(repeats)
        
        # locations = locations * self.L_terms# + self.phase_shift
        # if(self.n_dims == 2):
        #     locations[..., 0::4] = torch.sin(locations[..., 0::4])
        #     locations[..., 1::4] = torch.sin(locations[..., 1::4])
        #     locations[..., 2::4] = torch.cos(locations[..., 2::4])
        #     locations[..., 3::4] = torch.cos(locations[..., 3::4])
        # else:
        #     locations[..., 0::6] = torch.sin(locations[..., 0::6])
        #     locations[..., 1::6] = torch.sin(locations[..., 1::6])
        #     locations[..., 2::6] = torch.sin(locations[..., 2::6])
        #     locations[..., 3::6] = torch.cos(locations[..., 3::6])
        #     locations[..., 4::6] = torch.cos(locations[..., 4::6])
        #     locations[..., 5::6] = torch.cos(locations[..., 5::6])
        # return locations


class fVSRN(nn.Module):
    def __init__(self, opt):
        super().__init__()
        
        # print(opt['requires_padded_feats'])
        # exit()
        self.requires_padded_feats : bool = opt['requires_padded_feats']
        self.padding_size : int = 0
        self.full_shape = opt['n_dims']
        feature_grid_shape = [eval(i) for i in opt['feature_grid_shape'].split(",")]

        if(opt['requires_padded_feats']):
            self.padding_size : int = \
                16*int(math.ceil(max(1, (opt['n_features'] +opt['num_positional_encoding_terms']*opt['n_dims']*2)/16))) - \
                    (opt['n_features'] +opt['num_positional_encoding_terms']*opt['n_dims']*2)
            
        self.pe = PositionalEncoding(opt['num_positional_encoding_terms'], opt['n_dims'])        
        feat_shape : List[int] = [1, opt['n_features'] ] + feature_grid_shape
        self.full_shape = opt['n_dims']
        self.feature_grid = torch.rand(feat_shape, 
            dtype=torch.float32)
        self.feature_grid = torch.nn.Parameter(self.feature_grid, 
            requires_grad=True)

        
        def init_decoder_pytorch():
            decoder = nn.ModuleList()
            
            input_size:int = opt['n_features'] +opt['num_positional_encoding_terms']*opt['n_dims']*2
            if(opt['requires_padded_feats']):
                input_size = opt['n_features'] +opt['num_positional_encoding_terms']*opt['n_dims']*2 + self.padding_size
                                    
            layer = ReLULayer(input_size, 
                opt['nodes_per_layer'], bias=False)
            decoder.append(layer)
            
            for i in range(opt['n_layers'] ):
                if i == opt['n_layers']  - 1:
                    layer = nn.Linear(opt['nodes_per_layer'], opt['n_outputs'], bias=False)
                    decoder.append(layer)
                else:
                    layer = ReLULayer(opt['nodes_per_layer'], opt['nodes_per_layer'], bias=False)
                    decoder.append(layer)
            decoder = torch.nn.Sequential(*decoder)
            return decoder


        self.decoder = init_decoder_pytorch()

        self.register_buffer(
            "volume_min",
            torch.tensor([opt['data_min']], requires_grad=False, dtype=torch.float32),
            persistent=False
        )
        self.register_buffer(
            "volume_max",
            torch.tensor([opt['data_max']], requires_grad=False, dtype=torch.float32),
            persistent=False
        )
    
    def min(self):
        return self.volume_min

    def max(self):
        return self.volume_max
    
    def get_volume_feats(self): #* not used
        if self.feats is None:
            raise ValueError("No intermediate features available yet.")
        return self.feats
    
    def get_volume_extents(self):
        return self.full_shape
                       
    def forward(self, x, scale_output=False):     
        
        feats = F.grid_sample(self.feature_grid,
                x.reshape(([1]*x.shape[-1]) + list(x.shape)),
                mode='bilinear', align_corners=True) 
        pe = self.pe(x)  
        
        feats = feats.flatten(0, -2).permute(1, 0)
        
        feats = torch.cat([pe, feats], dim=1)
        
        if(self.requires_padded_feats):
            feats = F.pad(feats, (0, self.padding_size), value=1.0) 
        
        y = self.decoder(feats).float()
        if scale_output:
            y = y * (self.volume_max - self.volume_min) + self.volume_min
        return y

        