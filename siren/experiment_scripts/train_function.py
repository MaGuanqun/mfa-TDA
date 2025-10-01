# Enable import from parent package
import sys
import os
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )

import dataio, meta_modules, utils, training, loss_functions, modules

from torch.utils.data import DataLoader
import configargparse
from functools import partial
import torch

import function

p = configargparse.ArgumentParser()
p.add('-c', '--config_filepath', required=False, is_config_file=False, help='Path to config file.')

p.add_argument('--logging_root', type=str, default='./logs', help='root for logging')
p.add_argument('--experiment_name', type=str, required=True,
               help='Name of subdirectory in logging_root where summaries and checkpoints will be saved.')

# General training options
p.add_argument('--batch_size', type=int, default=1)
p.add_argument('--lr', type=float, default=1e-4, help='learning rate. default=1e-4')
p.add_argument('--num_epochs', type=int, default=5000,
               help='Number of epochs to train for.')

p.add_argument('--epochs_til_ckpt', type=int, default=2000,
               help='Time interval in seconds until checkpoint is saved.')
p.add_argument('--steps_til_summary', type=int, default=4999,
               help='Time interval in seconds until tensorboard summary is saved.')

p.add_argument('--model_type', type=str, default='sine',
               help='Options currently are "sine" (all sine activations), "relu" (all relu activations,'
                    '"nerf" (relu activations and positional encoding as in NeRF), "rbf" (input rbf layer, rest relu),'
                    'and in the future: "mixed" (first layer sine, other layers tanh)')
p.add_argument('--checkpoint_path', default=None, help='Checkpoint to trained model.')
p.add_argument('--function_name', default="expotential", help='Set different functions.')
p.add_argument('--hidden_features', type=int, default=16, help='Number of hidden features.')
p.add_argument('--num_hidden_layers', type=int, default=12, help='Number of hidden layers.')
p.add_argument('--mode', type=str, default='mlp', help='mode of modeling(mlp, nerf, rbf).')

opt = p.parse_args()


opt.steps_til_summary = opt.num_epochs-1



# ------------------------
# Dataset producing z(x,y)
# ------------------------
class XYDataset(torch.utils.data.Dataset):
    """
    Provides a scalar field z over a *normalized* [-1,1]^2 grid (SIREN-friendly).
    The underlying function is evaluated on a *physical domain* specified by `domain`.
    """
    def __init__(self, sidelength, domain,function_name):
        """
        sidelength: (H, W)
        domain: [x_min, x_max, y_min, y_max] for the *physical* coords where z = f(x,y) is evaluated.
        """
        super().__init__()
        if isinstance(sidelength, int):
            sidelength = (sidelength, sidelength)
        self.sidelength = sidelength
        self.domain = domain  # [xmin, xmax, ymin, ymax]

        # Canonical grid in [-1,1]^2; shape [H*W, 2]
        self.mgrid = dataio.get_mgrid(self.sidelength)

        # Precompute z once (dataset length is 1)
        x_phys, y_phys = self._denorm_to_physical(self.mgrid)  # map to requested range
        with torch.no_grad():
            z = function.function_2D(x_phys, y_phys,function_name).float().unsqueeze(1)  # [N,1]

        # Min–max to [-1,1] for SIREN stability
        z_min, z_max = z.min(), z.max()
        # Guard against degenerate ranges
        if (z_max - z_min) > 0:
            z = 2.0 * (z - z_min) / (z_max - z_min) - 1.0
        else:
            z = torch.zeros_like(z)

        self.func = z  # [N,1]

    def _denorm_to_physical(self, uv):
        """
        Map uv in [-1,1] to x,y in [xmin,xmax]x[ymin,ymax].
        uv: [N,2]
        """
        xmin, xmax, ymin, ymax = self.domain
        u = uv[:, 0]
        v = uv[:, 1]
        x = 0.5 * (u + 1.0) * (xmax - xmin) + xmin
        y = 0.5 * (v + 1.0) * (ymax - ymin) + ymin
        return x, y

    def __len__(self):
        return 1

    def __getitem__(self, index):
        # DO NOT mutate self.mgrid here.
        return {'coords': self.mgrid}, {'func': self.func}

sidelength=(401,401)
range=function.range(opt.function_name)
    
dataset = XYDataset(sidelength, range,opt.function_name)


coord_dataset = dataio.Implicit2DFuncWrapper(dataset, sidelength=sidelength)

# print(len(coord_dataset))
# print(coord_dataset[0])

dataloader = DataLoader(coord_dataset, shuffle=True, batch_size=opt.batch_size, pin_memory=True, num_workers=0)

# for (model_input, gt) in dataloader:
#     print(gt)


# Define the model.
if opt.model_type in ['sine', 'relu', 'tanh', 'selu', 'elu', 'softplus'] and opt.mode in ['mlp','nerf','rbf']:
    model = modules.SingleBVPNet(type=opt.model_type, mode=opt.mode, sidelength=sidelength,hidden_features = opt.hidden_features, num_hidden_layers=opt.num_hidden_layers)
# elif opt.model_type in ['rbf' , 'nerf']:
#     model = modules.SingleBVPNet(type='relu', mode=opt.model_type, sidelength=sidelength)
# elif opt.model_type in ['rbf', 'nerf']:
    
#     model = modules.SingleBVPNet(type=opt.model_type if opt.model_type not in ['rbf','nerf'] else 'tanh',
#                                  mode=opt.model_type, sidelength=sidelength)

else:
    raise NotImplementedError
model.cuda()

root_path = os.path.join(opt.logging_root, opt.experiment_name)

# Define the loss

loss_fn = partial(loss_functions.function_mse)

func_info = [sidelength[0],sidelength[1],range[0],range[1],range[2],range[3]]


summary_fn = partial(utils.write_function_summary, func_info)

training.train(
    model=model,
    train_dataloader=dataloader,
    epochs=opt.num_epochs,
    lr=opt.lr,
    steps_til_summary=opt.steps_til_summary,
    epochs_til_checkpoint=opt.epochs_til_ckpt,
    model_dir=root_path,
    loss_fn=loss_fn,
    summary_fn=summary_fn)
