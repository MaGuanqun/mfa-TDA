import sys
import os
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )

import modules

import torch
import numpy as np
import trimesh
import function
import dataio

import configargparse


p = configargparse.ArgumentParser()
p.add_argument('--hidden_features', type=int, default=16, help='Number of hidden features.')
p.add_argument('--num_hidden_layers', type=int, default=12, help='Number of hidden layers.')

p.add_argument('--model_type', type=str, default='sine', help='Activation type used during training.')
p.add_argument('--checkpoint', type=str, default='logs/expotential/checkpoints/model_final.pth')
p.add_argument('--out_root', type=str, default='logs/expotential')
p.add_argument('--function_name', default="expotential", help='Set different functions.')
p.add_argument('--mode', type=str, default='mlp', help='mode of modeling(mlp, nerf, rbf).')

opt = p.parse_args()

# --------------------
# Paths
# --------------------
model_path = opt.checkpoint
out_root = opt.out_root
os.makedirs(out_root, exist_ok=True)

out_f_siren = os.path.join(out_root, 'expotential_f.ply')
out_h_siren = os.path.join(out_root, 'expotential_h.ply')
out_f_true  = os.path.join(out_root, 'expotential_ori.ply')  # normalized like training

# --------------------
# Grid setup
# --------------------
sidelength = (401, 401)
# Physical domain used when generating the dataset: [xmin, xmax, ymin, ymax]
domain = function.range(opt.function_name)  

# UV grid in [-1,1]^2 (model input)
uv = dataio.get_mgrid(sidelength)                   # [N,2] in [-1,1]
uv_batched = uv.unsqueeze(0).requires_grad_(True)   # [1,N,2]

# Map uv -> physical xy (for true function eval & for exported vertices)
def uv_to_xy(uv_):
    u = uv_[:, 0]
    v = uv_[:, 1]
    xmin, xmax, ymin, ymax = domain
    x = 0.5 * (u + 1.0) * (xmax - xmin) + xmin
    y = 0.5 * (v + 1.0) * (ymax - ymin) + ymin
    return torch.stack([x, y], dim=-1)

xy = uv_to_xy(uv)                    # [N,2] physical coords
xy_np = uv.detach().numpy()

# --------------------
# Triangulate the image grid
# --------------------
def triangulate_grid(slen):
    H, W = slen
    faces = []
    for i in range(H - 1):
        for j in range(W - 1):
            p0 = i * W + j
            p1 = p0 + 1
            p2 = (i + 1) * W + j
            p3 = p2 + 1
            faces.append([p0, p1, p3])
            faces.append([p0, p3, p2])
    return np.array(faces, dtype=np.int64)

faces = triangulate_grid(sidelength)

# --------------------
# Evaluate ORIGINAL function (and normalize like training)
# --------------------
x_phys = xy[:, 0]
y_phys = xy[:, 1]
f_true = function.function_2D(x_phys, y_phys, opt.function_name).float().unsqueeze(1)

zmin, zmax = f_true.min(), f_true.max()

print(f"True function {opt.function_name} in [{zmin.item():.3f}, {zmax.item():.3f}]")

if (zmax - zmin) > 0:
    f_true_n = 2.0 * (f_true - zmin) / (zmax - zmin) - 1.0
else:
    f_true_n = torch.zeros_like(f_true)

print(f"Normalized to [{f_true_n.min().item():.3f}, {f_true_n.max().item():.3f}]")
# --------------------
# Load SIREN and forward on UV (no torch.no_grad so we can take derivatives)
# --------------------
model = modules.SingleBVPNet(
    type=opt.model_type,
    mode=opt.mode,
    sidelength=sidelength,
    hidden_features=opt.hidden_features,
    num_hidden_layers=opt.num_hidden_layers
)

state = torch.load(model_path)#, map_location='cpu'
model.load_state_dict(state)
model.eval()

# forward pass (track grad for derivatives)
out = model({'coords': uv_batched})
model_out = out['model_out']     # [1,N,1]
model_in  = out['model_in']      # [1,N,2] (== uv_batched, requires_grad=True)

# --------------------
# Derivatives & H for SIREN (with respect to UV)
# --------------------
def compute_first_and_hessian(model_out_, model_in_):
    grad1 = torch.autograd.grad(
        outputs=model_out_,
        inputs=model_in_,
        grad_outputs=torch.ones_like(model_out_),
        create_graph=True,
        allow_unused=False
    )[0]  # [1,N,2]
    H = torch.zeros((model_in_.shape[1], 2, 2), dtype=model_in_.dtype)
    for i in range(2):
        g = torch.autograd.grad(
            outputs=grad1[:, :, i],
            inputs=model_in_,
            grad_outputs=torch.ones_like(grad1[:, :, i]),
            create_graph=True,
            allow_unused=False
        )[0]  # [1,N,2]
        H[:, :, i] = g[0]
    return grad1[0], H  # grad: [N,2], H: [N,2,2]

def compute_h(fx, fy, fxx, fyy, fxy):
    # H(f) = 2 * ((fx^2 - fy^2) * fxy + fx*fy*(fyy - fxx))
    return 2.0 * ((fx**2 - fy**2) * fxy + fx * fy * (fyy - fxx))

grad_siren, H_siren = compute_first_and_hessian(model_out, model_in)
h_siren = compute_h(
    grad_siren[:, 0], grad_siren[:, 1],
    H_siren[:, 0, 0], H_siren[:, 1, 1], H_siren[:, 0, 1]
).unsqueeze(1)  # [N,1]

# --------------------
# Export helpers
# --------------------
def export_mesh(z_tensor, path):
    verts = np.concatenate([xy_np, z_tensor.detach().numpy()], axis=1)  # [N,3]
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    mesh.export(path)

# Write meshes
export_mesh(f_true_n,   out_f_true)          # original f (normalized)
export_mesh(model_out[0], out_f_siren)       # siren f
export_mesh(h_siren,    out_h_siren)         # siren h

print("Wrote:")
print(f"  True f (normalized): {out_f_true}")
print(f"  SIREN f:             {out_f_siren}")
print(f"  SIREN h:             {out_h_siren}")