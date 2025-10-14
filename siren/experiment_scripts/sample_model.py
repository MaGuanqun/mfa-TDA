##sample discrete data from a trained siren model

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import modules

import torch
import numpy as np
import trimesh
import function
import dataio
import configargparse


# Map uv -> physical xy (for true function eval & for exported vertices)
def uv_to_xy(uv_,domain):
    u = uv_[:, 0]
    v = uv_[:, 1]
    if len(domain)==4:
        xmin, xmax, ymin, ymax = domain
    else:
        xmin, xmax, ymin, ymax, zmin, zmax = domain
        
    x = 0.5 * (u + 1.0) * (xmax - xmin) + xmin
    y = 0.5 * (v + 1.0) * (ymax - ymin) + ymin
    
    if len(domain) ==6:
        zc = uv_[:, 2]
        z = 0.5 * (zc + 1.0) * (zmax - zmin) + zmin
        return torch.stack([x, y, z], dim=-1)
    return torch.stack([x, y], dim=-1)


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


# --------------------
# Derivatives & Hessian for SIREN (with respect to UV), batched
# --------------------
def compute_function_batched(model, uv, batch_size=10000):
    """
    Args:
        model: modules.SingleBVPNet
        uv: [N, 2] tensor in [-1, 1]
        batch_size: int
    Returns:
        f_all:     [N, 1]
        grad_all:  [N, 2]
        H_all:     [N, 2, 2]
        h_all:     [N, 1]  (computed from grad & Hessian)
    """
    f_list = []

    N = uv.shape[0]
    for i in range(0, N, batch_size):
        uv_chunk = uv[i:i + batch_size].unsqueeze(0).requires_grad_(True)  # [1, B, 2]

        # Forward
        out_chunk = model({'coords': uv_chunk})
        model_out_chunk = out_chunk['model_out']   # [1, B, 1]     
        # Accumulate
        f_list.append(model_out_chunk[0].detach())   # [B, 1]
    f_all = torch.cat(f_list, dim=0)
    return f_all


# --------------------
# Export helpers
# --------------------
def export_mesh(xy_np, faces, z_tensor, path):
    verts = np.concatenate([xy_np, z_tensor.detach().cpu().numpy()], axis=1)  # [N, 3]
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    mesh.export(path)

def write_volume_zyx(coords_N3: torch.Tensor,
                      values_N1: torch.Tensor,
                      sidelength, 
                      out_path: str,
                      dtype='<f8'):
    """
    coords_N3: [N,3] the SAME coords you fed to the model (columns: x,y,z)
    values_N1: [N,1] or [N] values evaluated at those coords
    sidelength: (nx, ny, nz)
    Writes binary matching converter's reshape((nz,ny,nx)) in C-order.
    """
    nx, ny, nz = map(int, sidelength)
    N = nx * ny * nz

    # Make flat & to CPU/NumPy
    coords = coords_N3.detach().cpu().numpy()
    vals   = values_N1.reshape(-1).detach().cpu().numpy()

    if coords.shape[0] != N or vals.size != N:
        raise RuntimeError(f"Count mismatch: coords={coords.shape[0]}, vals={vals.size}, expected={N}")

    # Sort by (z, y, x) so x varies fastest inside each-y, inside each-z
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]
    order = np.lexsort((x, y, z))  # z primary, then y, then x

    vals_sorted = vals[order].reshape(nz, ny, nx)      # (z, y, x)
    vals_sorted.astype(dtype).ravel(order='C').tofile(out_path)
    

p = configargparse.ArgumentParser()
p.add_argument('--hidden_features', type=int, default=16, help='Number of hidden features.')
p.add_argument('--num_hidden_layers', type=int, default=12, help='Number of hidden layers.')
p.add_argument('--model_type', type=str, default='sine', help='Activation type used during training.')
p.add_argument('--checkpoint', type=str, default='logs/expotential/checkpoints/model_final.pth')
p.add_argument('--out_root', type=str, default='logs/expotential')
p.add_argument('--function_name', default="expotential", help='Set different functions.')
p.add_argument('--mode', type=str, default='mlp', help='mode of modeling(mlp, nerf, rbf).')
p.add_argument('--omega', type=float, default=20, help='omega for sine')
p.add_argument('--batch_size', type=int, default=20000, help='Batch size for UV derivs.')
opt = p.parse_args()


# --------------------
# Paths
# --------------------
model_path = opt.checkpoint
out_root = opt.out_root

os.makedirs(out_root, exist_ok=True)

out_f_siren = os.path.join(out_root, 'f.ply')
out_h_siren = os.path.join(out_root, 'h.ply')
out_f_true  = os.path.join(out_root, 'ori.ply')  # normalized like training

out_f_3d_bin   = os.path.join(out_root, 'f.bin')
out_f_3d_true_bin = os.path.join(out_root, 'f_true.bin')
# --------------------
# Grid setup
# --------------------
sidelength = function.sample_size(opt.function_name)
# Physical domain used when generating the dataset: [xmin, xmax, ymin, ymax]
domain = function.range(opt.function_name)

# UV grid in [-1,1]^2 (model input)
uv = dataio.get_mgrid(sidelength,dim=len(sidelength))                   # [N, 2] in [-1,1]

xy = uv_to_xy(uv,domain)           

xy_np = uv.detach().cpu().numpy()

if uv.shape[1]==2:
    faces = triangulate_grid(sidelength)
# --------------------
# Evaluate ORIGINAL function (and normalize like training)

# --------------------
    x_phys = xy[:, 0]
    y_phys = xy[:, 1]
    f_true = function.function_2D(x_phys, y_phys, opt.function_name).float().unsqueeze(1)
else:
    x_phys = xy[:, 0]
    y_phys = xy[:, 1]
    z_phys = xy[:, 2]
    f_true = function.function_3D(x_phys, y_phys, z_phys, opt.function_name).float().unsqueeze(1)
    print(x_phys)
    print(z_phys)
    
zmin, zmax = f_true.min(), f_true.max()
print(f"True function {opt.function_name} in [{zmin.item():.3f}, {zmax.item():.3f}]")

if (zmax - zmin) > 0:
    f_true_n = 2.0 * (f_true - zmin) / (zmax - zmin) - 1.0
else:
    f_true_n = torch.zeros_like(f_true)

print(f"Normalized to [{f_true_n.min().item():.3f}, {f_true_n.max().item():.3f}]")


# # --------------------
# # Load SIREN and evaluate in batches (keep graph for autograd)
# # --------------------
# model = modules.SingleBVPNet(
#     type=opt.model_type,
#     mode=opt.mode,
#     sidelength=sidelength,
#     hidden_features=opt.hidden_features,
#     num_hidden_layers=opt.num_hidden_layers,
#     omega=opt.omega,
#     in_features=len(sidelength)
# )

# state = torch.load(model_path)  # add map_location='cpu' if needed
# model.load_state_dict(state)
# model.eval()

# # Batched derivatives & H functional
# f_siren = compute_function_batched(
#     model, uv, batch_size=opt.batch_size
# )

f_siren=None

# --------------------
# Write meshes
# --------------------
if uv.shape[1]==2:
    export_mesh(xy_np, faces, f_true_n,   out_f_true)   # original f (normalized)
    export_mesh(xy_np, faces, f_siren,    out_f_siren)  # siren f
    print("Wrote:")
    print(f"  True f (normalized): {out_f_true}")
    print(f"  SIREN f:             {out_f_siren}")
    print(f"  SIREN h:             {out_h_siren}")

else:
    # 3D write block — make the memory layout match the converter's expectation
    N = int(np.prod(sidelength))
    # if f_siren.numel() != N:
    #     raise RuntimeError(f"Value count mismatch: got {f_siren.numel()}, expected {N} from sidelength {sidelength}")

    nx, ny, nz = map(int, sidelength)  # expected (x, y, z)

    # # SIREN output: (N,1) -> (nx, ny, nz) then ravel in Fortran order (x fastest)
    # arr_siren = f_siren.detach().cpu().numpy().reshape(nx, ny, nz, order='C')
    # vals_siren = np.asarray(arr_siren, dtype='<f8').ravel(order='F')
    # with open(out_f_3d_bin, 'wb') as fh:
    #     vals_siren.tofile(fh)

    # True function: same treatment
    arr_true = f_true.detach().cpu().numpy().reshape(nx, ny, nz, order='C')
    vals_true = np.asarray(arr_true, dtype='<f8').ravel(order='F')
    with open(out_f_3d_true_bin, 'wb') as ft:
        vals_true.tofile(ft)

    print("Wrote:")
    print(f"  SIREN f:             {out_f_3d_bin}")
    print(f"  True f:              {out_f_3d_true_bin}")
