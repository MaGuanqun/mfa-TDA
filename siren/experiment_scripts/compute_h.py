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
def uv_to_xy(uv_):
    u = uv_[:, 0]
    v = uv_[:, 1]
    xmin, xmax, ymin, ymax = domain
    x = 0.5 * (u + 1.0) * (xmax - xmin) + xmin
    y = 0.5 * (v + 1.0) * (ymax - ymin) + ymin
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
# H functional for a scalar field f (2D): uses fx, fy, fxx, fyy, fxy
# --------------------
def compute_h(fx, fy, fxx, fyy, fxy):
    # H(f) = 2 * ((fx^2 - fy^2) * fxy + fx*fy*(fyy - fxx))
    return 2.0 * ((fx**2 - fy**2) * fxy + fx * fy * (fyy - fxx))


# --------------------
# Derivatives & Hessian for SIREN (with respect to UV), batched
# --------------------
def compute_first_and_hessian_batched(model, uv, batch_size=10000):
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
    grads_list = []
    hess_list = []
    h_list = []

    N = uv.shape[0]
    for i in range(0, N, batch_size):
        uv_chunk = uv[i:i + batch_size].unsqueeze(0).requires_grad_(True)  # [1, B, 2]

        # Forward
        out_chunk = model({'coords': uv_chunk})
        model_out_chunk = out_chunk['model_out']   # [1, B, 1]
        model_in_chunk  = out_chunk['model_in']    # [1, B, 2]

        # First derivatives: grad wrt inputs
        grad1 = torch.autograd.grad(
            outputs=model_out_chunk,
            inputs=model_in_chunk,
            grad_outputs=torch.ones_like(model_out_chunk),
            create_graph=True,
            allow_unused=False
        )[0]  # [1, B, 2]

        # Hessian via second derivatives
        B = model_in_chunk.shape[1]
        H_chunk = torch.zeros((B, 2, 2),
                              dtype=model_in_chunk.dtype,
                              device=model_in_chunk.device)
        for j in range(1):
            g = torch.autograd.grad(
                outputs=grad1[:, :, j],
                inputs=model_in_chunk,
                grad_outputs=torch.ones_like(grad1[:, :, j]),
                create_graph=True,
                allow_unused=False
            )[0]  # [1, B, 2]
            H_chunk[:, :, j] = g[0]

        j=1
        g = torch.autograd.grad(
                outputs=grad1[:, :, j],
                inputs=model_in_chunk,
                grad_outputs=torch.ones_like(grad1[:, :, j]),
                create_graph=False,
                allow_unused=False
            )[0]  # [1, B, 2]
        H_chunk[:, :, j] = g[0]


        # Compute H functional per chunk to keep memory low
        fx, fy = grad1[0][:, 0], grad1[0][:, 1]
        fxx, fyy, fxy = H_chunk[:, 0, 0], H_chunk[:, 1, 1], H_chunk[:, 0, 1]
        h_chunk = compute_h(fx, fy, fxx, fyy, fxy).unsqueeze(1)  # [B, 1]

        # Accumulate
        f_list.append(model_out_chunk[0].detach())   # [B, 1]
        grads_list.append(grad1[0].detach())         # [B, 2]
        # hess_list.append(H_chunk.detach())           # [B, 2, 2]
        h_list.append(h_chunk.detach())              # [B, 1]

        # (Optional) free graph references in loop iteration
        # Nothing explicit needed; Python refcounts will drop each iter.

    f_all = torch.cat(f_list, dim=0)
    grad_all = torch.cat(grads_list, dim=0)
    # H_all = torch.cat(hess_list, dim=0)
    h_all = torch.cat(h_list, dim=0)
    return f_all, grad_all, h_all


# --------------------
# Export helpers
# --------------------
def export_mesh(xy_np, faces, z_tensor, path):
    verts = np.concatenate([xy_np, z_tensor.detach().cpu().numpy()], axis=1)  # [N, 3]
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
    mesh.export(path)


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
uv = dataio.get_mgrid(sidelength)                   # [N, 2] in [-1,1]

xy = uv_to_xy(uv)                    # [N,2] physical coords
xy_np = uv.detach().cpu().numpy()
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
# Load SIREN and evaluate in batches (keep graph for autograd)
# --------------------
model = modules.SingleBVPNet(
    type=opt.model_type,
    mode=opt.mode,
    sidelength=sidelength,
    hidden_features=opt.hidden_features,
    num_hidden_layers=opt.num_hidden_layers,
    omega=opt.omega
)

state = torch.load(model_path)  # add map_location='cpu' if needed
model.load_state_dict(state)
model.eval()

# Batched derivatives & H functional
f_siren, grad_siren, h_siren = compute_first_and_hessian_batched(
    model, uv, batch_size=opt.batch_size
)

# --------------------
# Write meshes
# --------------------
export_mesh(xy_np, faces, f_true_n,   out_f_true)   # original f (normalized)
export_mesh(xy_np, faces, f_siren,    out_f_siren)  # siren f
export_mesh(xy_np, faces, h_siren,    out_h_siren)  # siren h

print("Wrote:")
print(f"  True f (normalized): {out_f_true}")
print(f"  SIREN f:             {out_f_siren}")
print(f"  SIREN h:             {out_h_siren}")
