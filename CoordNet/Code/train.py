import torch
from torch import nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
import os
import numpy as np
import torch.optim as optim
import time
from model import *
from utils import *
from skimage.io import imsave
from skimage.io import imread
from skimage import data,img_as_float,img_as_int
import lpips

def trainNet(model,args,dataset):
        # ----- logging file naming follows original conventions -----
    if args.application in ['spatial', 'super-spatial']:
        loss_path = args.model_path + args.dataset + '/' + \
            f'loss-{args.application}-{args.scale}-{args.init}-{args.factor}.txt'
    elif args.application == 'temporal':
        loss_path = args.model_path + args.dataset + '/' + \
            f'loss-{args.application}-{args.interval}-{args.init}-{args.factor}-{args.active}.txt'
    else:
        # includes 'super-spatial-temporal' and any other fallbacks
        loss_path = args.model_path + args.dataset + '/' + \
            f'loss-{args.application}-{args.init}-{args.num_res}.txt'

    loss = open(loss_path, 'w')
    
    # if args.application == 'spatial' or args.application == 'super-spatial':
    #     loss = open(args.model_path+args.dataset+'/'+'loss-'+args.application+'-'+str(args.scale)+'-'+str(args.init)+'-'+str(args.factor)+'.txt','w')
    # elif args.application == 'temporal':
    #     loss = open(args.model_path+args.dataset+'/'+'loss-'+args.application+'-'+str(args.interval)+'-'+str(args.init)+'-'+str(args.factor)+'-'+str(args.active)+'.txt','w')
    # else:
    #     loss = open(args.model_path+args.dataset+'/'+'loss-'+args.application+'-'+str(args.factor)+'.txt','w')

    optimizer = optim.Adam(model.parameters(), lr=args.lr,betas=(0.9,0.999),weight_decay=1e-6)
    criterion = nn.MSELoss()

    t = 0
    for itera in range(1,args.num_epochs+1):
        torch.cuda.empty_cache()
        train_loader = dataset.GetTrainingData()
        x = time.time()

        print('======='+str(itera)+'========')
        loss_mse = 0
        loss_grad = 0
        
        for batch_idx, (coord,v) in enumerate(train_loader):
            t1 = time.time()
            if args.cuda:
                coord = coord.cuda()
                v = v.cuda()
            optimizer.zero_grad()
            v_pred = model(coord)
            mse = criterion(v_pred.view(-1),v.view(-1))
            mse.backward()
            loss_mse += mse.mean().item()
            optimizer.step()
            #print(time.time()-t1)
        
        y = time.time()
        t += y-x
        print(y-x)
        print("Epochs "+str(itera)+": loss = "+str(loss_mse))
        loss.write("Epochs "+str(itera)+": loss = "+str(loss_mse))
        loss.write('\n')

        # if itera % args.checkpoint == 0 or itera == 1:
        if itera == args.num_epochs:
            if args.application in ['spatial', 'super-spatial']:
                ckpt = args.model_path + args.dataset + '/' + \
                    f'{args.application}-{args.scale}-{args.init}-{args.factor}-{itera}.pth'
            elif args.application == 'temporal':
                ckpt = args.model_path + args.dataset + '/' + \
                    f'{args.application}-{args.interval}-{args.init}-{args.factor}-{args.active}-{itera}.pth'
            else:
                # includes 'super-spatial-temporal'
                ckpt = args.model_path + args.dataset + '/' + \
                    f'{args.application}-{args.init}-{args.num_res}.pth'
            torch.save(model.state_dict(), ckpt)
    loss.write("Time = "+str(t))
    loss.write('\n')
    loss.close()
    
     # ==== Save full TorchScript model in float64 ====
    print("Converting model to float64 TorchScript .pt ...")

    # Move to CPU (so you can load in LibTorch without GPU dependency)
    model = model.cpu().to(torch.float64)
    model.eval()

    # Infer input dimension (3 for super-spatial-temporal, else 4)
    in_dim = 3 if args.application == 'super-spatial-temporal' else 4

    # Create example input in float64
    example = torch.zeros(1, in_dim, dtype=torch.float64)

    # Trace or script the model
    try:
        ts_model = torch.jit.trace(model, example)
    except Exception as e:
        print(f"[warn] trace failed ({e}), falling back to script()")
        ts_model = torch.jit.script(model)

    ts_model = torch.jit.freeze(ts_model)

    # Save TorchScript file
    ts_path = os.path.join(
        args.model_path, args.dataset,
        f"{args.application}-{args.init}-{args.num_res}-float64.pt"
    )

    ts_model.save(ts_path)
    print(f"[train] Saved full TorchScript model to: {ts_path}")
    


    

def adjust_lr(args, optimizer, epoch):
    lr = args.lr * (0.5 ** (epoch // 50))
    for param_group in optimizer.param_groups:
        param_group['lr'] = lr



def inf2(dataset,args):

    if args.application != 'viewsynthesis':
        in_dim = 3 if args.application == 'super-spatial-temporal' else 4
        if args.active == 'sine':
            model =  CoordNet(in_dim,1,args.omega_0,args.init,args.num_res)
    model.cuda()

    if args.application == 'super-spatial-temporal':
        # Build model if not built above (should be already built)
        # Load checkpoint using the same naming as training "else" branch
        model.load_state_dict(torch.load(
            args.model_path + args.dataset + '/' +
            f'{args.application}-{args.init}-{args.num_res}.pth'
        ))
        model.eval()

        out_dtype=np.float32        
        # Coords: [T*H*W, 3] with (t,y,x) in [-1,1], provided by ScalarDataSet.GetTestingData()
        coords = dataset.GetTestingData()
        T = dataset.total_samples
        H, W = dataset.dim

        loader = DataLoader(dataset=torch.FloatTensor(coords), batch_size=args.batch_size, shuffle=False)
        preds = []
        for batch in loader:
            with torch.no_grad():
                v_pred = model(batch.cuda())
            preds.append(v_pred.view(-1).detach().cpu().numpy())
        preds = np.concatenate(preds, axis=0).astype('<f')  # length T*H*W

        # Existing context: preds is length T*H*W in blocks of H*W per t,
        # where within each block values are in Fortran 'F' order (y-fastest).

        # Rebuild a [T, H, W] volume where A[t] is the HxW slice at time t.
        A = np.empty((T, H, W), dtype=np.float64)
        per_t = H * W
        for t in range(T):
            v = preds[t * per_t:(t + 1) * per_t]      # length H*W, y-fastest inside
            A[t] = v.reshape(W, H, order='F').transpose().astype(np.float64)         # restore the 2D slice

        # Now ravel in C-order so x (last axis) is fastest, then y, then t (slowest).
        onefile_path = os.path.join('../Result', args.dataset,
                                    f'{args.application}-{args.init}-{args.num_res}.dat')
        A.ravel(order='C').tofile(onefile_path)  # float64
        print(f"Saved single 3D file with x-fastest layout to: {onefile_path}")

        # If you prefer float64:
        # A.astype('<f8', copy=False).ravel(order='C').tofile(onefile_path.replace('.dat','-f64.dat'))
        
        ts_path = args.model_path + args.dataset + '/' + f'{args.application}-{args.init}-{args.num_res}.pt'


        
        model.eval()
        dev = next(model.parameters()).device
        example = torch.zeros(1, in_dim, dtype=torch.float32,device=dev)
        try:
            ts_model = torch.jit.trace(model, example)
        except Exception as e:
            print(f"[warn] trace failed ({e}), using script()")
            ts_model = torch.jit.script(model)
        ts_model = torch.jit.freeze(ts_model)
        ts_model.save(ts_path)


def _build_model_from_args(args):
    """Recreate the Python nn.Module with the right topology."""
    in_dim = 3 if args.application == 'super-spatial-temporal' else 4
    if getattr(args, 'active', 'sine') == 'sine':
        return CoordNet(in_dim, 1, args.omega_0, args.init, args.num_res)
    elif getattr(args, 'active', 'sine') == 'relu':
        return CoordNetReLU(in_dim, 1, args.init, args.num_res)
    else:
        # default to sine if unknown
        return CoordNet(in_dim, 1, args.omega_0, args.init, args.num_res)

def inf(dataset, args):
    """
    Inference with a float64 TorchScript (.pt) model.
    - Loads the full scripted/traced model (no state_dict needed)
    - Evaluates in float64 on CPU or CUDA (if available and requested)
    - Builds a fixed-z slice on a uniform (y,x) grid
    - Saves outputs as float64 in x-fastest layout
    # """
    # import os
    # import numpy as np
    # import torch
    # from torch.utils.data import DataLoader

    # ---- device selection
    use_cuda = getattr(args, 'cuda', False) and torch.cuda.is_available()
    device = torch.device('cuda') if use_cuda else torch.device('cpu')
    print(f"[inf] device = {device}")

    # ---- resolve model path (prefer '-float64.pt', fallback to '.pt')
    base = f"{args.application}-{args.init}-{args.num_res}"
    pt64_path = os.path.join(args.model_path, args.dataset, base + "-float64.pt")
    pt32_path = os.path.join(args.model_path, args.dataset, base + ".pt")
    
    ckpt_pth=os.path.join(
            args.model_path, args.dataset,
            f'{args.application}-{args.init}-{args.num_res}-{args.num_epochs}.pth'
        )

    ts_to_use = None
    
    if os.path.exists(ckpt_pth):
        print(f"[inf] Found checkpoint: {ckpt_pth} — exporting TorchScript float64...")
        py_model = _build_model_from_args(args)
        # load weights
        sd = torch.load(ckpt_pth, map_location="cpu")
        py_model.load_state_dict(sd, strict=True)
        py_model = py_model.cpu().to(torch.float64).eval()

        # export TS float64
        in_dim = 3 if args.application == 'super-spatial-temporal' else 4
        example = torch.zeros(1, in_dim, dtype=torch.float64)
        try:
            ts_model = torch.jit.trace(py_model, example)
        except Exception as e:
            print(f"[warn] trace failed ({e}), falling back to script()")
            ts_model = torch.jit.script(py_model)
        ts_model = torch.jit.freeze(ts_model)
        ts_model.save(pt64_path)
        print(f"[inf] Exported TorchScript float64 to: {pt64_path}")
        ts_to_use = pt64_path
    
    
    else: # ---- if only 32-bit .pt exists, convert to 64-bit and save
        if os.path.exists(pt64_path):
            print(f"[inf] Using existing TorchScript: {pt64_path}")
            ts_to_use = pt64_path
        elif os.path.exists(pt32_path):
            print(f"[inf] {pt64_path} not found; converting {pt32_path} -> float64 .pt ...")
            # Load on CPU to avoid CUDA dependency in saved file
            m32 = torch.jit.load(pt32_path, map_location="cpu")
            m32.eval()
            # Cast internal tensors to float64 if supported
            try:
                m32 = m32.to(dtype=torch.float64)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to cast TorchScript module to float64; "
                    f"re-save the model as float64 at training time. Details: {e}"
                )
            # Freeze and save float64 module
            m32 = torch.jit.freeze(m32)
            m32.save(pt64_path)
            ts_to_use = pt64_path
            print(f"[inf] Saved float64 TorchScript to: {pt64_path}")
        else:
            raise FileNotFoundError(f"No model found: {ckpt_pth} / {pt64_path} / {pt32_path}")

    # ---- load the float64 TorchScript for inference
    print(f"[inf] Loading TorchScript: {pt64_path}")
    model = torch.jit.load(pt64_path, map_location=device)
    model.eval()
    try:
        model = model.to(device)
    except Exception:
        pass  # Some TS wrappers don't expose .to; inputs will still be float64

    
    # Many TS models already carry dtype; we still feed float64 inputs to ensure double pipeline
    # If your TS is truly float32, sending float64 inputs may upcast internally (or you can
    # re-save as -float64.pt as shown in train()).

    # ---- build a fixed-z slice grid in float64 (normalized coords in [-1,1])
    # You can expose these via args if you like:
    # z_fixed = getattr(args, "z_fixed", 0.0)   # normalized z in [-1, 1]
    # H = getattr(args, "H", 3000)
    # W = getattr(args, "W", 3000)

    # y_lin = np.linspace(-1.0, 1.0, H, dtype=np.float64)
    # x_lin = np.linspace(-1.0, 1.0, W, dtype=np.float64)
    # yy, xx = np.meshgrid(y_lin, x_lin, indexing='ij')       # yy: HxW, xx: HxW
    # zz = np.full_like(yy, fill_value=float(z_fixed))        # fixed z everywhere
    # coords = np.stack([zz, yy, xx], axis=-1).reshape(-1, 3).astype(np.float64)
    
    # Coord order: (z, y, x)  -> shape [H*W, 3]
    coords = dataset.GetTestingData(up_sample_ratio=args.up_sample_ratio,type=np.float64)

    print(f"[inf] Total coords = {coords.shape[0]}, dim = {coords.shape[1]}")
    
    # ---- robust conversion to torch.DoubleTensor
    if torch.is_tensor(coords):
        coords_tensor = coords if coords.dtype == torch.float64 else coords.to(torch.float64)
    else:
        try:
            coords_np = np.asarray(coords, dtype=np.float64)
            coords_tensor = torch.from_numpy(coords_np)
        except Exception:
            coords_tensor = torch.tensor(coords, dtype=torch.float64)


    if coords_tensor.dtype != torch.float64:
        coords_tensor = coords_tensor.to(torch.float64)
    
    in_dim = int(coords_tensor.shape[1])
    preferred_dtype = None
    with torch.no_grad():
        # Try float64 first
        try:
            probe64 = torch.zeros(1, in_dim, dtype=torch.float64, device=device)
            _ = model(probe64)
            preferred_dtype = torch.float64
            print("[inf] model accepted float64 inputs")
        except Exception as e64:
            # Try float32
            try:
                probe32 = torch.zeros(1, in_dim, dtype=torch.float32, device=device)
                _ = model(probe32)
                preferred_dtype = torch.float32
                print("[inf] model accepted float32 inputs")
            except Exception as e32:
                raise RuntimeError(
                    "TorchScript model rejected both float64 and float32 inputs. "
                    "Please export a valid .pt (consider re-tracing without freeze) "
                    f"\nfloat64 error: {e64}\nfloat32 error: {e32}"
                )
        
    # ---- run inference in batches (float64 end-to-end)
    pin_mem = device.type == 'cuda'
    loader = DataLoader(
        dataset=coords_tensor,  # already float64
        batch_size=args.batch_size,
        shuffle=False,
        pin_memory=pin_mem
    )

    preds = []
    with torch.no_grad():
        for bi, batch in enumerate(loader, 1):
            if device.type == 'cuda':
                torch.cuda.empty_cache()
            # cast input to model dtype
            batch = batch.to(device, non_blocking=pin_mem).to(preferred_dtype)
            out = model(batch)  # [B,1] or [B]
            preds.append(out.view(-1).detach().cpu().to(torch.float64).numpy())
            # if bi % 10 == 0:
                # print(f"[inf] processed batch {bi}")
        if device.type == 'cuda':
            torch.cuda.synchronize()

    preds = np.concatenate(preds, axis=0)  # float64 output

    # ---- save outputs (float64, x-fastest when raveled in C-order)
    result_dir = os.path.join('../Result', args.dataset)
    os.makedirs(result_dir, exist_ok=True)
    out_path = os.path.join(result_dir, f"{base}-{args.up_sample_ratio}.dat")

    # reshape to [H,W] with x as last axis so C-order ravel is x-fastest
    # A = preds.reshape(H, W)                   # rows=y (slow), cols=x (fast)
    
    
    # T = dataset.total_samples
    T = getattr(dataset, 'total_samples', None)
    dim = getattr(dataset, 'dim', None)

    if T is not None and isinstance(dim, (tuple, list)) and len(dim) == 2:
        shape = dataset.span_num()
        shape = args.up_sample_ratio * shape

        H, W = shape[0], shape[1]
        T = shape[2]

        per_t = H * W
        if T * per_t == preds.size:
            # Rebuild a [T, H, W] volume; incoming per-slice stream is y-fastest inside
            A = np.empty((T, H, W), dtype=np.float64)
            for t in range(T):
                v = preds[t * per_t:(t + 1) * per_t]
                # reshape with Fortran order then transpose -> (H, W) with x-fastest when raveled in C
                A[t] = v.reshape(W, H, order='F').T
            A.ravel(order='C').astype('<f8').tofile(out_path)
            print(f"[inf] Saved 3D field (float64, x-fastest) to: {out_path}")
            return
        else:
            print(f"[inf][warn] Size mismatch: got {preds.size}, expected {T*per_t}. Saving flat.")

    # Fallback: save flat array
    preds.astype('<f8').tofile(out_path)
    print(f"[inf] Saved flat array (float64) to: {out_path}")