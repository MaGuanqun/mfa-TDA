"""Recover incomplete side components before assembling caps from shared boundary edges."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.optimize import brentq
from skimage.measure import marching_cubes
import torch

from join_isabel_refinements import read_ply, write_ply
from sample_coordnet_grid import run_parameters


class Query:
    def __init__(self, config, workers):
        self.lo, self.hi = np.array(config["box_min"]), np.array(config["box_max"])
        self.origin = np.array(config["domain_min"])
        self.scale = 2 / (config["domain_max"] - self.origin)
        self.triangle = np.array(config["triangle"])
        tangent = np.roll(self.triangle, -1, axis=0) - self.triangle
        self.normal = tangent[:, ::-1] * [-1, 1]
        self.normal /= np.linalg.norm(self.normal, axis=1)[:, None]
        self.offset = (self.normal * self.triangle).sum(1)
        self.models = []
        for name, digest in zip(config["models"], config["model_digests"]):
            h = 14695981039346656037
            for byte in Path(name).read_bytes():
                h = ((h ^ byte) * 1099511628211) & 0xffffffffffffffff
            if f"{h:016x}" != digest:
                raise ValueError(f"Model checksum differs from staged run: {name}")
            self.models.append(torch.jit.freeze(torch.jit.load(name, map_location="cpu").eval()))
        torch.set_num_threads(1)
        self.pool = ThreadPoolExecutor(workers)
        self.values(self.lo[None])  # initialize TorchScript before concurrent inference

    @torch.inference_mode()
    def values(self, p):
        x = torch.from_numpy(((p-self.origin)*self.scale-1)[:, ::-1].copy())
        y = np.column_stack([model(x).numpy().reshape(-1) for model in self.models])
        if y.shape != (len(p),2) or not np.isfinite(y).all():
            raise ValueError("Expected two finite CoordNet outputs per domain point")
        return y @ self.normal.T - self.offset

    def batch(self, p):
        return np.concatenate(list(self.pool.map(self.values, np.array_split(p, max(1, (len(p)+1023)//1024)))))


def sample(query, dims, path):
    shape = tuple(dims[::-1])
    field = np.lib.format.open_memmap(path, mode="w+", dtype="f8", shape=(3,*shape))
    axes = [np.linspace(a, b, n) for a, b, n in zip(query.lo, query.hi, dims)]
    xy = np.stack(np.meshgrid(axes[0], axes[1]), axis=-1).reshape(-1, 2)
    for k, z in enumerate(axes[2]):
        xyz = np.column_stack([xy, np.full(len(xy), z)])
        field[:,k] = query.batch(xyz).T.reshape(3,*shape[1:])
        if k % 20 == 0:
            print("Sampled plane", k+1, "of", shape[0], flush=True)
    field.flush()


def component_ids(p, faces):
    left = faces.ravel()
    right = np.roll(faces, -1, axis=1).ravel()
    graph = coo_matrix((np.ones(len(left), np.int8), (left, right)), shape=(len(p), len(p)))
    return connected_components(graph.tocsr(), directed=False)[1]


def rim(faces, n):
    edges = np.vstack([faces[:, [i,(i+1)%3]] for i in range(3)])
    keys = edges.min(1).astype(np.int64)*n + edges.max(1)
    _, first, counts = np.unique(keys, return_index=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("Nonmanifold side mesh")
    return edges[first[counts == 1]]


def planes(p, query):
    return np.stack([np.abs(p[:,d]-value) < 1e-9 for d in range(3)
                     for value in (query.lo[d], query.hi[d])], axis=1)


def close_excursions(p, f, labels, query, step):
    edges = rim(f, len(p))
    on = planes(p, query)
    off = edges[~(on[edges[:,0]] & on[edges[:,1]]).any(1)]
    incoming = {b:a for a,b in off}
    outgoing = {a:b for a,b in off}
    moved = merged = 0
    for tip in incoming.keys() & outgoing.keys():
        a, b = incoming[tip], outgoing[tip]
        face = np.flatnonzero(on[a] & on[b])
        if len(face) != 1 or on[tip].any() or a == b:
            continue
        axis = face[0]//2
        depth = abs(p[tip,axis]-p[a,axis])
        gap = np.linalg.norm(p[a]-p[b])
        # Two coincident ends of a slit can be identified if the only common
        # neighbor is the tip. This makes its existing triangle fan interior.
        merge = gap < .002 * step and np.array_equal(on[a], on[b])
        if merge:
            if set(f[(f==a).any(1)].ravel()) & set(f[(f==b).any(1)].ravel()) != {tip}:
                continue
            vertex, target = b, p[a].copy()
        elif depth < .002 * step and gap > .1*step:
            vertex, target = tip, p[tip].copy()
            target[axis] = p[a,axis]
            side = labels[np.flatnonzero((f==tip).any(1))[0]]-1
            free = np.flatnonzero(np.arange(3) != axis)
            delta = np.eye(3)[free] * step*1e-4
            for _ in range(6):
                val = query.values(target[None])[0,side]
                if abs(val) < 1e-12:
                    break
                grad = (query.values(target+delta)[:,side]-query.values(target-delta)[:,side])/(2*step*1e-4)
                if grad@grad < 1e-18:
                    break
                target[free] -= grad*val/(grad@grad)
            if (np.linalg.norm(target-p[tip]) > .02*step or
                    np.abs(query.values(target[None])[0,side]) > 1e-9 or
                    np.any(target < query.lo) or np.any(target > query.hi)):
                continue
        else:
            continue
        rows = (f==vertex).any(1)
        triangles = p[f[rows]]
        old = np.cross(triangles[:,1]-triangles[:,0], triangles[:,2]-triangles[:,0])
        triangles[f[rows]==vertex] = target
        new = np.cross(triangles[:,1]-triangles[:,0], triangles[:,2]-triangles[:,0])
        if np.any((old*new).sum(1) <= .5*(old*old).sum(1)):
            continue
        if merge:
            f[f==vertex] = a
            merged += 1
        else:
            p[vertex] = target
            moved += 1
    print("Closed native two-edge excursions:", merged, "slits welded and", moved,
          "tips fitted on the box", flush=True)


def fitted_sides(query, fields):
    points, faces, labels = [], [], []
    offset = 0
    for side, field in enumerate(fields):
        if not field.min() < 0 < field.max():
            continue
        coords, g, _, _ = marching_cubes(field.astype(np.float32), level=0, allow_degenerate=False)
        grid = coords[:, ::-1].astype(float)
        spacing = (query.hi-query.lo)/(np.array(field.shape[::-1])-1)
        q = query.lo+grid*spacing
        at_max = grid == np.array(field.shape[::-1])-1
        q[at_max] = np.broadcast_to(query.hi, q.shape)[at_max]
        # Fit within the sign-changing carrier edge. All other coordinates,
        # including every shared box coordinate, stay fixed.
        axis = np.abs(grid-np.rint(grid)).argmax(1)
        integer = np.rint(grid).astype(int)
        row = np.arange(len(q))
        integer[row,axis] = np.minimum(np.floor(grid[row,axis]).astype(int),
                                       np.array(field.shape[::-1])[axis]-2)
        end = integer.copy()
        end[row,axis] += 1
        left = field[tuple(integer[:,::-1].T)].copy()
        right = field[tuple(end[:,::-1].T)].copy()
        fractional = grid != np.rint(grid)
        bracket = ((left >= 0) != (right >= 0)) & (fractional.sum(1)==1)
        lo, hi = integer[row,axis].astype(float), end[row,axis].astype(float)
        values = query.batch(q)
        for iteration in range(48):
            take = np.flatnonzero(bracket & (np.abs(values[:,side]) > 5e-11))
            if not len(take):
                break
            x = lo[take]+(hi[take]-lo[take])*left[take]/(left[take]-right[take])
            if iteration >= 8:
                x = (lo[take]+hi[take])/2
            ax = axis[take]
            q[take,ax] = query.lo[ax]+x*spacing[ax]
            values[take] = query.batch(q[take])
            val = values[take,side]
            same = (val >= 0) == (left[take] >= 0)
            lo[take[same]], left[take[same]] = x[same], val[same]
            hi[take[~same]], right[take[~same]] = x[~same], val[~same]
        # Lewiner's ambiguous cube cases can add an interior vertex. Float32
        # output can also round a crossing to a lattice node. Neither has a
        # uniquely identifiable carrier edge. Fit locally on the same box face.
        original = q.copy()
        free = (q != query.lo) & (q != query.hi)
        for _ in range(10):
            take = np.flatnonzero(~bracket & (np.abs(values[:,side]) > 5e-11))
            if not len(take):
                break
            x = q[take]
            grad = np.zeros_like(x)
            for d in range(3):
                active = free[take,d]
                if active.any():
                    delta = np.eye(3)[d]*spacing[d]*1e-4
                    grad[active,d] = (query.batch(x[active]+delta)[:,side]-
                                      query.batch(x[active]-delta)[:,side])/(2*delta[d])
            norm = np.einsum("ij,ij->i",grad,grad)
            if np.any(norm < 1e-20):
                raise ValueError("Singular local surface fit")
            q[take] -= grad*(values[take,side]/norm)[:,None]
            if (np.any(np.abs(q[take]-original[take])>spacing) or
                    np.any(q[take]<query.lo) or np.any(q[take]>query.hi)):
                raise ValueError("An interior cube vertex needs nonlocal recovery")
            values[take] = query.batch(q[take])
        ids = component_ids(q, g)
        other = np.delete(values, side, axis=1).min(1)
        keep = []
        for group in np.unique(ids):
            region = other[ids == group]
            if region.min() < 0 <= region.max():
                raise ValueError(f"Range edge {side+1} crosses a corner on component {group}")
            if region.min() >= 0:
                keep.append(group)
        g = g[np.isin(ids[g[:,0]],keep)]
        if not len(g):
            continue
        used, local = np.unique(g, return_inverse=True)
        g, q = local.reshape(-1,3), q[used]
        worst = used[np.abs(values[used,side]).argmax()]
        if abs(values[worst,side]) > 1e-9:
            raise ValueError(f"Range edge {side+1} residual {values[worst,side]:g} at {grid[worst]} (bracket {bracket[worst]})")
        # All triangles from the same scalar marching-cubes call have the
        # same sign convention. Evaluate a well-conditioned face to orient it.
        tri = q[g]
        normals = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
        i = np.argmax(np.linalg.norm(normals, axis=1))
        center = tri[i].mean(0)
        delta = normals[i]/np.linalg.norm(normals[i])*1e-6
        signs = query.values(np.stack([center+delta, center-delta]))[:,side]
        if signs[0] > signs[1]:
            g = g[:,::-1]
        faces.append(g+offset)
        points.append(q)
        labels.append(np.full(len(q), side+1, dtype=np.int8))
        offset += len(q)
        print("Fitted range edge",side+1,"components",len(keep),"vertices",len(q),
              "max residual",np.abs(values[used,side]).max(),flush=True)
    return np.vstack(points), np.vstack(faces), np.concatenate(labels)


def finish_box_path(p, faces, labels, query, step):
    edges = rim(faces,len(p))
    masks = planes(p,query)
    off = edges[~masks[edges].all(1).any(1)]
    if not len(off):
        return p,faces,labels
    starts, ends = set(off[:,0])-set(off[:,1]), set(off[:,1])-set(off[:,0])
    if len(starts)!=1 or len(ends)!=1 or len(off)>2:
        return None
    a,b = starts.pop(),ends.pop()
    side = int(labels[0])-1
    if (masks[a]&masks[b]).any() and len(off)==2:
        tip = off[off[:,0]==a,1][0]
        extra = np.array([[tip,a,b]])
    else:
        first,last = np.flatnonzero(masks[a]),np.flatnonzero(masks[b])
        if len(first)!=1 or len(last)!=1 or first[0]//2==last[0]//2:
            return None
        x = (p[a]+p[b])*.5
        d1,d2 = first[0]//2,last[0]//2
        free = 3-d1-d2
        x[d1],x[d2] = p[a,d1],p[b,d2]
        def value(t):
            y=x.copy()
            y[free]=t
            return query.values(y[None])[0,side]
        radius = 3*step
        ticks = np.linspace(max(query.lo[free],x[free]-radius),min(query.hi[free],x[free]+radius),91)
        values = [value(t) for t in ticks]
        roots = [brentq(value,u,v) for u,v,s,t in zip(ticks,ticks[1:],values,values[1:]) if s*t<0]
        if not roots:
            return None
        x[free] = min(roots,key=lambda t:abs(t-x[free]))
        if max(np.linalg.norm(p[a]-x),np.linalg.norm(p[b]-x))>radius:
            return None
        extra = np.array([[z,u,len(p)] for u,z in off])
        p = np.vstack([p,x])
    tri = p[extra]
    n = np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    sites = np.concatenate([tri,(tri+np.roll(tri,1,axis=1))*.5,tri.mean(1)[:,None]],axis=1)
    flat = sites.reshape(-1,3)
    # Use one-sided differences on the box so all model queries are in-domain.
    gradients = []
    for delta in np.eye(3)*step*2e-4:
        u,v = np.maximum(query.lo,flat-delta),np.minimum(query.hi,flat+delta)
        gradients.append((query.batch(v)[:,side]-query.batch(u)[:,side])/(v-u).sum(1))
    grad = np.stack(gradients,axis=1).reshape(sites.shape)
    dot = -(grad*n[:,None]).sum(2)
    if not (dot>.09*np.linalg.norm(grad,axis=2)*np.linalg.norm(n,axis=1)[:,None]).all():
        return None
    faces = np.vstack([faces,extra])
    edges = rim(faces,len(p))
    if not planes(p,query)[edges].all(1).any(1).all():
        return None
    return p,faces,np.r_[labels,np.repeat(side+1,len(extra))]


def prefer_saved_components(query, run, fine, repair_bin, q, g, qside):
    config, donor = run_parameters(run), run_parameters(fine)
    if any(config[k]!=donor[k] for k in ("model_digests","triangle","box_min","box_max","domain_min","domain_max")):
        raise ValueError("Saved finer run has different models, range query or box")
    p,f,l,corners = read_ply(fine/"volume_boundary.ply")
    if any(corners):
        raise ValueError("Finer run contains corner fibers")
    f,l = f[l<=3],l[l<=3]
    close_excursions(p,f,l,query,donor["step"])
    ids,groups = component_ids(p,f),component_ids(q,g)
    taken = set()
    parts = []
    for side in np.unique(qside):
        ix = np.unique(f[l==side])
        if not len(ix):
            continue
        tree = cKDTree(p[ix])
        for group in np.unique(groups[qside==side]):
            rows = groups[g[:,0]]==group
            patch = q[np.unique(g[rows])]
            nearest = tree.query(patch)[1]
            match = np.bincount(ids[ix[nearest]]).argmax()
            if match in taken:
                continue
            keep = ids[f[:,0]]==match
            used,local = np.unique(f[keep],return_inverse=True)
            candidate = p[used]
            def agree(a,b,tol):
                return (len(a)==len(b)==0 or len(a)>0 and len(b)>0 and
                        max(cKDTree(a).query(b)[0].max(),cKDTree(b).query(a)[0].max())<tol)
            if not agree(candidate,patch,1.4*donor["step"]):
                continue
            pm,cm = planes(patch,query),planes(candidate,query)
            if not all(agree(candidate[cm[:,k]],patch[pm[:,k]],1.2*donor["step"]) for k in range(6)):
                continue
            newfaces = local.reshape(-1,3)
            newlabels = l[keep]
            if not cm[rim(newfaces,len(candidate))].all(1).any(1).all():
                with TemporaryDirectory() as tmp:
                    src,dest = Path(tmp)/"side.ply",Path(tmp)/"filled"
                    write_ply(src,candidate,newfaces,newlabels,corners)
                    result = subprocess.run([repair_bin,run/"stages",src,dest,"--repair"],check=False)
                    if result.returncode:
                        continue
                    candidate,newfaces,newlabels,_ = read_ply(dest.with_suffix(".ply"))
            restored = finish_box_path(candidate,newfaces,newlabels,query,donor["step"])
            if restored is not None:
                taken.add(match)
                parts.append((group,restored))
    retained = ~np.isin(groups[g[:,0]],[i for i,_ in parts])
    points,faces,labels = [q],[g[retained]],[qside[g[retained,0]]]
    count = len(q)
    for _,(p,f,l) in parts:
        points.append(p);faces.append(f+count);labels.append(l)
        count += len(p)
    print("Reused finer native components",len(parts),"triangles",sum(len(v[1]) for _,v in parts),flush=True)
    return np.vstack(points),np.vstack(faces),np.concatenate(labels)


def recover(query, run, field, output, fine=None, repair_bin=None):
    if field.ndim != 4 or field.shape[0] != 3 or min(field.shape[1:]) < 2:
        raise ValueError("Expected three signed range-edge sample grids")
    dims = np.array(field.shape[:0:-1])
    probe = np.random.default_rng(77411).integers(0,dims,(256,3))
    xyz = query.lo+(query.hi-query.lo)*probe/(dims-1)
    if not np.allclose(field[(slice(None),*probe[:,::-1].T)].T,query.values(xyz),rtol=1e-12,atol=1e-12):
        raise ValueError("Grid samples disagree with the staged field and query")
    p, f, label, corners = read_ply(run / "volume_boundary.ply")
    if np.any(label > 3) or any(corners):
        raise ValueError("Expected uncapped sides and disjoint range-edge surfaces")
    step = run_parameters(run)["step"]
    close_excursions(p, f, label, query, step)
    ids = component_ids(p, f)
    native_side = np.zeros(len(p), int)
    native_side[f.ravel()] = np.repeat(label, 3)
    edges = rim(f, len(p))
    on = planes(p, query)
    common = on[edges[:,0]] & on[edges[:,1]]
    bad = set(ids[edges[~common.any(1),0]])
    q, g, qside = fitted_sides(query, field)
    gids = component_ids(q, g)
    owner = {}
    replace = set()
    for side in np.unique(qside):
        a = np.flatnonzero(native_side == side)
        b = np.flatnonzero(qside == side)
        if not len(a):
            for group in np.unique(gids[b]):
                owner[group] = set()
                replace.add(group)
            continue
        dist, nn = cKDTree(p[a]).query(q[b], workers=4)
        for group in np.unique(gids[b]):
            take = gids[b] == group
            if group in owner:
                raise ValueError("A grid component crosses a range corner")
            # Match in a neighborhood scaled by the marching step. Missing nearby
            # components require additional contact and shape checks below.
            matches, counts = np.unique(ids[a[nn[take & (dist < 1.2*step)]]], return_counts=True)
            owner[group] = set(matches[counts >= max(1, int(take.sum()*.012))])
            if not owner[group] or dist[take].max() > 1.8*step:
                replace.add(group)
    ge = rim(g, len(q))
    qplane = planes(q, query)
    gcommon = qplane[ge[:,0]] & qplane[ge[:,1]]
    if not gcommon.any(1).all():
        raise ValueError("Grid recovery has interior boundary edges")
    for d in range(6):
        side_ids = np.unique(edges[common[:,d]])
        donor_ids = np.unique(ge[gcommon[:,d]])
        for group in np.unique(gids[donor_ids]):
            b = donor_ids[gids[donor_ids] == group]
            a = side_ids[np.isin(ids[side_ids], list(owner[group]))]
            if not len(a) or cKDTree(p[a]).query(q[b])[0].max() > 1.8*step:
                replace.add(group)
    changed = True
    while changed:
        changed = False
        for group, matches in owner.items():
            if group in replace or matches & bad:
                if group not in replace or not matches.issubset(bad):
                    changed = True
                replace.add(group)
                bad.update(matches)
    unmatched = bad - set().union(*(owner[j] for j in replace))
    if unmatched:
        raise ValueError(f"No field component matched broken native components {unmatched}")
    retained = ~np.isin(ids[f[:,0]], list(bad))
    selected = np.isin(gids[g[:,0]], list(replace))
    g = g[selected]
    used, local = np.unique(g, return_inverse=True)
    g = local.reshape(-1, 3)
    q, qside = q[used], qside[used]
    print("Retained native triangles", int(retained.sum()), "of", len(f),
          "recovered triangles", len(g), "donor components", len(replace), flush=True)
    extra_labels = qside[g[:,0]]
    if fine and len(g):
        q,g,extra_labels = prefer_saved_components(query,run,fine,repair_bin,q,g,qside)
    write_ply(output, np.vstack([p,q]), np.vstack([f[retained],g+len(p)]),
              np.r_[label[retained], extra_labels], corners)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("grid", type=Path, help="signed range-edge distances in side,t,y,x order")
    parser.add_argument("--sample", nargs=3, type=int, metavar=("NX","NY","NT"))
    parser.add_argument("--output", type=Path, help="repaired uncapped side PLY")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--fine-run", type=Path, help="reuse matching native components from another saved march")
    parser.add_argument("--repair-bin", type=Path, default=Path("build/src/fiber_surface/fiber_repair_boundary"))
    args = parser.parse_args()
    if not (args.sample or args.output):
        parser.error("Specify --sample, --output, or both")
    if args.fine_run and not args.output:
        parser.error("--fine-run needs a repaired mesh output")
    inputs = [args.grid,args.run/"volume_boundary.ply"]
    if args.fine_run:
        inputs.append(args.fine_run/"volume_boundary.ply")
    if args.output and args.output.resolve() in {path.resolve() for path in inputs}:
        parser.error("The repaired mesh needs a separate output path")
    config = run_parameters(args.run)
    query = Query(config, args.workers)
    if args.sample:
        if args.grid.exists() or min(args.sample) < 2:
            parser.error("Sampling needs a new grid path and at least two nodes per axis")
        with TemporaryDirectory(dir=args.grid.parent) as temp:
            staging = Path(temp)/"grid.npy"
            sample(query, args.sample, staging)
            staging.replace(args.grid)
    if args.output:
        recover(query, args.run, np.load(args.grid, mmap_mode="r"), args.output,args.fine_run,args.repair_bin)
