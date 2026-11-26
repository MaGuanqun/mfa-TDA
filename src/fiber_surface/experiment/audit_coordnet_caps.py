"""Compare triangulated box caps to direct CoordNet membership at face cell centers."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from repair_coordnet_caps import Query, read_ply, rim, run_parameters


def cover(triangles, axes):
    counts = np.zeros((len(axes[0]),len(axes[1])), dtype=np.int32)
    for tri in triangles:
        slices = [slice(np.searchsorted(x, tri[:,i].min()),
                        np.searchsorted(x, tri[:,i].max(), side="right")) for i,x in enumerate(axes)]
        u = axes[0][slices[0]][:,None]
        v = axes[1][slices[1]]
        if not u.size or not v.size:
            continue
        inside = np.ones((len(u),len(v)), bool)
        for a,b in zip(tri, np.roll(tri,-1,axis=0)):
            inside &= (b[0]-a[0])*(v-a[1]) >= (b[1]-a[1])*(u-a[0])
        counts[tuple(slices)] += inside
    return counts


def distances(points, segments):
    if not len(points):
        return np.zeros(0)
    if not len(segments):
        raise ValueError("Disagreement on a face without cap contours")
    mid = segments.mean(1)
    half = np.linalg.norm(segments[:,0]-mid, axis=1)
    tree = cKDTree(mid)
    nearest, _ = tree.query(points)
    # A segment's midpoint must lie within the current upper bound + half length.
    result = []
    for x, ids in zip(points,tree.query_ball_point(points,nearest+half.max()+1e-12)):
        a,b = segments[ids,0],segments[ids,1]
        delta = b-a
        t = ((x-a)*delta).sum(1)/(delta*delta).sum(1)
        result.append(np.linalg.norm(x-(a+t.clip(0,1)[:,None]*delta),axis=1).min())
    return np.asarray(result)


def audit(run, mesh, spacing, workers):
    config = run_parameters(run)
    query = Query(config,workers)
    p,f,labels,_ = read_ply(mesh)
    edges = rim(f[labels<=3],len(p))
    result = {}
    for face,name in enumerate(("xmin","xmax","ymin","ymax","tmin","tmax")):
        d = face//2
        u,v = (d+1)%3,(d+2)%3
        plane = (query.hi if face%2 else query.lo)[d]
        tris = p[f[labels == 4+face]]
        axes = []
        for i in (u,v):
            n = max(1, int(np.ceil((query.hi[i]-query.lo[i])/spacing)))
            axes.append(query.lo[i]+(np.arange(n)+.5)*(query.hi[i]-query.lo[i])/n)
        uv_tris = tris[:,:,[u,v]]
        if not face%2:
            uv_tris = uv_tris[:,::-1]
        counts = cover(uv_tris, axes)
        xy = np.stack(np.meshgrid(*axes,indexing="ij"), axis=-1).reshape(-1,2)
        x = np.full((len(xy),3),plane)
        x[:,[u,v]] = xy
        truth = query.batch(x).min(1) >= 0
        predicted = counts.ravel()>0
        contour = edges[np.max(np.abs(p[edges,d]-plane),axis=1)<1e-9]
        error = distances(xy[truth != predicted], p[contour][:,:,[u,v]])
        result[name] = dict(triangles=len(tris), samples=len(x), selected=int(truth.sum()),
                            false_positive=int((predicted & ~truth).sum()),
                            false_negative=int((~predicted & truth).sum()),
                            max_mismatch_distance_to_rim=float(error.max(initial=0)),
                            mismatches_beyond_h=int((error>config["step"]).sum()))
        print(name, result[name], flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("--spacing", type=float, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if not args.spacing > 0:
        parser.error("Spacing must be positive")
    result = audit(args.run,args.mesh,args.spacing,args.workers)
    args.report.write_text(json.dumps(result, indent=2)+"\n")
