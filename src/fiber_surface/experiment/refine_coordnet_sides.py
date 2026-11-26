"""Locate and repair local side-facet errors against a staged CoordNet field."""
import argparse
from hashlib import sha256
from pathlib import Path
import shlex

import numpy as np
import torch

from repair_coordnet_caps import Query, component_ids, read_ply, run_parameters, write_ply
from fit_coordnet_tube import Tube


def settings(run):
    cfg=run_parameters(run)
    records=dict((row[0],row[1:]) for line in (run/'stages/run.conf').read_text().splitlines()
                 if (row:=shlex.split(line)))
    options=records['options']
    surface=float(records.get('solver_tolerances',[0,0])[1])
    cfg.update(spatial=float(options[3]),derivative=float(options[5]),
               solve=min(float(options[4]),surface if surface>0 else 1e-10),
               angle=float(options[7]),iterations=int(options[9]))
    return cfg


def samples(tri):
    a,b,c = tri.transpose(1,0,2)
    return np.stack([a,b,c,a+(b-a)/2,b+(c-b)/2,c+(a-c)/2,a+((b-a)+(c-a))/3],axis=1)


class SideField(Query):
    def __init__(self,config,workers):
        super().__init__(config,workers)
        self.spatial=config.get('spatial',config['step'])
        self.solve=config.get('solve',1e-10)
        self.derivative=config.get('derivative',1e-10)
        self.iterations=config.get('iterations',40)

    def differential(self,p):
        x = torch.from_numpy(((p-self.origin)*self.scale-1)[:,::-1].copy()).requires_grad_(True)
        values, gradients = [], []
        for model in self.models:
            y = model(x).reshape(-1)
            values.append(y.detach().numpy())
            gradients.append(torch.autograd.grad(y.sum(),x)[0].numpy()[:,::-1]*self.scale)
        phi = np.stack(values,axis=1)@self.normal.T-self.offset
        jac = np.einsum('ki,nij->nkj', self.normal, np.stack(gradients,axis=1))
        return phi,jac

    def derivatives(self,p):
        blocks = np.array_split(p,max(1,(len(p)+1023)//1024))
        vals,grads = zip(*self.pool.map(self.differential,blocks))
        return np.concatenate(vals),np.concatenate(grads)

    def project(self,p,side,phi,jac,spatial,step,*,follow_wall=False):
        """Check Newton roots. Mesh edits also have to pass with the original, interior-only steps."""
        x = p.copy()
        fixed = (np.abs(x-self.lo)<=step*1e-9) | (np.abs(x-self.hi)<=step*1e-9)
        x = np.where(fixed,np.where(np.abs(x-self.lo)<np.abs(x-self.hi),self.lo,self.hi),x)
        active = np.zeros(phi.shape,bool)
        active[np.arange(len(x)),side] = True
        todo = np.arange(len(x))
        ok = np.zeros(len(x),bool)
        for _ in range(self.iterations):
            active[todo] |= phi < -self.solve
            done = (~active[todo] | (np.abs(phi)<=self.solve)).all(1)
            ok[todo[done]] = True
            todo,phi,jac = todo[~done],phi[~done],jac[~done]
            if not len(todo):
                break
            j = jac*(~fixed[todo])[:,None,:]
            if follow_wall:
                start=x[todo].copy()
            single = active[todo].sum(1)==1
            eligible = np.ones(len(todo),bool)
            if single.any():
                ids = np.flatnonzero(single)
                g = j[ids,side[todo[ids]]]
                norm = np.einsum('ij,ij->i',g,g)
                valid = norm>1e-24
                eligible[ids[~valid]] = False
                take = ids[valid]
                x[todo[take]] -= g[valid]*(phi[take,side[todo[take]]]/norm[valid])[:,None]
            for k in np.flatnonzero(~single):
                rows = active[todo[k]]
                delta,_,rank,_ = np.linalg.lstsq(j[k,rows],-phi[k,rows],rcond=self.derivative)
                if rank<rows.sum():
                    eligible[k] = False
                else:
                    x[todo[k]] += delta
            bounded = ((x[todo]>=self.lo)&(x[todo]<=self.hi)).all(1)
            if follow_wall:
                for k in np.flatnonzero(eligible & ~bounded):
                    idx=todo[k]
                    delta=x[idx]-start[k]
                    lower=x[idx]<self.lo
                    hit=lower | (x[idx]>self.hi)
                    wall=np.where(lower,self.lo,self.hi)
                    times=np.full(3,np.inf)
                    times[hit]=(wall[hit]-start[k,hit])/delta[hit]
                    fraction=times.min()
                    contact=times==fraction
                    x[idx]=np.clip(start[k]+fraction*delta,self.lo,self.hi)
                    x[idx,contact]=wall[contact]
                    fixed[idx,contact]=True
                    bounded[k]=True
            bounded &= np.linalg.norm(x[todo]-p[todo],axis=1)<=spatial
            todo = todo[eligible & bounded]
            if not len(todo):
                break
            phi,jac = self.derivatives(x[todo])
        return x,ok


def inspect(field,p,f,labels,rows,step,spatial):
    tri=p[f[rows]]
    sites=samples(tri)
    n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    n/=np.maximum(np.linalg.norm(n,axis=1)[:,None],1e-300)
    side=labels[rows].astype(int)-1
    cos=np.empty((len(rows),7))
    valid=np.zeros((len(rows),7),bool)
    for k in np.unique(side):
        group=np.flatnonzero(side==k)
        for ix in np.array_split(group,max(1,(len(group)+15999)//16000)):
            xyz,back=np.unique(sites[ix].reshape(-1,3),axis=0,return_inverse=True)
            phi,jac=field.derivatives(xyz)
            gradient=-jac[:,k]
            gradient/=np.maximum(np.linalg.norm(gradient,axis=1)[:,None],1e-300)
            cos[ix]=np.einsum("ij,ikj->ik",n[ix],gradient[back.reshape(-1,7)])
            _,ok=field.project(xyz,np.full(len(xyz),k),phi,jac,spatial,step)
            valid[ix]=ok[back.reshape(-1,7)]
    return cos,valid


def diagnose(run,mesh,output,workers):
    config=settings(run)
    query=SideField(config,workers)
    p,f,l,_=read_ply(mesh)
    rows=np.flatnonzero(l<=3)
    cos,ok=inspect(query,p,f,l,rows,config['step'],query.spatial)
    with output.open('wb') as stream:
        np.savez(stream,rows=rows,cos=cos,ok=ok,mesh_hash=sha256(mesh.read_bytes()).hexdigest(),
                 config_hash=config['config_sha256'])
    ids=component_ids(p,f[rows])
    tri=p[f[rows]]
    sites=samples(tri)
    depth=np.minimum(sites-query.lo,query.hi-sites).min(2)
    for name,mask in [('failed projection',~ok),('reversed',cos<=0)]:
        distances=depth[mask]
        print(name,len(distances),'sites on box',int((distances<1e-9).sum()),
              'within h',int((distances<=config['step']).sum()),
              'max distance',distances.max(initial=0),flush=True)
    for group in np.unique(ids[f[rows[(cos.min(1)<=0)|~ok.all(1)],0]]):
        choose=ids[f[rows,0]]==group
        vertices=np.unique(f[rows[choose]])
        print(group,'faces',int(choose.sum()),'label',np.unique(l[rows[choose]]).tolist(),
              'failed sites',int((~ok[choose]).sum()),'reversed centers',int((cos[choose,6]<=0).sum()),
              'reversed faces',int((cos[choose].min(1)<=0).sum()),
              'box',p[vertices].min(0).round(5),p[vertices].max(0).round(5),flush=True)


def edge_points(field,p,edges,sides,step):
    a,b=p[edges[:,0]],p[edges[:,1]]
    mid=a+(b-a)/2
    x=mid.copy()
    fixed=((a==field.lo)&(b==field.lo)) | ((a==field.hi)&(b==field.hi))
    direction=(b-a)*(~fixed)
    direction/=np.maximum(np.linalg.norm(direction,axis=1)[:,None],1e-300)
    phi,jac=field.derivatives(x)
    todo=np.arange(len(x))
    limit=np.linalg.norm(a-b,axis=1)*.8
    for _ in range(25):
        r=phi[np.arange(len(todo)),sides[todo]]
        use=np.abs(r)>3e-11
        if not use.any():
            break
        g=jac[np.arange(len(todo)),sides[todo]]*(~fixed[todo])
        todo,r,g=todo[use],r[use],g[use]
        phase=direction[todo]
        g-=phase*np.einsum('ij,ij->i',g,phase)[:,None]
        dx=-g*(r/np.maximum(np.einsum('ij,ij->i',g,g),1e-30))[:,None]
        # Stay inside the box during the solve. Edges close to a tangential
        # contact can require a different in-domain point on the same side.
        raw=x[todo]+dx
        x[todo]=np.clip(raw,field.lo,field.hi)
        phi,jac=field.derivatives(x[todo])
    phi=field.batch(x)
    ok=(np.abs(phi[np.arange(len(x)),sides])<1e-9)&(phi.min(1)>-1e-9)
    ok&=np.linalg.norm(x-mid,axis=1)<np.maximum(limit,step*1e-6)
    # Only original rim edges may introduce a new vertex on a box plane.
    ok&=~(((x==field.lo)|(x==field.hi))&~fixed).any(1)
    return x,ok


def tighten_rims(run,mesh,output,workers,tolerance,angle):
    """Split curved rim chords while checking both new side facets on the fitted field."""
    cfg=settings(run)
    field=SideField(cfg,workers)
    p,f,l,c=read_ply(mesh)
    side=l<=3
    f,l=f[side],l[side]
    step=cfg['step']
    for iteration in range(12):
        edges=np.stack([f,f[:,[1,2,0]]],axis=2)
        keys=edges.min(2).astype(np.int64)*len(p)+edges.max(2)
        _,first,count=np.unique(keys.ravel(),return_index=True,return_counts=True)
        faces,local=np.divmod(first[count==1],3)
        a,b,cid=f[faces,local],f[faces,(local+1)%3],f[faces,(local+2)%3]
        on=((p[a]==field.lo)&(p[b]==field.lo)) | ((p[a]==field.hi)&(p[b]==field.hi))
        if not on.any(1).all():
            raise ValueError("Expected complete side rims on box planes")
        x,fit=edge_points(field,p,np.column_stack([a,b]),l[faces].astype(int)-1,step)
        error=np.linalg.norm(x-(p[a]+(p[b]-p[a])/2),axis=1)
        candidates=np.flatnonzero(fit & (error>tolerance))
        # If a triangle touches two rims, refine the larger deviation first.
        candidates=candidates[np.argsort(-error[candidates],kind='stable')]
        _,order=np.unique(faces[candidates],return_index=True)
        k=candidates[order]
        if not len(k):
            break
        new=np.arange(len(k))+len(p)
        tris=np.stack([np.column_stack([a[k],new,cid[k]]),np.column_stack([new,b[k],cid[k]])],axis=1)
        d,valid=inspect(field,np.vstack([p,x[k]]),tris.reshape(-1,3),np.repeat(l[faces[k]],2),
                        np.arange(2*len(k)),step,field.spatial)
        accepted=valid.reshape(-1,14).all(1)&(d.reshape(-1,14).min(1)>=np.cos(np.deg2rad(angle)))
        print('Rim pass',iteration,'accepted',accepted.sum(),'of',len(k),flush=True)
        if not accepted.any():
            break
        k=k[accepted]
        new=len(p)+np.arange(len(k))
        p=np.vstack([p,x[k]])
        f[faces[k]]=np.column_stack([a[k],new,cid[k]])
        f=np.vstack([f,np.column_stack([new,b[k],cid[k]])])
        l=np.r_[l,l[faces[k]]]
    write_ply(output,p,f,l,c)


def move_vertices(field,p,f,l,cos,ok,step,threshold):
    """Accept field-fitted vertex moves only after checking their complete triangle fans."""
    bad=(cos.min(1)<threshold)|~ok.all(1)
    vertices=np.unique(f[bad])
    if not len(vertices):
        return 0
    flat=f.ravel()
    order=flat.argsort()
    offsets=np.searchsorted(flat[order],np.arange(len(p)+1))
    targets=[]; orig=[]; stars=[]
    for v in vertices:
        fan=np.unique(order[offsets[v]:offsets[v+1]]//3)
        nb=np.setdiff1d(np.unique(f[fan]),[v])
        boundary=(p[v]==field.lo)|(p[v]==field.hi)
        if boundary.any():
            # Smooth the contour in its own plane, using its two rim neighbors.
            nb=nb[((p[nb][:,boundary]==p[v,boundary]).all(1))]
        if len(nb)<2:
            continue
        delta=p[nb].mean(0)-p[v]
        proposals=[delta*.4,delta*.8]
        depth=np.minimum(p[v]-field.lo,field.hi-p[v])
        d=depth.argmin()
        if not boundary.any() and depth[d]<step*.3:
            inward=np.zeros(3)
            inward[d]=1 if p[v,d]-field.lo[d]<field.hi[d]-p[v,d] else -1
            scale=np.median(np.linalg.norm(p[nb]-p[v],axis=1))
            proposals.extend(inward*scale*factor for factor in (.04,.12,.3))
        for delta in proposals:
            x=p[v]+delta
            x[boundary]=p[v,boundary]
            if ((x>=field.lo)&(x<=field.hi)).all():
                targets.append(x);orig.append(v);stars.append(fan)
    if not targets:
        return 0
    x=np.asarray(targets); orig=np.asarray(orig)
    side=np.array([l[fan[0]]-1 for fan in stars],int)
    fixed=(p[orig]==field.lo)|(p[orig]==field.hi)
    valid=np.ones(len(x),bool)
    for _ in range(16):
        ix=np.flatnonzero(valid)
        if not len(ix):
            break
        phi,j=field.derivatives(x[ix])
        val=phi[np.arange(len(ix)),side[ix]]
        take=np.abs(val)>2e-11
        if not take.any():
            break
        g=j[np.arange(len(ix)),side[ix]]*(~fixed[ix])
        ids=ix[take]
        dx=g[take]*(val[take]/np.maximum((g[take]**2).sum(1),1e-25))[:,None]
        x[ids]-=dx
        valid[ids]=((x[ids]>=field.lo)&(x[ids]<=field.hi)).all(1)&(np.linalg.norm(x[ids]-p[orig[ids]],axis=1)<step)
    ix=np.flatnonzero(valid)
    if not len(ix):
        return 0
    phi=field.batch(x[ix])
    valid[ix]=(np.abs(phi[np.arange(len(ix)),side[ix]])<1e-10)&(phi.min(1)>=-1e-10)
    ix=np.flatnonzero(valid)
    if not len(ix):
        return 0
    triangles=[];labels=[]
    for k in ix:
        faces=f[stars[k]]
        a=p[faces]
        a[faces==orig[k]]=x[k]
        triangles.extend(a);labels.extend(l[stars[k]])
    pts=np.asarray(triangles).reshape(-1,3)
    faces=np.arange(len(pts)).reshape(-1,3)
    d,success=inspect(field,pts,faces,np.asarray(labels),np.arange(len(faces)),step,field.spatial)
    def cost(normals,valid):
        return (10000*(~valid).sum()+100*(normals<=0).sum()
                +np.maximum(0,threshold-normals.min(1)).sum())
    changes=[]
    pos=0
    for k in ix:
        fan=stars[k]; sl=slice(pos,pos+len(fan));pos+=len(fan)
        if ((d[sl]<=0).sum()>(cos[fan]<=0).sum() or
                (~success[sl]).sum()>(~ok[fan]).sum() or
                d[sl].min()<min(.02,cos[fan].min())-1e-9):
            continue
        gain=cost(cos[fan],ok[fan])-cost(d[sl],success[sl])
        if gain>1e-6:
            changes.append((gain,k,sl))
    occupied=set();count=0
    for _,k,sl in sorted(changes,key=lambda a:-a[0]):
        fan=stars[k]
        if any(i in occupied for i in fan):
            continue
        p[orig[k]]=x[k];cos[fan]=d[sl];ok[fan]=success[sl]
        occupied.update(fan);count+=1
    return count


def repair(run,mesh,output,workers,angle,diagnostic):
    cfg=settings(run)
    field=SideField(cfg,workers)
    p,f,l,c=read_ply(mesh)
    if any(c):
        raise ValueError("This component repair requires disjoint range sides")
    rows=np.flatnonzero(l<=3)
    if diagnostic:
        with np.load(diagnostic) as data:
            if ('mesh_hash' not in data or data['mesh_hash']!=sha256(mesh.read_bytes()).hexdigest()
                    or data['config_hash']!=cfg['config_sha256'] or not np.array_equal(rows,data['rows'])):
                raise ValueError("Diagnostic does not match this input mesh and staged run")
            cos,ok=data['cos'],data['ok']
    else:
        cos,ok=inspect(field,p,f,l,rows,cfg['step'],field.spatial)
    ids=component_ids(p,f[rows])
    replaced=np.unique(ids[f[rows[cos[:,6]<=0],0]])
    groups=ids[f[rows,0]]
    keep=~np.isin(groups,replaced)
    replacements=[]
    for group in replaced:
        old=rows[groups==group]
        label=np.unique(l[old])
        if len(label)!=1:
            raise ValueError("Expected one range edge per tube")
        seed=p[np.unique(f[old])]
        tube=Tube(field,int(label[0])-1,seed,cfg['step'])
        vertices,faces=tube.mesh(seed,cfg['step'],32)
        replacements.append((vertices,faces,np.full(len(faces),label[0],dtype=l.dtype)))
    f,l,cos,ok=f[rows[keep]],l[rows[keep]],cos[keep],ok[keep]
    threshold=np.cos(np.deg2rad(angle))
    for k in range(22):
        bad=(cos.min(1)<threshold)|~ok.all(1)
        print('Pass',k,'flagged',bad.sum(),'failed',int((~ok).sum()),
              'angle',np.rad2deg(np.arccos(cos.min().clip(-1,1))),flush=True)
        if not bad.any() or not move_vertices(field,p,f,l,cos,ok,cfg['step'],threshold):
            break
    points,faces,labels=[p],[f],[l]
    offset=len(p)
    for vertices,triangles,patch in replacements:
        points.append(vertices)
        faces.append(triangles+offset)
        labels.append(patch)
        offset+=len(vertices)
    write_ply(output,np.vstack(points),np.vstack(faces),np.concatenate(labels),c)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    parser.add_argument("mesh",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--workers",type=int,default=8)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--diagnose",action="store_true",help="write a per-site diagnostic NPZ")
    mode.add_argument("--repair",action="store_true",help="write fitted side triangles for cap-sides")
    mode.add_argument("--rim-tolerance",type=float,help="refine side rim chords above this fitted midpoint deviation")
    parser.add_argument("--use-diagnostic",type=Path,help="reuse the diagnostic NPZ for this mesh and run")
    parser.add_argument("--angle",type=float,help="normal angle target in degrees, defaults to the staged limit")
    args=parser.parse_args()
    if args.mesh.resolve()==args.output.resolve() or (args.use_diagnostic and args.output.resolve()==args.use_diagnostic.resolve()):
        parser.error("Choose a separate output path")
    if args.use_diagnostic and not args.repair:
        parser.error("--use-diagnostic requires --repair")
    angle=args.angle if args.angle is not None else settings(args.run)['angle']
    if not 0<angle<90 or args.workers<1:
        parser.error("Expected a positive worker count and an angle between 0 and 90 degrees")
    if args.rim_tolerance is not None:
        if args.rim_tolerance<=0:
            parser.error("Rim tolerance must be positive")
        tighten_rims(args.run,args.mesh,args.output,args.workers,args.rim_tolerance,angle)
    elif args.repair:
        repair(args.run,args.mesh,args.output,args.workers,angle,args.use_diagnostic)
    else:
        diagnose(args.run,args.mesh,args.output,args.workers)
