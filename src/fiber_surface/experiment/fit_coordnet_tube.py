"""Fit time slices of a thin, locally star-shaped CoordNet level-set component."""

import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree


class Tube:
    def __init__(self,field,side,points,step):
        self.field,self.side=field,side
        self.lo=np.maximum(field.lo[:2],points[:,:2].min(0)-2*step)
        self.hi=np.minimum(field.hi[:2],points[:,:2].max(0)+2*step)
        t=np.median(points[:,2])
        corners=np.array([[a,b,t] for a in (self.lo[0],self.hi[0]) for b in (self.lo[1],self.hi[1])])
        outer=field.batch(corners)[:,side]
        if not ((outer>0).all() or (outer<0).all()):
            raise ValueError("The cross section is not isolated inside the local search box")
        self.sign=np.sign(outer[0])
        self.times=np.linspace(field.lo[2],field.hi[2],1+int(np.ceil((field.hi[2]-field.lo[2])/(step/8))))
        axes=[np.linspace(a,b,33) for a,b in zip(self.lo,self.hi)]
        xy=np.stack(np.meshgrid(*axes),axis=-1).reshape(-1,2)
        self.centers=np.empty((len(self.times),2))
        for group in np.array_split(np.arange(len(self.times)),max(1,len(self.times)//32)):
            x=np.column_stack([np.tile(xy,(len(group),1)),np.repeat(self.times[group],len(xy))])
            v=self.value(x).reshape(len(group),len(xy))
            self.centers[group]=xy[v.argmin(1)]
        self.centers,self.minima=self.minimize(self.times,self.centers)

    def value(self,x):
        return self.sign*self.field.batch(x)[:,self.side]

    def minimize(self,t,xy):
        x=np.column_stack([xy,t])
        for _ in range(12):
            val,j=self.field.derivatives(x)
            grad=j[:,self.side,:2]*self.sign
            h=[]
            for delta in np.eye(3)[:2]*1e-6:
                plus=self.field.derivatives(x+delta)[1][:,self.side,:2]
                minus=self.field.derivatives(x-delta)[1][:,self.side,:2]
                h.append(self.sign*(plus-minus)/2e-6)
            h=np.stack(h,axis=-1)
            eig,vec=np.linalg.eigh(h)
            eig=np.maximum(eig,1.)
            dx=np.einsum('nij,nj->ni',vec,np.einsum('nji,nj->ni',vec,grad)/eig)
            dx*=np.minimum(1,.002/np.maximum(1e-16,np.linalg.norm(dx,axis=1)))[:,None]
            proposal=np.clip(x[:,:2]-dx,self.lo,self.hi)
            if np.max(np.abs(proposal-x[:,:2]))<2e-12:
                break
            new=np.column_stack([proposal,t])
            better=self.value(new)<=self.sign*val[:,self.side]+1e-15
            x[better]=new[better]
            if not better.any():
                break
        return x[:,:2],self.value(x)

    def centers_at(self,t):
        xy=np.column_stack([np.interp(t,self.times,x) for x in self.centers.T])
        return self.minimize(t,xy)

    def pole(self,i,j):
        def fun(t):
            return self.centers_at(np.array([t]))[1][0]
        t=brentq(fun,self.times[i],self.times[j],xtol=2e-15)
        xy,v=self.centers_at(np.array([t]))
        if abs(v[0])>1e-10:
            raise ValueError("Terminal point misses the fitted field")
        return np.r_[xy[0],t]

    def mesh(self,points,step,n):
        active=self.minima<0
        starts=np.flatnonzero(active & ~np.r_[False,active[:-1]])
        ends=np.flatnonzero(active & ~np.r_[active[1:],False])
        angles=2*np.pi*np.arange(n)/n
        rays=np.column_stack([np.cos(angles),np.sin(angles)])
        meshes=[]
        for i,j in zip(starts,ends):
            if self.times[j]<points[:,2].min()-step or self.times[i]>points[:,2].max()+step:
                continue
            lower=self.pole(i-1,i) if i else None
            upper=self.pole(j,j+1) if j<len(active)-1 else None
            ta=self.times[i] if lower is None else lower[2]
            tb=self.times[j] if upper is None else upper[2]
            t=np.linspace(ta,tb,max(3,1+int(np.ceil((tb-ta)/(step/4)))))
            dt=t[1]-t[0]
            # Geometric layers resolve the change from a horizontal normal on
            # the wall to a vertical normal at a smooth terminal point.
            if lower is not None:
                t=np.r_[ta+dt*.5**np.arange(15,0,-1),t[1:]]
            if upper is not None:
                t=np.r_[t[:-1],tb-dt*.5**np.arange(1,16)]
            xy,v=self.centers_at(t)
            if (v>=0).any():
                raise ValueError("An intermediate cross section disappears")
            center=np.column_stack([np.repeat(xy,n,axis=0),np.repeat(t,n)])
            direction=np.column_stack([np.tile(rays,(len(t),1)),np.zeros(len(center))])
            # Select the first outward crossing. Guard against a contact with
            # the edge of the local search region or a second disconnected sheet.
            with np.errstate(divide='ignore',invalid='ignore'):
                limits=np.where(direction[:,:2]>0,(self.hi-center[:,:2])/direction[:,:2],
                                (self.lo-center[:,:2])/direction[:,:2])
            radius=np.where(np.abs(direction[:,:2])>1e-12,limits,np.inf).min(1)
            if (radius<=0).any():
                raise ValueError("Extremum reached the outside of the local search box")
            left=np.zeros(len(center))
            right=radius.copy()
            searching=np.ones(len(center),bool)
            for fraction in np.linspace(1/24,1,24):
                ids=np.flatnonzero(searching)
                if not len(ids):
                    break
                s=radius[ids]*fraction
                positive=self.value(center[ids]+direction[ids]*s[:,None])>0
                right[ids[positive]]=s[positive]
                left[ids[~positive]]=s[~positive]
                searching[ids[positive]]=False
            if searching.any():
                raise ValueError("Some rays have no outward crossing")
            for _ in range(34):
                mid=(left+right)/2
                positive=self.value(center+direction*mid[:,None])>0
                right=np.where(positive,mid,right)
                left=np.where(positive,left,mid)
            p=center+direction*((left+right)/2)[:,None]
            values,jac=self.field.derivatives(p)
            radial=self.sign*(jac[:,self.side]*direction).sum(1)
            if (radial<=0).any() or np.abs(values[:,self.side]).max()>1e-9 or values.min() < -1e-9:
                raise ValueError("Cross section fails its continuous-field checks")
            ring=np.arange(len(p)).reshape(-1,n)
            a,b=ring[:-1],np.roll(ring[:-1],-1,axis=1)
            c,d=ring[1:],np.roll(ring[1:],-1,axis=1)
            f=np.vstack([np.stack([a,b,d],axis=-1).reshape(-1,3),np.stack([a,d,c],axis=-1).reshape(-1,3)])
            if lower is not None:
                f=np.vstack([f,np.column_stack([ring[0],np.full(n,len(p)),np.roll(ring[0],-1)])])
                p=np.vstack([p,lower])
            if upper is not None:
                f=np.vstack([f,np.column_stack([ring[-1],np.roll(ring[-1],-1),np.full(n,len(p))])])
                p=np.vstack([p,upper])
            if self.sign>0:
                f=f[:,::-1]
            print('Fitted interval',ta,tb,'rings',len(ring),'triangles',len(f),flush=True)
            meshes.append((p,f))
        if not meshes:
            raise ValueError("No matching time interval")
        p=np.vstack([x for x,_ in meshes])
        f=np.vstack([f+offset for (_,f),offset in zip(meshes,np.cumsum([0]+[len(x) for x,_ in meshes[:-1]]))])
        missed=cKDTree(p).query(points)[0].max()
        if missed>step:
            raise ValueError(f"Refined surface misses original samples by {missed:g}")
        return p,f
