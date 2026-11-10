"""Check both model formats against an independent exact field, without Python Torch."""
import argparse
import json
from pathlib import Path
import struct
import subprocess


def run(command):
    result = subprocess.run([str(x) for x in command],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(result.stdout)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler",required=True)
    parser.add_argument("--fixtures",required=True)
    parser.add_argument("--work",required=True)
    args=parser.parse_args()
    work=Path(args.work)
    run([args.fixtures,work])
    for name,lo,hi in [("boussinesq_3d",[-.5,-.5,0],[.5,2.5,1.5]),
                       ("vortex_street_3d",[-.5,-.5,13.5],[7.5,.5,15])]:
        region=[-.23,.17,-.18,.12,lo[2]+.37,lo[2]+.61]
        grids=[]
        for backend in ("mfa","inr"):
            prefix=work/(name+"_"+backend)
            model=work/(name+".mfa" if backend=="mfa" else "linear.pt")
            run([args.sampler,"--backend",backend,"--dataset",name,"--input",model,
                 "--output",prefix,"--region",*region,"--ratios",3,5,7,
                 "--batch-size",23,"--threads",2])
            info=json.loads(Path(str(prefix)+".json").read_text())
            assert info["sample_min"]==region[::2] and info["sample_max"]==region[1::2]
            nx,ny,nt=info["dimensions"]
            assert [nx,ny,nt]==[13,16,18],info["dimensions"]
            raw=Path(str(prefix)+".bin").read_bytes()
            assert len(raw)==nx*ny*nt*8
            values=struct.unpack("<"+str(nx*ny*nt)+"d",raw)
            max_error=0
            for k in range(nt):
                for j in range(ny):
                    for i in range(nx):
                        p=[region[2*a]+index*(region[2*a+1]-region[2*a])/(n-1)
                           for a,index,n in zip(range(3),(i,j,k),(nx,ny,nt))]
                        expected=sum(w*(2*(p[a]-lo[a])/(hi[a]-lo[a])-1)
                                     for a,w in enumerate((5,3,2)))
                        max_error=max(max_error,abs(values[(k*ny+j)*nx+i]-expected))
            assert max_error<1e-11,(name,backend,max_error)
            grids.append(values)
            print(name,backend,"max error",max_error)
        assert max(abs(a-b) for a,b in zip(*grids))<1e-11
        # Each fixture has one span per axis, despite the dataset's different
        # legacy base counts. This catches using hard-coded ratios as MFA spans.
        prefix=work/(name+"_spans")
        command=[args.sampler,"--backend","mfa","--dataset",name,
                 "--input",work/(name+".mfa"),"--output",prefix,
                 "--region",*region,"--points-per-span",17,"--batch-size",23]
        run(command)
        info=json.loads(Path(str(prefix)+".json").read_text())
        assert info["points_per_span"]==17 and info["mfa_span_counts"]==[1,1,1]
        assert info["sampling_mode"]=="mfa-knot-spans" and info["ratios"] is None
        for a in range(3):
            h=(hi[a]-lo[a])/17
            assert abs(info["spacing"][a]-h)<1e-14
            assert info["sample_min"][a]<=region[2*a]+1e-14
            assert info["sample_max"][a]>=region[2*a+1]-1e-14
            assert region[2*a]-info["sample_min"][a]<h+1e-14
            assert info["sample_max"][a]-region[2*a+1]<h+1e-14
        nx,ny,nt=info["dimensions"]
        values=struct.unpack("<"+str(nx*ny*nt)+"d",Path(str(prefix)+".bin").read_bytes())
        for k in range(nt):
            for j in range(ny):
                for i in range(nx):
                    p=[info["sample_min"][a]+idx*info["spacing"][a] for a,idx in enumerate((i,j,k))]
                    expected=sum(w*(2*(p[a]-lo[a])/(hi[a]-lo[a])-1) for a,w in enumerate((5,3,2)))
                    assert abs(values[(k*ny+j)*nx+i]-expected)<1e-11
        dry=subprocess.run([str(x) for x in command]+["--dry-run"],check=True,text=True,stdout=subprocess.PIPE)
        assert json.loads(dry.stdout[dry.stdout.index('{'):])==info
        print(name,"MFA span grid and exact linear field passed")


if __name__=="__main__":
    main()
