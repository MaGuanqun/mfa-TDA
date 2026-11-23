"""Replay the interior RK4 and phase-correction steps at reported wake trace failures."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from wake_reference import LO, HI, WakeField


def value_jacobian(field, xyz):
    p = torch.tensor(xyz, dtype=torch.float64, requires_grad=True)
    normalized = ((p-torch.from_numpy(LO))*torch.from_numpy(2/(HI-LO))-1).flip([0])[None, :]
    values, jac = [], []
    for model in field.models:
        out = model(normalized).reshape(())
        jac.append(torch.autograd.grad(out, p, retain_graph=True)[0].numpy())
        values.append(out.item())
    return np.array(values), np.stack(jac)


def replay(field, p, corner, h, sign):
    def tangent(x):
        _, j = value_jacobian(field, x)
        t = np.cross(j[0], j[1])
        return sign*t/np.linalg.norm(t)

    v1 = tangent(p)
    v2 = tangent(p+h*v1/2)
    v3 = tangent(p+h*v2/2)
    v4 = tangent(p+h*v3)
    anchor = p+h*(v1+2*v2+2*v3+v4)/6
    x = anchor.copy()
    for it in range(40):
        displacement = float(np.linalg.norm(x-anchor))
        if (x < LO).any() or (x > HI).any():
            return dict(outcome="left domain", iterations=it, displacement=displacement)
        if displacement > h:
            return dict(outcome="Newton displacement exceeded h", iterations=it, displacement=displacement)
        value, j = value_jacobian(field,x)
        phase = v1@(x-anchor)
        if np.linalg.norm(value-corner)<=1e-9 and abs(phase)<=1e-6*h:
            return dict(outcome="accepted", iterations=it, displacement=displacement)
        a = np.vstack((j, v1))
        if np.linalg.matrix_rank(a, tol=1e-10*np.linalg.norm(a,2))<3:
            return dict(outcome="singular phase system", iterations=it, displacement=displacement)
        x += np.linalg.solve(a, np.r_[corner-value,-phase])
    return dict(outcome="Newton iteration limit", iterations=40, displacement=float(np.linalg.norm(x-anchor)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    torch.set_num_threads(1)
    field = WakeField()
    conf = (args.trace/"run.conf").read_text().splitlines()
    h = float(next(line for line in conf if line.startswith("options ")).split()[1])
    results = []
    for i in range(3):
        report = json.loads((args.trace/f"segment_{i}/endpoint_fibers_report.json").read_text())
        for event in report["unresolved"]:
            p = np.array(event["position"])
            value, j = value_jacobian(field,p)
            corner = int(np.argmin(np.linalg.norm(field.query-value, axis=1)))
            clearance = float(np.min(np.minimum(HI-p,p-LO)))
            row = dict(position=p.tolist(), corner="ABC"[corner], box_clearance=clearance, nominal_h=h,
                       singular_values=np.linalg.svd(j,compute_uv=False).tolist(), trials=[])
            if clearance > h:
                for sign in (1, -1):
                    for step in (h, h/2, h/4):
                        row["trials"].append(dict(sign=sign,h=step,**replay(field,p,field.query[corner],step,sign)))
            results.append(row)
    args.output.write_text(json.dumps(results,indent=2)+"\n")


if __name__ == "__main__":
    main()
