"""Measure a full-domain upward-flow run, including partial marches and all original traces."""

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
import torch
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from experiment.evaluate_curved_prism import components, incidence, read_mesh, save_polydata
from experiment.ttk_geometry import compare, original_curves
from wake_reference import WAKE, WakeField, LO, HI


def read_reference(folder, lines):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(folder / ("fibers.vtp" if lines else "sides.vtp")))
    reader.Update()
    m = reader.GetOutput()
    p = vtk_to_numpy(m.GetPoints().GetData()).astype(float)
    a = m.GetLines() if lines else m.GetPolys()
    cells = vtk_to_numpy(a.GetConnectivityArray()).reshape(-1, 2 if lines else 3)
    labels = vtk_to_numpy(m.GetCellData().GetScalars())
    return [p[cells[labels == i+int(not lines)]] for i in range(3)]


def fiber_summary(group, reference, i, field, sample_data=None):
    lines = np.concatenate([np.stack((xyz[:-1], xyz[1:]), axis=1) for xyz, _, _ in group])
    vertices = np.concatenate([xyz[:-1] if closed else xyz for xyz, closed, _ in group])
    residual = np.linalg.norm(field(vertices)-field.query[i], axis=1)
    midpoints = lines.mean(axis=1)
    chord_residual = np.linalg.norm(field(midpoints)-field.query[i], axis=1)
    if sample_data is not None:
        sample_data["fiber_"+"ABC"[i]] = np.column_stack((vertices, residual))
        sample_data["fiber_"+"ABC"[i]+"_midpoints"] = np.column_stack((midpoints, chord_residual))
    return dict(corner="ABC"[i], traces=len(group), complete=sum(v[2] for v in group),
                closed=sum(v[1] for v in group), vertices=len(vertices), segments=len(lines),
                x_extent=[float(vertices[:, 0].min()),float(vertices[:, 0].max())],
                max_vertex_residual=float(residual.max()), p95_vertex_residual=float(np.quantile(residual,.95)),
                max_chord_midpoint_residual=float(chord_residual.max()),
                p95_chord_midpoint_residual=float(np.quantile(chord_residual,.95)),
                **compare(lines,reference,810+2*i,LO,HI))


def seam_counts(curves, lookup, owned, patch, corner):
    required = sorted((1+(corner-1)%3, 1+corner))
    counts = Counter()
    for curve, _, _ in curves:
        for a, b in zip(curve, curve[1:]):
            options = [owned.get(tuple(sorted((v, w))), []) for v in lookup.get(tuple(a), ())
                       for w in lookup.get(tuple(b), ())]
            matched = [sorted(patch[faces].tolist()) for faces in options if faces]
            if not matched:
                counts["no_incident_triangle"] += 1
            elif matched == [required]:
                counts["two_expected_patches"] += 1
            elif len(matched) == 1 and len(matched[0]) == 1 and matched[0][0] in required:
                counts["one_expected_patch"] += 1
            else:
                counts["other_or_ambiguous"] += 1
    return dict(counts)


def patch_metrics(field, points, faces, i, sample_data):
    a = field.query[i]
    d = field.query[(i+1)%3]-a
    edge = np.array(list(incidence(faces, 2)))
    sample_sets = dict(vertices=points[np.unique(faces)], edge_midpoints=points[edge].mean(axis=1),
                       centroids=points[faces].mean(axis=1))
    result = {}
    for name, xyz in sample_sets.items():
        delta = field(xyz)-a
        s = delta@d/(d@d)
        residual = abs(delta@[-d[1], d[0]])/np.linalg.norm(d)
        sample_data[("AB","BC","CA")[i]+"_"+name] = np.column_stack((xyz, residual, s))
        outside = np.maximum(np.maximum(-s, s-1), 0)
        result[name] = dict(n=len(s), line_max=float(residual.max()), line_p95=float(np.quantile(residual,.95)),
                            projection_outside_strict=int((outside>0).sum()),
                            projection_outside=int((outside>1e-12).sum()),
                            projection_outside_by_more_than_fiber_tolerance=int((outside*np.linalg.norm(d)>1e-9).sum()),
                            max_overshoot_range_units=float(outside.max()*np.linalg.norm(d)))
    return result


def evaluate(folder, trace, ref, dest):
    report = json.loads((folder / "volume_report.json").read_text())
    field = WakeField()
    if not np.allclose(report["triangle"], field.query, atol=1e-14, rtol=0):
        raise ValueError("Unexpected range triangle")
    dest.mkdir(exist_ok=True, parents=True)
    curves = original_curves(trace)
    config = (trace/"run.conf").read_text().splitlines()
    opts = next(row.split() for row in config if row.startswith("options "))
    if float(opts[1]) != report["step"] or list(map(len, curves)) != report["corner_fiber_components"]:
        raise ValueError("Trace and boundary runs do not match")
    source = "volume_tet_boundary.vtk" if report["complete"] else "volume_boundary.vtk"
    p, cells, patch = read_mesh(folder/source)
    owners = incidence(cells, 2)
    lookup = defaultdict(list)
    for j, pt in enumerate(p):
        lookup[tuple(pt)].append(j)
    open_edges = np.array([e for e, own in owners.items() if len(own)==1], dtype=int).reshape(-1, 2)
    on_box = np.zeros(len(open_edges), dtype=bool)
    for axis in range(3):
        for bound in (LO[axis], HI[axis]):
            on_box |= (np.abs(p[open_edges, axis]-bound)<1e-9).all(axis=1)
    reasons = Counter(event["reason"] for event in report["unresolved"])
    result = dict(step=report["step"], spatial_tolerance=report["spatial_tolerance"], complete=report["complete"],
                  tetrahedra=report["tetrahedra"], boundary_triangles=len(cells), open_edges=len(open_edges),
                  open_edges_on_box=int(on_box.sum()), open_edges_off_box=int((~on_box).sum()),
                  unresolved=dict(reasons), projection_epsilon=1e-12, fibers=[], patches=[])
    refcurves = read_reference(ref, True)
    refsides = read_reference(ref, False)
    all_curves, fiber_ids, curve_offset, indexed = [], [], 0, []
    sample_data = {}
    for i, group in enumerate(curves):
        for xyz, closed, _ in group:
            used = xyz[:-1] if closed else xyz
            all_curves.append(used)
            ids = curve_offset + np.arange(len(used))
            indexed.append(np.r_[ids, ids[0]] if closed else ids)
            curve_offset += len(used)
            fiber_ids.append(i)
        result["fibers"].append(dict(**fiber_summary(group,refcurves[i],i,field,sample_data),
                                     seam_incidence=seam_counts(group,lookup,owners,patch,i)))
        faces = cells[patch == i+1]
        result["patches"].append(dict(patch=("AB","BC","CA")[i], triangles=len(faces),
            edge_connected_components=components(faces,2), vertex_connected_components=components(faces,1),
            samples=patch_metrics(field,p,faces,i,sample_data),
            **compare(p[faces],refsides[i],910+2*i,LO,HI)))
    save_polydata(dest/"corner_fibers.vtp", np.vstack(all_curves), indexed, fiber_ids, lines=True)
    sides = patch <= 3
    save_polydata(dest/"side_patches.vtp", p, cells[sides], patch[sides])
    if len(open_edges):
        save_polydata(dest/"open_edges.vtp", p, open_edges, [int(patch[owners[tuple(e)][0]]) for e in open_edges],
                      lines=True, scalar_name="incident_patch")
    np.savez_compressed(dest/"geometry_samples.npz", **sample_data)
    (dest/"geometry.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--trace", type=Path)
    parser.add_argument("--reference", type=Path, default=WAKE/"ttk/grid0025")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trace-only", action="store_true", help="Measure the traces in a seed-density trial")
    args = parser.parse_args()
    torch.set_num_threads(1)
    trace = args.trace or args.folder/"trace"
    dest = args.output or args.folder
    if args.trace_only:
        field = WakeField()
        ref = read_reference(args.reference, True)
        rows = [fiber_summary(g,ref[i],i,field) for i,g in enumerate(original_curves(trace))]
        (dest/"trace_geometry.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n")
    else:
        evaluate(args.folder, trace, args.reference, dest)
