"""Assemble the paper's supplemental table from the saved experiment reports."""
import csv
import json
from pathlib import Path
import re

import numpy as np
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy

from join_isabel_refinements import read_ply
from ttk_geometry import compare, original_curves, sizes

ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "build"
EXPERIMENT = Path(__file__).resolve().parent


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def load(path):
    return json.loads(path.read_text())


def vtk_sides(path):
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    mesh = reader.GetOutput()
    points = vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)
    faces = vtk_to_numpy(mesh.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    labels = vtk_to_numpy(mesh.GetCellData().GetArray("patch"))
    return [points[faces[labels==i]] for i in (1,2,3)]


def main():
    base = EXPERIMENT / "model_metrics.csv"
    with base.open(newline="") as stream:
        head, *records = csv.reader(stream)
    cases = head[4:]
    table = {tuple(row[:4]): dict(zip(cases, row[4:])) for row in records}
    if len(table) != len(records):
        raise ValueError("Duplicate metric identifiers in starting table")

    def add(case, key, value):
        if case not in cases:
            cases.append(case)
        if value is not None and value != "":
            table.setdefault(key, {})[case] = value

    def add_report(case, metrics, group, path, exclude=()):
        add(case, ("Provenance","Measurement file",group,"repository-relative path"), str(path.relative_to(ROOT)))
        for key, value in metrics.items():
            if key in exclude:
                continue
            if key in aliases:
                add(case, aliases[key], value)
            else:
                add(case, (group, key.replace("_"," "), "", ""), value)

    def results(key):
        return ("Results",key,"","")

    topology = lambda name: ("Topology",name,"Candidate boundary","count")
    aliases = {}
    for names, key in [
            (("triangles",), results("Boundary triangles")),
            (("open_edges",), results("Open boundary edges")),
            (("nonmanifold_vertices",), topology("Nonmanifold vertices")),
            (("intersections","intersecting_pairs"), topology("Intersections")),
            (("source_iou","integer_node_surface_source_iou"), results("Mesh vs original simulation IoU")),
            (("fitted_iou","integer_node_surface_fitted_iou"), results("Mesh vs fitted model IoU at grid points")),
            (("shells","edge_connected_shells"), ("Surface","Edge-connected components","including open components","count")),
            (("projection_p95","projection_p95_20000_samples"), ("Surface","Projection displacement, P95","20,000 area-weighted samples","grid units")),
            (("surface_area","area"), ("Surface","Total side area","Isabel","square grid units")),
            (("seam_failures","seam_test_failures","corner_failures"), ("Surface","Corner vertex tests failed","","count")),
            (("seam_tests","corner_vertices_checked"), ("Surface","Corner vertices tested","","count"))]:
        aliases.update(dict.fromkeys(names, key))

    normal_study = BUILD/"isabel_normal_step/comparison.csv"
    for r in rows(normal_study):
        angle = (", " + r["normal_limit_deg"] + " deg") if r["normal_limit_deg"] else ""
        name = f"Isabel outer ring, {r['domain']}, {r['method']}, spacing {r['step']}{angle}"
        if (r["domain"],r["method"],r["step"],r["normal_limit_deg"]) == ("whole","marcher","1","55"):
            name = "Isabel outer ring (500k, 1h verified arc recovery)"
        add_report(name, r, "Isabel step and grid study", normal_study, ("domain","method"))

    half = BUILD/"isabel_half_study"
    same_run = {"Uniform h=0.5, assembled":"Isabel outer ring, Uniform h=0.5",
                "Both native h=0.125 refinements joined":"Isabel outer ring, h=0.5, two sites re-extracted at h=0.125",
                "Uniform h=0.8, previous native bounded repair":
                    "Isabel outer ring, whole, marcher plus bounded repair, spacing 0.8, 85 deg"}
    for r in rows(half/"comparison.csv"):
        add_report(same_run.get(r["stage"],"Isabel outer ring, "+r["stage"]), r,
                   "Isabel local re-extraction", half/"comparison.csv", ("stage",))
    for r in rows(half/"local_reference.csv"):
        name = f"Isabel local reference, x={r['site']}, spacing {r['spacing']}"
        add_report(name, r, "Isabel local reference", half/"local_reference.csv", ("site","spacing"))

    two = BUILD/"isabel_two_methods"
    names = dict(h05_raw="Uniform h=0.5", h025_raw="Uniform h=0.25",
                 h05_final="h=0.5, reconnect only (4h)", h05_fitted="h=0.5, reconnect and fit (4h)",
                 h05_local_h0125="h=0.5, two sites re-extracted at h=0.125")
    for r in rows(two/"comparison.csv"):
        case = "Isabel outer ring, "+names[r["method"]]
        add_report(case,r,"Isabel final comparison",two/"comparison.csv", ("method",))
    for name, path in [(names["h025_raw"],two/"quarter/audit.json"),
                       (names["h05_fitted"],two/"half_fitted/final_audit.json"),
                       ("h=0.5, reconnect and fit (2h)",two/"half_strict/final_audit.json")]:
        audit = load(path)
        case = "Isabel outer ring, "+name
        add_report(case,audit,"Whole-mesh audit",path,("topology",))
        for field, value in re.findall(r"(\w+)=(-?\d+)", audit["topology"]):
            if field=="boundary_edges":
                add(case,results("Open boundary edges"),value)
            else:
                label = "Euler" if field=="Euler" else field.replace("_"," ").capitalize()
                key = ("Topology",label,"Candidate boundary","Euler characteristic") if label=="Euler" else topology(label)
                add(case,key,value)

    ttk_path = next((BUILD/"fiber-volume-results").glob("*/trace_ablation/arc_recovery/ttk/reference.json"))
    ttk = load(ttk_path)
    case = "Isabel outer ring, TTK spacing 1"
    add_report(case,{"triangles":sum(ttk["side_triangles"]), **ttk["welded_side_topology"]}, "TTK side extraction",ttk_path)
    for label,n in zip(("AB","BC","CA"),ttk["side_triangles"]):
        add(case,("TTK side extraction","Triangles",label,"count"),n)
    case = "Isabel outer ring, "+names["h05_fitted"]
    path = two/"half_fitted/final.ply"
    x,f,labels,_ = read_ply(path)
    counts = original_curves(next(half.glob("hurricane*/stages")))
    for corner,curves in zip("ABC",counts):
        add(case,("Corner traces","Closed traces",corner,"count"),sum(closed for _,closed,_ in curves))
    reference = vtk_sides(ttk_path.with_name("sides.vtp"))
    for i,patch in enumerate(("AB","BC","CA")):
        measures = compare(x[f[labels==i+1]],reference[i],910+2*i,np.zeros(3),np.array([499,499,99]))
        for direction in ("extracted_to_ttk","ttk_to_extracted"):
            for metric in ("median","p95","sampled_max"):
                add(case,("Surface vs TTK",direction.replace("_"," ")+" "+metric,patch,
                          "unit cube, 10,000 area-weighted samples"),measures[direction][metric])
    add(case,("Provenance","Surface mesh","Repaired surface vs TTK","repository-relative path"),str(path.relative_to(ROOT)))

    wake = BUILD/"fiber-volume-results/vortex_street_up"
    domain = load(wake/"ttk/membership.json")["regions"]["domain"]
    reference_area = float(rows(wake/"step_comparison.csv")[0]["ttk_physical_side_area"])
    for grid, h in [("grid005",0.05),("grid0025",0.025)]:
        case = f"Upward wake, TTK spacing {h}"
        path = wake/"ttk"/grid/"reference.json"
        ref = load(path)
        add_report(case,{"triangles":sum(ref["side_triangles"]),"implicit_tetrahedra":ref["implicit_tetrahedra"]},
                   "TTK side extraction",path)
        for label,n in zip(("AB","BC","CA"),ref["side_triangles"]):
            add(case,("TTK side extraction","Triangles",label,"count"),n)
        area = sum(sizes(side).sum() for side in vtk_sides(path.with_name("sides.vtp")))
        add(case,("TTK side extraction","Side area","","model domain units"),area)
        add(case,("TTK side extraction","Side area / finer TTK area","","ratio"),area/reference_area)
        # This membership test evaluates the PL interpolant, independently of the exported triangles.
        for k,v in domain["pl_"+grid].items():
            add(case,("Field membership",k,"INR vs piecewise-linear field","100,000 uniform samples"),v)
    path = BUILD/"full-domain-experiments/wake-samples.json"
    for r in load(path)["results"]:
        for which in ("up","down"):
            measure = "volume" if r["time"] is None else "area"
            where = "whole domain" if r["time"] is None else f"t={r['time']}"
            add("Wake "+which+", direct field sampling",("Field membership",measure,where,"model domain units"),
                r[which+"_"+measure])
    path = wake/"h0025_100k_retry/fiber_validation.json"
    corners = load(path)["corners"]
    add("Vortex up 100k triangle retry",("Geometry totals","Completed traces","",""),sum(c["complete"] for c in corners))
    down = original_curves(BUILD/"full-domain-experiments/wake-down/corners")
    add("Vortex down h=0.05",("Geometry totals","Completed traces","",""),sum(c[2] for group in down for c in group))

    dest = ROOT/"paper/data/model-results.csv"
    dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open("w",newline="") as stream:
        w = csv.writer(stream, lineterminator="\n")
        w.writerow(head[:4]+cases)
        for key,values in table.items():
            w.writerow([*key, *(values.get(case,"") for case in cases)])


if __name__=="__main__":
    main()
