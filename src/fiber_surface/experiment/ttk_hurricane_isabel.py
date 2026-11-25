"""Compare Isabel MFA boundary patches and corner seams against a TTK grid reference."""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from evaluate_curved_prism import incidence, read_mesh, save_polydata
from ttk_geometry import compare, extract_triangle, original_curves, sizes
from evaluate_hurricane_isabel import compare_original_grid, plot_comparison_section, plot_ttk_footprints


FIELDS = (("Pf30", "pressure"), ("TCf30", "temperature"))
NAMES = tuple(name for _, name in FIELDS)
DIMS = (500, 500, 100)
LO = np.zeros(3)
HI = np.asarray(DIMS) - 1


def sampled_image(folder):
    image = vtk.vtkImageData()
    image.SetDimensions(*DIMS)
    image.SetOrigin(*LO)
    image.SetSpacing(1, 1, 1)
    image.mapped_fields = []  # Keep the zero-copy buffers alive as long as the VTK image.
    models = {}
    for stem, name in FIELDS:
        meta = json.loads((folder / f"{stem}.json").read_text())
        if (meta["backend"] != "mfa" or meta["dataset"] != "hurricane_isabel_3d"
                or meta["dimensions"] != list(DIMS) or meta["sample_min"] != LO.tolist()
                or meta["sample_max"] != HI.tolist() or meta["spacing"] != [1, 1, 1]
                or meta["layout"] != "t,y,x" or meta["dtype"] != "<f8"):
            raise ValueError(f"Unexpected MFA sample grid in {stem}.json")
        data = np.memmap(folder / f"{stem}.bin", mode="r", dtype=meta["dtype"])
        if data.size != image.GetNumberOfPoints():
            raise ValueError(f"Wrong sample count in {stem}.bin")
        # The sampler writes x fastest. VTK can use the mapped array directly.
        array = numpy_to_vtk(data)
        array.SetName(name)
        image.GetPointData().AddArray(array)
        image.mapped_fields.append(data)
        models[name] = meta["model"]
    return image, models


def compare_group(extracted, reference, seed):
    result = dict(mfa_cells=len(extracted), ttk_cells=len(reference))
    if len(extracted) and len(reference):
        result.update(compare(extracted, reference, seed, LO, HI))
    else:
        result.update(physical_measure=float(sizes(extracted).sum()),
                      ttk_physical_measure=float(sizes(reference).sum()),
                      extracted_to_ttk=None, ttk_to_extracted=None)
    return result


def edge_incidence(faces):
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    _, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    return dict(open_edges=int(np.count_nonzero(counts == 1)),
                nonmanifold_edges=int(np.count_nonzero(counts > 2)))


def compare_corner_traces(run, dest, reference):
    stages = run / "stages"
    if not (stages / "run.conf").exists():
        return
    xyz, _, _, edges, corners = reference
    results = {}
    for i, curves in enumerate(original_curves(stages)):
        segments = [np.stack((c[:-1], c[1:]), axis=1) for c, _, _ in curves if len(c) > 1]
        lines = np.concatenate(segments) if segments else np.empty((0, 2, 3))
        result = compare_group(lines, xyz[edges[corners == i]], 810 + 2*i)
        result.update(traces=len(curves), closed=sum(closed for _, closed, _ in curves),
                      completed=sum(done for _, _, done in curves))
        results["ABC"[i]] = result
    (dest / "corner_traces.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")


def write_reference(dest, image, models, query, threads):
    print(f"Running ttkFiberSurface on {DIMS}, triangle {query.tolist()}", flush=True)
    start = perf_counter()
    xyz, faces, patches, edges, corners, collapsed = extract_triangle(
        image, NAMES, query, threads, range_octree=False)
    seconds = perf_counter() - start
    dest.mkdir(parents=True, exist_ok=True)
    save_polydata(dest / "sides.vtp", xyz, faces, patches)
    if len(edges):
        save_polydata(dest / "fibers.vtp", xyz, edges, corners, lines=True)

    metrics = dict(models=models, triangle=query.tolist(),
                   dimensions=list(DIMS), origin=LO.tolist(), spacing=list(image.GetSpacing()),
                   grid_points=image.GetNumberOfPoints(), implicit_tetrahedra=int(6 * np.prod(HI)),
                   vtk_version=vtk.vtkVersion.GetVTKVersion(), threads=threads, range_octree=False,
                   extraction_and_weld_seconds=seconds, collapsed_triangles=collapsed,
                   welded_side_topology=edge_incidence(faces),
                   distance_units="Euclidean distance after mapping the full domain to the unit cube",
                   side_triangles=np.bincount(patches, minlength=4)[1:4].tolist(),
                   corner_segments=np.bincount(corners, minlength=3).tolist())
    (dest / "reference.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    plot_ttk_footprints(dest, xyz, faces, patches, DIMS)
    print(f"{dest}: {len(faces):,} TTK triangles and {len(edges):,} corner segments", flush=True)
    return (xyz, faces, patches, edges, corners), metrics


def compare_boundary(run, dest, reference, reference_metrics, image, section_z=None):
    report = json.loads((run / "volume_report.json").read_text())
    if (report["box_min"] != LO.tolist() or report["box_max"] != HI.tolist()
            or report["triangle"] != reference_metrics["triangle"]):
        raise ValueError(f"MFA domain or query differs from TTK in {run}")
    filled = report["complete"] and report["boundary_valid"]
    path = run / ("volume_tet_boundary.vtk" if filled else "volume_boundary.vtk")
    if filled or report["boundary_triangles"]:
        mfa_xyz, mfa_faces, mfa_patches = read_mesh(path)
    else:
        mfa_xyz = np.empty((0, 3))
        mfa_faces = np.empty((0, 3), dtype=int)
        mfa_patches = np.empty(0, dtype=int)
    corner_patches = ((1, 3), (1, 2), (2, 3))
    mfa_seams = {pair: [] for pair in corner_patches}
    for edge, owners in incidence(mfa_faces, 2).items():
        pair = tuple(sorted(mfa_patches[owners]))
        if pair in mfa_seams:
            mfa_seams[pair].append(edge)
    xyz, faces, patches, edges, corners = reference
    metrics = dict(reference_metrics, mfa_boundary=str(path.resolve()), mfa_complete=filled,
                   mfa_boundary_valid=report["boundary_valid"], mfa_failed_stage=report["failed_stage"],
                   patches={}, corner_seams={})
    for i, (name, pair) in enumerate(zip(("AB", "BC", "CA"), corner_patches)):
        metrics["patches"][name] = compare_group(mfa_xyz[mfa_faces[mfa_patches == i+1]],
                                                 xyz[faces[patches == i+1]], 910 + 2*i)
        seam = np.asarray(mfa_seams[pair], dtype=int).reshape(-1, 2)
        metrics["corner_seams"]["ABC"[i]] = compare_group(mfa_xyz[seam], xyz[edges[corners == i]], 810 + 2*i)
    metrics["section_z"] = plot_comparison_section(
        dest, path, filled, [a.reshape(DIMS[::-1]) for a in image.mapped_fields], report["triangle"], section_z)
    (dest / "comparison.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    print(f"Compared TTK with {'filled' if filled else 'unfilled'} MFA boundary in {run}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", type=Path, help="folder with the two sampler .bin/.json pairs")
    parser.add_argument("runs", nargs="+", type=Path, help="Isabel run directories with a volume report")
    parser.add_argument("--threads", type=int, default=3)
    parser.add_argument("--source-grid", type=Path, help="time-30 simulation grids for fitted/source IoU")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    runs = [run.resolve() for run in args.runs]
    if len({run.name for run in runs}) != len(runs):
        parser.error("Run directory names must be distinct within one comparison")
    image, models = sampled_image(args.samples)
    for run in runs:
        dest = args.samples / run.name
        report = json.loads((run / "volume_report.json").read_text())
        query = np.asarray(report["triangle"], dtype=float)
        reference, metrics = write_reference(dest, image, models, query, args.threads)
        compare_corner_traces(run, dest, reference)
        section_z = None
        if args.source_grid:
            originals = [np.memmap(args.source_grid / f"{stem}_3d.bin", mode="r", dtype="<f4", shape=DIMS[::-1])
                         for stem, _ in FIELDS]
            fitted = [data.reshape(DIMS[::-1]) for data in image.mapped_fields]
            agreement = compare_original_grid(originals, fitted, query)
            section_z = agreement["fitted_peak_z"]
            (dest / "fit_vs_source.json").write_text(json.dumps(agreement, indent=2) + "\n")
        compare_boundary(run, dest, reference, metrics, image, section_z)


if __name__ == "__main__":
    main()
