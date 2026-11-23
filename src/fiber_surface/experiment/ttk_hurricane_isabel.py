"""Compare Isabel MFA boundary patches and corner seams against a TTK grid reference."""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

from evaluate_curved_prism import incidence, read_mesh, save_polydata
from ttk_geometry import compare, extract_triangle, sizes


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


def compare_run(run, dest, image, models, threads):
    report = json.loads((run / "volume_report.json").read_text())
    if (not report["complete"] or not report["boundary_valid"] or report["box_min"] != LO.tolist()
            or report["box_max"] != HI.tolist()):
        raise ValueError(f"Expected a completed full-domain Isabel volume in {run}")
    query = np.asarray(report["triangle"], dtype=float)
    path = run / "volume_tet_boundary.vtk"
    mfa_xyz, mfa_faces, mfa_patches = read_mesh(path)
    # The seam for A, B or C is shared by the two side patches incident on it.
    corner_patches = ((1, 3), (1, 2), (2, 3))  # A, B, C
    mfa_seams = {pair: [] for pair in corner_patches}
    for edge, faces in incidence(mfa_faces, 2).items():
        pair = tuple(sorted(mfa_patches[faces]))
        if pair in mfa_seams:
            mfa_seams[pair].append(edge)

    print(f"Running ttkFiberSurface on {DIMS} for {run.name}, triangle {query.tolist()}", flush=True)
    start = perf_counter()
    xyz, faces, patches, edges, corners, collapsed = extract_triangle(
        image, NAMES, query, threads, range_octree=False)
    seconds = perf_counter() - start
    dest.mkdir(parents=True, exist_ok=True)
    save_polydata(dest / "sides.vtp", xyz, faces, patches)
    if len(edges):
        save_polydata(dest / "fibers.vtp", xyz, edges, corners, lines=True)

    metrics = dict(mfa_boundary=str(path.resolve()), models=models, triangle=query.tolist(),
                   dimensions=list(DIMS), origin=LO.tolist(), spacing=list(image.GetSpacing()),
                   grid_points=image.GetNumberOfPoints(), implicit_tetrahedra=int(6 * np.prod(HI)),
                   vtk_version=vtk.vtkVersion.GetVTKVersion(), threads=threads, range_octree=False,
                   extraction_and_weld_seconds=seconds, collapsed_triangles=collapsed,
                   distance_units="Euclidean distance after mapping the full domain to the unit cube",
                   patches={}, corner_seams={})
    for i, (name, pair) in enumerate(zip(("AB", "BC", "CA"), corner_patches)):
        metrics["patches"][name] = compare_group(mfa_xyz[mfa_faces[mfa_patches == i+1]],
                                                 xyz[faces[patches == i+1]], 910 + 2*i)
        seam = np.asarray(mfa_seams[pair], dtype=int).reshape(-1, 2)
        metrics["corner_seams"]["ABC"[i]] = compare_group(mfa_xyz[seam], xyz[edges[corners == i]], 810 + 2*i)
    (dest / "comparison.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    print(f"{dest}: {len(faces):,} TTK triangles and {len(edges):,} corner segments", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples", type=Path, help="folder with the two sampler .bin/.json pairs")
    parser.add_argument("runs", nargs="+", type=Path, help="completed Isabel MFA run directories")
    parser.add_argument("--threads", type=int, default=3)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    runs = [run.resolve() for run in args.runs]
    if len({run.name for run in runs}) != len(runs):
        parser.error("Run directory names must be distinct within one comparison")
    image, models = sampled_image(args.samples)
    for run in runs:
        compare_run(run, args.samples / run.name, image, models, args.threads)


if __name__ == "__main__":
    main()
