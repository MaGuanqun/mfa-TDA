"""Compare the two full-domain Hurricane meshes with the original simulation grid."""

import argparse
import json
from pathlib import Path

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator
from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader


def grid_selection(pressure, temperature, triangle):
    """Membership in all three closed half-planes of a CCW range triangle."""
    p, t = pressure.astype(np.float64), temperature.astype(np.float64)
    selected = np.ones(p.shape, dtype=bool)
    for (a, b), (c, d) in zip(triangle, np.roll(triangle, -1, axis=0)):
        selected &= (c-a)*(t-b) >= (d-b)*(p-a)
    return selected


def mesh_locator(folder, report):
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(str(folder / "volume_volume.vtk"))
    reader.Update()
    mesh = reader.GetOutput()
    if mesh.GetNumberOfCells() != report["tetrahedra"]:
        raise ValueError(f"Cell count differs from the extraction report in {folder}")
    locator = vtkStaticCellLocator()
    locator.SetDataSet(mesh)
    locator.BuildLocator()
    return locator


def measure_grid(folder, report, pressure, temperature):
    if report["box_min"] != [0, 0, 0] or report["box_max"] != [499, 499, 99]:
        raise ValueError(f"Expected a full-domain run in {folder}")
    triangle = np.asarray(report["triangle"])
    if not report["complete"]:
        locator = None
    else:
        locator = mesh_locator(folder, report)
        bounds = np.asarray(locator.GetDataSet().GetBounds()).reshape(3, 2)
        xlo, ylo, zlo = np.ceil(bounds[:, 0]).astype(int)
        xhi, yhi, zhi = np.floor(bounds[:, 1]).astype(int)

    # All points outside the mesh bounds have known mesh membership. Evaluate
    # source-field membership on every grid point, including those outside.
    footprint = np.zeros(pressure.shape[1:], dtype=int)
    layers = []
    both = missed = excess = 0
    for z in range(pressure.shape[0]):
        selected = grid_selection(pressure[z], temperature[z], triangle)
        footprint += selected
        layers.append(int(selected.sum()))
        if locator is None:
            continue
        in_mesh = np.zeros_like(selected)
        if zlo <= z <= zhi:
            for y in range(ylo, yhi + 1):
                for x in range(xlo, xhi + 1):
                    in_mesh[y, x] = locator.FindCell((x, y, z)) >= 0
        both += int(np.count_nonzero(selected & in_mesh))
        missed += int(np.count_nonzero(selected & ~in_mesh))
        excess += int(np.count_nonzero(~selected & in_mesh))
    stats = {"triangle": triangle.tolist(), "source_grid_points": pressure.size, "selected_source_points": sum(layers),
             "selected_xy_columns": int(np.count_nonzero(footprint)),
             "total_xy_columns": footprint.size, "selected_by_z": layers}
    if triangle[:2].tolist() == [[-1100, -27], [-1100, -22]] and triangle[2, 1] == -22:
        stats["third_pressure"] = float(triangle[2, 0])
    if locator is not None:
        union = both + missed + excess
        stats.update(both=both, source_only=missed, mesh_only=excess,
                     neither=pressure.size - union, iou=both / union if union else None)
    return stats, footprint


def compare(folder, third_pressure, pressure, temperature):
    report = json.loads((folder / "volume_report.json").read_text())
    triangle = [[-1100, -27], [-1100, -22], [third_pressure, -22]]
    if not report["complete"] or report["triangle"] != triangle:
        raise ValueError(f"Expected a completed mesh for {triangle} in {folder}")
    stats, _ = measure_grid(folder, report, pressure, temperature)
    return {"third_pressure": third_pressure, **{key: stats[key] for key in
            ("source_grid_points", "both", "source_only", "mesh_only", "neither", "iou")}}


if __name__ == "__main__":
    build = Path(__file__).resolve().parents[2] / "build"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_dir", type=Path, nargs="?",
                        default=build / "fiber-volume-results/hurricane_isabel_3d_full")
    args = parser.parse_args()
    data = build / "examples/hurricane_isabel"
    pressure, temperature = (np.memmap(data / name, mode="r", dtype="<f4", shape=(100, 500, 500))
                             for name in ("Pf30_3d.bin", "TCf30_3d.bin"))
    results = [compare(args.mesh_dir / str(-p), p, pressure, temperature) for p in (-1200, -1300)]
    print(json.dumps(results, indent=2))
