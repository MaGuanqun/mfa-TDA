#!/usr/bin/env bash
# Filter the final colored vortex-street INR trajectories for all three regions.
set -eo pipefail
if [[ $# -gt 1 || ${1:-} == --help || ${1:-} == -h ]]; then
    printf 'Usage: bash critical_point_tracking_local_remove_obstacle.sh [RESULTS_ROOT]\n'
    printf 'Default: build/local_tracking/vortex_street_3d_INR/ratio_80\n'
    exit 0
fi
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
region2_reference="$repo_root/build/src/vortex_street_3d_INR/vortex_street_3d_8_ttk_tracking_remove_obstacle_colorid.vtp"
output_dir=${1:-"$repo_root/build/local_tracking/vortex_street_3d_INR/ratio_80"}
output_dir=$(realpath -e -- "$output_dir")
pvpython_bin=${PVPYTHON:-pvpython}
if [[ $pvpython_bin == */* ]]; then pvpython_bin=$(realpath -e -- "$pvpython_bin"); fi
command -v "$pvpython_bin" >/dev/null
cd "$repo_root"
test -s "$region2_reference"
for region in 1 2 3; do
    prefix="$output_dir/region_$region/inr/samples_tracking_all_critical_points"
    test -s "${prefix}_regionid_colorid.vtp"
    test -s "${prefix}_points.csv"
done
for region in 1 2 3; do
    folder="$output_dir/region_$region/inr"
    prefix="$folder/samples_tracking_all_critical_points"
    command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
             "$repo_root/src/python/remove_obstacle.py" --data vortex_street_3d
             -i "${prefix}_regionid_colorid.vtp"
             -o "${prefix}_regionid_colorid_remove_obstacle.vtp"
             -c "${prefix}_points.csv" -d "${prefix}_points_remove_obstacle.csv")
    printf '%q ' "${command[@]}" > "$folder/remove_obstacle_commands.txt"
    printf '\n' >> "$folder/remove_obstacle_commands.txt"
    "${command[@]}" 2>&1 | tee "$folder/remove_obstacle.log"
    if [[ $region -eq 2 ]]; then
        # Match the surviving paths, since removed portions can invalidate
        # the original reference match. Preserve geometry and fallback colors.
        command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
                 "$repo_root/src/python/recolor_filtered_trajectories.py"
                 --input-vtp "${prefix}_regionid_colorid_remove_obstacle.vtp"
                 --output-vtp "${prefix}_regionid_colorid_remove_obstacle.vtp"
                 --output-csv "${prefix}_points_remove_obstacle.csv"
                 --tracking-json "${prefix}.json"
                 --reference-vtp "$region2_reference"
                 --matches-csv "${prefix}_remove_obstacle_color_matches.csv"
                 --report-json "${prefix}_remove_obstacle_colors.json"
                 --color-override 342 39)
        printf '%q ' "${command[@]}" >> "$folder/remove_obstacle_commands.txt"
        printf '\n' >> "$folder/remove_obstacle_commands.txt"
        "${command[@]}" 2>&1 | tee -a "$folder/remove_obstacle.log"
    fi
done

"$pvpython_bin" --no-mpi --force-offscreen-rendering - "$output_dir" <<'PY'
import csv
from pathlib import Path
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()

root = Path(sys.argv[1])
rows = []
for region in (1, 2, 3):
    folder = root / f'region_{region}' / 'inr'
    stem = 'samples_tracking_all_critical_points'
    source = read(folder / (stem + '_regionid_colorid.vtp'))
    vtp = folder / (stem + '_regionid_colorid_remove_obstacle.vtp')
    result = read(vtp)
    xyz = vtk_to_numpy(source.GetPoints().GetData())
    outside = np.sum(xyz[:, :2]**2, axis=1) > .12**2
    kept_xyz = vtk_to_numpy(result.GetPoints().GetData())
    np.testing.assert_array_equal(kept_xyz, xyz[outside])
    for i in range(source.GetPointData().GetNumberOfArrays()):
        before = source.GetPointData().GetArray(i)
        name = before.GetName()
        if region == 2 and name in ('ColorId', 'ReferenceColorMatched'):
            name += 'BeforeObstacleRematch'
        after = result.GetPointData().GetArray(name)
        np.testing.assert_array_equal(vtk_to_numpy(after), vtk_to_numpy(before)[outside])
    before_lines = vtk_to_numpy(source.GetLines().GetConnectivityArray()).reshape(-1, 2)
    after_lines = vtk_to_numpy(result.GetLines().GetConnectivityArray()).reshape(-1, 2)
    keep_lines = outside[before_lines].all(axis=1)
    np.testing.assert_array_equal(kept_xyz[after_lines], xyz[before_lines[keep_lines]])
    for i in range(source.GetCellData().GetNumberOfArrays()):
        before = source.GetCellData().GetArray(i)
        name = before.GetName()
        if region == 2 and name in ('ColorId', 'ReferenceColorMatched'):
            name += 'BeforeObstacleRematch'
        after = result.GetCellData().GetArray(name)
        np.testing.assert_array_equal(vtk_to_numpy(after)[result.GetNumberOfVerts():],
                                      vtk_to_numpy(before)[source.GetNumberOfVerts():][keep_lines])
    used = set()
    for i in range(result.GetNumberOfCells()):
        cell = result.GetCell(i)
        ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
        used.update(ids)
        for name in ('ColorId', 'RegionId', 'ReferenceColorMatched'):
            point_values = vtk_to_numpy(result.GetPointData().GetArray(name))
            cell_value = result.GetCellData().GetArray(name).GetTuple1(i)
            assert np.all(point_values[ids] == cell_value), (region, name, i)
    assert used == set(range(result.GetNumberOfPoints()))
    csv_path = folder / (stem + '_points_remove_obstacle.csv')
    with csv_path.open() as stream:
        csv_rows = list(csv.DictReader(stream))
    csv_xyz = np.array([[float(row[k]) for k in ('PositionX', 'PositionY', 'PositionZ')]
                        for row in csv_rows]).reshape(-1, 3)
    np.testing.assert_array_equal(csv_xyz, kept_xyz)
    for name in ('ColorId', 'RegionId', 'ReferenceColorMatched'):
        np.testing.assert_array_equal([int(row[name]) for row in csv_rows],
                                      vtk_to_numpy(result.GetPointData().GetArray(name)))
    if region == 2:
        import json
        settings = json.loads((folder / (stem + '.json')).read_text())
        color_report = json.loads((folder / (stem + '_remove_obstacle_colors.json')).read_text())
        reference = read(color_report['reference'])
        ref_ids = vtk_to_numpy(reference.GetPointData().GetArray('RegionId'))
        ref_colors = vtk_to_numpy(reference.GetPointData().GetArray('ColorId'))
        with (folder / (stem + '_remove_obstacle_color_matches.csv')).open() as stream:
            matches = list(csv.DictReader(stream))
        for data in (result.GetPointData(), result.GetCellData()):
            components = vtk_to_numpy(data.GetArray('ColorMatchComponentId'))
            colors = vtk_to_numpy(data.GetArray('ColorId'))
            automatic_colors = vtk_to_numpy(data.GetArray('AutomaticColorId'))
            overrides = vtk_to_numpy(data.GetArray('ManualColorOverride'))
            flags = vtk_to_numpy(data.GetArray('ReferenceColorMatched'))
            for row in matches:
                mask = components == int(row['ComponentId'])
                assert np.all(colors[mask] == int(row['ColorId']))
                assert np.all(automatic_colors[mask] == int(row['AutomaticColorId']))
                assert np.all(overrides[mask] == int(row['ManualColorOverride']))
                assert np.all(flags[mask] == int(row['MatchedReference']))
                if int(row['ManualColorOverride']):
                    assert int(row['OriginalRegionId']) == 342 and int(row['ColorId']) == 39
                else:
                    assert int(row['ColorId']) == int(row['AutomaticColorId'])
                if int(row['MatchedReference']):
                    assert float(row['Distance']) <= settings['color_match_threshold']
                    expected = ref_colors[ref_ids == int(row['ReferenceRegionId'])]
                    assert len(expected) and np.all(expected == int(row['AutomaticColorId']))
                else:
                    assert int(row['AutomaticColorId']) == int(row['PreviousColorId'])
    rows.append([region, source.GetNumberOfPoints(), int((~outside).sum()),
                 result.GetNumberOfPoints(), result.GetNumberOfLines(),
                 result.GetNumberOfVerts(), str(vtp), str(csv_path)])
with (root / 'remove_obstacle_summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['region', 'input_points', 'removed_points', 'kept_points',
                     'kept_segments', 'singleton_vertices', 'colored_vtp', 'points_csv'])
    writer.writerows(rows)
for row in rows:
    print(f'Region {row[0]}: removed {row[2]} obstacle points; kept {row[3]} points, '
          f'{row[4]} segments and {row[5]} singleton vertices; attributes/colors verified.')
print(f'Summary: {root / "remove_obstacle_summary.csv"}')
PY
