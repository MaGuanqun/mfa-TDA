#!/usr/bin/env bash
# Local MFA/CoordNet sampling -> critical-point trajectories -> reference colors.
set -eo pipefail

# EDIT THESE SETTINGS, then run: bash critical_point_tracking_local.sh
# z is physical time. Each row is: XMIN XMAX YMIN YMAX ZMIN ZMAX.
dataset="boussinesq_3d_INR"  # Also accepts boussinesq_3d or vortex_street_3d for MFA.
backend="auto"              # _INR selects inr; otherwise mfa. Override with --backend.
points_per_span=160          # Sampling density: spacing = span_width / this number.
inr_sampling_ratio=80        # Intervals per INR block axis; 0.1/80 = 0.00125.
tracker="auto"              # INR: all critical points (including saddles); MFA: legacy TTK.
threads=3
vortex_street_regions=(
    "0.5 1.0 -0.5 0.5 14.5 15.0"
    "0.36 0.76 -0.1 0.42 13.6 14.0"
    "0.46 0.51 -0.1 0.1 14.3 14.6"
)
boussinesq_regions=(
    "-0.05 0.10 0.05 0.22 0.81 1.20"
    "0.25 0.30 1.16 1.26 0.28 0.66"
    "0.32 0.40 1.62 1.73 1.07 1.30"
)

usage() {
    cat <<'EOF'
Usage:
  bash critical_point_tracking_local.sh [options]

NAME: boussinesq_3d or vortex_street_3d (both are 2D space + time).
Append _INR to select CoordNet, e.g. boussinesq_3d_INR (the default).
All region bounds and padding are in ORIGINAL domain coordinates.

  --dataset NAME            Default: editable dataset setting at top of script
  --backend auto|mfa|inr|both Default: auto (_INR selects inr; otherwise mfa)
  --tracker auto|all|legacy  Default: auto (INR uses all; MFA/both uses legacy)
  --region X0 X1 Y0 Y1 T0 T1 Replace configured boxes; repeat for multiple boxes
  --points-per-span N       MFA span sampling (default: 160 intervals per span)
  --ratio R                 Physical spacing 0.1/R; INR defaults to 80 per block axis
  --ratios RX RY RT          Separate spatial/time ratios
  --padding PX PY PT        Expand the region in physical units (default: 0 0 0)
  --output-dir DIR          Default: build/local_tracking/NAME/points_per_span_N (MFA)
                            or build/local_tracking/NAME_INR/ratio_80 (INR)
                            Results go in DIR/region_N/BACKEND
  --mfa-model FILE          Default: build/src/NAME/NAME.mfa
  --inr-model FILE          Default: CoordNet/Exp/NAME/super-spatial-temporal-64-5-float64.pt
  --batch-size N            Model evaluation batch size (default: 65536)
  --threads N               Sampling and TTK thread limit (default: 3)
  --device cpu|cuda         CoordNet device (default: cpu)
  --pvpython PATH           Default: $PVPYTHON or pvpython on PATH
  --ttk-plugin PATH         Default: $TTK_PLUGIN or automatic plugin discovery
  --persistence P           TTK persistence threshold (default: 0)
  --spatial-weight W        Override TTK X/Y weights (default: original plugin defaults)
  --destruction-cost C      Optional; only supported by newer TTK builds
  --tetrahedralize          Original explicit triangulation (default)
  --implicit-grid           Opt into implicit triangulation; results can differ
  --skip-vti                Legacy tracker: skip optional sampled VTI copy
  --sample-only             Stop after writing samples.bin and samples.json
  --tracking-only           Reuse samples.bin/.json and rerun the selected tracker
  --track-all-critical-points  Reuse saved per-slice points, track all types and color
  --tracking-max-displacement D  Maximum XY motion per slice (default: spatial-step-size)
  --critical-points-only    Reuse samples.bin/.json; extract critical points at every time slice
  --critical-points-chunk-size N  Slices per critical-point batch (default: 8)
  --postprocess-only        Legacy: recolor saved trajectories; all: retrack/recolor saved points
  --skip-postprocess        Legacy tracker only: stop after TTK without coloring
  --reference-vtp FILE      Default: build/src/NAME/NAME_16_regionid_colorid.vtp
                            Vortex street uses NAME_16_remove_obstacle_regionid_colorid.vtp
                            Boussinesq INR: build/src/boussinesq_3d_INR/
                              boussinesq_3d_8_remove_obstacle_regionid_colorid.vtp
  --spatial-step-size H     Color X/Y scale (Boussinesq INR: 0.0125, MFA: 0.00694444;
                            vortex MFA: 0.00625)
  --t-step-size H           Color time scale (Boussinesq INR: 0.0125, MFA: 0.0133929;
                            vortex MFA: 0.00520833)
  --match-distance-threshold D  Mean nearest-point distance in step units (default: 1)
  --dry-run                 Print all grids (reads knots for MFA); no field evaluation
  --no-build                Reuse the existing sample_local_field executable
  -h, --help                Show this help

MFA uses each axis's actual uniform knot-span width/N.
Bounds are snapped outward to the same domain-anchored grid in every region,
so partial spans get proportional point counts with identical per-axis spacing.
Unequal/nonuniform spans are rejected in this mode. With --ratio/--ratios,
the physical-spacing mode is used, including the exact box endpoints.
CoordNet uses its configured block counts: Boussinesq (10,30,15), vortex (80,10,15).
The default INR ratio is 80 intervals per block axis; partial blocks receive
proportional counts. All three Boussinesq boxes align exactly at spacing 0.00125.
CoordNet input is normalized against the full dataset domain and ordered (t,y,x).
Padding is clipped to the dataset domain. Sampling uses bounded batches;
all-point extraction processes eight time slices per chunk by default.
Legacy trajectories retain TTK's native schema and need the local time series
in memory. Legacy postprocessing removes vertices on all six sampled-box faces (X/Y/Z, including
the first/last time slices) and their incident segments, joins exact endpoints,
then runs
seperate_diffferent_branches.py and seperate_branch_assign_color2.py with seed 0.
Matching uses (x/Hxy, y/Hxy, z/Ht), where z is physical time. Accepted branches
copy the reference ColorId; unmatched branches keep the matcher's seeded random
fallback (1..40). The legacy final file is samples_tracking_regionid_colorid.vtp.
Use the same ColorId lookup table/range as the reference in ParaView.
Critical-point references use the existing explicit TTK per-slice extractor.
samples_critical_points.csv/.vtp exclude XY boundary points, retaining every
time slice. The _interior versions also exclude the first/last time slices;
the _all versions retain all boundary points for inspection.
--tracker all runs sampling, per-slice extraction, per-type adjacent-frame assignment
and coloring in one invocation. --track-all-critical-points reuses saved detections.
Both use per-type adjacent-frame assignment on those
references, with no persistence/track-length filtering. Unlinked points remain
vertex cells. It writes samples_tracking_all_critical_points*.vtp and a coverage
report. The _regionid_colorid.vtp version excludes all six boundary faces.
The raw trajectory retains boundary points too. Color matching uses Hxy=Ht=0.0125
for the Boussinesq INR step-8 reference, independently of the sampling density.
Legacy TTK options (--persistence, --implicit-grid, etc.) require --tracker legacy.

Examples:
  bash critical_point_tracking_local.sh --dry-run
  bash critical_point_tracking_local.sh
  bash critical_point_tracking_local.sh --postprocess-only
  bash critical_point_tracking_local.sh --critical-points-only
  bash critical_point_tracking_local.sh --track-all-critical-points
  bash critical_point_tracking_local.sh --dataset vortex_street_3d
  bash critical_point_tracking_local.sh --dataset boussinesq_3d \
    --region -0.2 0.2 0.5 0.9 0.4 0.6 --output-dir build/local_tracking/my_box
EOF
}

fail() { printf '%s\n' "$*" >&2; exit 1; }
need() { [[ $remaining -ge $1 ]] || fail "Missing arguments after $option"; }

output_dir=""
mfa_model=""
inr_model=""
pvpython_bin="${PVPYTHON:-pvpython}"
build=1
sample_only=0
tracking_only=0
postprocess_only=0
critical_points_only=0
track_all_critical_points=0
tracking_max_displacement=""
critical_points_chunk_size=8
postprocess=1
reference_vtp=""
spatial_step_size=""
t_step_size=""
match_distance_threshold=1
dry_run=0
regions=()
sampling_mode=spans
explicit_spans=0
sampling_ratios=()
sample_args=()
tracking_args=()
original_args=("$@")
invocation_dir=$PWD

while [[ $# -gt 0 ]]; do
    option=$1
    remaining=$#
    case "$option" in
        --dataset) need 2; dataset=$2; shift 2 ;;
        --backend) need 2; backend=$2; shift 2 ;;
        --tracker) need 2; tracker=$2; shift 2 ;;
        --region) need 7; regions+=("$2 $3 $4 $5 $6 $7"); shift 7 ;;
        --points-per-span) need 2; points_per_span=$2; explicit_spans=1; shift 2 ;;
        --ratio) need 2; sampling_mode=ratio; sampling_ratios=("$2" "$2" "$2"); shift 2 ;;
        --ratios) need 4; sampling_mode=ratio; sampling_ratios=("$2" "$3" "$4"); shift 4 ;;
        --batch-size|--device) need 2; sample_args+=("$1" "$2"); shift 2 ;;
        --padding) need 4; sample_args+=("${@:1:4}"); shift 4 ;;
        --threads) need 2; threads=$2; shift 2 ;;
        --output-dir) need 2; output_dir=$2; shift 2 ;;
        --mfa-model) need 2; mfa_model=$2; shift 2 ;;
        --inr-model) need 2; inr_model=$2; shift 2 ;;
        --pvpython) need 2; pvpython_bin=$2; shift 2 ;;
        --ttk-plugin|--persistence|--destruction-cost|--spatial-weight) need 2; tracking_args+=("$1" "$2"); shift 2 ;;
        --tetrahedralize|--implicit-grid|--skip-vti) tracking_args+=("$1"); shift ;;
        --no-build) build=0; shift ;;
        --sample-only) sample_only=1; shift ;;
        --tracking-only) tracking_only=1; build=0; shift ;;
        --postprocess-only) postprocess_only=1; build=0; shift ;;
        --critical-points-only) critical_points_only=1; build=0; shift ;;
        --track-all-critical-points) track_all_critical_points=1; build=0; shift ;;
        --tracking-max-displacement) need 2; tracking_max_displacement=$2; shift 2 ;;
        --critical-points-chunk-size) need 2; critical_points_chunk_size=$2; shift 2 ;;
        --skip-postprocess) postprocess=0; shift ;;
        --reference-vtp) need 2; reference_vtp=$2; shift 2 ;;
        --spatial-step-size) need 2; spatial_step_size=$2; shift 2 ;;
        --t-step-size) need 2; t_step_size=$2; shift 2 ;;
        --match-distance-threshold) need 2; match_distance_threshold=$2; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) fail "Unknown option: $1 (see --help)" ;;
    esac
done
if [[ $dataset == *_INR ]]; then
    [[ $backend == auto || $backend == inr ]] || fail "A dataset ending in _INR requires --backend inr (or auto)"
    dataset=${dataset%_INR}
    backend=inr
fi
if [[ $backend == auto ]]; then backend=mfa; fi
if [[ $tracker == auto ]]; then
    tracker=legacy
    if [[ $backend == inr ]]; then tracker=all; fi
fi
[[ "$dataset" == boussinesq_3d || "$dataset" == vortex_street_3d ]] || fail "Choose --dataset boussinesq_3d or vortex_street_3d"
[[ "$backend" == mfa || "$backend" == inr || "$backend" == both ]] || fail "Choose --backend mfa, inr, or both"
[[ $tracker == all || $tracker == legacy ]] || fail "Choose --tracker all or legacy"
if [[ $tracker == all ]]; then
    [[ $postprocess -eq 1 ]] || fail "--skip-postprocess requires --tracker legacy; all-critical-point tracking includes coloring"
    for option in "${tracking_args[@]}"; do
        case "$option" in
            --persistence|--destruction-cost|--spatial-weight|--implicit-grid|--skip-vti)
                fail "$option requires --tracker legacy" ;;
        esac
    done
fi
if [[ ${#regions[@]} -eq 0 ]]; then
    if [[ "$dataset" == vortex_street_3d ]]; then
        regions=("${vortex_street_regions[@]}")
    else
        regions=("${boussinesq_regions[@]}")
    fi
fi
[[ ${#regions[@]} -gt 0 ]] || fail "Configure regions at the top of this script, or supply repeated --region arguments"
[[ $sample_only -eq 0 || $tracking_only -eq 0 ]] || fail "--sample-only and --tracking-only are mutually exclusive"
[[ $dry_run -eq 0 || $tracking_only -eq 0 ]] || fail "--dry-run previews sampling; omit --tracking-only"
if [[ $postprocess_only -eq 1 ]]; then
    [[ $sample_only -eq 0 && $tracking_only -eq 0 && $dry_run -eq 0 && $postprocess -eq 1 ]] || \
        fail "--postprocess-only cannot be combined with --sample-only, --tracking-only, --dry-run or --skip-postprocess"
fi
if [[ $critical_points_only -eq 1 ]]; then
    [[ $sample_only -eq 0 && $tracking_only -eq 0 && $postprocess_only -eq 0 && $dry_run -eq 0 ]] || \
        fail "--critical-points-only cannot be combined with --sample-only, --tracking-only, --postprocess-only or --dry-run"
fi
if [[ $track_all_critical_points -eq 1 ]]; then
    [[ $sample_only -eq 0 && $tracking_only -eq 0 && $postprocess_only -eq 0 && \
       $critical_points_only -eq 0 && $dry_run -eq 0 && $postprocess -eq 1 ]] || \
        fail "--track-all-critical-points is a separate mode; omit other --*-only modes, --dry-run and --skip-postprocess"
fi
[[ $critical_points_chunk_size =~ ^[0-9]+$ && $critical_points_chunk_size -ge 1 ]] || \
    fail "--critical-points-chunk-size must be an integer >= 1"
[[ $explicit_spans -eq 0 || $sampling_mode == spans ]] || fail "Choose --points-per-span or --ratio/--ratios, not both"
if [[ $backend == inr && $sampling_mode == spans && $explicit_spans -eq 0 ]]; then
    sampling_mode=ratio
    sampling_ratios=("$inr_sampling_ratio" "$inr_sampling_ratio" "$inr_sampling_ratio")
fi
if [[ $sampling_mode == spans ]]; then
    [[ $backend == mfa ]] || fail "Span sampling requires --backend mfa; for INR/both specify --ratio or --ratios"
    [[ $points_per_span =~ ^[0-9]+$ && $points_per_span -ge 2 ]] || fail "--points-per-span must be an integer >= 2"
    sample_args+=(--points-per-span "$points_per_span")
    resolution="points_per_span_$points_per_span"
else
    sample_args+=(--ratios "${sampling_ratios[@]}")
    resolution="ratio_${sampling_ratios[0]}_${sampling_ratios[1]}_${sampling_ratios[2]}"
    if [[ ${sampling_ratios[0]} == "${sampling_ratios[1]}" && ${sampling_ratios[0]} == "${sampling_ratios[2]}" ]]; then
        resolution="ratio_${sampling_ratios[0]}"
    fi
fi

# Resolve user paths relative to the invoking directory, before changing to root.
[[ -z "$output_dir" ]] || output_dir=$(realpath -m -- "$output_dir")
[[ -z "$mfa_model" ]] || mfa_model=$(realpath -m -- "$mfa_model")
[[ -z "$inr_model" ]] || inr_model=$(realpath -m -- "$inr_model")
[[ -z "$reference_vtp" ]] || reference_vtp=$(realpath -m -- "$reference_vtp")
if [[ "$pvpython_bin" == */* ]]; then pvpython_bin=$(realpath -m -- "$pvpython_bin"); fi
for ((i=0; i<${#tracking_args[@]}; ++i)); do
    if [[ "${tracking_args[i]}" == --ttk-plugin ]]; then
        tracking_args[i+1]=$(realpath -m -- "${tracking_args[i+1]}")
    fi
done
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$repo_root"
output_dataset=$dataset
if [[ $backend == inr ]]; then output_dataset="${dataset}_INR"; fi
output_dir=${output_dir:-"$repo_root/build/local_tracking/$output_dataset/$resolution"}
mfa_model=${mfa_model:-"$repo_root/build/src/$dataset/$dataset.mfa"}
inr_model=${inr_model:-"$repo_root/CoordNet/Exp/$dataset/super-spatial-temporal-64-5-float64.pt"}
# Match reference step sizes independently of the dense sampling grid.
if [[ $dataset == boussinesq_3d && $backend == inr ]]; then
    reference_vtp=${reference_vtp:-"$repo_root/build/src/boussinesq_3d_INR/boussinesq_3d_8_remove_obstacle_regionid_colorid.vtp"}
    spatial_step_size=${spatial_step_size:-0.0125}
    t_step_size=${t_step_size:-0.0125}
elif [[ $dataset == vortex_street_3d ]]; then
    reference_vtp=${reference_vtp:-"$repo_root/build/src/$dataset/${dataset}_16_remove_obstacle_regionid_colorid.vtp"}
    spatial_step_size=${spatial_step_size:-0.00625}
    t_step_size=${t_step_size:-0.00520833}
else
    reference_vtp=${reference_vtp:-"$repo_root/build/src/$dataset/${dataset}_16_regionid_colorid.vtp"}
    spatial_step_size=${spatial_step_size:-0.00694444}
    t_step_size=${t_step_size:-0.0133929}
fi
backends=("$backend")
tracking_max_displacement=${tracking_max_displacement:-$spatial_step_size}
if [[ $track_all_critical_points -eq 1 || $tracker == all ]]; then
    python3 - "$tracking_max_displacement" <<'PY'
import math
import sys
displacement = float(sys.argv[1])
if not math.isfinite(displacement) or displacement <= 0:
    raise SystemExit('Tracking maximum displacement must be finite and positive')
PY
fi
if [[ "$backend" == both ]]; then backends=(mfa inr); fi
if [[ $postprocess -eq 1 && $sample_only -eq 0 && $dry_run -eq 0 && $critical_points_only -eq 0 ]]; then
    [[ -f "$reference_vtp" ]] || fail "Color reference missing: $reference_vtp (set --reference-vtp)"
    python3 - "$spatial_step_size" "$t_step_size" "$match_distance_threshold" <<'PY'
import math
import sys
spatial, time, threshold = map(float, sys.argv[1:])
if not all(math.isfinite(x) and x > 0 for x in (spatial, time)):
    raise SystemExit('Color matching step sizes must be finite and positive')
if math.isnan(threshold) or threshold < 0:
    raise SystemExit('Color matching threshold must be nonnegative')
PY
fi

# load-env.sh requires Spack shell support and must be sourced in this shell.
if [[ $(type -t spack || true) != function ]]; then
    setup=${SPACK_SETUP_ENV:-"${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"}
    [[ -f "$setup" ]] || fail "Initialize Spack shell support or set SPACK_SETUP_ENV to share/spack/setup-env.sh"
    source "$setup"
fi
# The existing environment script may deactivate an inactive environment.
set +e
source ./load-env.sh
env_status=$?
set -e
[[ $env_status -eq 0 && "${SPACK_ENV##*/}" == mfa-TDA-env && -d "$MFA_PATH" && -d "$TORCH_PATH" ]] || fail "Could not activate mfa-TDA-env"

postprocess_trajectories() {
    local region_index kind result_dir prefix
    local -a command
    command -v "$pvpython_bin" >/dev/null || fail "pvpython not found; set --pvpython or PVPYTHON"
    # Check every input before producing any postprocessed files.
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            [[ -f "$output_dir/region_$region_index/$kind/samples_tracking.vtp" ]] || \
                fail "Missing trajectory: $output_dir/region_$region_index/$kind/samples_tracking.vtp"
            [[ -f "$output_dir/region_$region_index/$kind/samples.json" ]] || \
                fail "Missing sample bounds: $output_dir/region_$region_index/$kind/samples.json"
        done
    done
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            result_dir="$output_dir/region_$region_index/$kind"
            prefix="$result_dir/samples_tracking"
            printf '\nPostprocessing region %s (%s)\n' "$region_index" "$kind"
            {
                printf 'cd %q\n' "$invocation_dir"
                printf 'bash %q ' "$repo_root/critical_point_tracking_local.sh"
                printf '%q ' "${original_args[@]}"
                printf '\n'
            } >> "$result_dir/postprocess_invocation.txt"
            # Remove all boundary vertices and incident segments before joining
            # TTK's duplicate endpoints. Never bridge across a removed vertex.
            "$pvpython_bin" --no-mpi --force-offscreen-rendering - "$prefix" "$repo_root/src/python" "$result_dir/samples.json" <<'PY' 2>&1 | tee "$result_dir/postprocess.log"
import json
from pathlib import Path
import sys
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
sys.path.insert(0, sys.argv[2])
from seperate_diffferent_branches import vtk_read_polydata, vtk_write_polydata
prefix = sys.argv[1]
source = vtk_read_polydata(prefix + '.vtp')
metadata = json.loads(Path(sys.argv[3]).read_text())
lo, hi = np.asarray(metadata['sample_min']), np.asarray(metadata['sample_max'])
spacing = np.asarray(metadata['spacing'])
if (lo.shape != (3,) or hi.shape != (3,) or spacing.shape != (3,)
        or not np.isfinite([lo, hi, spacing]).all() or np.any(hi <= lo) or np.any(spacing <= 0)):
    raise ValueError('Invalid sampled-box bounds/spacing')
# TTK exports float32 coordinates. Allow rounding at the faces without removing
# the adjacent sampling layer, even for a very fine grid.
tolerance = np.minimum(8*np.finfo(np.float32).eps*np.maximum(1., np.maximum(abs(lo), abs(hi))),
                       spacing/4)
xyz = (vtk_to_numpy(source.GetPoints().GetData()) if source.GetNumberOfPoints()
       else np.empty((0, 3), dtype=np.float32))
if not np.isfinite(xyz).all():
    raise ValueError('Nonfinite trajectory coordinates')
on_boundary = np.any((np.abs(xyz-lo) <= tolerance) | (np.abs(xyz-hi) <= tolerance), axis=1)
if (source.GetNumberOfCells() != source.GetNumberOfLines()
        or np.any(np.diff(vtk_to_numpy(source.GetLines().GetOffsetsArray())) != 2)):
    raise ValueError('Expected native TTK trajectories with two-point line segments')
segments = vtk_to_numpy(source.GetLines().GetConnectivityArray()).reshape(-1, 2)
keep_cells = np.flatnonzero(~on_boundary[segments].any(axis=1))
keep_points = np.unique(segments[keep_cells])
point_map = np.full(len(xyz), -1, dtype=np.int64)
point_map[keep_points] = np.arange(len(keep_points))
interior = vtk.vtkPolyData()
points = vtk.vtkPoints()
points.SetData(numpy_to_vtk(xyz[keep_points], deep=True))
interior.SetPoints(points)
interior.GetPointData().CopyAllocate(source.GetPointData(), len(keep_points))
for new_id, old_id in enumerate(keep_points):
    interior.GetPointData().CopyData(source.GetPointData(), int(old_id), new_id)
lines = vtk.vtkCellArray()
interior.GetCellData().CopyAllocate(source.GetCellData(), len(keep_cells))
for new_id, old_id in enumerate(keep_cells):
    lines.InsertNextCell(2)
    for point in point_map[segments[old_id]]:
        lines.InsertCellPoint(int(point))
    interior.GetCellData().CopyData(source.GetCellData(), int(old_id), new_id)
interior.SetLines(lines)
interior.GetFieldData().ShallowCopy(source.GetFieldData())
vtk_write_polydata(interior, prefix + '_interior.vtp', 'binary')
stats = {
    'sample_min': lo.tolist(), 'sample_max': hi.tolist(), 'tolerance': tolerance.tolist(),
    'boundary_axes': ['x', 'y', 'z'],
    'input_points': len(xyz), 'boundary_vertices_removed': int(on_boundary.sum()),
    'unique_boundary_vertices_removed': len(np.unique(xyz[on_boundary], axis=0)),
    'unused_interior_vertices_removed': int((~on_boundary).sum())-len(keep_points),
    'input_segments': len(segments), 'boundary_segments_removed': len(segments)-len(keep_cells),
    'interior_points': len(keep_points), 'interior_segments': len(keep_cells),
}
Path(prefix + '_boundary_filter.json').write_text(json.dumps(stats, indent=2) + '\n')
print(f"Removed {stats['boundary_vertices_removed']} boundary endpoint copies "
      f"({stats['unique_boundary_vertices_removed']} unique vertices) and "
      f"{stats['boundary_segments_removed']} incident segments on X/Y/Z faces")
# Join exact duplicates before splitting, never after branch separation.
clean = vtk.vtkCleanPolyData()
clean.SetInputData(interior)
clean.PointMergingOn()
clean.ToleranceIsAbsoluteOn()
clean.SetAbsoluteTolerance(0.)
clean.ConvertLinesToPointsOff()
clean.ConvertPolysToLinesOff()
clean.ConvertStripsToPolysOff()
clean.Update()
# vtkCleanPolyData drops Points and attribute schemas on empty input.
joined = clean.GetOutput() if interior.GetNumberOfCells() else interior
if joined.GetNumberOfCells() != interior.GetNumberOfCells():
    raise RuntimeError('Joining endpoints changed the number of trajectory cells')
vtk_write_polydata(joined, prefix + '_connected.vtp', 'binary')
print(f'Joined exact endpoints: {interior.GetNumberOfPoints()} -> {joined.GetNumberOfPoints()} points; '
      f'{joined.GetNumberOfCells()} cells')
PY
            command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
                     "$repo_root/src/python/seperate_diffferent_branches.py"
                     --input-vtp "${prefix}_connected.vtp" --output-vtp "${prefix}_regionid.vtp" --seed 0)
            printf '%q ' "${command[@]}" > "$result_dir/postprocess_commands.txt"
            printf '\n' >> "$result_dir/postprocess_commands.txt"
            "${command[@]}" 2>&1 | tee -a "$result_dir/postprocess.log"
            command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
                     "$repo_root/src/python/seperate_branch_assign_color2.py"
                     --input-vtp "${prefix}_regionid.vtp" --reference-vtp "$reference_vtp"
                     --output-vtp "${prefix}_regionid_colorid.vtp" --seed 0
                     --distance-dim 3d --distance-mode mean
                     --spatial-step-size "$spatial_step_size" --t-step-size "$t_step_size"
                     --match-distance-threshold "$match_distance_threshold"
                     --matches-csv "${prefix}_color_matches.csv")
            printf '%q ' "${command[@]}" >> "$result_dir/postprocess_commands.txt"
            printf '\n' >> "$result_dir/postprocess_commands.txt"
            "${command[@]}" 2>&1 | tee -a "$result_dir/postprocess.log"
        done
    done
    "$pvpython_bin" --no-mpi --force-offscreen-rendering - "$output_dir" "${#regions[@]}" \
        "$repo_root/src/python" "$reference_vtp" "${backends[@]}" <<'PY'
import csv
import json
from pathlib import Path
import sys
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
root, count, scripts, reference, *backends = sys.argv[1:]
sys.path.insert(0, scripts)
from seperate_diffferent_branches import vtk_read_polydata
root = Path(root)
ref = vtk_read_polydata(reference)
ref_regions = vtk_to_numpy(ref.GetPointData().GetArray('RegionId'))
ref_colors = vtk_to_numpy(ref.GetPointData().GetArray('ColorId'))
rows = []
for region in range(1, int(count)+1):
    for backend in backends:
        folder = root / f'region_{region}' / backend
        final = folder / 'samples_tracking_regionid_colorid.vtp'
        poly = vtk_read_polydata(str(final))
        source = vtk_read_polydata(str(folder / 'samples_tracking.vtp'))
        stats = json.loads((folder / 'samples_tracking_boundary_filter.json').read_text())
        metadata = json.loads((folder / 'samples.json').read_text())
        lo, hi = np.asarray(metadata['sample_min']), np.asarray(metadata['sample_max'])
        tolerance = np.asarray(stats['tolerance'])
        before_xyz = (vtk_to_numpy(source.GetPoints().GetData()) if source.GetNumberOfPoints()
                      else np.empty((0, 3)))
        after_xyz = (vtk_to_numpy(poly.GetPoints().GetData()) if poly.GetNumberOfPoints()
                     else np.empty((0, 3)))
        if np.any((abs(after_xyz-lo) <= tolerance) | (abs(after_xyz-hi) <= tolerance)):
            raise RuntimeError(f'Boundary vertices remain: {final}')
        boundary = np.any((abs(before_xyz-lo) <= tolerance) | (abs(before_xyz-hi) <= tolerance), axis=1)
        before_ids = vtk_to_numpy(source.GetLines().GetConnectivityArray()).reshape(-1, 2)
        keep = ~boundary[before_ids].any(axis=1)
        after_ids = vtk_to_numpy(poly.GetLines().GetConnectivityArray()).reshape(-1, 2)
        if poly.GetNumberOfCells() != int(keep.sum()):
            raise RuntimeError(f'Interior trajectory segment count changed: {final}')
        if not np.array_equal(before_xyz[before_ids[keep]], after_xyz[after_ids]):
            raise RuntimeError(f'Interior trajectory segment coordinates changed: {final}')
        for index in range(source.GetCellData().GetNumberOfArrays()):
            before = source.GetCellData().GetArray(index)
            after = poly.GetCellData().GetArray(before.GetName())
            if after is None or not np.array_equal(vtk_to_numpy(before)[keep], vtk_to_numpy(after), equal_nan=True):
                raise RuntimeError(f'Native cell array changed: {before.GetName()} in {final}')
        with (folder / 'samples_tracking_color_matches.csv').open() as stream:
            matches = list(csv.DictReader(stream))
        for data in (poly.GetPointData(), poly.GetCellData()):
            regions = vtk_to_numpy(data.GetArray('RegionId'))
            colors = vtk_to_numpy(data.GetArray('ColorId'))
            if set(regions) != {int(row['RegionId']) for row in matches}:
                raise RuntimeError(f'Branch labels disagree with match report: {final}')
            for row in matches:
                if not np.all(colors[regions == int(row['RegionId'])] == int(row['ColorId'])):
                    raise RuntimeError(f'Inconsistent branch color: {final}')
        for row in matches:
            if int(row['MatchedReference']):
                expected = ref_colors[ref_regions == int(row['ReferenceRegionId'])]
                if not len(expected) or not np.all(expected == int(row['ColorId'])):
                    raise RuntimeError(f'Color differs from reference: {final}')
        matched = sum(int(row['MatchedReference']) for row in matches)
        rows.append([region, backend, poly.GetNumberOfPoints(), poly.GetNumberOfCells(),
                     len(matches), matched, len(matches)-matched,
                     stats['unique_boundary_vertices_removed'], stats['boundary_segments_removed'], str(final)])
with (root / 'postprocess_summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['region', 'backend', 'points', 'segments', 'branches',
                     'matched_reference', 'fallback_colors', 'unique_boundary_vertices_removed',
                     'boundary_segments_removed', 'colored_trajectory_vtp'])
    writer.writerows(rows)
print(f'Validated boundary removal, interior segments and branch colors. '
      f'Postprocessing summary: {root / "postprocess_summary.csv"}')
PY
}

track_all_trajectories() {
    local region_index kind prefix result_dir
    local -a command
    command -v "$pvpython_bin" >/dev/null || fail "pvpython not found; set --pvpython or PVPYTHON"
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            prefix="$output_dir/region_$region_index/$kind/samples"
            [[ -f "${prefix}_critical_points_all.vtp" && -f "${prefix}_critical_points.json" && -f "$prefix.json" ]] || \
                fail "Missing critical-point reference for $prefix; run --critical-points-only first"
        done
    done
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            result_dir="$output_dir/region_$region_index/$kind"
            command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
                     "$repo_root/src/critical_point_tracking/local_reference/track_all_critical_points.py"
                     --sample "$result_dir/samples" --reference-vtp "$reference_vtp"
                     --max-displacement "$tracking_max_displacement"
                     --spatial-step-size "$spatial_step_size" --t-step-size "$t_step_size"
                     --match-distance-threshold "$match_distance_threshold")
            printf '%q ' "${command[@]}" > "$result_dir/tracking_all_critical_points_commands.txt"
            printf '\n' >> "$result_dir/tracking_all_critical_points_commands.txt"
            "${command[@]}" 2>&1 | tee "$result_dir/tracking_all_critical_points.log"
        done
    done
    python3 - "$output_dir" "${#regions[@]}" "${backends[@]}" <<'PY'
import csv
import json
from pathlib import Path
import sys
root, count, *backends = sys.argv[1:]
root = Path(root)
with (root / 'tracking_all_critical_points_summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['region', 'backend', 'slices', 'all_points', 'legacy_ttk_points',
                     'interior_points', 'interior_linked_points', 'interior_singletons',
                     'interior_segments', 'interior_tracks', 'colored_vtp'])
    for region in range(1, int(count)+1):
        for backend in backends:
            path = root / f'region_{region}' / backend / 'samples_tracking_all_critical_points.json'
            info = json.loads(path.read_text())
            interior = info['interior']
            writer.writerow([region, backend, info['slices'], info['raw']['points'],
                             info['legacy_ttk']['covered_points'], interior['points'],
                             interior['linked_points'], interior['singleton_points'],
                             interior['segments'], interior['tracks'], info['colored_vtp']])
print(f'All-critical-point tracking summary: {root / "tracking_all_critical_points_summary.csv"}')
PY
}

extract_critical_points() {
    local region_index kind prefix result_dir i
    local -a command critical_points_args
    command -v "$pvpython_bin" >/dev/null || fail "pvpython not found; set --pvpython or PVPYTHON"
    critical_points_args=()
    for ((i=0; i<${#tracking_args[@]}; ++i)); do
        if [[ ${tracking_args[i]} == --ttk-plugin ]]; then
            critical_points_args+=(--ttk-plugin "${tracking_args[i+1]}")
        fi
    done
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            prefix="$output_dir/region_$region_index/$kind/samples"
            [[ -f "$prefix.bin" && -f "$prefix.json" ]] || fail "Missing samples: $prefix.bin/.json"
        done
    done
    for ((region_index=1; region_index<=${#regions[@]}; ++region_index)); do
        for kind in "${backends[@]}"; do
            result_dir="$output_dir/region_$region_index/$kind"
            printf '\nCritical-point reference: region %s (%s)\n' "$region_index" "$kind"
            command=("$pvpython_bin" --no-mpi --force-offscreen-rendering
                     "$repo_root/src/critical_point_tracking/local_reference/critical_points_local.py"
                     --sample "$result_dir/samples" --threads "$threads"
                     --chunk-size "$critical_points_chunk_size" "${critical_points_args[@]}")
            printf '%q ' "${command[@]}" > "$result_dir/critical_points_commands.txt"
            printf '\n' >> "$result_dir/critical_points_commands.txt"
            "${command[@]}" 2>&1 | tee "$result_dir/critical_points.log"
        done
    done
    python3 - "$output_dir" "${#regions[@]}" "${backends[@]}" <<'PY'
import csv
import json
from pathlib import Path
import sys
root, count, *backends = sys.argv[1:]
root = Path(root)
with (root / 'critical_points_summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['region', 'backend', 'slices', 'raw_points', 'reference_points',
                     'interior_points', 'reference_csv', 'reference_vtp'])
    for region in range(1, int(count)+1):
        for backend in backends:
            path = root / f'region_{region}' / backend / 'samples_critical_points.json'
            info = json.loads(path.read_text())
            writer.writerow([region, backend, info['slices'], info['raw_points'],
                             info['reference_points'], info['interior_points'],
                             info['reference_csv'], info['reference_vtp']])
print(f'Critical-point summary: {root / "critical_points_summary.csv"}')
PY
}

if [[ $track_all_critical_points -eq 1 ]]; then
    track_all_trajectories
    exit 0
fi

if [[ $critical_points_only -eq 1 ]]; then
    extract_critical_points
    exit 0
fi

if [[ $postprocess_only -eq 1 ]]; then
    if [[ $tracker == all ]]; then track_all_trajectories; else postprocess_trajectories; fi
    exit 0
fi

if [[ $build -eq 1 ]]; then
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
    cmake --build build --target sample_local_field -j 3
fi
sampler="$repo_root/build/src/critical_point_tracking/local_reference/sample_local_field"
[[ $tracking_only -eq 1 || -x "$sampler" ]] || fail "Sampler is missing; rerun without --no-build"
if [[ $sample_only -eq 0 && $dry_run -eq 0 ]]; then
    command -v "$pvpython_bin" >/dev/null || fail "pvpython not found; set --pvpython or PVPYTHON"
fi

# Validate every region before starting an expensive run.
for region_text in "${regions[@]}"; do
    read -r -a region <<< "$region_text"
    [[ ${#region[@]} -eq 6 ]] || fail "Each region needs six numbers: $region_text"
    if [[ $tracking_only -eq 0 ]]; then
        for kind in "${backends[@]}"; do
            model=$mfa_model
            if [[ $kind == inr ]]; then model=$inr_model; fi
            "$sampler" --backend "$kind" --dataset "$dataset" --input "$model" \
                --output "$output_dir/unused" --region "${region[@]}" --threads "$threads" \
                "${sample_args[@]}" --dry-run
        done
    fi
done
[[ $dry_run -eq 0 ]] || exit 0
region_index=0
for region_text in "${regions[@]}"; do
    region_index=$((region_index+1))
    read -r -a region <<< "$region_text"
for kind in "${backends[@]}"; do
    model=$mfa_model
    if [[ "$kind" == inr ]]; then model=$inr_model; fi
    result_dir="$output_dir/region_$region_index/$kind"
    prefix="$result_dir/samples"
    printf '\nRegion %s (%s): %s\n' "$region_index" "$kind" "$region_text"
    command=("$sampler" --backend "$kind" --dataset "$dataset" --input "$model"
             --output "$prefix" --region "${region[@]}" --threads "$threads" "${sample_args[@]}")
    mkdir -p "$result_dir"
    # Record a shell-escaped invocation and commands alongside the results.
    {
        printf 'cd %q\n' "$invocation_dir"
        printf 'bash %q ' "$repo_root/critical_point_tracking_local.sh"
        printf '%q ' "${original_args[@]}"
        printf '\n'
    } >> "$result_dir/invocation.txt"
    if [[ $tracking_only -eq 0 ]]; then
        printf '%q ' "${command[@]}" > "$result_dir/commands.txt"
        printf '\n' >> "$result_dir/commands.txt"
        "${command[@]}" 2>&1 | tee "$result_dir/sampling.log"
    else
        # A resume must not accidentally track samples from a different box/model.
        python3 - "$prefix" "$dataset" "$kind" "$model" "$sampling_mode" "$points_per_span" \
            "${region[@]}" "${sampling_ratios[@]}" <<'PY'
import json
from pathlib import Path
import sys
prefix, dataset, backend, model, mode, density, *numbers = sys.argv[1:]
info = json.loads(Path(prefix + '.json').read_text())
bounds = list(map(float, numbers[:6]))
ratios = list(map(float, numbers[6:])) if mode == 'ratio' else None
if (info['dataset'] != dataset or info['backend'] != backend
        or Path(info['model']).resolve() != Path(model).resolve()
        or info.get('points_per_span') != (int(density) if mode == 'spans' else None)
        or info.get('ratios') != ratios
        or info['region_min'] != bounds[::2] or info['region_max'] != bounds[1::2]):
    raise SystemExit('Existing samples do not match dataset/backend/model/region/sampling density')
PY
    fi
    if [[ $sample_only -eq 0 && $tracker == legacy ]]; then
        command=("$pvpython_bin" --force-offscreen-rendering
                 "$repo_root/src/critical_point_tracking/local_reference/track_local_region.py"
                 --sample "$prefix" --threads "$threads" "${tracking_args[@]}")
        printf '%q ' "${command[@]}" >> "$result_dir/commands.txt"
        printf '\n' >> "$result_dir/commands.txt"
        "${command[@]}" 2>&1 | tee "$result_dir/tracking.log"
    fi
    printf '\nResults: %s\n' "$result_dir"
done
done

if [[ $sample_only -eq 0 && $tracker == all ]]; then
    extract_critical_points
    track_all_trajectories
fi

# A compact index of the grids and resulting trajectory files for this run.
python3 - "$output_dir" "${#regions[@]}" "$sample_only" "$tracker" "${backends[@]}" <<'PY'
import csv
import json
from pathlib import Path
import sys
root, count, sample_only, tracker, *backends = sys.argv[1:]
root = Path(root)
with (root / 'summary.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['region', 'backend', 'nx', 'ny', 'nt', 'dx', 'dy', 'dt',
                     'points', 'segments', 'trajectory_vtp', 'sample_metadata'])
    for region in range(1, int(count)+1):
        for backend in backends:
            folder = root / f'region_{region}' / backend
            metadata = folder / 'samples.json'
            info = json.loads(metadata.read_text())
            tracking = {}
            trajectory = ''
            if not int(sample_only):
                stem = 'samples_tracking_all_critical_points' if tracker == 'all' else 'samples_tracking'
                report = json.loads((folder / (stem + '.json')).read_text())
                tracking = report['raw'] if tracker == 'all' else report
                trajectory = str(folder / (stem + '.vtp'))
            writer.writerow([region, backend, *info['dimensions'], *info['spacing'],
                             tracking.get('points', ''), tracking.get('segments', ''),
                             trajectory, str(metadata)])
print(f'Run summary: {root / "summary.csv"}')
PY

if [[ $postprocess -eq 1 && $sample_only -eq 0 && $tracker == legacy ]]; then
    postprocess_trajectories
fi
