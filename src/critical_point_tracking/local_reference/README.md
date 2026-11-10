# Local sampled tracking reference

Run the complete three-region Boussinesq CoordNet workflow from the repository root:

```bash
bash critical_point_tracking_local.sh --dry-run
bash critical_point_tracking_local.sh
```

The driver activates the existing Spack environment, builds `sample_local_field`
through the root CMake configuration, validates every region, samples the INR,
extracts TTK critical points at every time slice, tracks every detection (including
saddles and singletons), and assigns reference colors. Default output is
`build/local_tracking/boussinesq_3d_INR/ratio_80/region_N/inr/`.
The editable configuration at the **top of `critical_point_tracking_local.sh`**
contains `dataset`, `backend`, `points_per_span`, `inr_sampling_ratio`, `tracker`,
`threads`, and one region array per dataset. `_INR` selects the INR backend;
unsuffixed dataset names select MFA. Explicit `--backend inr` also works with
unsuffixed names. `--tracker all` runs the full all-point workflow for either
backend; `--tracker legacy` selects the original TTK trajectory pipeline.

Configured Boussinesq boxes (original coordinates, z = time):

| Region | x | y | z = time | Grid (nx, ny, nt) at ratio 80 |
| --- | --- | --- | --- | --- |
| 1 | [-0.05, 0.10] | [0.05, 0.22] | [0.81, 1.20] | 121 × 137 × 313 |
| 2 | [0.25, 0.30] | [1.16, 1.26] | [0.28, 0.66] | 41 × 81 × 305 |
| 3 | [0.32, 0.40] | [1.62, 1.73] | [1.07, 1.30] | 65 × 89 × 185 |

`src/INRModel.h` and `CoordNet/Code/sample_compute_cp_by_slice.py` define
Boussinesq's (10,30,15) blocks over a (1,3,1.5) domain extent. Each block is
0.1 units wide. Ratio 80 means 80 intervals per block axis, giving
**dx = dy = dt = 0.00125**. Counts scale with each requested region's extent,
including both endpoints. These three boxes align exactly with this spacing.
This is a sampling density, independent of the root-finder's initial seed count.

The model is `CoordNet/Exp/boussinesq_3d/super-spatial-temporal-64-5-float64.pt`.
The final colored file is
`samples_tracking_all_critical_points_regionid_colorid.vtp` in each region's
directory. Colors refer to
`build/src/boussinesq_3d_INR/boussinesq_3d_8_remove_obstacle_regionid_colorid.vtp`,
using **Hxy = Ht = 0.0125** (block width / 8). The default motion cutoff is
also 0.0125; `--tracking-max-displacement` changes it independently of coloring.
The raw `samples_tracking_all_critical_points.vtp` retains every detection,
including boundary points; the colored file applies the established six-face
boundary exclusion. Per-slice coverage is checked and written to
`samples_tracking_all_critical_points_coverage.csv`.

The ratio-80 run retained all **8,435 raw detections** and all **2,437 interior
detections**, with exact per-slice coverage:

| Region | Raw points | Interior points | Interior paths/singletons | Reference colors matched | Fallback colors |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 3,506 | 356 | 23 | 23 | 0 |
| 2 | 3,383 | 1,389 | 152 | 148 | 4 |
| 3 | 1,546 | 692 | 46 | 46 | 0 |

The four unmatched region-2 paths retain seeded fallback colors; they are marked
`ReferenceColorMatched=0`. Use the reference's `ColorId` lookup table and fixed
range **0..47** in ParaView. Isolated points are vertex cells, so use Surface or
Points representation to display them alongside the trajectory lines.

To run the previous vortex-street MFA workflow, pass
`--dataset vortex_street_3d`; its output remains
`build/local_tracking/vortex_street_3d/points_per_span_160/region_N/mfa/`.

Both supported datasets have coordinates **(x, y, time)**; your z bounds are time:

| Dataset | x | y | time |
| --- | --- | --- | --- |
| boussinesq_3d | [-0.5, 0.5] | [-0.5, 2.5] | [0, 1.5] |
| vortex_street_3d | [-0.5, 7.5] | [-0.5, 0.5] | [13.5, 15] |

Configured vortex-street boxes (original coordinates):

| Region | x | y | z = time |
| --- | --- | --- | --- |
| 1 | [0.5, 1.0] | [-0.5, 0.5] | [14.5, 15.0] |
| 2 | [0.36, 0.76] | [-0.1, 0.42] | [13.6, 14.0] |
| 3 | [0.46, 0.51] | [-0.1, 0.1] | [14.3, 14.6] |

`--region XMIN XMAX YMIN YMAX TMIN TMAX` replaces the configured boxes;
repeat it to supply multiple regions. `--output-dir` sets their parent directory.
Each region always gets a `region_N/BACKEND` subdirectory, including single-box runs.
`--no-build` reuses the executable. See `--help` for all options.

## MFA points per span

`--points-per-span 160` reads the scalar MFA model's actual nonzero knot spans.
For each axis, with uniform physical span width W, use **h = W / 160**.
The editable `points_per_span=160` is the sampling density: N equal intervals
per span, so S complete spans have `N*S + 1` grid points including both endpoints.
Change that setting at the top of the script, or pass `--points-per-span N`,
to use another density. The grid spacing is always exactly **W / N**.
Spacing is constant along each axis and identical between regions; different
axes can have different spacings if their span widths differ.

The current vortex-street model has **(80, 10, 18)** spans, giving
`dx = dy = 0.000625` and `dt = 0.0005208333333333333` at density 160:

| Region | nx | ny | nt | Scalar storage (GiB) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 801 | 1601 | 961 | 9.182 |
| 2 | 641 | 833 | 769 | 3.059 |
| 3 | 81 | 321 | 577 | 0.112 |

These values are derived from the model, not hard-coded into the sampler.

Let D be the domain minimum and [L,U] the requested (optionally padded) bounds.
Choose `i0 = floor((L-D)/h)` and `i1 = ceil((U-D)/h)`, then sample
`D + i*h` for `i0 <= i <= i1`. The point count is **i1-i0+1**, so partial
spans receive proportional counts without changing h. Bounds snap **outward**
by less than one grid step and stay within the full domain. Floating-point
roundoff near integer grid indices is removed. Very small boxes are expanded
to at least two spatial points and three time steps, as required by legacy TTK.

`--padding PX PY PT` first expands the requested box in physical units, clipped
to the domain; there is no extra padding by default. Metadata distinguishes
requested and sampled bounds. CSV `InsideRequestedRegion` identifies points
inside the original box; native VTK trajectories retain the small snapped margins.

This mode supports scalar, single-block, single-tensor models with uniform knots.
Nonuniform/adaptive knots fail explicitly because a uniformly spaced VTI cannot
have the same number of points in unequal spans. The model's core domain must
match the selected dataset. `--dry-run` loads the MFA knot metadata and reports
dimensions, spacing, sampled bounds and scalar storage, without evaluating the field.

The default model is `build/src/DATASET/DATASET.mfa`; override using `--mfa-model`.
MFA evaluation uses the existing tracker's core-domain normalization and the
same `VolPt` evaluation as `DecodeVar`, reusing decoder buffers within each worker.
Sampling uses bounded batches (`--batch-size`, default 65536); `--threads` defaults
to 3. TTK still loads the local time series and additional ordering/matching data
into memory. Regions run sequentially; temporal-window stitching is not implemented.

## Physical-spacing mode and CoordNet

Explicit `--ratio R` or `--ratios RX RY RT` selects the existing physical-spacing
mode instead of MFA span sampling. Its nominal spacing is `0.1/ratio` on each
axis for both backends. It includes the exact requested/padded endpoints with
`n = ceil(length/h)+1`, so actual spacing may be slightly smaller and may differ
between boxes. This mode is independent of MFA knot counts.

INR defaults to `--ratio 80`; override with `--ratio`/`--ratios`.
`--backend both` requires an explicit `--ratio`/`--ratios`.
Physical-ratio output directories are named `ratio_R` or `ratio_RX_RY_RT`;
use `--output-dir` to resume older `physical_ratio` directories.
CoordNet defaults to `CoordNet/Exp/DATASET/super-spatial-temporal-64-5-float64.pt`;
override with `--inr-model`. It uses full-domain [-1,1] normalization and reversed
(t,y,x) input order as in `INRModel.h`; `--device cuda` accelerates CoordNet only.

```bash
bash critical_point_tracking_local.sh --dataset boussinesq_3d --backend both \
  --region -0.2 0.2 0.5 0.9 0.4 0.6 --ratios 16 16 32 \
  --output-dir build/local_tracking/boussinesq_comparison
```

## Outputs

Each region/backend writes a separate directory beneath `--output-dir`.
The following trajectory outputs describe `--tracker legacy`; the default INR
all-point outputs are listed under “Tracking every per-slice critical point.”

- `samples.bin`: little-endian float64 values, layout (nt,ny,nx), x fastest.
- `samples.json`: requested and sampled bounds, dimensions, spacing, sampling mode, point/span counts,
  model path, and dataset. These are the authoritative conversion metadata.
- `samples.vti`: one XY image with zero-padded scalar arrays for time steps and
  physical time values in FieldData. Its geometric z is the sampled **tmin**, so
  the original script's `ZTranslation = dt` produces original-domain time
  coordinates. Omit this optional copy using `--skip-vti`.
- `samples_tracking.vtu` and `.vtp`: connected TTK trajectories, with **z equal
  to physical time** (at TTK's native coordinate precision). These files preserve
  TTK's native array names, types, point/cell associations, IDs and connectivity;
  no diagnostic arrays, point merging or boundary filtering are applied. The VTU
  is the native tracking mesh; VTP is its line-preserving surface conversion.
- `samples_tracking_points.csv`: coordinates first as `PositionX,PositionY,PositionZ`,
  followed by `PointId`, native point attributes and diagnostic fields. This
  column order works with the existing `matching_ratio.py` CSV reader.
  Connectivity and IDs stored on cells remain in the VTK files; CSV is not a
  replacement for those files.
- `samples_tracking.json`: settings actually used, output counts and time range.
- `sampling.log`, `tracking.log`, `commands.txt`, `invocation.txt`.

The parent directory also contains `summary.csv` with all region dimensions,
spacings, trajectory counts and output paths.

CSV-only fields `TimeIndex` and `Time` identify the frame and physical time.
`InsideRequestedRegion` marks points in the unpadded box;
`OnSampleSpatialBoundary` marks points on the sampled XY boundary. Native output keeps
all tracks, including padding/boundary tracks, for inspection. Postprocessing
removes vertices on all six sampled-box faces as described below. In figures, clip
to the requested box. Do not count a trajectory's entry/exit as a birth/death.

Use distinct output directories for different boxes/resolutions; rerunning a
directory replaces its generated results. `--sample-only` stops after sampling.
Use `bash critical_point_tracking_local.sh --tracking-only` with the same dataset,
boxes and output settings to rerun the selected tracker for all regions, or run
one box directly with legacy TTK:

```bash
source ./load-env.sh
pvpython --force-offscreen-rendering \
  src/critical_point_tracking/local_reference/track_local_region.py \
  --sample build/local_tracking/vortex_street_3d/points_per_span_160/region_1/mfa/samples
```

## Branch separation and reference colors

The driver now postprocesses trajectories after TTK. To process existing results
without resampling or rerunning TTK:

```bash
bash critical_point_tracking_local.sh --dataset vortex_street_3d --postprocess-only
```

It first removes every vertex on the sampled box's X/Y/Z min/max faces, including
the first and last time slices, and all segments incident to those vertices.
Interior portions of trajectories are retained without reconnecting across removed
vertices; unused points are dropped. Bounds come from `samples.json` (`sample_min`
and `sample_max`). The face tolerance allows float32 coordinate rounding and is
capped at one quarter of the sampling spacing per axis.
It then joins exactly coincident segment endpoints, runs
`src/python/seperate_diffferent_branches.py`, then runs
`src/python/seperate_branch_assign_color2.py`. Native `samples_tracking.vtp`
and the reference are left unchanged. The default vortex-street reference is
`build/src/vortex_street_3d/vortex_street_3d_16_remove_obstacle_regionid_colorid.vtp`.
Boussinesq MFA uses `build/src/boussinesq_3d/boussinesq_3d_16_regionid_colorid.vtp`; override with
`--reference-vtp`. Use `--skip-postprocess` to stop after TTK.

Matching uses the step-size-16 reference spacing, independently of sampling density:
vortex street uses `Hxy=0.00625`, `Ht=0.00520833`; Boussinesq uses
`Hxy=0.00694444`, `Ht=0.0133929`. Override with `--spatial-step-size` and
`--t-step-size`. The existing mean nearest-point matcher operates on
`(x/Hxy, y/Hxy, z/Ht)`, where z is physical time, with a threshold of 1 step
(`--match-distance-threshold`). It samples up to 48 points per branch.
Accepted branches copy the reference's `ColorId`; unmatched branches use the
existing random 1..40 fallback with seed 0. A fallback is not a reference match.

Each region writes `samples_tracking_interior.vtp` (boundary vertices and incident
segments removed), `samples_tracking_boundary_filter.json` (bounds, tolerance,
and removal counts), `samples_tracking_connected.vtp`,
`samples_tracking_regionid.vtp` (intermediate branch labels/random colors), and
**`samples_tracking_regionid_colorid.vtp`** (final reference-matched colors).
`samples_tracking_color_matches.csv` records each branch's nearest reference
branch, reference/output colors, distance in steps, and `MatchedReference` status.
Commands and output are saved in `postprocess_commands.txt` and `postprocess.log`.
The parent `postprocess_summary.csv` lists branch and match counts and final paths.
The driver validates that no boundary vertices remain, every retained interior
segment preserves its geometry and cell data, and each branch has a constant
color on both points and cells.

In ParaView, color the final files by `ColorId` using the same lookup table as
the reference. For the supplied vortex reference, keep the color range fixed
at **0..47** across all files; independent rescaling changes displayed colors.

## Critical-point references at every time slice

To extract critical points independently from the saved scalar samples in every
configured region, run:

```bash
bash critical_point_tracking_local.sh --critical-points-only
```

This reuses `mfa_bin_critical_point.py`'s existing pipeline:
`TTKArrayPreconditioning -> Tetrahedralize -> TTKScalarFieldCriticalPoints`.
It reads `samples.bin` and `samples.json` directly, without resampling or rerunning
tracking. Eight slices are loaded per batch by default; change this with
`--critical-points-chunk-size`. `--threads` and `--ttk-plugin` also apply.
Critical points are extracted from each XY scalar field; the Z coordinate in
the exports is its physical time, not a third spatial derivative direction.

Each region/backend directory receives:

- `samples_critical_points.csv` and `.vtp`: the reference at all time slices,
  excluding XY boundary critical points as in `extract_all_critical_points.py`.
- `samples_critical_points_interior.csv` and `.vtp`: also exclude the first and
  last time slices, for comparison with the trajectory files filtered on all
  six X/Y/Z faces.
- `samples_critical_points_all.csv` and `.vtp`: all extracted critical points,
  including boundary points, for inspection.
- `samples_critical_points_slices.csv`: one row for **every** time index, even if
  it has no interior points, with counts and critical-type breakdowns.
- `samples_critical_points.json`: sample path, time range, counts and settings.
- `critical_points.log` and `critical_points_commands.txt`: progress and command.

CSV coordinates come first as `PositionX,PositionY,PositionZ`, followed by
`TimeIndex,Time,CriticalType,IsOnBoundary,OnTimeBoundary,VertexId,Scalar`.
The VTPs contain vertex cells with the same attributes. `VertexId` indexes the
sampled XY grid (X fastest); scalar values and XY positions are checked against
that grid on every slice. Time is stored at the precision of `samples.json`.
`CriticalType` retains TTK's codes (0 minimum, 1/2 saddle, 3 maximum).
The parent directory contains `critical_points_summary.csv` for all regions.

These are discrete critical points on the sampled grid, useful as a reference
for the trajectories. They are not continuous MFA/INR gradient roots.

## Tracking every per-slice critical point

The installed legacy `TTKTrackingFromFields` already uses persistence threshold
zero. Comparing its native output against the explicit per-slice references for
the three sampled vortex regions found **no saddle detections in its tracking
mesh**, and some extrema were absent too. It covered 11,070 / 5,609 / 1,347 unique
detections out of 25,668 / 13,730 / 3,181 respectively, before boundary removal.
The legacy persistence-pair tracking mesh is therefore not a complete export of
scalar-field critical points. Setting the persistence threshold to zero does not
fix that limitation. The driver also explicitly sets first frame 0, last frame
-1 (all frames), and time sampling 1 when the plugin supports those settings.

To track the complete saved per-slice reference, run:

```bash
bash critical_point_tracking_local.sh --track-all-critical-points
```

Run `--critical-points-only` first if the references have not been generated.
The default INR run (or `--tracker all`) automatically performs extraction,
tracking and coloring; no separate commands are required. `--tracking-only`
reuses scalar samples and reruns extraction plus tracking. With `--tracker all`,
`--postprocess-only` reuses saved detections and reruns assignment plus coloring.

This workflow uses **direct assignment of TTK-extracted critical points**,
implemented with SciPy's `linear_sum_assignment`; it does not invoke
`TTKTrackingFromFields`. Minima, saddles, and maxima are matched separately
between adjacent frames. Match costs are squared XY distances divided by the
maximum displacement squared; unmatched birth/death costs are 0.5 each.
The default maximum XY displacement is `--spatial-step-size` (0.00625 for vortex
street); override it with `--tracking-max-displacement`. This controls links,
not which detections survive. It does not use scalar persistence, simplify the
field, skip time slices, filter short tracks, or interpolate across missing frames.
Unmatched and single-frame detections remain explicit vertex cells.

The raw tracking output includes every input detection. For the colored interior
output, all six X/Y/Z faces are excluded **before** matching, so boundary detections
cannot compete for interior matches. Each connected component is already a path
with strictly increasing time (or a singleton), so the branch splitter is not
needed. Existing connectivity labeling and reference-color matching are reused,
including the previously specified XY/time distance scales and reference file.
All original TTK outputs and per-slice references remain unchanged.

Each region/backend directory receives:

- `samples_tracking_all_critical_points.vtp`: all detections, with links and singleton vertices.
- `samples_tracking_all_critical_points_interior.vtp`: independently matched interior detections.
- **`samples_tracking_all_critical_points_regionid_colorid.vtp`**: colored interior result.
- `samples_tracking_all_critical_points_points.csv`: coordinates and point attributes.
- `samples_tracking_all_critical_points_coverage.csv`: old/new/reference coverage for every time slice.
- `samples_tracking_all_critical_points_color_matches.csv`: reference color assignments and fallback status.
- `samples_tracking_all_critical_points.json`: method, settings, counts, and legacy coverage audit.
- `tracking_all_critical_points.log` and `tracking_all_critical_points_commands.txt`.

`TimeIndex`, `Time`, `CriticalType`, `VertexId`, `Scalar`, and coordinates are
preserved exactly. `RegionId` identifies each path, `TrackLength` counts its
detections, `IsSingleton` identifies isolated detections, and `HasPrevious` /
`HasNext` identify its links. `MatchDistance` records the XY distance per segment.
`ReferenceColorMatched` is 1 when the reference color passed the distance threshold,
otherwise 0; fallback colors are not established reference correspondences.
In ParaView, use Surface or Points to see singleton vertices; rendering only
lines hides them. Keep the reference's fixed `ColorId` range (0..47 for vortex MFA).

For the configured three vortex-street MFA regions, all 42,579 raw detections and all 18,476 interior
detections are retained. Interior counts are 11,283 / 6,390 / 803, including
453 / 364 / 94 singletons. Region 2 has 12 paths/singletons using fallback colors;
all other paths match the reference within the configured color distance.
The parent `tracking_all_critical_points_summary.csv` records counts and paths.
Complete detection coverage does not establish unique temporal correspondence:
close features, grid-scale noise, births/deaths, and the motion cutoff can still
produce short or ambiguous paths.

## TTK compatibility and interpretation

The driver reuses `../tracking_script.py` with the original defaults: array
preconditioning, explicit triangulation, persistence threshold zero, unchanged
plugin matching weights, and `ForceZtranslation=1`, `ZTranslation=dt`. The input
plane origin supplies the physical time offset. Output is exported without
rewriting TTK's coordinates or arrays. The old script requested output port 1,
which does not exist in the installed plugin; its native mesh is on port 0.

`--implicit-grid` opts into implicit triangulation to reduce memory, but can
change the result. `--tetrahedralize` explicitly selects the default. For an unusual
installation set `--pvpython` and `--ttk-plugin` (or `PVPYTHON` / `TTK_PLUGIN`).
Set `SPACK_SETUP_ENV` if Spack's shell initialization script is not under
`$SPACK_ROOT/share/spack/setup-env.sh` or `$HOME/spack/share/spack/setup-env.sh`.

Persistence defaults to zero; `--persistence` changes that setting. Matching
weights retain their plugin defaults, as in the original script. Use
`--spatial-weight W` only to explicitly change X/Y matching weights.
The installed TTK version determines which tracking backend and critical types are supported.
Older builds use persistence matching and do not expose the newer critical-point
matching/destruction-cost options. `--destruction-cost` fails clearly when the
plugin does not support it. The run summary records the exposed matching settings.

The compatibility test independently reproduces the original pipeline and
compares coordinates, coordinate precision, connectivity, array names/types/
associations and array values against the local runner. It also reads the outputs
with the existing VTP and CSV analysis readers.

Crop boundaries can change persistence pairings and temporal assignments even
inside the requested box. Check stability under both grid refinement and larger
padding before treating a result as a local numerical reference. Padding does
not establish equivalence to full-domain tracking.

## Verification

```bash
source ./load-env.sh
cmake --build build --target sample_local_field local_sampling_grid_tests local_sampling_fixtures -j 3
ctest --test-dir build -R '^local_sampling_' --output-on-failure
pvpython --force-offscreen-rendering \
  src/critical_point_tracking/local_reference/test_tracking.py
pvpython --no-mpi --force-offscreen-rendering \
  src/critical_point_tracking/local_reference/test_critical_points_local.py
pvpython --no-mpi --force-offscreen-rendering \
  src/critical_point_tracking/local_reference/test_track_all_critical_points.py
```
