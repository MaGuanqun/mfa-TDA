# Continuous-model fiber surfaces

To extract a **filled volume** from the vortex-street CoordNet models, run the
script below from the repository root. It applies the upward-flow triangle to
the full space–time domain. Use `./fiber_volume.sh coordnet all down` for the
downward-flow triangle or `./fiber_volume.sh mfa all 1200` and
`./fiber_volume.sh mfa all 1300` for the two full-domain Hurricane Isabel runs.
The boundary stage builds the fiber-surface patches and caps and saves them
as a binary PLY. If validation succeeds, a second process fills the boundary
with TetGen. Both Hurricane queries passed on the full domain. The vortex-street
preset allows 100,000 triangles per range side. The earlier run at 4,000
triangles per side reached that limit with open edges.

```bash
./fiber_volume.sh
```

Use `./fiber_volume.sh coordnet boundary` to stop and inspect `volume_boundary.ply`.
Use `./fiber_volume.sh coordnet fill` to resume from that file. Both stages also work
with `analytical` and `mfa`. Open `volume_volume.vtk` in ParaView to inspect the tetrahedra.
The default CoordNet results are under `build/fiber-volume-results/vortex_street_3d_2_variables_full/up/`.
See [the volume guide](VOLUME.md#run-the-two-stage-examples) for input paths,
the PLY layout, and validation.

Set `--step` once for both RK4 fiber tracing and surface marching. When a fiber
step fails, the tracer retries it as two corrected half-steps, down to
`--min-step` (default `h/8`). It saves the endpoint of the nominal step or
the first verified box contact. The marcher adapts each triangle height to
its front and curvature. The legacy
`--fiber-step` spelling is an alias for `--step`. To resume a staged run
saved with different steps, rerun the seeds stage with the shared value.

The analytical, MFA, and CoordNet extractors accept a directed segment or a simple open
or closed polyline in range space. Each segment produces a surface patch.
Adjacent patches share a single trace of the fiber at their common range
vertex. The extractor evaluates the supplied continuous field directly.

For `F=(F1,F2)` and `C=[a,b]`, the extractor uses
`d=b-a`, `n=(-d2,d1)/|d|`, `w=d/|d|²`,
`phi=n·(F-a)=0`, and `0 <= s=w·(F-a) <= 1`.

## Build and test

Run from the repository root, with Spack shell support initialized:

```bash
source ./load-env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
cmake --build build -j 3
ctest --test-dir build -R '^fiber_surface_' --output-on-failure
```

If `spack` is unavailable, first source your installed Spack's
`share/spack/setup-env.sh`. Source `load-env.sh` in every new shell that builds
or runs C++ executables. The targets are registered through `src/CMakeLists.txt`;
enabling the isosurface executables is unnecessary.

The surface tests cover derivatives and constrained correction, shared-vertex
intersection predicates, the truncated and
closed spheres, planar and disconnected strips, current-format MFA loading and
physical derivatives, proportional counts for both field types, parallel seed ordering,
pipeline restarts, and CLI rejection/partial-result behavior. The MFA test writes
exact polynomial fixtures with unequal physical extents, degrees, and span counts
to check physical derivative scaling and concurrent evaluation. No fitting is performed during
extraction. Run MFA executables directly with one MPI rank; MPI needs permission
to initialize local communication sockets.

## Staged runs and dataset folders

The repository-root scripts run the five stages and save every output directly
in one folder per dataset:

```bash
./fiber_surface_analytical.sh all
./fiber_surface_mfa.sh
./fiber_surface_coordnet.sh
```

The scripts write under `build/fiber-surface-results/` in
`analytical_<example>/` (with underscores replacing hyphens),
`hurricane_isabel_3d_full/`, and `vortex_street_3d_2_variables_full/`, respectively.
`FIBER_RESULTS_ROOT` overrides the results root, and `FIBER_DATASET` overrides
the MFA dataset folder name. There are no step or log subfolders.

Each single-segment dataset folder contains:

- `run.conf`, `surface_seeds.dat`, `endpoint_seeds.dat`, both seed VTK files,
  and `seeds_report.json`.
- `fibers.state`, `marching.state`, and `stitched.state` for restarting stages.
- Intermediate `*_mesh.vtk`, `*_fibers.vtk`, `*_fronts.vtk`, and `*_report.json`
  files prefixed by `endpoint_fibers`, `marching`, or `stitching`.
- Final `surface.obj`, `surface.vtk`, `surface_fibers.vtk`,
  `surface_report.json`, `report.json`, and `evaluation.conf`.
- `<stage>.log` and `<stage>_command.sh` from the scripts, for stages `seeds`,
  `fibers`, `march`, `stitch`, and `evaluate`.

Rerun one stage with, for example, `./fiber_surface_analytical.sh sphere stitch`
or `./fiber_surface_mfa.sh evaluate`. For a single segment, direct stage executables use this same flat
layout in `--run-dir`. Runs saved in the older numbered-folder layout need to
be rerun from `seeds` or have their files moved to the names above before resuming.

## Analytical example

```bash
source ./load-env.sh
mkdir -p build/fiber-surface-results
build/src/fiber_surface/fiber_surface_explicit \
  --example sphere --output build/fiber-surface-results/sphere
```

The default is the summary's `F=(x²+y²+z²,z)`, `C=[(1,-2),(1,0.5)]`,
`Omega=[-2,2]^3`. The mesh retains the circular opening at `z=0.5`.
Other `--example` values are `closed-sphere`, `strip`, and `disconnected`.
`--a u,v --b u,v` overrides the segment; `--min x,y,z --max x,y,z` selects a
sub-box. Use `--help` for solver, sampling, and triangle-budget settings.

## Connected segments and triangle loops

Pass ordered range vertices with `--polyline`. Repeat the first vertex to
close the curve. These examples use `F(x,y,z)=(x,y)` in the strip example's
box, where each range vertex lifts to a vertical line.

```bash
source ./load-env.sh
build/src/fiber_surface/fiber_surface_explicit --example strip \
  --polyline '1,0;0,0;0,1' --output build/elbow
build/src/fiber_surface/fiber_surface_explicit --example strip \
  --polyline '0,0;1,0;0,1;0,0' --output build/triangle_sides
```

The elbow has two rectangular patches meeting along the fiber at `(0,0)`.
The triangle has three lateral faces with open rims where they meet the box.
For the filled prism, use the volume executable, which takes the same sides,
orients them outward, and adds the two box caps.

The surface accepts either traversal direction. It removes consecutive
duplicate vertices and collinear subdivisions along the same direction.
Retracing, overlapping segments, and self-crossings are rejected.
`--max-triangles` applies to each segment. In ParaView, color the mesh by
the cell scalar `range_segment` to see the patches. The fiber file's
`endpoint` label indexes the normalized range vertices listed in the JSON
report. The point scalar `curve_parameter` is the zero-based segment index
plus its local parameter, with the closing corner mapped to the start.

`fiber_surface_mfa` and the staged `fiber_find_seeds` entry point accept the
same option. For multi-segment runs, the staged executables store intermediate
artifacts in `segment_0/`, `segment_1/`, and so on under `--run-dir`. Final
mesh, fiber, and report files are written in the run directory.

## Vortex street CoordNet example

Run `./fiber_surface_coordnet.sh` with no arguments. The full-domain surface
and volume experiments share the model paths, bounds, range edge and mesh
settings in [model_experiments.sh](model_experiments.sh). The script activates Spack, builds
the stage executables, and runs all five stages. Open `surface.vtk` for the
mesh or `surface_fibers.vtk` for the endpoint fibers.

The default pair in `CoordNet/Exp/vortex_street_3d_2_variables/` predicts
independently normalized velocity components. The experiment selects
`v = 0` and `0.60 ≤ u ≤ 0.90` in physical velocity units over
`[-0.5,7.5] × [-0.5,0.5] × [13.5,15]`. Range endpoints passed to the CLI and the residual tolerance use the
network outputs directly. To select a velocity in original data units, convert
each component with its own training extrema using
`q = 2 * (velocity - min) / (max - min) - 1`.

For other pairs, supply two differentiable scalar TorchScript modules that take
a float64 tensor of shape `[1,3]` on CPU. The seed stage accepts
`--model coordnet --first U.pt --second V.pt`. Set `--domain-min` and
`--domain-max` to the full physical training bounds in `(x,y,t)` order.
The adapter follows critical-point tracking's convention of normalizing each
axis to `[-1,1]` and passing `(t,y,x)` to the networks. It differentiates
through this mapping to compute both Jacobian rows and both full Hessians.
Use `--min` and `--max` to crop the mesh. Later stages recover both boxes and
the two model file digests from `run.conf`.

## Hurricane Isabel MFA example

`--first` supplies F1 and `--second` supplies F2. The supplied files are separate
scalar pressure and temperature models with the same physical coordinate box
`[0,499] × [0,499] × [0,99]`. Inspect their sampled range first:

```bash
source ./load-env.sh
build/src/fiber_surface/fiber_surface_mfa \
  --first build/src/hurricane_isabel_3d/Pf30_3d.mfa \
  --second build/src/hurricane_isabel_3d/TCf30_3d.mfa --inspect
```

This also compares the center Jacobian with finite differences. Run
`./fiber_surface_mfa.sh` to extract the shared pressure–temperature edge over
this full domain. The following command reproduces the older local demo:

```bash
source ./load-env.sh
mkdir -p build/fiber-surface-results
build/src/fiber_surface/fiber_surface_mfa \
  --first build/src/hurricane_isabel_3d/Pf30_3d.mfa \
  --second build/src/hurricane_isabel_3d/TCf30_3d.mfa \
  --a -1100,-27 --b -1100,-22 \
  --min 240,240,45 --max 260,260,55 \
  --step 1 --spatial-tolerance 0.08 \
  --residual-tolerance 0.005 --seed-grid 7 --max-triangles 8000 \
  --output build/fiber-surface-results/isabel_tolerance
```

This is a **local sub-box demonstration**, not an extraction of the entire
hurricane volume. Pressure `-1100` is close to the sampled center pressure
`-1107.6`; the temperature interval surrounds the center value `-24.48`.
The choice is a geometric demonstration, not a meteorological feature criterion.

The local run produces 287 triangles. Each of the 44 traced edges belongs to
exactly one triangle.
The following results predate the exact segment interval and
automatic proportional seeding. Those runs reported zero unresolved events,
intersections, and nonmanifold links.

| Case | Vertices | Triangles | Components | Boundary edges |
| --- | ---: | ---: | ---: | ---: |
| Analytical truncated sphere | 437 | 836 | 1 | 36 |
| MFA polynomial truncated sphere | 437 | 836 | 1 | 36 |
| Closed sphere | 632 | 1260 | 1 | 0 |
| Planar strip | 261 | 447 | 1 | 73 |
| Disconnected strips | 292 | 452 | 2 | 128 |
| Hurricane local patch | 371 | 649 | 1 | 91 |

The truncated sphere has area `9.37445` versus `3*pi` (about 0.53% error),
maximum vertex residual below `1e-9`, and minimum angle about 19.3 degrees.
The Hurricane patch has maximum vertex residual `0.00428`, endpoint residual
`1.12e-9`, and minimum angle about 3.06 degrees. Interior chord samples have
separate residual statistics in each report; their raw residuals are not
expected to match the vertex solver tolerance.

Both MFA inputs must use the same physical coordinate system. Different grids
and extents are allowed; extraction uses the intersection of their core boxes.
The adapter checks Cartesian affine geometry, scalar variables, unweighted
single-tensor splines, and C2 continuity. Unsupported models are rejected.
For two scalar variables in one file, omit `--second` and use `--variable1 0
--variable2 1`. Use `--block1` and `--block2` for explicit block selection;
cross-block stitching is not implemented.

## Seed counts and parallel discovery

All three field types use `--seed-grid N` to set the seed count on the
shortest physical axis (default `N=7`). Each other axis uses
`ceil(N * length[d] / shortest_length)`, based on the selected spatial box.
For example, lengths `1:2:4` and `N=7` produce `7 x 14 x 28` starts. MFA degree
and span counts do not affect this grid. `--seed-grid 0` is retained as an alias
for the default count of 7. Positive counts must be at least 2.

Seeds are uniformly spaced and include both box endpoints. The total grid is
limited to 1,000,000 starts by default. Volume extractors accept an explicit
`--max-seed-starts N` allocation limit for larger lattices.
These are starting points for correction, so the number of accepted seeds can
be smaller. Rerun from the seeds stage to apply a change.

TBB processes the independent endpoint and surface seed corrections in parallel.
Each worker has separate MFA decoder buffers, and accepted seeds are collected
in grid order for reproducibility. The tracer skips an endpoint seed within
`--fiber-seed-radius` of any polyline segment already traced for the same
range point. This radius defaults to `--spatial-tolerance`. `seeds_report.json` and the final report include the
actual `discovery_grid_shape`; seed files save this shape for restart validation.

The CoordNet marcher also uses TBB to predict several front edges concurrently.
It checks each prediction against the current front and inserts triangles in
queue order. For a Slurm run, allocate multiple CPUs per segment, launch each
segment with `srun --cpu-bind=cores`, and set `OMP_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`, and `OPENBLAS_NUM_THREADS=1` for the scalar network queries.

When the corrected advance is shorter than `0.1 * step`, tracing saves the
last accepted point. A second pass tries later corrected seeds within one
nominal step of that tip, including seeds skipped by the duplicate check.
It traces back toward the tip using spacing `max(min_step, step/8)`. A match
requires agreement over **one full nominal step of arclength**. The two paths
must cross the same plane at the same corrected fiber point, have compatible
tangents, and stay within one quarter of the recovery spacing over the overlap.
Recovery publishes a component after verifying its return to the seed or both
domain contacts and accounting for the original seed and every stalled tip.
Ambiguous or incomplete matches retain their diagnostic events.

## Outputs and interpretation

For prefix `P`, each run writes:

- `P.obj`: the triangle mesh, including an incomplete mesh when applicable.
- `P.vtk`: ParaView-readable POLYDATA with segment parameter, endpoint index,
  fiber-component ID, box-face bitmask, and connected mesh-component ID.
- `P_fibers.vtk`: the traced endpoint edges using the same vertex IDs.
- `P_report.json`: selection, domain, tolerances, discovery counts, raw residuals
  at vertices/edge midpoints/centroids, minimum-angle statistics, area,
  connectivity/intersections, extraction time including discovery, and unresolved
  locations/reasons. Loading, fitting, output, and final evaluation are outside
  the reported extraction time.

Open both VTK files in ParaView and display the fibers as colored lines.
An optional static preview requires NumPy and Matplotlib:

```bash
python3 src/fiber_surface/preview.py build/fiber-surface-results/sphere.obj \
  --title 'Analytical truncated unit sphere'
python3 src/fiber_surface/preview.py build/fiber-surface-results/isabel_tolerance.obj \
  --title 'Hurricane Isabel: pressure -1100, temperature [-27, -22]'
```

Exit status **0** means the discovered mesh passes validation with no unresolved
events; **1** indicates invalid input/runtime failure; **2** means an incomplete
or invalid partial result was exported. A successful finite-grid discovery does
not certify that every component exists in the output. An endpoint with no seed
hits is not proved empty. Increase `--seed-grid` to examine discovery sensitivity.

To keep a candidate triangle, the marcher requires `0 <= tau <= 1` at its
centroid and at each vertex that does not belong to a traced fiber, including
untraced box-boundary vertices. At a traced vertex it checks the distance
between the field value and the assigned range endpoint, with tolerance
`1e-9`. The report records raw interval violations at all vertices and edge
midpoints. A straight chord of a curved fiber can have an out-of-range
midpoint even when both ends map to the assigned endpoint. This midpoint
value is diagnostic. The residual tolerance controls surface root accuracy.
The spatial tolerance controls the curvature-based surface predictor. All field backends use an
absolute field tolerance of `1e-9` by default for fiber seeds, traced samples, and box
contacts. Interior fiber corrections also require the distance to the correction
plane to be at most `1e-6 * step`.
Frozen fiber samples are never moved during surface attachment.

The reported samples describe the field at discrete locations on the mesh.
Discovery, geometry comparisons, and cavity repair use simple
serial scans; this is a correctness prototype for small/moderate meshes.
Rank loss, failed transverse contacts, insufficient boundary sampling, exhausted
limits, and failed attachments remain unresolved. Use the report before trusting
an application result. The core does not support surfaces coincident with box
faces or ambiguous tangential contacts.

## Implementation map

- `field.h`, `analytical_field.h`, `mfa_field.h`: physical-coordinate field
  interface, segment derivatives including Hessians, examples, and MFA adapter.
- `coordnet_field.cpp` implements the same field interface with LibTorch autograd.
- `corrector.h`: minimum-norm scalar/fiber discovery, phase-plane fiber Newton,
  and simultaneous box face/edge/corner correction, with domain guards.
- `trace.h`: RK4/Newton tracing with local retries in both directions, loop/component
  matching, and fiber-box contacts before any triangle growth.
- `extract.h`: inward lifting of fiber edges, independent surface seeds,
  neighboring-edge step adaptation and curvature limits, triangle validation,
  front joining, atomic box clipping, contour stitching, and local cavity repair.
  New ordinary triangles use the Delaunay surface test. Boundary/stitching
  triangles use constraint, orientation, overlap, and intersection checks.
  `--min-step` is the weighted recovery threshold, not a hard minimum edge length.
- `curve_surface.cpp` assembles directed patches with shared corner fibers.
- `output.h`: independent reconstruction of topology and vertex links, raw
  model-agreement measurements, component labeling, OBJ/VTK/JSON export.

The mesh predicates, contour traversal, and MFA decoder come from `../isosurface`.
The fiber constraints and boundary lifecycle are implemented here; the scalar
isosurface marcher is not used as a substitute for finite-segment extraction.
