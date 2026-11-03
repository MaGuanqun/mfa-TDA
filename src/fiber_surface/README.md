# Continuous-model fiber surfaces

This implements the **main single-segment method** in `summary_fiber.tex` for
analytical fields and MFA. It evaluates the supplied continuous model directly;
the discovery grid is not a volumetric interpolant. Polyline creases, range
crossings, cone tips, and neural-field adapters are outside this implementation.

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

The ten tests cover derivatives and constrained correction, the truncated and
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
```

The default folders are `build/fiber-surface-results/analytical_<example>/`
(with underscores replacing hyphens) and
`build/fiber-surface-results/hurricane_isabel_3d/`.
`FIBER_RESULTS_ROOT` overrides the results root, and `FIBER_DATASET` overrides
the MFA dataset folder name. There are no step or log subfolders.

Each dataset folder contains:

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
or `./fiber_surface_mfa.sh evaluate`. Direct stage executables use this same flat
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

This also compares the center Jacobian with finite differences. A reproducible
local extraction around the data center is:

```bash
source ./load-env.sh
mkdir -p build/fiber-surface-results
build/src/fiber_surface/fiber_surface_mfa \
  --first build/src/hurricane_isabel_3d/Pf30_3d.mfa \
  --second build/src/hurricane_isabel_3d/TCf30_3d.mfa \
  --a -1100,-27 --b -1100,-22 \
  --min 240,240,45 --max 260,260,55 \
  --step 1 --fiber-step 0.6 --spatial-tolerance 0.08 \
  --residual-tolerance 0.005 --seed-grid 7 --max-triangles 8000 \
  --output build/fiber-surface-results/isabel_tolerance
```

This is a **local sub-box demonstration**, not an extraction of the entire
hurricane volume. Pressure `-1100` is close to the sampled center pressure
`-1107.6`; the temperature interval surrounds the center value `-24.48`.
The choice is a geometric demonstration, not a meteorological feature criterion.

Historical results before automatic/proportional seeding (all have no unresolved events, detected
intersections, or nonmanifold links):

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

Both analytical and MFA fields use `--seed-grid N` to set the seed count on the
shortest physical axis (default `N=7`). Each other axis uses
`ceil(N * length[d] / shortest_length)`, based on the selected spatial box.
For example, lengths `1:2:4` and `N=7` produce `7 x 14 x 28` starts. MFA degree
and span counts do not affect this grid. `--seed-grid 0` is retained as an alias
for the default count of 7. Positive counts must be at least 2.

Seeds are uniformly spaced and include both box endpoints. The total grid is
limited to 1,000,000 starts.
These are starting points for correction, so the number of accepted seeds can
be smaller. Rerun from the seeds stage to apply a change.

TBB processes the independent endpoint and surface seed corrections in parallel.
Each worker has separate MFA decoder buffers, and accepted seeds are collected
in grid order for reproducibility. Fiber tracing and mesh growth still use their
ordered coverage checks. `seeds_report.json` and the final report include the
actual `discovery_grid_shape`; seed files save this shape for restart validation.

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

The range tolerance is in field units; `epsilon_s=tau/|b-a|`. It is separate from
spatial tolerance. Curved boundary chords can violate a very tight range bound
even though their endpoints satisfy `F=q`. Reduce `--fiber-step` or choose a
range tolerance appropriate to the input scale. Fiber roots are solved more
accurately than the allowed chord slack to stabilize component matching.
Frozen fiber samples are never moved during surface attachment.

The implementation checks triangle edge midpoints and centroids, not all points
of a curved preimage. Its checks are numerical evidence, not a global geometric
error certificate. Discovery, geometry comparisons, and cavity repair use simple
serial scans; this is a correctness prototype for small/moderate meshes.
Rank loss, failed transverse contacts, insufficient boundary sampling, exhausted
limits, and failed attachments remain unresolved. Use the report before trusting
an application result. The core does not support surfaces coincident with box
faces or ambiguous tangential contacts.

## Implementation map

- `field.h`, `analytical_field.h`, `mfa_field.h`: physical-coordinate field
  interface, segment derivatives including Hessians, examples, and MFA adapter.
- `corrector.h`: minimum-norm scalar/fiber discovery, phase-plane fiber Newton,
  and simultaneous box face/edge/corner correction, with domain guards.
- `trace.h`: fixed-step RK4/Newton tracing in both directions, loop/component
  matching, and fiber-box contacts before any triangle growth.
- `extract.h`: inward lifting of fiber edges, independent surface seeds,
  neighboring-edge step adaptation and curvature limits, triangle validation,
  front joining, atomic box clipping, contour stitching, and local cavity repair.
  New ordinary triangles use the Delaunay surface test. Boundary/stitching
  triangles use constraint, orientation, overlap, and intersection checks.
  `--min-step` is the weighted recovery threshold, not a hard minimum edge length.
- `output.h`: independent reconstruction of topology and vertex links, raw
  model-agreement measurements, component labeling, OBJ/VTK/JSON export.

The mesh predicates, contour traversal, and MFA decoder come from `../isosurface`.
The fiber constraints and boundary lifecycle are implemented here; the scalar
isosurface marcher is not used as a substitute for finite-segment extraction.
