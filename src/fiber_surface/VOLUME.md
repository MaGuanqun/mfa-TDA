# Filled range triangles

This implements the four regular-volume stages in
`fiber-surface-tex/summary_fiber.tex`, “Extension to filled range triangles:
volume extraction,” for analytical, MFA, and CoordNet fields.

For counterclockwise range vertices `q0,q1,q2`, the selected set is
`{x in box : n_i dot (F(x)-q_i) >= 0 for i=0,1,2}`, where `n_i` is the inward
unit normal of directed range edge `(q_i,q_(i+1))`. Clockwise, collinear, and
repeated vertices are rejected. The triangle selects a spatial region, which
can have several components, holes through caps, or enclosed cavities.

## Build

From the repository root, initialize Spack shell support if necessary, then:

```bash
source ./load-env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
cmake --build build -j 3
ctest --test-dir build -R '^fiber_(surface|volume)_' --output-on-failure
```

The volume targets include `fiber_volume_explicit`, `fiber_volume_mfa`, `fiber_volume_coordnet`,
`fiber_volume_tests`, `fiber_volume_review_tests`, and `fiber_volume_mfa_tests`. Run MFA executables with one
MPI rank; MPI must be allowed to initialize local sockets.

Volume extraction uses [TetGen 1.6.0](https://github.com/TetGen/TetGen/tree/v1.6.0)
(AGPLv3 or a separate commercial license) and
[earcut.hpp 2.2.4](https://github.com/mapbox/earcut.hpp/tree/v2.2.4) (ISC).
Both dependencies are pinned in the repository's `env.yaml`, with package
recipes in `spack-repo`. To add them to an existing Spack environment, or install
them on a new machine after registering MFA and LibTorch, run from the repository
root:

```bash
source ./create-env.sh
source ./load-env.sh
```

`create-env.sh` preserves an existing environment and its installed packages.
`load-env.sh` exports `TETGEN_PATH`, `EARCUT_PATH`, and their CMake search paths.
CMake uses those installed headers/libraries, so configuring and building no
longer requires downloading either dependency. The recipes install upstream
licenses under each package's `share/licenses` directory. See the
[root setup instructions](../../README.md#setting-up-spack-environment) for
installation on another computer.

For builds outside that environment, CMake still falls back to hash-pinned
source archives in `build/_deps`. It never provisions the Spack environment.
An offline fallback configure can use existing source checkouts after activating
your build dependencies:

```bash
cmake -S . -B build \
  -DFETCHCONTENT_SOURCE_DIR_FIBER_TETGEN=/path/to/tetgen-source \
  -DFETCHCONTENT_SOURCE_DIR_FIBER_EARCUT=/path/to/earcut-source
```

`-DFIBER_ENABLE_VOLUME=OFF` disables only the new volume targets and their
external dependencies.

## Run the two-stage examples

Run both stages for the vortex-street CoordNet pair from the repository root.
Set `default_model` at the top of `fiber_volume.sh` to `analytical` or `mfa`
to change the no-argument example.

```bash
./fiber_volume.sh
```

An explicit model argument overrides the script setting for one run. To
inspect the MFA boundary before filling it, run the two stages separately.

```bash
./fiber_volume.sh mfa boundary
./fiber_volume.sh mfa fill
```

The script activates the Spack environment and builds the selected executable.
The analytical case uses the planar prism. The MFA case uses
`build/src/hurricane_isabel_3d/{Pf30_3d,TCf30_3d}.mfa`. Set `FIRST_MFA` and
`SECOND_MFA` if the pressure and temperature files are elsewhere.
Each model writes under `build/fiber-volume-results/` in its own directory
(`vortex_street_3d_2_variables/`, `analytical_prism/`, or `hurricane_isabel_3d/`). Set `FIBER_RESULTS_ROOT`
to use a different root, and `FIBER_SKIP_BUILD=1` to reuse a build.

The **boundary** stage calls the shared fiber-surface extractor on the three
range edges, adds the necessary box-face caps, validates the assembled shell,
and writes `volume_boundary.ply`. Its JSON report records `boundary_valid`.
The tetrahedral output has zero cells at this point. The **fill** stage reads
that exact PLY in a fresh process, checks it against the same field, range
triangle, and box, reconstructs cavity seeds, runs TetGen, and checks the
result. Open `volume_volume.vtk` for the tetrahedra and read `complete` and
`volume` in `volume_report.json` for the final result. Reusing the prefix
replaces the report and VTK files on each stage.

The PLY stores three little-endian binary64 coordinates per vertex. Each face
has a `list uchar int vertex_indices` with three zero-based indices in outward
winding order and a `uchar patch` (1–3 for range edges, 4–9 for box faces).
The reader expects this layout and the exporter's `fiber_volume_boundary`
header comment, which records the validation flag, three corner-fiber counts,
and the box-contact triangle count. Keep the original PLY for a restart since
a viewer may drop the patch property or comment when resaving. The companion
boundary VTK exposes the patch and shell labels for coloring in ParaView.

All three volume executables accept `--stage boundary --output PREFIX` and
`--stage fill --output PREFIX`. Fill defaults to `PREFIX_boundary.ply`.
Use `--boundary-input FILE` to read a saved PLY from another path. Supply
the same model and selection to both invocations. A direct executable run
with the default `--stage all` passes the boundary to TetGen in memory and
also exports the PLY.

## CoordNet example

The script loads the two float64 TorchScript checkpoints in
`CoordNet/Exp/vortex_street_3d_2_variables/`. Edit `q0`, `q1`, and `q2`
in counterclockwise order to select a filled triangle in the models'
independently normalized velocity components. The output measures volume in
the physical `(x,y,time)` domain.

The default triangle has vertices `(0.402,-0.344)`, `(0.429,-0.320)`, and
`(0.402,-0.296)`. At the default marching step of `0.006`, all three corner
fibers are traced. Two fronts near the lower time boundary remain open,
which prevents the cap from closing. The boundary stage exports 164 triangles
with 14 open edges and stops before tetrahedralization. Reducing only
`--step` to `0.005` closes the shell with 292 boundary triangles. TetGen
then fills it with 363 tetrahedra.

For another scalar model pair, change `first_coordnet`, `second_coordnet`,
and the **full training box** (`domain_min`, `domain_max`). The shared surface
and volume adapter uses these bounds to normalize coordinates and pass them to
the models in `(time,y,x)` order. Change `box_min` and `box_max` to crop the
extraction. Keep the same models, boxes, and range triangle when resuming a fill.

```bash
./fiber_volume.sh coordnet boundary
./fiber_volume.sh coordnet fill
```

To sample the output range when choosing a triangle, run `fiber_volume_coordnet`
with `--inspect`, the model paths (`--first`, `--second`), the training bounds
(`--domain-min`, `--domain-max`), and the desired crop (`--min`, `--max`).

## Analytical examples

```bash
source ./load-env.sh
build/src/fiber_surface/fiber_volume_explicit --example prism \
  --output build/fiber-volume-results/analytical_prism/volume
build/src/fiber_surface/fiber_volume_explicit --example curved-prism \
  --output build/fiber-volume-results/analytical_curved_prism/volume
```

The planar example is exactly the paper's `F=(x,y)`,
`box=[-1,2]^2 x [0,1]`, and triangle `(0,0),(1,0),(0,1)`. Its exact volume is
`0.5`. The curved benchmark applies that same triangle to the paper's
`F=(x+z²/2,y+z²/2)`, `box=[-2,2]^2 x [-1,1]`. Each horizontal section is a
translated unit right triangle, so its exact volume is `1`.

Both default to step `0.2`, fiber step `0.15`, and spatial tolerance `0.015`.
As in the MFA and CoordNet examples, corner fibers use an absolute field
tolerance of `1e-9` and an interior plane tolerance of `1e-6 * fiber_step`.
On the curved prism, a chord along the fiber of `(0,0)` maps at its midpoint
to `(-Δz²/8,-Δz²/8)`, where `Δz` is the difference in endpoint heights.
This midpoint falls outside the segment interval on both axis-aligned sides.
The side marcher records that value and tests each traced vertex against
its assigned range corner. The two side patches can therefore share the chord.
Volume validation separately checks a local correction at the chord midpoint.
Use `--step`, `--fiber-step`,
`--spatial-tolerance`, `--residual-tolerance`, and `--seed-grid` for resolution
studies. Smaller steps do not guarantee successful front closure on every field;
failed fronts are reported, and filling is blocked.

Additional regression examples include `hollow-sphere`, `hollow-cylinder`
(use `--step 0.25`), and the previously implemented `sphere` and `disconnected`
fields. Arbitrary triangles use `--q0 u,v --q1 u,v --q2 u,v`; a spatial sub-box
uses `--min x,y,z --max x,y,z`. All three range vertices must be supplied together.

## MFA examples

The MFA adapter is shared with the segment extractor. It supports separate
scalar files or two scalar variables in one file, physical derivatives, explicit
block selection, and the intersection of the two model domains. No fitting is
performed during volume extraction.

The regression test writes exact polynomial MFA representations of both prism
fields using the current serializer, checks values/Jacobians/Hessians against
the analytical functions, then extracts and validates their volumes:

```bash
source ./load-env.sh
ctest --test-dir build -R '^fiber_volume_mfa_adapter$' --output-on-failure
build/src/fiber_surface/fiber_volume_mfa \
  --first build/src/fiber_surface/fiber_volume_mfa_fixture/curved-prism.mfa \
  --q0 0,0 --q1 1,0 --q2 0,1 \
  --step 0.2 --fiber-step 0.15 --spatial-tolerance 0.015 \
  --residual-tolerance 0.005 \
  --output build/fiber-volume-results/mfa_curved_prism/volume
```

The default variables in a single file are `--variable1 0 --variable2 1`.
For separate scalar files the defaults are both zero. `--block1` and `--block2`
select blocks. Use `--inspect` to sample the field range without extracting.

The existing Hurricane Isabel pressure/temperature models were also tested on
the same local box as the surface example. This triangle selects the side with
pressure at least `-1100`; the other two range edges lie outside the sampled
local range:

```bash
source ./load-env.sh
build/src/fiber_surface/fiber_volume_mfa \
  --first build/src/hurricane_isabel_3d/Pf30_3d.mfa \
  --second build/src/hurricane_isabel_3d/TCf30_3d.mfa \
  --q0 -1100,-40 --q1 10000,-40 --q2 -1100,0 \
  --min 240,240,45 --max 260,260,55 \
  --step 1 --fiber-step 0.6 --spatial-tolerance 0.08 \
  --residual-tolerance 0.005 --seed-grid 7 --max-triangles 8000 \
  --output build/fiber-volume-results/hurricane_isabel_3d/volume
```

This is a local geometry demonstration, not a whole-dataset extraction or a
meteorological feature criterion. The run produced 1,050 boundary triangles and
1,052 tetrahedra, volume `2313.0340112143`, maximum sampled correction distance
`0.02358`, and no reported topology, coverage, or filling failures. There is no
known analytical reference volume for this dataset; these are numerical
consistency checks against the MFA model.

## Measured reference results

The analytical and exact polynomial MFA prism runs pass closure, model checks,
and tetrahedral conformity. Cell counts can differ under floating-point tie
breaking; both representations agree on the enclosed volume.

| Case | Step | Exact volume | Extracted volume | Relative error |
| --- | ---: | ---: | ---: | ---: |
| Analytical planar prism | 0.2 | 0.5 | 0.5000000000 | roundoff |
| MFA planar prism | 0.2 | 0.5 | 0.5000000000 | roundoff |
| Analytical curved prism | 0.2 | 1 | 1.0009565295 | 0.09565% |
| MFA curved prism | 0.2 | 1 | 1.0009565295 | 0.09565% |
| Analytical curved prism, finer boundary | 0.1 | 1 | 1.0003013949 | 0.03014% |
| Original sphere field, filled selection | 0.2 | 4.4284415639 | 4.4119714288 | 0.37192% |
| Disconnected volume, two components | 0.2 | 0.4705142440 | 0.4697524883 | 0.16190% |

For the sphere row use `--example sphere`. Its triangle is
`(0.5,-1),(2,-1),(0.5,1)`: an outer spherical region clipped at `z=-1` with a
complete inner spherical cavity. The test checks automatic cavity exclusion.
The suite also checks an annular cap, a wholly selected box, rejected invalid
boundaries/cells, CLI errors, and replacement of stale output after failure.

Additional validation regressions cover invalid patch IDs, vertices outside the
model domain, imported boundaries with missing adjacency caches, and finite
statistics when every output cell is invalid. Cap correction enforces its
assigned plane and every other active box face, including edges and corners.
Coordinate tests assess a valid mesh rescaled by `1e-8` and `1e-55`, and extract
a wholly selected box translated to coordinates near `1e9`. Local coordinates
for cap triangulation, stable sample averages, and direct orientation-sign
comparisons prevent scale/translation-dependent false failures in these cases.

## Outputs and validation

For output prefix `P`:

- `P_boundary.ply` is the binary64 boundary written by `all` and `boundary`
  and consumed by `fill`.
- `P_boundary.obj`, `P_boundary.vtk`: assembled triangle boundary.
- `P_boundary.poly`: TetGen PLC, containing only boundary vertices, facets,
  and verified cavity seeds. Its comment refers to the validation report;
  an unresolved PLC is diagnostic and must not be filled independently.
- `P_volume.vtk`: tetrahedral mesh (VTK cell type 10).
- `P_tet_boundary.vtk`: actual output boundary faces, including subdivisions.
- `P_report.json`: selection, tolerances, topology, shell signed volumes,
  cavity seeds, model residuals, correction distances, volume, cell/face
  validation, intersection counts, stage times, and unresolved events.

Patch IDs are `1..3` for range edges and `4..9` for
`xmin,xmax,ymin,ymax,zmin,zmax`. VTK point data `incident_patch_mask` preserves
all incident labels at creases. Each PLC facet has a unique parent marker;
its semantic patch is recorded in a comment and in the boundary VTK. Output
facets retain parent IDs and patch labels, even after TetGen subdivisions.

The implementation:

1. Extracts the closed range polyline through the three corners with the shared
   fiber-surface core. Each corner fiber is traced once. The two incident
   surfaces share its mesh edges and vertices. On box faces it
   assembles oriented contact/perimeter graphs, checks for missed contacts by
   constrained face searches, and triangulates selected regions with holes.
   Collinear graph samples are restored after cap triangulation. A short front
   between two distinct box faces can be completed by solving their joint box
   edge contact and validating the added triangle; frozen fibers never move.
2. Validates closure, vertex links, intersections, active-constraint independence,
   model agreement, and retained-side orientation. Classifies connected shells
   by containment and signed volume. Cavity seeds are nearest-hit ray midpoints
   checked against all shells, the field, and boundary clearance. Independent
   grid membership checks assess component coverage.
3. Invokes TetGen as a library with exactly `pV`. No interior points, background
   tetrahedra, quality/size switches, or convex-hull retention are requested.
   TetGen may insert Steiner vertices as needed.
4. Reconstructs cell-face incidence; verifies retained-side orientation with
   adaptive orientation predicates, interior adjacency, duplicate/degenerate
   cells, and geometric overlaps. Verifies every output boundary face against
   its parent and checks coverage of every input facet. Checks output boundary
   topology and field agreement, including new vertices, then compares summed
   cell volume with signed boundary volume.

Raw residuals are in field units. Curved chord samples need not satisfy the
vertex root tolerance: validation corrects temporary copies against their
assigned boundary and violated range constraints, preserving box constraints,
and requires displacement within the spatial tolerance. Raw violations remain
visible in the report. The target inequalities themselves are unchanged.

Exit **0** means the requested stage passed. Exit **1** means invalid input or
a runtime error, and **2** means a partial or invalid result. Filling is blocked after failed
boundary validation. Reusing an output prefix replaces old tetrahedral outputs,
including writing an empty volume if the new boundary fails validation.

Coverage and field-distance checks are sampled, not certified. A missing side
with no discovered seed is treated as empirically inactive and checked through
independent boundary/volume coverage searches; this is not a proof of emptiness.
An entirely empty selection is reported unresolved. Singular contacts, dependent
active constraints, unresolved fronts, and ambiguous face graphs are unsupported.
The serial geometry checks target small to moderate meshes. Curve-only selections end at their contacts with the spatial box.

The default `disconnected` example now passes validation. Its earlier two
reported boundary contacts were numerical false positives: recomputing a shared
vertex's plane distance introduced roundoff that was amplified along an almost
parallel edge. The intersection checker now uses the exact zero plane distance
implied by shared vertex IDs. Regression tests cover both original triangle
pairs under every vertex ordering, genuine intersections, and the complete
two-component volume. Its exact volume is `17/3 - 3*sqrt(3)`.

## Implementation map

`volume_boundary.cpp`, `volume_validation.cpp`, `volume_tetgen.cpp`, and
`volume_io.cpp` implement the four stages and outputs. `volume_examples.h`
contains the analytical benchmarks. `volume_cli.h` shares CLI parsing across
the three entry points. `volume_tests.cpp`, `volume_review_tests.cpp`,
`volume_mfa_tests.cpp`, and `volume_cli_test.py` cover reference volumes,
cap holes/cavities, coordinate transformations, invalid input, missing/duplicate
cells, incorrect facet markers, and partial-result handling.
