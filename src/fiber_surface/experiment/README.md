# Fiber volume experiments

To fill the paper's **Ours** and **TTK** rows for a new CoordNet model, use the
[CoordNet table runbook](COORDNET_TABLE.md). It includes a known-answer test,
the required model inputs, and the commands that write both rows to CSV.

## Gray–Scott Bubbles over time

The two CoordNets fit concentrations A and B from test trajectory 1
(zero-based) of [The Well's Bubbles experiment](https://polymathic-ai.org/the_well/datasets/gray_scott_reaction_diffusion/).
The window covers frames 2 through 129 inclusive at the native 128 × 128
spatial resolution. The last mesh coordinate is simulation time, from 20 to
1290. In ParaView a Z scale of `2/1270` (about `0.001575`) gives the time
axis the same displayed length as either spatial axis.

With the project's CUDA-enabled `siren` environment, run from the repository root

```bash
archive=gray_scott_reaction_diffusion_bubbles_F_0.098_k_0.057.hdf5
mkdir -p CoordNet/Data/gray_scott_bubbles
curl -fL -o "CoordNet/Data/gray_scott_bubbles/$archive" \
  "https://huggingface.co/datasets/polymathic-ai/gray_scott_reaction_diffusion/resolve/main/data/test/$archive"
python CoordNet/Code/prepare_gray_scott_bubbles.py "CoordNet/Data/gray_scott_bubbles/$archive"
bash CoordNet/Code/fit_vortex_street_components.sh gray_scott_bubbles_A gray_scott_bubbles_B
python src/fiber_surface/experiment/gray_scott_bubbles.py sample --device cuda
```

The sampler reads the float64 TorchScript exports, converts the predictions
back to concentration units and compares triangle membership on all 2,097,152
original grid points. The counterclockwise triangle in (A, B) has vertices
(0.32, 0.35), (0.77, 0.13) and (0.32, 0.56).

Run the following with a Python environment containing conda-forge
`topologytoolkit`, `vtk`, `scipy`, `numpy` and `meshoptimizer` (shared library).

```bash
out=build/gray-scott-bubbles
python src/fiber_surface/experiment/gray_scott_bubbles.py extract --size 128 --times 128 --output "$out/native"
python src/fiber_surface/experiment/reduce_gray_scott_surface.py "$out/native" "$out/budget80k" --budget 80000
```

The native output contains TTK side patches and endpoint fibers. The reduced
output keeps the patch IDs and measures the distance to the native mesh on
each patch in both directions. To obtain a smaller mesh directly from TTK,
use fewer spatial samples and time slices with `extract`. It selects time
slices according to the concentration interpolation error, includes both
endpoints and records all sampled indices in `surface.json`.

In this run the native TTK surface had 1,322,222 triangles. The reduced
mesh had 79,995 triangles, all three patch labels, zero nonmanifold edges
and 16 open edges in the interior of the domain. Only the (0.77, 0.13)
range corner produced an endpoint fiber on the sampled models. Native-grid
source versus CoordNet membership IoU was 0.8818. For the direct TTK
alternative, `--size 61 --times 13` produced 79,644 triangles. The
full-resolution reduction tracked the native surface more closely on all
three sides in the sampled 95th-percentile distance check.

### Direct TTK output for a 200³ grid with time from 0 to 2

The run under `build/gray-scott-bubbles/t0-2-grid200` uses the same 128
training frames and concentration triangle. The new time coordinate is
`(simulation_time - 20) / 635`. The trainer normalizes each axis to
[-1, 1], so frame `i` in the window has network time `2*i/127 - 1`.
The fresh 300-epoch checkpoints are under the run's `models/` directory.

From the repository root, use the training Python environment to evaluate
the models, then use the TTK Python environment to extract the surface.

```bash
run=build/gray-scott-bubbles/t0-2-grid200
experiment=src/fiber_surface/experiment/gray_scott_uniform_ttk.py
python "$experiment" "$run" sample A --device cuda
python "$experiment" "$run" sample B --device cuda
# In the TTK environment
python "$experiment" "$run" extract
```

`sampled_coordnet.vti` contains the two concentration fields at eight
million uniformly spaced coordinates on [-1, 1] × [-1, 1] × [0, 2].
The XML writer takes its input directly from the `ttkFiberSurface` output
port. The extraction reads back `fiber_surface.vtp` and compares its
points, cell connectivity and attributes with the in-memory output.

## Hurricane Isabel and vortex-street metrics

For the Hurricane Isabel outer ring, both vertex metrics use the continuous pressure and temperature MFA fields evaluated at the vertices of each mesh.

| Method | Side vertices checked | Max vertex residual (model-output units) | Max vertex overshoot (fraction of side length) |
| --- | ---: | ---: | ---: |
| Ours, reconnect and fit | 325,275 | $9.99\times10^{-10}$ | $1.15\times10^{-10}$ |
| TTK, native point merging | 1,321,085 | 1.66 | 0.00509 |

Residual measures the distance from the field value to its labeled range segment or shared corner. Overshoot measures how far the projection along an incident segment extends past an endpoint. TTK's largest overshoot was 0.509% of that segment's length. The runbook gives the formulas and a worked example for Isabel. [Vertex residual and overshoot](COORDNET_TABLE.md#how-vertex-residual-and-overshoot-are-measured)

[Open the comparison table (CSV)](model_metrics.csv). Each column is a run.
The first rows give closure, mesh size, volume, and membership IoU. Filter
the metric group for corner fibers, side patches, TTK distances in both
directions, sampled cap agreement, or timings. Blank cells indicate an
unmeasured or inapplicable value. A measured zero appears as 0. The final
columns label the cropped validation examples separately. Isabel volumes use
cubic grid units. Vortex-street volumes use the native `x × y × t` coordinates.

For the Isabel `narrow` and `wedge` triangles, the original data, fitted
pressure and temperature, and tetrahedral mesh were classified at the same
25 million integer grid points. The fitted fields are the float64 samples
saved for TTK, evaluated at grid nodes before piecewise-linear interpolation.
IoU is the number selected by both sets divided by the number selected by
either set. For `narrow`, the mesh and source agree on 626 selected points,
with 124 source-only and 126 mesh-only (IoU 0.7146). The mesh and fitted
field disagree at only 13 points (IoU 0.9830). For `wedge`, the corresponding
IoUs are 0.7451 and 0.9884. These are exhaustive grid counts. The Wilson
intervals in the table apply to the separately labeled random-sample tests.

The full-domain vortex-street boundaries remained open. Their surface
measurements are available, while their tetrahedral-volume and mesh-IoU
cells are blank. The vortex piecewise-linear-reference IoUs use fitted INR
membership as the target. Those reference samples are shared by the three
upward-flow step sizes.

The manuscript's [combined measurement table](../../../paper/data/model-results.csv)
also includes the later Isabel step-size and repair studies. The
[model-results figure](../../../paper/figure/model-results.pdf) shows the repaired
Isabel surface, a local shell before and after repair, and a common time slice
of the saved upward-wake and TTK surfaces. To regenerate both from the saved
experiments under `build/`, run from the repository root with NumPy, SciPy,
Matplotlib and VTK available:

```bash
python3 src/fiber_surface/experiment/export_model_results_table.py
python3 src/fiber_surface/experiment/plot_model_results.py
```

## Vortex street INR with the Isabel outer-ring settings

Place the trained `u-super-spatial-temporal-64-5-float64.pt` and
`v-super-spatial-temporal-64-5-float64.pt` files under
`CoordNet/Exp/vortex_street_3d_2_variables/`. With the project Spack
environment installed, run from the repository root

```bash
./src/fiber_surface/experiment/run_vortex_street.sh up
```

The default upward-flow triangle has velocity vertices `(0.60,0)`, `(0.90,0)`,
and `(0.75,0.35)`. Use `down` for the triangle with its tip at `(0.75,-0.35)`.
The two float64 models predict normalized `u` and `v` over the full
`(x,y,t)` box from `(-0.5,-0.5,13.5)` to `(7.5,0.5,15)`. The nominal
vortex step is one two-hundredth of its unit y extent. The reference is
the Isabel outer-ring run with a nominal step of 0.5 grid units.

| Extraction setting | Ratio to nominal step | Isabel | Vortex |
| --- | ---: | ---: | ---: |
| Nominal tracing and marching step | 1 | 0.5 | 0.005 |
| Minimum retry step | 1/8 | 0.0625 | 0.000625 |
| Endpoint seed duplicate radius | 1 | 0.5 | 0.005 |
| Spatial tolerance for the surface predictor | 1 | 0.5 | 0.005 |
| Fiber corrector plane tolerance | 1e-6 | 5e-7 | 5e-9 |
| Surface range residual tolerance | 0.01 | 0.005 | 0.00005 |

`VORTEX_STEP` updates these settings together. The residual ratio is
numerical because Isabel uses pressure and temperature in its range,
while the vortex model uses normalized velocity.

At every step size, the vortex run uses the Isabel corner-fiber field residual
tolerance of `1e-9`, surface-root and audit projection residual tolerance
of `1e-10`, and derivative and rank threshold of `1e-10`. Both experiments use
a 0.1° minimum triangle angle, an 85° face-to-field normal limit,
40 Newton iterations, 100,000 steps per trace direction and 500,000
triangles per side. The shortest-axis seed count is 40 and the seed
allocation limit is 1,632,160. The vortex domain gives a 320 × 40 × 60
lattice in `(x,y,t)`, or 768,000 starting locations.

The launcher prints its run directory under `build/fiber-volume-results/`.
Open `volume_boundary.ply` for the assembled surface and `audit.json` for
topology and field agreement. The directory also contains the stage commands
in `run.txt` and per-side checkpoints under `stages/`.

## Hurricane Isabel at time step 30

Run `./src/fiber_surface/experiment/run_hurricane_isabel.sh` for the broad
outer ring (preset `outer`, candidate 2 from the larger-ring screen). Its
pressure–temperature triangle is `A=(-25,-19)`, `B=(165,-22)`,
`C=(165,-15)`. The run uses a 0.5-grid-unit step, 202 × 202 × 40 seed
locations, 500,000 triangles per side, and up to 100,000 steps per trace
direction. The outer-ring preset sets both the endpoint seed duplicate-check
radius (`--fiber-seed-radius`) and the surface spatial tolerance to one step
(`1.0 × h`). Set `ISABEL_OUTER_STEP` to override all three distances together. The minimum retry step follows `h/8` and the
range residual tolerance is 0.005.
The outer-ring preset uses an 85-degree limit for the angle between an
oriented triangle normal and the field normal. The marcher checks the three
vertices, three edge midpoints and centroid. The final repair uses the limit
saved with the run. `ISABEL_OUTER_MAX_NORMAL_ANGLE` accepts angles strictly
between 0 and 90 degrees.

### Half-step outer ring and local repairs

For the time-30 candidate-2 run, keep the 85-degree normal limit explicit. The
commands below use the full run directory printed by the first script and a
fresh output directory for the repair. Set `PYTHON` to an interpreter with
NumPy and SciPy for the repair script.

```bash
ISABEL_OUTER_MAX_NORMAL_ANGLE=85 ./src/fiber_surface/experiment/run_hurricane_isabel.sh outer
./src/fiber_surface/experiment/repair_hurricane_isabel_outer.sh "$full_run" "$new_output"
```

| Full-domain native surface | Faces | Open edges | Nonmanifold vertices | Intersecting pairs |
| --- | ---: | ---: | ---: | ---: |
| Uniform `h=0.5`, assembled | 650,258 | 20 | 3 | 0 |
| Two holes filled at `h=0.5` | 650,264 | 12 | 3 | 0 |
| Isolated shell replaced at `h=0.125` | 650,522 | 8 | 2 | 0 |
| Both local refinements joined | 651,410 | 0 | 0 | 0 |

The fixed-step repair closes holes near `(383.33,309.08,42.18)` and
`(425.19,304.53,37.38)`. At the second site, two front vertices move
on their own MFA sides before the fill. This experiment then uses local re-extraction at two sites.
At `(325.06,353.24,42.67)`, two BC faces share one vertex with a small
shell. Its crop contains the whole shell. At `(279.09,350.09,42.96)`,
two vertices each have disconnected incident triangle fans. The crop here
crosses the main shell. The experiment replaces its interior and connects
corresponding AB and BC rims using 146 transition faces. The two crop boxes
and joining contours are specific to this dataset and triangle. Both native
crop runs use `h=0.125` and seven seed locations along their shortest axis.

The script exports `adaptive.ply` with per-face patch labels. It checks
intersections on the joined mesh, all shared seam vertices, and seven field
sites on every triangle. The final surface has 37 edge-connected shells.
All 4,559,870 face-site projections succeeded, with a largest correction of
0.1590 grid units. The largest face-normal error was 82.23 degrees. On the
transition faces alone these maxima were 0.0549 grid units and 28.97 degrees.
At the 25 million original grid locations, surface-interior membership gave
IoU 0.99785 against the fitted fields and 0.90458 against the source samples.
These results describe the boundary and its enclosed grid points.

### Uniform quarter-step and final repair of a saved half-step mesh

These two experiments use the same candidate-2 range triangle, 202 × 202 × 40
seed lattice, 85-degree normal limit and 100,000 steps per trace direction.
The duplicate radius and the surface spatial tolerance both equal the nominal
step. The uniform quarter-step run needs more than 500,000 triangles on each
of its three sides. From the repository root, run

```bash
ISABEL_OUTER_STEP=0.25 ISABEL_OUTER_MAX_NORMAL_ANGLE=85 \
ISABEL_OUTER_MAX_TRIANGLES=2000000 ISABEL_OUTER_STAGE=surface \
./src/fiber_surface/experiment/run_hurricane_isabel.sh outer
```

The script writes the assembled, labeled `volume_boundary.ply`, the audit in
`audit.json`, and the stage commands in `run.txt`. The audit checks open
edges, vertex fans, triangle intersections, corner seams, and field agreement
at seven sites on every triangle.

For the second experiment, set `full_half_run` to the saved native
`h=0.5` run with `volume_boundary.ply` and `stages/`. Select a new
`output_dir` and a Python with NumPy, SciPy and VTK. From the repo root,
with Spack initialized, run

```bash
source ./load-env.sh
cmake --build build --target fiber_repair_boundary fiber_refine_saved_patch fiber_grid_recovery -j 3
python3 src/fiber_surface/experiment/repair_saved_fronts.py "$full_half_run" "$output_dir"
```

The first pass fills simple open rims. For a remaining pinch near a small,
completed corner trace, the script removes its incident faces, reconstructs
the two labeled sides using the saved trace points, and fits new midpoints
to the corresponding MFA sides or shared corner fiber. The exposed patch rim
stays fixed during subdivision. For a saved loop that crosses itself in its
best-fitting plane, it also tests a cyclic order around the loop. Proposed
patches undergo topology, intersection, and field checks with their neighbors
and again in the assembled mesh. The process uses the original half-step
traces throughout. A localized native run at `h=0.125` belongs to the
preceding experiment.

| Surface | Triangles | Open edges | Nonmanifold vertices | Intersecting pairs |
| --- | ---: | ---: | ---: | ---: |
| Uniform `h=0.5`, assembled | 650,258 | 20 | 3 | 0 |
| Uniform `h=0.25`, assembled | 2,480,827 | 15 | 3 | 0 |
| `h=0.5` with final repair | 650,490 | 0 | 0 | 0 |

Halving the nominal step lowered the 95th-percentile distance from 20,000
area-weighted surface samples to their fitted sides from 0.00904 to 0.00233
grid units. The quarter-step run closed the native gaps near x=383 and
x=425. It retained four open edges near x=279 and three near x=325, and
had eight near x=400 where the half-step surface was closed. All 17,365,789
quarter-step face-site projections succeeded. Its 39 edge-connected components
include open patches, so whole-domain interior IoU is left blank.

The `--max-edge-factor 2` control left eight open edges and two pinched
vertices. With the default factor of four, the longest proposed edge was
1.6434 grid units (3.29 times the nominal step). The 85-degree normal
limit and 0.5-grid-unit tolerance stayed in effect. Subdividing the replacement
faces improved the small shell near `(325,353,43)` from 45.8% to 94.7% of
the volume measured on a fitted-field reference with grid spacing 1/32.
The whole-domain membership IoU on the original grid was unchanged by that
subdivision. Set `--fit-levels 0` for the connectivity-only control.

`final.ply` contains the resulting surface, `reconstructed_faces.ply`
contains the replacement faces from the corner-loop repairs, and `repair.json`
records the sites, decisions and closure status. Exit code 2 means the script
saved a partially repaired mesh with remaining topology defects. The resulting
surface with the default settings has 37 closed components. At the 25 million
simulation grid points its interior has IoU 0.99785 with the fitted fields and
0.90458 with the original simulation.

### Recorded normal-check comparison on the half-step outer ring

The earlier ablation used identical seeds and corner traces for both extractions.
The table compares the assembled marching checkpoints and the result of
applying `repair_saved_fronts.py` directly to each checkpoint mesh. This
comparison bypasses the separate per-side cavity search in both runs.
A face counts once in the final column if its normal points against the
outward field normal at any of its three vertices, three edge midpoints or
centroid.

| Face-to-field angle check | Stage | Triangles | Open edges | Pinched vertices | Intersecting pairs | Faces above 90° |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 85° | After marching | 650,248 | 50 | 4 | 0 | 0 |
| Disabled | After marching | 633,384 | 15,728 | 2,100 | 15 | 7,749 |
| 85° | After direct final repair | 650,490 | 0 | 0 | 0 | 0 |
| Disabled | After direct final repair | 640,817 | 6,199 | 1,688 | 7 | 7,004 |

The 85-degree repaired surface had 37 closed components. The unchecked
repaired mesh had 284 edge-connected components including open patches.
All seven projection tests on every face succeeded in both runs.

Starting with the **same** native 85-degree boundary (20 open edges),
both final repairs closed all 37 components, with zero intersecting triangle
pairs. Disabling the comparison only in repair moved two vertices near
`(425,305,37)`, affecting nine faces while preserving every vertex index and
triangle connection. On one CA face, the third vertex landed on the other
side of its fixed edge in the original triangle plane. Its largest angle
at the seven audit sites increased from 2.11 to 179.58 degrees. In both
meshes the maximum fitted-field correction on this face was 0.000634 grid
units. Both repairs gave fitted-field IoU 0.99785 and source IoU 0.90458.

### Isabel comparison and volume cutaway

To regenerate the paper figure from the half-step repaired surface, run the
native MFA volume fill from the repository root in the active build environment.

```bash
mkdir -p build/isabel-presentation/native
source ./load-env.sh
build/src/fiber_surface/fiber_volume_mfa \
    --first build/src/hurricane_isabel_3d/Pf30_3d.mfa \
    --second build/src/hurricane_isabel_3d/TCf30_3d.mfa \
    --variable1 0 --variable2 0 --q0 '-25,-19' --q1 '165,-22' --q2 '165,-15' \
    --step 0.5 --fiber-seed-radius 0.5 --spatial-tolerance 0.5 \
    --residual-tolerance 0.005 --max-normal-angle 85 \
    --seed-grid 40 --max-seed-starts 1632160 --max-triangles 500000 --max-trace-steps 100000 \
    --stage fill-repaired --boundary-input build/isabel_two_methods/half_fitted/final.ply \
    --output build/isabel-presentation/native/isabel
python src/fiber_surface/experiment/plot_isabel_comparison.py
```

The plotter uses NumPy, Matplotlib and VTK. It matches the native boundary to
the displayed surface, slices the tetrahedra at `x=197` and checks both section
areas against the surface contours. The smaller section cuts 606 of the 996,518
tetrahedra. The overviews and cross section show the vertical coordinate at
twice its horizontal scale. The matched surface details use equal spatial scales.

### Earlier outer-ring trials

The 50,000-triangle trial reached its limit on AB. The larger run saves shared
seeds and corner traces under `stages/`, marches and stitches the three sides
independently, then assembles and validates the full boundary. Follow
`AB.log`, `BC.log` and `CA.log` inside the results directory printed by the script.
Before stalled-trace recovery, at `h=1` with both spatial thresholds equal to the step, candidate 2 gave
1, 58 and 40 traces at A, B and C, with two incomplete at B. All three sides
finished below budget. The assembled boundary has 176,888 triangles, 406 open
edges, 37 nonmanifold vertices and 21 intersections. Of the open edges, 382
were already open on a marching front, 17 lost a neighbor during assembly,
and seven have a corner trace with only one attached patch. Validation stopped
before tetrahedralization.

Rerunning candidate 2 from the same saved seeds with a required shared arc
of one grid unit recovered the two incomplete B traces as one closed loop
of length 21.747. The resulting corner counts were A = 1, B = 57 and C = 40,
all closed. Replaying the native MFA tracer separately from each of the three
original stalls at a step of 1/256 gave a first return after 21.757 grid units.
TTK extracted two local loops on the unit grid and at spacing 1/16. On local
grids of spacing 1/32 and 1/64 it extracted one loop, with lengths 21.771
and 21.760. Both adjacent patches attached to all 175 edges of the recovered
loop. The assembled boundary had 177,328 triangles and 350 open edges,
including 334 unfinished patch-front edges, nine exposed by rejected faces
and seven corner edges with only one patch attached. Its 27 nonmanifold
vertices and 18 intersecting triangle pairs also failed boundary validation.

The measurements in the next paragraphs come from the earlier run with a 0.35-grid-unit
surface tolerance and duplicate-check radius. Its corner stage found 1 trace at A, 66 at B and 43 at C. Two of the
66 B traces remained incomplete, with three insufficient-progress events
near `(400.3,161.2,41.9)`. The volume validator requires all events to be resolved.
The original grid selects 98,788 points across levels 36–44. Of 291 groups
connected through grid faces, the largest contains 96,108 points. The main ring's
horizontal footprint has area-equivalent inner and outer radii of 84.0 and
131.4 grid units. At the integer grid points, the fitted fields select
98,537 points and agree with the original fields on 93,711 selected points
(fitted/source IoU 0.9044).
At the 500,000 limit, AB, BC and CA finished marching with 74,748, 44,226
and 50,930 triangles. Assembly rejected 336 BC triangles at inconsistent
corner incidences. The resulting boundary has 169,568 triangles, 1,772
open edges, 139 nonmanifold vertices and 1,262 intersections. Validation
stopped before tetrahedralization. Its total side area is 75,749.2 square
grid units, compared with 75,860.1 for TTK. The sewn B seam has 76.7% of
the TTK B-fiber length.
With 50,000 triangles per side, AB supplied 66.36% of the TTK AB area,
and the assembled surface had 5,766 open edges. With the higher limit, AB's
assembled area reached 100.004% of TTK.
Across the three sides, the MFA/TTK area ratio rose from 83.62% to 99.85%.

Run `./src/fiber_surface/experiment/run_hurricane_isabel.sh narrow` for the
earlier narrow experiment. Its pressure–temperature triangle is
`A=(-1000,-22)`, `B=(-900,-22)`, `C=(-950,-18)`. It uses a 1-grid-unit step
for tracing and marching, a 202 × 202 × 40 seed lattice, a budget of 50,000
triangles per side, and up to 100,000 steps per trace direction. The range
residual tolerance is 0.005 and the spatial tolerance is 0.35 grid units.
At spatial tolerances of 0.08 and 0.20 the surface retained 146 and 30 open
edges, respectively. At 0.35 it closed with 5,234 boundary triangles and
7,317 tetrahedra. The maximum sampled spatial correction was 0.2104 grid units.
All tetrahedra belong to one face-connected volume. Counting face neighbors
on the original grid gives 67 groups, with 310, 269, and 79 samples in the
three largest. Sampling the validated mesh on those same integer grid points
also splits it into 60 face-connected groups. The pressure and temperature
panels show the same selection at `z=46`, where it contains the most source points.

The `wedge` comparison uses `A=(-1000,-21)`, `B=(-850,-21)`,
`C=(-850,-15)`, the same step, seed grid, and budgets, and a spatial
tolerance of 0.20. It gives a wider band with 9,424 boundary triangles and
13,655 tetrahedra. Its pressure and temperature sections are at `z=45`.

To reproduce the validated storm-core volume, run
`./src/fiber_surface/experiment/run_hurricane_isabel.sh core`. The triangle
has vertices `A=(-3300,-35)`, `B=(-550,-55)`, and `C=(-550,19)`. The constant
pressure edge BC traces the outer wall of the selected core. Portions of all
three range edges contribute surface patches. Corner fibers were found at B and C.
The experiment sets the tracing and marching step to 5 grid units, the shortest-axis
seed count to 40, the spatial tolerance to 1.6 grid units, the range residual
tolerance to 0.005, and the triangle budget to 15,000 per side. The run validates
the resulting mesh before comparing its membership against all 25 million source
points. Run outputs start with `hurricane_isabel_t30_core_h5_tol1p6_seed40.`.

For the broader, diagnostic pressure–temperature selection, run
`./src/fiber_surface/experiment/run_hurricane_isabel.sh storm`. This preset
uses `A=(-1500,-36)`, `B=(100,-36)`, `C=(150,10)` on the full domain with
40 seeds on the shortest axis. The shared tracing and marching step is 5
grid units, the minimum retry step is 0.625, the spatial tolerance is 0.4,
the range residual tolerance is 0.005, and the diagnostic limit is 3,500
triangles per side. At this step size, corner tracing reported stalled and
step-limited traces before marching. The script retains the partial mesh.
Output directories start with `hurricane_isabel_t30_storm_h5_budget3500_seed40.`.

| Preset | Shared step | Source-grid samples inside the triangle | Horizontal columns selected (of 250,000) | Validated mesh volume (cubic grid units) | Mesh/source IoU |
| --- | ---: | ---: | ---: | ---: | ---: |
| `outer` | 1 | 98,788 | 32,834 (13.13%) | — | — |
| `narrow` | 1 | 750 | 522 (0.21%) | 747.69 | 0.7146 |
| `wedge` | 1 | 1,774 | 1,124 (0.45%) | 1,759.94 | 0.7451 |
| `core` | 5 | 263,230 | 9,514 (3.81%) | 262,213.01 | 0.9746 |
| `storm` | 5 | 929,974 | 42,879 (17.15%) | — | — |
| `loop` | 1 | 480 | 302 (0.12%) | 475.79 | 0.8402 |

The validated core extends from approximately `(198.5,171.7,3.5)` to
`(309.3,281.3,69.3)` in grid coordinates. Its maximum sampled spatial correction
was 1.09 grid units. The wider `storm` query shows more of the horizontal storm
structure on the source grid, but its extracted boundary remained open.

`isabel_footprint.png` and `.pdf` show the number of selected source-grid
levels in each horizontal column, followed by pressure and temperature maps
at the level with the most selected points. Both field maps show the same
source-grid selection outline and MFA mesh section. Each field has its own
color scale, and all three panels use the full horizontal domain.
`isabel_triangle.png` and `.pdf` show the range triangle and the MFA side patches.
If boundary validation fails, the script plots the partial surface and the
source footprint, saves diagnostics, and exits with status 2. Mesh volume and
source-grid IoU require a validated fill.

### TTK comparison on a 500 × 500 × 100 grid

From the repository root, compare the outer ring with TTK using

```bash
./src/fiber_surface/experiment/run_hurricane_isabel_ttk.sh
```

Pass `both` to compare the earlier narrow and wedge cases together, or pass
either name individually. For each named preset, the script uses a saved
completed MFA fill when available, then a saved partial boundary. If neither
exists it runs the MFA experiment first. To compare particular meshes, pass
their run directories as arguments.
For example, after running the `1200` and `1300` cases from `fiber_volume.sh`,
use `./src/fiber_surface/experiment/run_hurricane_isabel_ttk.sh build/fiber-volume-results/hurricane_isabel_3d_full/{1200,1300}`.
Each supplied mesh must come from the same time-30 pressure and temperature models.
Its `volume_report.json` supplies the three range vertices in order.

The runner samples `Pf30_3d.mfa` and `TCf30_3d.mfa` once at all 25,000,000
integer locations in `[0,499] × [0,499] × [0,99]`. The two float64 sample
files take 400 MB combined. All triangles in an invocation share these samples. TTK processes
the full image with an implicit triangulation of 147,906,594 tetrahedra.
The reference tests how well the continuous-model extraction agrees with the
piecewise-linear interpolation of the **fitted MFA fields**. The source-grid
IoU above uses the original simulation values.

Each invocation creates a fresh `hurricane_isabel_t30_ttk.*` folder under
`FIBER_RESULTS_ROOT` (default `build/fiber-volume-results/`). Its per-run
folders contain the following results.

| File | Contents |
| --- | --- |
| `sides.vtp` | TTK triangles labeled `patch` 1 = AB, 2 = BC, 3 = CA. Overlay with the MFA boundary identified in `comparison.json`. |
| `fibers.vtp` | TTK line segments labeled `corner_id` 0 = A, 1 = B, 2 = C, when corner fibers exist. |
| `reference.json` | TTK query, grid, counts by patch and corner, edge incidence after welding, and extraction and welding time. |
| `corner_traces.json` | For checkpointed runs, trace counts, lengths and distances to TTK for each corner before marching. |
| `fit_vs_source.json` | Membership counts and fitted/source IoU over all 25 million grid points when original simulation fields are available. |
| `comparison.json` | MFA fill and boundary status, compared mesh path, physical patch areas and seam lengths, and closest-point distances in both directions. On an open mesh, the seam measurement covers edges shared by the two labeled patches. |
| `section.png`, `section.pdf` | Side-by-side cross sections of our mesh and TTK, colored by range edge, over the selected fitted-grid points on the level with the most hits. |
| `side_footprints.png` | Separate horizontal projections of all TTK triangles on AB, BC and CA across all heights. |

Distances use 10,000 seeded samples per patch or corner seam in each direction.
Each domain axis is scaled to `[0,1]` before sampling by triangle area or line
length. The report gives the median, 95th percentile, and sampled maximum.
Physical areas and lengths use the original grid coordinates. The MFA
experiment's reports supply the volume and membership measurements.

| Case | MFA side triangles | TTK side triangles | TTK corner segments |
| --- | ---: | ---: | ---: |
| `outer` (50,000 per side, open boundary) | 143,894 | 675,857 | 14,245 |
| `outer` (500,000 per side, open boundary) | 169,568 | 675,857 | 14,245 |
| `narrow` | 5,234 | 23,632 | 2,001 |
| `wedge` | 9,424 | 40,761 | 2,160 |

The sampler builds in the repository's Spack environment. The TTK stage uses
`TTK_PYTHON`, the active Python 3, or the Conda environment `mfa`, with TTK
Python bindings, VTK, NumPy, SciPy and Matplotlib. `TTK_THREADS` sets the thread count
for sampling and extraction (default 3). To reuse an existing sample pair,
run `python ttk_hurricane_isabel.py SAMPLES_DIR MFA_RUN_DIR [MFA_RUN_DIR ...]`
from this folder with that Python environment active.

### Shared-grid surface recovery

`fiber_grid_recovery` reads a triangle and model from an existing `stages/`
directory. It evaluates the fitted fields on a lattice, splits each spatial
cell into six tetrahedra, and clips their affine range maps against the three
query sides. Neighboring tetrahedra reuse the same crossing on each grid edge.
Adjacent range sides reuse the same corner-fiber vertex on each grid face.
The result contains the side triangles, with patch labels 1 (AB), 2 (BC), and
3 (CA). A cropped domain may have open edges on the crop faces.

From the repository root, with the Spack environment active, the following
commands reproduce the half-grid reconstruction and the B-loop refinement.
`STAGES` identifies the candidate-2 run and `OUT` an existing output directory.

```bash
cmake --build build -j 3 --target fiber_grid_recovery
bin=build/src/fiber_surface/fiber_grid_recovery
"$bin" grid "$STAGES" 0.5 "$OUT/grid_half"
"$bin" grid "$STAGES" 0.5 "$OUT/grid_refined" --refine 397.5,156.5,41.5 404.5,166.5,42.5 0.0625
"$bin" recover-grid "$STAGES" 0.5 "$OUT/model_refined" --refine 397.5,156.5,41.5 404.5,166.5,42.5 0.0625
```

Refinement subdivides coordinate intervals on each axis. Those coordinate
planes extend across the domain, so neighboring cells agree on every face.
The grid spacing controls field sampling independently of the marcher's
`--step` and triangle budget. A query that passes exactly through a grid
vertex or has a corner exactly on a grid edge can require special handling.
The experimental extractor reports those coincidences.

`grid` retains the affine crossing positions. `carrier-grid` fits side points
on their grid edges and corner points in their grid faces. `recover-grid`
also tries a free corner solve when the fitted fiber misses that grid face.
It restores linear positions around inverted, degenerate or intersecting
triangles and rechecks the mesh. `fit-grid` is a comparison that fits both
range coordinates independently at each vertex.

The binary PLY contains double coordinates and a per-face patch label. The
JSON reports topology, intersections, seam residuals against the continuous
field, and field-projection distances from 20,000 area-weighted side samples
in grid units. The reconstruction exports a surface. The production volume
validator has its own orientation, corner and coverage tests.

### Recovering box caps on the Gray–Scott CoordNet mesh

For a saved uncapped run whose range sides have no shared corner fibers,
`repair_coordnet_caps.py` closes two-edge excursions and compares the remaining
side components and their box edges with independently sampled field contours.
It substitutes fitted components where the native mesh has interior openings
or misses a contact. The cap assembler triangulates the resulting boundary
loops, preserving their vertex IDs across seams.

Run from the repository root on a compute node. Set `RUN` to the saved run,
`OUT` to an existing output directory and `PYTHON` to the environment with
PyTorch, NumPy, SciPy and scikit-image. The example uses the Gray–Scott box
and a 0.005 recovery lattice for the 0.01 marching result.

```bash
source ./load-env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
scripts=src/fiber_surface/experiment
workers=${SLURM_CPUS_PER_TASK:-8}
"$PYTHON" "$scripts/repair_coordnet_caps.py" "$RUN" "$OUT/sides-grid.npy" \
  --sample 401 401 201 --workers "$workers" --output "$OUT/sides.ply"
bin=build/src/fiber_surface/fiber_grid_recovery
"$bin" cap-sides "$RUN/stages" "$OUT/sides.ply" "$OUT/capped"
"$bin" audit "$RUN/stages" "$OUT/capped.ply" "$OUT/audit" --all-face-samples
"$PYTHON" "$scripts/audit_coordnet_caps.py" "$RUN" "$OUT/capped.ply" \
  "$OUT/coverage.json" --spacing .0025 --workers "$workers"
```

Omit `--sample` to reuse a completed grid. The audit reports edge incidence,
triangle intersections and bounded projections for seven sites per triangle.
The face comparison measures the distance of each membership disagreement
from the side rim. The PLY keeps its volume-validation flag clear pending the
separate volume validator.

Add `--fine-run "$FINE_RUN"` to prefer matching components from a saved finer
march. The command checks the model digests, query, spatial coverage and all
six face contacts before substituting a component. It uses
`fiber_repair_boundary` for interior holes and fits missing corner contacts
on the shared box edge. Both the Python environment and the Spack runtime
must be active for this option.

### Repairing side geometry on the capped Gray–Scott mesh

`refine_coordnet_sides.py` checks the three vertices, three edge midpoints and
centroid of every side triangle against the saved field. For this run, five
thin components had reversed normals at their triangle centers. The repair
fits those components in time slices, locating their terminal points and
fitting each ring to its range edge. This method requires an isolated,
star-shaped cross section around a tracked extremum. On the other components,
each vertex adjustment is checked against every incident side triangle.

Use the same compute-node environment and variables as above. Set `INPUT` to
the capped PLY from the first repair. The normal-angle target defaults to the
value in `stages/run.conf`. The optional rim pass subdivides chords according to
the displacement between their midpoint and the fitted box contour.

The volume sample corrector follows an outgoing Newton step up to the first
box face it reaches and continues the solve in that plane. Acceptance requires
an in-domain point that meets the original range and distance tolerances.
To replay particular samples, put one `patch x y z` record per line in a
text file and run `fiber_grid_recovery project-samples STAGES SAMPLES.txt PREFIX`.
The CSV records the corrected position, solve status and displacement for
each input point.

```bash
"$PYTHON" "$scripts/refine_coordnet_sides.py" "$RUN" "$INPUT" "$OUT/side-sites.npz" \
  --diagnose --workers "$workers"
"$PYTHON" "$scripts/refine_coordnet_sides.py" "$RUN" "$INPUT" "$OUT/refitted-sides.ply" \
  --repair --use-diagnostic "$OUT/side-sites.npz" --workers "$workers"
"$PYTHON" "$scripts/refine_coordnet_sides.py" "$RUN" "$OUT/refitted-sides.ply" \
  "$OUT/refined-rims.ply" --rim-tolerance 0.00001 --workers "$workers"
"$bin" cap-sides "$RUN/stages" "$OUT/refined-rims.ply" "$OUT/refitted-boundary"
"$bin" audit "$RUN/stages" "$OUT/refitted-boundary.ply" "$OUT/refitted-audit" --all-face-samples
"$PYTHON" "$scripts/audit_coordnet_caps.py" "$RUN" "$OUT/refitted-boundary.ply" \
  "$OUT/refitted-coverage.json" --spacing .0025 --workers "$workers"
```

For the whole-domain Isabel candidate 2, the radius used to suppress duplicate
seeds and the surface tolerance were both `h` in every native trial.
Open length is the sum of one-sided mesh-edge lengths in physical grid units.
The extra repair pass tried bridges between open edges, local changes of
diagonal, and splits of a rim triangle. It allowed overlap in tangent-plane
projection for repairs of a bounded three-edge hole. Each accepted edit
passed the patch intersection tests and improved open-edge or nonmanifold
vertex counts on the assembled surface.

| Extraction | March step | Face normal limit | Triangles | Open edges | Open length | Nonmanifold vertices | Intersecting pairs |
|---|---:|---:|---:|---:|---:|---:|---:|
| Native march | 1 | 55° | 177,328 | 350 | 291.53 | 27 | 18 |
| Native march | 1 | 85° | 177,402 | 180 | 138.13 | 21 | 18 |
| Native march and bounded-gap repair | 1 | 85° | 177,408 | 174 | 133.15 | 19 | 18 |
| Native march | 0.8 | 85° | 269,913 | 137 | 94.25 | 8 | 18 |
| Native march and bounded-gap repair | 0.8 | 85° | 269,921 | 113 | 75.88 | 8 | 18 |
| Shared grid with guarded fits, 0.5 spacing and 1/16 through B | | | 4,596,808 | 0 | 0 | 0 | 0 |

The recovered surface had 36 edge-connected shells. All 20,000 sampled
points projected onto the continuous fitted sides, with a 95th-percentile
displacement of 0.00468 grid units. Of 41,064 seam vertices, 2,059 exceeded
the production corner residual tolerance of `1e-9`. The output is a surface
mesh with measured field error. It has no validated tetrahedral fill.

The small B loop distinguishes field accuracy from mesh closure. The uniform
0.5 grid split it into two loops. Refinement followed by guarded vertex fits
recovered one loop, with a 95th-percentile distance of 0.00019 grid units
from the independent 1/64-spacing TTK loop to the assembled seam.
At `h=1`, the native tracer represented this loop with 175 chords (median
length 0.125), and both sides attached to every chord. At `h=0.8`, it closed
the loop with 27 chords (median length 0.799). Both sides attached to 26 and
the remaining chord was absent from the assembled mesh. The smaller nominal
step improved the global defect counts while this loop received coarser
local sampling.

### Original small-loop selection

Run `./src/fiber_surface/experiment/run_hurricane_isabel.sh loop` from the repository
root, or `./run_hurricane_isabel.sh loop` here. It loads the Spack environment and
builds `fiber_volume_mfa`. The analysis uses NumPy, VTK and Matplotlib from
`PYTHON`, the active Python 3, or the Conda environment `mfa_env`.

1. Read pressure from `Pf30_3d.mfa` and temperature from `TCf30_3d.mfa`
   in `build/src/hurricane_isabel_3d/`, using scalar variable 0 from each.
   Inspect their common box on the full seed lattice.
2. Select the counterclockwise triangle with vertices `A=(-1100,-27)`,
   `B=(-1100,-22)`, `C=(-1200,-22)` in the stored field units. It shares
   the constant-pressure edge used in the existing Isabel experiment.
3. Trace the three corner fibers and march the three sides at **h = 1.0 grid
   unit**. The box is `[0,499] × [0,499] × [0,99]`. Forty starts on its
   shortest axis give `202 × 202 × 40 = 1,632,160` locations per lattice.
   Set the volume CLI's `--max-seed-starts` to this exact count. The minimum
   retry step is 0.125, the spatial tolerance is 0.08 grid units, the range
   residual tolerance is 0.005, and the budget is 8,000 triangles per side.
4. Check boundary closure, orientation, constraints and sampled coverage.
   The extractor fills a validated boundary with TetGen and checks cell
   validity and boundary conformity.
5. Compare tetrahedral membership against the pressure–temperature selection
   at all 25 million original grid points in `build/examples/hurricane_isabel/`.
   This measures the combined effect of the MFA fit and the mesh extraction.

Every invocation creates a fresh `hurricane_isabel_t30_p1200_h1_seed40.*`
directory under `build/fiber-volume-results/`. Set `FIBER_RESULTS_ROOT`
to choose another parent. The directory contains the command (`run.txt`),
model inspection and extraction logs, labeled boundary and tetrahedral meshes,
`volume_report.json`, source-grid membership counts and IoU (`source_grid.json`),
and the range triangle and boundary figure (`isabel_triangle.png` and `.pdf`).
The figure shows the boundary from above and below. The orange patch
corresponds to the constant-temperature edge BC and is visible from below.

## Curved-prism experiment

From the repository root, run

```bash
./src/fiber_surface/experiment/run_curved_prism.sh
```

From this folder, use `./run_curved_prism.sh`. The script loads the repository's
Spack environment, builds `fiber_volume_explicit`, extracts the geometry,
evaluates it against ground truth, and renders the figures. It uses the active
Python 3 if it has NumPy, SciPy, VTK and Matplotlib, then tries the Conda
environment `mfa_env`. Set `PYTHON` to select a particular interpreter.

The field is `F(x,y,z) = (x + z²/2, y + z²/2)` on `[-2,2]² × [-1,1]`.
The range triangle has vertices `A=(0,0)`, `B=(1,0)`, and `C=(0,1)`.
Both fiber tracing and triangle marching use `h=0.1` in domain coordinates.
The experiment sets `--seed-grid 10`, placing 10 seeds along the shortest axis
(`z`). The `x` and `y` axes are twice as long, so the lattice is 20 × 20 × 10.
The exact volume is 1.

### Outputs

The script writes to `build/fiber-volume-results/curved_prism_h010_seed10/` and
creates `curved_prism_h010_seed10.zip` alongside it. Set `FIBER_RESULTS_ROOT` to
choose a different parent directory.

| File | Contents |
| --- | --- |
| `boundary_fibers.pdf`, `boundary_fibers.png` | The three side patches with the corner fibers highlighted and the range triangle labeled |
| `volume_cutaway.pdf`, `volume_cutaway.png` | The tetrahedral volume cut at `z=0`, with the exact section overlaid |
| `curved-prism-teaser.pdf`, `curved-prism-teaser.png` | Range triangle, fiber surface with corner fibers, volume cutaway, and tetrahedral section at `z=0` in a single row |
| `corner_fibers.vtp`, `side_patches.vtp` | Traced fibers labeled A=0, B=1, C=2 and final side patches labeled AB=1, BC=2, CA=3 |
| `volume_volume.vtk`, `volume_tet_boundary.vtk` | The tetrahedral mesh and its exterior, including the two caps |
| `metrics.json` | Maximum vertex-to-segment distance in each direction between traced and reference fibers in unit-cube coordinates, patch residuals and projections, connectivity, enclosed volume, and cap sampling |
| `fiber_vertices.csv`, `patch_samples.csv` | Per-sample coordinates and measurements |
| `volume_report.json` | Extractor settings, timings and validation results |
| `run.txt` | Extractor command with the seed count and output prefix, for use from the repository root with the Spack environment loaded |

To recompute measurements and figures from existing meshes, run these commands
in this folder with the plotting Python environment active. Both scripts
also accept the results directory as a positional argument.

```bash
python3 evaluate_curved_prism.py
python3 plot_curved_prism.py
```
