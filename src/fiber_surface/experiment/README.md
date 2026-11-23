# Fiber volume experiments

## Hurricane Isabel at time step 30

Run `./src/fiber_surface/experiment/run_hurricane_isabel.sh` for the default
`narrow` experiment. Its pressure–temperature triangle is
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

From the repository root, run both current Isabel cases with

```bash
./src/fiber_surface/experiment/run_hurricane_isabel_ttk.sh
```

Pass `narrow` or `wedge` to run one case. The script uses the newest completed
MFA run for each named preset and runs the existing MFA experiment if one is
missing. To compare particular meshes, pass their run directories as arguments.
For example, after running the `1200` and `1300` cases from `fiber_volume.sh`,
use `./src/fiber_surface/experiment/run_hurricane_isabel_ttk.sh build/fiber-volume-results/hurricane_isabel_3d_full/{1200,1300}`.
Each supplied mesh must come from the same time-30 pressure and temperature models.
Its `volume_report.json` supplies the three range vertices in order.

The runner samples `Pf30_3d.mfa` and `TCf30_3d.mfa` once at all 25,000,000
integer locations in `[0,499] × [0,499] × [0,99]`. The two float64 sample
files take 400 MB combined. Both cases use the same samples. TTK processes
the full image with an implicit triangulation of 147,906,594 tetrahedra.
The reference tests how well the continuous-model extraction agrees with the
piecewise-linear interpolation of the **fitted MFA fields**. The source-grid
IoU above uses the original simulation values.

Each invocation creates a fresh `hurricane_isabel_t30_ttk.*` folder under
`FIBER_RESULTS_ROOT` (default `build/fiber-volume-results/`). Its per-run
folders contain the following results.

| File | Contents |
| --- | --- |
| `sides.vtp` | TTK triangles labeled `patch` 1 = AB, 2 = BC, 3 = CA. Open alongside the MFA run's `volume_tet_boundary.vtk` to overlay matching patches. |
| `fibers.vtp` | TTK line segments labeled `corner_id` 0 = A, 1 = B, 2 = C, when corner fibers exist. |
| `comparison.json` | Query, sample grid, model and mesh paths, counts, physical patch areas and seam lengths, and closest-point distances in both directions. |

Distances use 10,000 seeded samples per patch or corner seam in each direction.
Each domain axis is scaled to `[0,1]` before sampling by triangle area or line
length. The report gives the median, 95th percentile, and sampled maximum.
Physical areas and lengths use the original grid coordinates. The MFA
experiment's reports supply the volume and membership measurements.

| Case | MFA side triangles | TTK side triangles | TTK corner segments |
| --- | ---: | ---: | ---: |
| `narrow` | 5,234 | 23,632 | 2,001 |
| `wedge` | 9,424 | 40,761 | 2,160 |

The sampler builds in the repository's Spack environment. The TTK stage uses
`TTK_PYTHON`, the active Python 3, or the Conda environment `mfa`, with TTK
Python bindings, VTK, NumPy and SciPy. `TTK_THREADS` sets the thread count
for sampling and extraction (default 3). To reuse an existing sample pair,
run `python ttk_hurricane_isabel.py SAMPLES_DIR MFA_RUN_DIR [MFA_RUN_DIR ...]`
from this folder with that Python environment active.

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
