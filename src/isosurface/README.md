# Adaptive marching triangles: explicit functions and MFA models

This directory implements the two-stage algorithm in Akkouche and Galin,
*Adaptive Implicit Surface Polygonization using Marching Triangles*, Computer
Graphics Forum 20(2), 2001, pp. 67–80. The explicit example and root finder build
with C++17 and Eigen; they do not require MPI, MFA, TBB, or Torch.

## Build and run

The isosurface executables build with the other enabled project folders. Load
the Spack environment in your current shell, then configure and compile from
the repository root:

```sh
source ./load-env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 3
ctest --test-dir build -R '^isosurface_' --output-on-failure
./isosurface_explicit.sh
```

This produces `root_finding_explicit` and `isosurface_explicit` in
`build/src/isosurface`, which is also the run script's default `binary_dir`.
To rebuild only these executables, use
`cmake --build build --target root_finding_explicit isosurface_explicit -j 3`.
`cmake --install build` installs them under `<install-prefix>/src/isosurface`,
matching the other project folders. Tests are enabled by default; configure
with `-DBUILD_TESTING=OFF` to omit them.

Use `source ./load-env.sh`, rather than executing it in a separate process, so
its exported dependency paths remain available to CMake and the executables.
Set the function name, isovalue, numerical parameters, executable paths, and
output location directly in `isosurface_explicit.sh`.

The script executes the compiled root finder and mesher without rebuilding.
Roots and meshes are saved in `build/src/isosurface/results/FUNCTION/` by default.
Uncomment entries in its `mesh_options` array to allow open surfaces, disable
crack closing, or supply singular points. The script can also be launched using
its full path.

An optional standalone build of the explicit code still uses only C++17 and
Eigen. Use a separate build directory:

```sh
cmake -S src/isosurface -B build/isosurface -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/isosurface -j 3
ctest --test-dir build/isosurface --output-on-failure

build/isosurface/root_finding_explicit -f ellipsoid -v 1 -s build/isosurface/roots.dat
build/isosurface/isosurface_explicit -f ellipsoid -v 1 \
    -s build/isosurface/roots.dat -o build/isosurface/ellipsoid_sheet
```

The last command writes `ellipsoid_sheet0.obj`, after validating it. To use this
optional build with the run script, set its `binary_dir` to
`${script_dir}/build/isosurface`. Standalone installation places executables in
`<install-prefix>/bin`.

A single seed can be supplied in a text file containing `1 0 0`. Binary root files
from the previous `utility::writeMatrixVector` format are supported: one matrix,
three columns, native `size_t` dimensions, column-major doubles. The root finder
now writes **exactly the `-s` path**. It no longer appends a resolution and a second
`.dat` extension. Text roots may contain xyz rows, OBJ `v` rows, and `#` comments.

The ellipsoid field is `x² + y²/4 + z²/9`. Its level **1**, not level 0, is a
surface. The explicit adapter computes a padded domain enclosing the requested
positive level; the previous `[-2,2]³` domain cut off the level-1 ellipsoid at its
z poles. These bounds are local to this example and do not change the shared
analytic functions used by other applications.

Useful options:

| Option | Meaning |
| --- | --- |
| `-g` | Resolution divisor; initial target edge length is `min(domain_max - domain_min) / (20 * g)`. |
| `-n` | `d_min / seed_edge_length`, default `0.1`. A value around `0.5` can prevent excessive shrinking on tubular surfaces. |
| `-x` | Absolute scalar residual tolerance, default `1e-8`. |
| `-m` | Projection iteration budget, default 80. |
| `-d`, `-c` | Optional singular-point stop set and its exclusion distance. |
| `--max-triangles` | Per-sheet limit; reaching it with an open front is reported. |
| `--no-crack-closing` | Debug the growing mesh before closing cracks. |
| `--quality-iterations` | Closed-surface quality improvement passes, default `8`; `0` disables it; maximum `100`. |
| `--allow-open` | Explicitly permit a valid mesh with boundary. All other validation remains enabled. |

The initial target edge length is one twentieth of the shortest dimension of
the final sampling domain, divided by `-g` (default `1`). Root separation uses
the same length. For the level-1 ellipsoid, the padded domain spans
`2.2 x 4.4 x 6.6`, giving a target seed edge of `0.11`; `-g 2` gives `0.055`.
Projection onto the surface can slightly change the actual seed edges. This
shortest-dimension convention is an implementation choice inspired by the paper's
relative seed-size experiment.

Crack closing is enabled by default. The old `-r` flag remains accepted. The old
Hessian threshold `-H` is accepted for compatibility but no longer determines
regularity. Non-default `-k` shrink ranges are rejected rather than silently
ignored. `quartic_potential_3d` has a four-dimensional domain and is rejected by
this three-dimensional triangle extractor. The other existing three-dimensional
quartic fields remain available, but their chosen level sets can be singular or
intersect the domain boundary; such surfaces need not close.

## Critical curves

The existing CMake target directory also builds three separate critical-curve
tools. These produce polylines; attaching triangles to the curves is future work.
The existing regular marching-triangles pipeline and critical-point extractor
are unchanged.

```sh
source ./load-env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --target critical_curves_explicit critical_values_mfa critical_curves_mfa \
    critical_curves_tests critical_curves_mfa_tests -j 2
ctest --test-dir build -R '^critical_curves_' --output-on-failure
```

### Explicit functions: choose a known value and trace

The examples live in `critical_explicit_field.h`. Their field offset `b` defaults
to **1.25** and is independent of `--iso`. The domain is `[-2,2]^3`.

| `--function` | Field | Critical crossing at value `b` |
| --- | --- | --- |
| `crossing_line` | `xy + b` | `x=y=0`, ending on the two z faces |
| `crossing_circle` (default) | `z(x²+y²−1) + b` | `z=0`, `x²+y²=1` |
| `parallel_lines` | `z(z−y²+1/4) + b` | Two lines: `z=0`, `y=±1/2` |

```sh
build/src/isosurface/critical_curves_explicit -f crossing_circle -v 1.25 \
    --samples 9 --step 0.05 -o build/src/isosurface/results/critical_curves/circle
build/src/isosurface/critical_curves_explicit -f crossing_line -v 1.25 \
    -o build/src/isosurface/results/critical_curves/line
```

No critical-value search is performed for these examples: seeds are found at the
requested level, then traced. Omitting `--iso` uses the known offset value.
For example, `--offset 2 --iso 2` moves the critical value to 2, while
`--offset 2 --iso 1.25` should produce no curves. The analytic curve positions
are used only for testing, not provided to the search as known answers.

Negative controls are `definite_line` (`x²+y²+b`), `isolated_rank_two`
(`x²−y²+z⁴+b`), and `tilted_line` (`xy+0.01z+b`). They respectively test rejection
of definite normal Hessians, isolated rank-two critical points, and nonzero
tangential gradients.

### MFA: discover values, select one, then trace

Run value discovery without providing an isovalue:

```sh
build/src/isosurface/critical_values_mfa \
    -f build/src/boussinesq_3d_iso/boussinesq_3d.mfa --samples 9 \
    -o build/src/isosurface/results/critical_curves/boussinesq_s9
```

Inspect `boussinesq_s9.values.csv`. Each row has a candidate ID, estimated value,
minimum/maximum seed values, and seed count. Its positions are saved separately
as `boussinesq_s9.value_ID.seeds.txt`. The catalog contains **candidate** critical
values: pointwise gradient/Hessian checks cannot prove a critical curve exists.
An empty catalog is a valid result, not proof that the model has no curves.
Use the catalog from the current run; seed files from older runs with the same
prefix may still exist.

Given an actual row, substitute its value and ID in the following command:

```sh
build/src/isosurface/critical_curves_mfa \
    -f build/src/boussinesq_3d_iso/boussinesq_3d.mfa -v VALUE \
    -s build/src/isosurface/results/critical_curves/boussinesq_s9.value_ID.seeds.txt \
    -o build/src/isosurface/results/critical_curves/boussinesq_selected
```

The tracer revalidates and refines supplied seeds against the current model and
selected value. Omit `-s` to perform a fresh grid search at that value. Different
disconnected curves can share one value; value grouping does not join them.

`critical_curves_mfa_tests` creates an exact polynomial MFA fixture at
`critical_curve_test_files/circle.mfa` relative to its working directory. Its
known curve is the unit circle at **0.37**, providing a positive MFA example even
when a real dataset has no critical curves. Tests verify the full
discovery → seed-file → tracing CLI workflow, analytic/mixed Hessians, and
physical derivative scaling on unequal axis extents.

The supplied `boussinesq_3d.mfa` (variable 0) was tested with the default residual
and Hessian tolerances. It has cubic degree in each axis, `13×33×13` control
points, and domain `[-0.5,0.5] × [-0.5,2.5] × [0,1.5]`. These searches found:

| Seed grid | Cell offset | Starts | Full-rank convergences | Failed solves | Curve-candidate seeds |
| --- | --- | --- | --- | --- | --- |
| `9³` | `0.5` | 729 | 199 | 530 | 0 |
| `17³` | `0.5` | 4,913 | 1,448 | 3,465 | 0 |
| `17³` | `0.25` | 4,913 | 1,427 | 3,486 | 0 |

Convergences count solver starts, not distinct critical points. All converged
starts had resolved full-rank Hessians, so none qualified as a crossing-curve
seed. There was no discovered value to pass to the tracer for this dataset.
The failed solves and finite sampling prevent a claim of global absence.
Reports are under `build/src/isosurface/results/critical_curves/boussinesq_*`.

### Method, tolerances, and output

`critical_field.h` contains seed discovery and Hessian classification.
`critical_curves.h` contains continuation. Both backends use the same routines.
Discovery solves the gradient residual with a damped augmented SVD
Levenberg–Marquardt step; it never rejects a point merely because the Hessian
determinant is small. When a value is supplied, the scalar residual `F-c` is
included as well. A candidate must have small gradient and exactly one unresolved
Hessian eigenvalue, with opposite signs on the two resolved eigenvalues.

Continuation predicts along the Hessian null direction and corrects in a fixed
normal plane using a `3x2` SVD solve involving **all three** gradient components.
It checks the selected value, full gradient, Hessian signature, correction
distance, and tangent change at every accepted sample. It traces both directions,
detects closed loops, and corrects transverse endpoints on box faces. Rank loss,
failed corrections, point limits, and ambiguous self-contacts are recorded.
Short failed traces are reported as unvalidated events rather than accepted
curve components. Partial longer curves remain explicitly labeled in the output.

All positions and `--step`/`--min-step` are in physical coordinates. Numerical
residuals are normalized using `L = longest box side` and a fixed field scale
`A`. By default, a fixed 5×5×5 grid estimates
`A = max(sampled value range, L * maximum sampled gradient norm)`; a constant
field uses 1. Set `--field-scale` to supply a reference scale. Consequently:

- `--value-tolerance` bounds `|F-c|/A` and the width of a candidate value group.
- `--gradient-tolerance` bounds `L*|gradient F|/A`.
- `--hessian-absolute` and `--hessian-relative` set the rank threshold
  `absolute + relative * spectral_norm(L² H/A)`.
- Output eigenvalues and gradient norms use these normalized units; field values
  use the original units. The report records both scales and all settings.

Use `--samples` and `--seed-offset` to repeat discovery at different densities and
offsets. Seeding, eigenvalue thresholds, seed deduplication, and branch coverage
are resolution dependent; none certifies global completeness. A fitted MFA may
also lose an exact crossing, even when scalar fitting errors are small.

The MFA adapter evaluates six analytic second partial derivatives, applying
`H_x(i,j) = H_u(i,j)/(extent[i]*extent[j])`. It requires the existing Cartesian,
single-tensor-product, unweighted scalar model and verifies at least C2 spline
continuity across interior knots in the selected block. Endpoint knot repetition
does not reduce interior smoothness. As with the regular MFA tools, run one MPI
rank and select a block explicitly for multiblock files. Junction merging,
cross-block stitching, higher-order contacts, and surface-to-curve stitching
are not implemented in this stage.

For a tracing output prefix `P`, the files are:

- `P.obj` and `P.vtk`: indexed curve polylines (VTK opens in ParaView).
- `P.points.csv`: coordinates, original values, scaled gradient norms, and Hessian eigenvalues.
- `P.curves.csv`: per-component endpoint/closure status and residual summaries.
- `P.events.csv`: unresolved locations and reasons.
- `P.report.txt`: configuration, search counts, evaluation counts, and tracing summary.

Discovery writes `.values.csv`, `.seeds.csv`, per-value `.seeds.txt` files, and
`.report.txt`. Exit code 0 means the search ran successfully, including empty
or partial results; inspect the reports for curve counts and unresolved events.
Invalid input or an I/O failure returns a nonzero exit code.

## MFA models

The main project build also creates `root_finding_mfa` and `isosurface_mfa` in
`build/src/isosurface`. These targets use MFA, DIY/MPI, Eigen, and the project's
threading libraries; they do not link Torch. The optional standalone build
continues to build only the explicit programs.

```sh
source ./load-env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 3
./isosurface_mfa.sh
```

Edit `model_file`, `isovalue`, `variable`, mesh/seed resolution, and the other settings in
`isosurface_mfa.sh`. Its defaults use
`build/src/boussinesq_3d_iso/boussinesq_3d.mfa`, variable 0, and level 0.5, with
`--allow-open` because this level set intersects the domain boundary. Output is saved in
`build/src/isosurface/results/boussinesq_3d/iso_0.5/`. Run these executables directly
with one MPI rank; each process writes the whole selected block.

The equivalent commands are:

```sh
build/src/isosurface/root_finding_mfa -f build/src/boussinesq_3d_iso/boussinesq_3d.mfa -v 0.5 \
    -g 2 --samples 20 -s roots.dat
build/src/isosurface/isosurface_mfa -f build/src/boussinesq_3d_iso/boussinesq_3d.mfa -v 0.5 \
    -g 1 -n 0.5 --allow-open -s roots.dat -o boussinesq_3d_sheet
```

The MFA root finder orders seeds by distance from the domain boundary, with
interior seeds first. This reduces tiny duplicate patches caused by starting at
a clipped edge. The example uses seed resolution 2 (separation 0.025) and mesh
resolution 1 (initial edge length 0.05); these can be adjusted independently.

The MFA programs retain the explicit programs' root-file format and numerical
options. `-f` selects the MFA file, `--variable` selects a zero-based scalar
science variable, and `--block` selects a zero-based block index when the input
contains multiple blocks. A multiblock input requires explicit selection;
stitching between blocks is not implemented. Input paths must not contain
spaces because of the repository's existing argument parser.

`mfa_model.h` uses the project's standard DIY `Block<double>::load` reader.
Only files written in the current installed MFA format are supported; there is
no old-format detection or conversion. Extraction evaluates the spline model
directly. Files use the native double-precision serialization used by this
project.

`mfa_field.h` adapts a model to the same `FieldSetup<double>` and
`SurfaceField<double>` interface used by `explicit_field.h`. It checks for
Cartesian affine geometry, scalar unweighted variables, and one tensor product
per model. Non-Cartesian, rational/weighted, and hierarchical models are rejected
explicitly. The geometry maps physical coordinates to parameter coordinates;
MFA derivatives are divided by the geometry extent to give physical gradients.
Out-of-domain queries return non-finite values for the safeguarded projector
to reject. Reusable decoder buffers make this adapter serial; use separate
adapters for independent threads.

Both paths share `find_initial_root.h`, `surface_field.h`,
`marching_triangles.h`, crack closing, mesh quality, and mesh validation.
`extract_isosurface.h` also shares the CLI extraction and output checks,
including scalar residuals and domain bounds. Open-front seed coverage checks
a short projected connection to the existing front, avoiding duplicate
patches from seeds just outside its triangle footprints.

### Boussinesq experiment

The example uses the current-format model
`build/src/boussinesq_3d_iso/boussinesq_3d.mfa` at level 0.5. Boundary sampling
finds values on both sides of 0.5 on the faces `x = 0.5`, `z = 0`, and `z = 1.5`,
so the example enables `--allow-open`. With the settings above, the experiment
produced:

| Check | Result |
| --- | --- |
| Seeds | 1,215, ordered from the domain interior toward its boundary |
| Components | 8 (one closed, seven open) |
| Vertices / triangles | 6,815 / 13,498 |
| Boundary edges | 134 |
| Maximum vertex scalar residual | 9.99925e-9 |
| VTK nonmanifold edges / nonadjacent triangle intersections | 0 / 0 |
| Gradient finite-difference relative error | 1.49318e-9 |

The VTK check includes intersections between different output files. With the
optional Python `vtk` package installed, reproduce it using:

```sh
python src/isosurface/tests/verify_obj_vtk.py --allow-open --components 8 \
    build/src/isosurface/results/boussinesq_3d/iso_0.5/boussinesq_3d_sheet*.obj
```

A separate VTK contour of MFA samples on a `65 x 193 x 97` grid also has eight
components. Its area is about 2.50834, compared with 2.47242 for the direct mesh.
This sampled comparison does not establish complete coverage or exact topology.

The mesher does not compute exact clipped boundary curves or cap them. Additional
gaps and thin triangles remain: this run has 907 triangles with an angle below
5 degrees, and a maximum face-normal error of about 89.7 degrees. Quality
optimization currently applies only to closed sheets. Increasing resolution can
make crack closing substantially slower; disabling it leaves many more gaps.
Seed ordering improves this example but does not guarantee duplicate-free
extraction for every model, so retain independent validation when changing inputs.

`ctest --test-dir build -R '^isosurface_' --output-on-failure` runs the explicit
suite plus a clipped-plane coverage regression and a current-format MFA fixture.
The MFA fixture checks analytic scalar values, anisotropic physical gradients,
finite differences, domain guards, and a closed ellipsoid extraction. To check
derivatives against finite differences on another model:

```sh
build/src/isosurface/isosurface_mfa_tests build/src/boussinesq_3d_iso/boussinesq_3d.mfa
```

## Paper-to-code correspondence

### Growth (sections 3.1 and 3.2.1; Figures 1–6)

1. Project a seed onto `f(x) = isovalue`, construct a tangent frame, and project a
   **60-degree**, initially equilateral triangle.
2. Store the face's three directed boundary halfedges. For a front edge `(a,b)`,
   the adjacent new face must use the opposite direction: `(a,p,b)`.
3. Project the edge midpoint onto the surface to obtain `x_s`.
4. Predict along `t = normalize((b-a) × grad f(x_s))` with
   `d = sqrt(3)/2 * (previous_length + edge_length + next_length)/3`.
5. If `d < d_min`, apply the paper's blend `d = 0.75*d + 0.25*d_min`.
6. Project the prediction using the paper's `alpha = 1.5` gradient march until a
   root is bracketed, then bisect. Backtracking, finite-value checks, displacement
   limits, and a final domain check safeguard the numerical iteration.
7. Apply the orientation-dependent, relaxed empty-sphere test globally, not just
   to triangles incident on the candidate's vertices. Disjoint triangle interiors
   are tested against the sphere. Shared edges necessarily intersect that sphere;
   adjacent faces use non-shared vertex tests and direct intersection checks.
8. If the prediction is rejected, try the predecessor and successor wedges. A
   candidate is committed only after all checks pass. Process newly created edges
   through a FIFO queue; an unsuccessful edge remains available for crack closing.

The spatial step and edge-length caps are practical guards, not additional
curvature formulas from the paper. In particular, `d_min` is a blended target,
not a hard lower bound on every edge. Very small values can produce many faces.

### Crack closing (section 3.2.2; Figures 7–8)

The implementation processes **all** remaining contours. It splits long contours
at nearby facing boundary vertices, triangulates small contours through validated
ears, and tries joining facing contours when they cannot be meshed separately.
It adds no new surface vertices during crack closing.

Boundary edges are derived from face incidence. Halfedge walks that revisit a
vertex after a split are decomposed into the smaller contours shown in Figure 8.
Failed insertions never remove contours or consume edges. Blind triangle fans
are not used for nonplanar or concave holes.

There are additional safeguards beyond the paper's outline:

- Every insertion checks undirected edge incidence and directed winding, rejects
  zero-area and duplicate faces (including reversed duplicates), and tests actual
  triangle intersections, including coplanar overlaps and shared-vertex cases.
- Local projection into a surface tangent plane prevents two chord patches from
  covering the same region without physically intersecting in 3D. Oppositely
  facing, nearby parts are permitted, as in Figure 2.
- Greedy splitting can leave a sliver that cannot be capped without an inverted
  face. In that case, bounded local cavities are retriangulated with different
  diagonals. A trial is retained only when it reduces the boundary-edge count and
  preserves all vertices. This local repair is an extension, not a step specified
  by the paper.
- The tangent-plane guard can be conservative for a final triangular hole. Its
  fixed-edge cap may use the paper's small-hole rule while still satisfying 3D
  intersection, winding, nondegeneracy, and surface-projection checks.

`mesh_validation.h` reconstructs connectivity from the final vertex/face arrays.
The CLI checks boundary edges, edge orientation, nonmanifold edges and vertex
links, repeated faces, degenerate faces, invalid indices, unused vertices, and
triangle intersections before writing any OBJ. An incomplete result produces a
nonzero exit code unless `--allow-open` was explicitly selected.

### Triangle quality after closure

Topology and nonzero area alone do not ensure useful triangles. In the ellipsoid
at the previous default seed edge length of `0.4`, greedy crack closing produced
extremely small altitudes, and a positive gradient dot product still allowed face
normals almost 90 degrees
away from the surface normal. Of 554 triangles with an angle below 5 degrees,
496 were introduced by crack closing or its local repairs.

`mesh_quality.h` now improves closed sheets with diagonal flips and tangential
vertex relaxation followed by projection onto the implicit surface. Each accepted
change improves the worst local score, combining scale-independent triangle shape
and alignment with the field gradient at the vertices and centroid. Changes must
preserve winding, pass intersection checks against the mesh, and satisfy the
existing local surface-projection and edge-length guards. Vertex and face counts
stay fixed. This optional postprocessing is an extension to the paper's two
stages. It is skipped for open sheets and can be disabled with
`--quality-iterations 0`.

The CLI reports minimum angle, number of faces below 5 degrees, minimum and mean
shape quality (1 is equilateral), maximum normal error, and reversed faces. It
rejects reversed faces and reports when very thin or poorly aligned faces remain.
A finite number of improvement passes does not guarantee a minimum angle for
arbitrary surfaces; inspect these diagnostics when changing functions or settings.

### Components and limits

A regular point is determined by its **gradient**, not the rank of its Hessian.
The mesher processes uncovered seeds instead of unconditionally selecting seed
zero. Coverage uses nearby triangle patches and surface projection, rather than
requiring a seed to coincide with a mesh vertex. Disconnected components still
need seeds. Uniform root sampling does not prove that every component was found.

The paper itself notes in section 3.2.2 that its small-contour/tubular-contour
heuristic may be topologically consistent without being topologically correct.
These checks and regression cases do not establish a universal topology guarantee
for arbitrary scalar fields, singularities, sharp features, or insufficient
resolution. A clipped surface is intentionally open; this implementation does
not add artificial domain caps. The incremental skeletal-model editor from
section 4 is not implemented here.

The numerical and topology core accepts `SurfaceField<T>` callbacks for scalar
values and gradients, plus `MeshingOptions<T>`. This replaces the old constructor
that directly depended on MFA/INR. No other repository caller used that
constructor. Both the explicit path and the MFA adapter above are implemented
and tested; an INR adapter is not provided here. The return value of
`extract_all_sheets` means at least one sheet was generated; library callers must
also inspect `diagnostic()` and run `validate_mesh` when closure is required.

The current global collision checks prioritize correctness and are quadratic in
mesh size; local repair adds further work. Large meshes would benefit from a
spatial acceleration structure. The explicit executables are serial: run one
process, not multiple MPI ranks writing the same output names.

## Bugs repaired and verification

The old implementation had independent problems in successor-wedge updates,
three-edge closure, triangle winding and duplicate keys, premature vertex
insertion, global collision detection, contour splitting/merging, crack-closing
failure handling, component iteration, numerical projection, default isovalue,
domain clipping, root filenames, and error reporting. Its packed triangle key
also assumed vertex IDs fit into 21 bits. The new edge and face keys do not have
that packing assumption.

A reproduction of the original mesher with a corrected enclosing ellipsoid box,
step `0.2`, crack closing enabled, and a 3,000-growth-triangle limit ended with
**1,988 boundary edges, 604 nonmanifold edges, and 309 unused vertices**. This is a
bounded failure reproduction, not a performance comparison.

The regression suite covers:

- A sphere and reversed field orientation.
- Ellipsoids with seed edge lengths `0.6`, `0.4`, and `0.25`.
- A torus (`Euler = 0`, preserving its handle).
- Two disconnected spheres and closely spaced concentric shells.
- Coordinate scales from radius `0.001` to `1000`.
- An unoptimized ellipsoid with slivers, followed by quality improvement without
  changing vertex/face counts or topology; cached edge incidence is checked against
  independently reconstructed connectivity.
- Shared-edge/vertex intersections, coplanar overlap, winding, projection failure,
  invalid parameters, binary root round trips, and truncated files.
- The root-finder → OBJ pipeline, unsupported 4D inputs, invalid isovalues,
  intentionally open output, triangle limits, and output I/O failures.

All closed-surface cases require zero boundary/nonmanifold edges, zero
nonmanifold vertices, consistent winding, no duplicates or unused vertices, no
detected triangle intersections, the expected Euler characteristic, and vertex
residuals within tolerance. These regular-surface fixtures also require minimum
angles above 5 degrees and maximum normal errors below 45 degrees. The dedicated
ellipsoid quality regression uses stronger bounds of 10 and 35 degrees.

For the generated-root ellipsoid example at seed edge length `0.4` (the previous
default), the checked output contains 3,417 vertices and 6,830 faces. Independent
VTK validation found one connected
component, zero boundary/nonmanifold edges, and zero intersections between
nonadjacent triangles. With the default eight quality passes, its volume is
25.0022 versus the analytic `8*pi = 25.1327` (about 0.52% low, as expected from
planar chord faces). The minimum angle improves from 0.0041 to 18.687 degrees,
the number of triangles below 5 degrees falls from 554 to zero, and the maximum
normal error falls from 89.934 to 21.229 degrees. Normal error compares the face
normal with the analytic gradient at all three vertices and the face centroid.

With the current shortest-domain-dimension default (`0.11` at level 1), the
generated-root ellipsoid has 11,397 vertices and 22,790 faces. It is closed with
no detected intersections or reversed faces, a minimum angle of 17.929 degrees,
zero faces below 5 degrees, and a maximum normal error of 12.986 degrees.
