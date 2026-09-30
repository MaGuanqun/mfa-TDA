# Adaptive marching triangles: explicit scalar fields

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
constructor. The explicit path is implemented and tested here; MFA/INR adapters
and their derivatives have not been validated by this work. The return value of
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
