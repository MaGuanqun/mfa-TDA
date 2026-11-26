# Reproduce the model-results table with a new CoordNet pair

Use one full-domain run to produce two rows, **Ours** and **TTK**, for the same pair of scalar fields and the same counterclockwise range triangle. The deliverable is `metrics/table.csv`. If you already have `stages/run.conf` and a labeled boundary PLY from your new model, start at **Measure both rows** after building the tools.

The saved Hurricane Isabel control helps identify the intended table and units.

| Model and method | Triangles | Open edges | Pinches | Intersections | $d_{95}$ | IoU |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Isabel MFA, Ours (reconnect and fit) | 650,490 | 0 | 0 | 0 | 0.00891 | 0.998 |
| Isabel MFA, TTK (native point merging) | 2,642,087 | 57 | 0 | 411 | 0.0118 | — |

On Isabel's `1000 × 1000 × 200` evaluation grid, the native mesh and the fitted field both selected 785,376 nodes. Another 728 were field-only and 693 were mesh-only, giving IoU $785376/786797$. For TTK, 87 degenerate faces were excluded from the intersection test. The 57 open edges are counted on the full mesh.

The 164,115 TTK triangles in the draft's vortex row came from the older upward-wake model, a `161 × 21 × 31` grid at spacing 0.05, and exact-coordinate welding. Recompute every cell for the new model. The procedure below uses the native TTK point merging used for the Isabel comparison. The older experiment measurements are in [model_metrics.csv](model_metrics.csv).

## Build and test the new machine

Run Bash from the repository root on branch `dev-no-batch`. Use the existing `mfa-TDA-env` Spack environment. For a first installation, follow the [repository setup](../../../README.md). Select Python executables available on this machine. One interpreter can serve all three roles.

```bash
# If Spack shell support is not initialized, first source your Spack's share/spack/setup-env.sh.
source ./load-env.sh &&
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb &&
cmake --build build -j 3 --target fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface \
    fiber_stitch_surface fiber_volume_from_patches fiber_grid_recovery fiber_repair_boundary \
    fiber_refine_saved_patch fiber_surface_coordnet_tests fiber_volume_coordnet &&
ctest --test-dir build -R 'fiber_(surface_(coordnet_adapter|grid_recovery|boundary_repair)|volume_coordnet_cli)$' --output-on-failure

export TORCH_PYTHON=python3    # PyTorch and NumPy
export TTK_PYTHON=python3      # topologytoolkit.ttkFiberSurface, VTK, SciPy and NumPy
export METRICS_PYTHON=python3  # full VTK with ImagingStencil, SciPy and NumPy
exp=src/fiber_surface/experiment
bash "$exp/test_coordnet_table.sh"
```

The smoke test creates two float64 TorchScript models in a fresh directory under `build/`. Their outputs are normalized $x$ and $y$, so the selected volume is a triangular prism. It checks that both surfaces project onto their planar sides, the native mesh is closed, and the field selects 16 nodes on each of 15 time slices. It also checks that the open TTK surface gets a reason for its missing IoU. On the tested VTK build the zero-tolerance stencil omitted one cap slice (224 of 240 selected nodes, IoU 0.93333). These are fixture values.

## Set up the model and extract Ours

Bring the two **float64 TorchScript exports** from training and the training metadata for each component. The scalar models accept a batch of normalized coordinates in `(t,y,x)` order and return one value per point. The extractor takes the training box in physical `(x,y,t)` order and applies the normalization and axis reversal. Use the actual exports' training box even if your source arrays use different time indices.

Specify three counterclockwise corners in **model-output units**, in the same order as the two files. If training independently mapped each component to $[-1,1]$, convert each physical range coordinate $v$ with that component's training extrema, $2(v-v_{min})/(v_{max}-v_{min})-1$. Models that already return physical values take the physical triangle directly. The min/max conversion applies to the query once.

Set the variables below for your run. The numeric settings shown are an example for the affine smoke model. The draft describes the Gray–Scott Bubbles box as `0,0,20` to `1,1,650`; confirm those bounds against your exports and choose the concentration triangle for your experiment. The historical vortex preset has different bounds, output scales and velocity corners.

```bash
MODEL_NAME='my CoordNet pair'
FIRST=/path/to/first-float64.pt
SECOND=/path/to/second-float64.pt
LO=0,0,0
HI=1,1,1
A=-0.42,-0.26
B=0.69,-0.18
C=-0.13,0.68
H=0.2                 # nominal step in extraction coordinates
SEEDS=3               # starts on the shortest axis
RANGE_TOL=0.00005     # tolerance in model-output units
DIMS=(11 13 15)       # TTK/evaluation grid nodes in x,y,t order, endpoints included
```

Choose the production `H`, seed count, triangle budget and comparison grid before recording results. Seed counts scale with the extraction-box aspect ratio. For the full Gray–Scott box, 3 starts on the shortest axis give `3 × 3 × 1890`, while 40 give 40,320,000 starts, above the staged extractor's ten-million limit. TTK grid dimensions are independent of seed counts. The lattice spacing on each axis is `(hi - lo)/(n - 1)`; choose dimensions that resolve the selection and record them with the rows.

```bash
source ./load-env.sh
exp=src/fiber_surface/experiment
source "$exp/run_triangle_stages.sh"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RUN=$(mktemp -d "$PWD/build/coordnet-table.XXXXXXXX")
echo "$RUN"
git rev-parse HEAD > "$RUN/revision.txt"
run_triangle_stages "$RUN" surface "$H" --model coordnet --first "$FIRST" --second "$SECOND" \
    --domain-min "$LO" --domain-max "$HI" --polyline "$A;$B;$C;$A" \
    --seed-grid "$SEEDS" --max-seed-starts 1632160 --max-triangles 500000 \
    --residual-tolerance "$RANGE_TOL" --fiber-field-tolerance 1e-9 \
    --surface-solve-tolerance 1e-10 --derivative-tolerance 1e-10 --max-normal-angle 85
MESH="$RUN/volume_boundary.ply"
```

The helper records the commands in `run.txt`, seed/trace/march/stitch checkpoints in `stages/`, per-side logs, and the assembled surface audit in `audit.json`. Check `extract.log` for unresolved events. A partial surface can still supply the geometry columns. For a row with the paper's final reconnect-and-fit repair, run the following against the original checkpoints and use the resulting candidate for **all six** columns.

```bash
source ./load-env.sh
"$METRICS_PYTHON" "$exp/repair_saved_fronts.py" "$RUN" "$RUN/repair" --max-edge-factor 4
# A saved final candidate may still be open (exit 2). Inspect repair.json and final_audit.json.
# Select this mesh only when final.ply was produced and this is the candidate being reported.
MESH="$RUN/repair/final.ply"
```

## Measure both rows

For an existing run, set `RUN` to its directory and `MESH` to the assembled or repaired PLY. Keep the checkpoints and the original model files at the paths recorded in `run.conf`. The sampler verifies their digests. Set `MODEL_NAME` and `DIMS` as above. Each command takes a **new output directory** so a different model, mesh or resolution gets its own results.

```bash
source ./load-env.sh
exp=src/fiber_surface/experiment
"$TORCH_PYTHON" "$exp/sample_coordnet_grid.py" "$RUN" "$RUN/grid" --dims "${DIMS[@]}"
"$METRICS_PYTHON" "$exp/evaluate_coordnet_table.py" "$RUN" "$RUN/grid" "$RUN/metrics" \
    --mesh "$MESH" --model-name "$MODEL_NAME" --ttk-python "$TTK_PYTHON"
cat "$RUN/metrics/table.csv"
```

The sampler writes two float64 fields (16 bytes per grid node) to `grid/fields.npy`, with $x$ varying fastest. TTK uses both fields on the same VTK image and receives AB, BC and CA together. It uses its native point merge at distance `1e-6` and preserves the resulting connectivity in `metrics/ttk.ply` (patch IDs 1, 2 and 3). `grid/grid.json` and `metrics/evaluation.json` record model paths and digests, training and extraction boxes, ordered triangle, grid dimensions and spacing, and software settings. To compare another native candidate against this grid, repeat only the evaluator with `--mesh` and a fresh output directory.

| Table column | Measurement and check |
| --- | --- |
| Triangles | `triangles` counts faces on the full exported candidate. The native boundary can include domain-face caps. TTK exports range sides. |
| Open edges | `open_edges` counts distinct undirected edges incident to exactly one face, including domain-boundary edges. |
| Pinches | `pinches` counts nonmanifold vertices from their incident triangle fans. Read `nonmanifold_vertices` in the full mesh audit. |
| Intersections | `intersections` counts conflicting triangle pairs. If degeneracies prevent the full intersection test, the evaluator reruns it on nondegenerate faces and records the excluded count. Use the **full** mesh for triangle, edge and pinch counts and projection. An unchecked intersection count is blank. |
| $d_{95}$ | `d95` is the 95th percentile of Euclidean displacement after projecting 20,000 points onto their corresponding **continuous CoordNet** sides. The audit uses seed 331 and samples across AB/BC/CA in proportion to physical area, excluding caps. Distances use the native $(x,y,t)$ coordinates. Isabel's units are grid units; state the coordinate units used for the new model. If any projection fails, the CSV leaves the all-sample percentile blank and records the failure count. |
| IoU | `iou` compares membership inside a closed surface (VTK stencil, tolerance 0) with membership in the closed range triangle at **every node of the TTK grid**. Read `both / (both + fitted_only + mesh_only)`. The four counts including `neither` must sum to `nx × ny × nt`. An open, nonmanifold, degenerate, or intersecting surface, an unchecked intersection test, or an empty union gets a blank and an `iou_note`. |

Use the same grid and query for both rows. A measured zero stays zero. Render an unavailable value as an em dash and keep the reason with the results. Per-patch nearest distances to TTK and random-point comparisons with a piecewise-linear field answer different questions from these two table columns.

## What travels with Git

The extraction stages, repair, TTK wrapper, grid sampler, evaluator, smoke test and this guide live in the code repository. A clean checkout can generate the two rows without the separate manuscript repository or any old `build/` tree. The new trained weights and their scientific choices (component order, units, training extrema, box and range triangle) come from your model run. Archive those with `revision.txt`, `run.txt`, `stages/run.conf`, the chosen PLY, `grid.json`, and the contents of `metrics/`. Keep full-precision CSV values for the paper and round only for display.
