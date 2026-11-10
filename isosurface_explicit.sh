#!/usr/bin/env bash

set -e

# Edit the parameters here, then run: ./isosurface_explicit.sh
function_name="ellipsoid"
isovalue="1.0"
# Initial target edge length = shortest domain dimension / (20 * resolution).
resolution="1.0"
min_step_ratio="0.5"
samples="20"
tolerance="1e-8"
iterations="80"
max_triangles="200000"
stop_distance="0.05"

# Optional mesher parameters: uncomment the options you want to use.
mesh_options=(
    # --allow-open
    # --no-crack-closing
    # -d /path/to/singular_points.dat
)

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Executables from the main project build.
binary_dir="${script_dir}/build/src/isosurface"
root_finder="${binary_dir}/root_finding_explicit"
mesher="${binary_dir}/isosurface_explicit"
output_dir="${binary_dir}/results/${function_name}"
root_file="roots.dat"
mesh_prefix="${function_name}_sheet"

for executable in "$root_finder" "$mesher"; do
    if [[ ! -x "$executable" ]]; then
        printf 'Executable not found: %s\nBuild the project first; see src/isosurface/README.md.\n' "$executable" >&2
        exit 1
    fi
done

mkdir -p "$output_dir"
# Use local filenames because the C++ argument parser rejects paths with spaces.
cd "$output_dir"

"$root_finder" \
    -f "$function_name" \
    -v "$isovalue" \
    -g "$resolution" \
    -x "$tolerance" \
    -m "$iterations" \
    --samples "$samples" \
    -s "$root_file"

"$mesher" \
    -f "$function_name" \
    -v "$isovalue" \
    -g "$resolution" \
    -n "$min_step_ratio" \
    -x "$tolerance" \
    -m "$iterations" \
    -c "$stop_distance" \
    --max-triangles "$max_triangles" \
    -s "$root_file" \
    -o "$mesh_prefix" \
    "${mesh_options[@]}"

printf '\nOutput directory: %s\n' "$output_dir"
