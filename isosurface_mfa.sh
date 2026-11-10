#!/usr/bin/env bash

set -e

# Edit the parameters here, then run: ./isosurface_mfa.sh
model_name="boussinesq_3d"
isovalue="0.5"
variable="0"
# Initial target edge length = shortest domain dimension / (20 * resolution).
resolution="1.0"
# Keep seed sampling denser than the initial mesh edges.
seed_resolution="2.0"
min_step_ratio="0.5"
samples="20"
tolerance="1e-8"
iterations="80"
max_triangles="200000"
stop_distance="0.05"

# Optional mesher parameters: uncomment the options you want to use.
mesh_options=(
    --allow-open  # Boussinesq at level 0.5 intersects the model domain boundary.
    # --no-crack-closing
    # -d /path/to/singular_points.dat
)

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Executables from the main project build.
binary_dir="${script_dir}/build/src/isosurface"
model_file="${script_dir}/build/src/boussinesq_3d_iso/boussinesq_3d.mfa"
root_finder="${binary_dir}/root_finding_mfa"
mesher="${binary_dir}/isosurface_mfa"
output_dir="${binary_dir}/results/${model_name}/iso_${isovalue}"
root_file="roots.dat"
mesh_prefix="${model_name}_sheet"

for executable in "$root_finder" "$mesher"; do
    if [[ ! -x "$executable" ]]; then
        printf 'Executable not found: %s\nBuild the project first; see src/isosurface/README.md.\n' "$executable" >&2
        exit 1
    fi
done

if [[ ! -f "$model_file" ]]; then
    printf 'MFA model not found: %s\n' "$model_file" >&2
    exit 1
fi

mkdir -p "$output_dir"
# Use local filenames because the C++ argument parser rejects paths with spaces.
cd "$output_dir"

"$root_finder" \
    -f "$model_file" \
    --variable "$variable" \
    -v "$isovalue" \
    -g "$seed_resolution" \
    -x "$tolerance" \
    -m "$iterations" \
    --samples "$samples" \
    -s "$root_file"

"$mesher" \
    -f "$model_file" \
    --variable "$variable" \
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
