#!/usr/bin/env bash
set -eo pipefail

if (( $# )); then
    printf 'Usage: %s\nFIBER_RESULTS_ROOT sets the output root. PYTHON selects the plotting interpreter.\n' "$0"
    [[ $# == 1 && ( $1 == --help || $1 == -h ) ]]
    exit
fi

experiment_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$experiment_dir/../../.." && pwd)
results_root=${FIBER_RESULTS_ROOT:-"$repo_root/build/fiber-volume-results"}
mkdir -p -- "$results_root"
results_root=$(cd -- "$results_root" && pwd)
seed_grid=10
name=curved_prism_h010_seed${seed_grid}
work=$(mktemp -d "$results_root/.curved_prism.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
run_dir="$work/$name"
mkdir "$run_dir"

# Keep Spack's build libraries separate from the plotting interpreter's VTK.
export MPLCONFIGDIR=${MPLCONFIGDIR:-"$work/matplotlib"}
python_cmd=("${PYTHON:-python3}")
imports='import numpy, scipy, vtk, matplotlib.pyplot'
if ! "${python_cmd[@]}" -c "$imports" 2>/dev/null; then
    if [[ -n ${PYTHON:-} ]] || ! command -v conda >/dev/null; then
        printf 'Set PYTHON to an interpreter with NumPy, SciPy, VTK and Matplotlib.\n' >&2
        exit 1
    fi
    python_cmd=(conda run --no-capture-output -n mfa_env python)
    "${python_cmd[@]}" -c "$imports"
fi

(
    cd "$repo_root"
    if [[ $(type -t spack) != function ]]; then
        source "${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"
    fi
    source ./load-env.sh || exit 1
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
    cmake --build build --target fiber_volume_explicit -j 3
    extractor=(build/src/fiber_surface/fiber_volume_explicit --example curved-prism
        --q0 0,0 --q1 1,0 --q2 0,1 --step 0.1
        --spatial-tolerance 0.015 --residual-tolerance 0.005
        --seed-grid "$seed_grid" --max-triangles 30000)
    "${extractor[@]}" --output "$run_dir/volume"
    printf '%q ' "${extractor[@]}" --output > "$run_dir/run.txt"
    printf '%q\n' "$results_root/$name/volume" >> "$run_dir/run.txt"
)

"${python_cmd[@]}" "$experiment_dir/evaluate_curved_prism.py" "$run_dir"
"${python_cmd[@]}" "$experiment_dir/plot_curved_prism.py" "$run_dir"
"${python_cmd[@]}" -m zipfile -c "$work/$name.zip" "$run_dir"

mkdir -p "$results_root/$name"
mv -- "$run_dir/"* "$results_root/$name/"
mv -- "$work/$name.zip" "$results_root/$name.zip"
printf 'Results and figures: %s\nArchive: %s\n' "$results_root/$name" "$results_root/$name.zip"
