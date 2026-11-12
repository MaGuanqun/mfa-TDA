#!/bin/bash
# Tracking ablation; initial seed finding still uses Newton.
# Usage: bash critical_point_tracking_no_newton.sh [dataset] [step_ratio] [output_dir]
# Datasets: quartic_potential_2 (default, ratio 4), boussinesq_3d (MFA, ratio 16).
# The original invocation with only a step ratio and output directory still works.
# Build the selected tracker through the root CMake project first.
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [[ $(type -t spack || true) != function ]]; then
    spack_setup=${SPACK_SETUP_ENV:-"${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"}
    source "${spack_setup}"
fi
# load-env.sh first deactivates Spack, which can return 1 if none is active.
set +e
source ./load-env.sh
env_status=$?
set -e
if [[ $env_status -ne 0 || ${SPACK_ENV##*/} != mfa-TDA-env || ! -d $MFA_PATH || ! -d $TORCH_PATH ]]; then
    printf 'Could not activate mfa-TDA-env\n' >&2
    exit 1
fi

dataset=quartic_potential_2
default_step_ratio=4
case "${1:-}" in
    quartic_potential_2) shift ;;
    boussinesq_3d) dataset=boussinesq_3d; default_step_ratio=16; shift ;;
esac
step_ratio="${1:-${default_step_ratio}}"
if [[ ! $step_ratio =~ ^[1-9][0-9]*$ ]]; then
    printf 'Step ratio must be a positive integer; supported datasets: quartic_potential_2, boussinesq_3d\n' >&2
    exit 1
fi
output_dir="${2:-./build/src/${dataset}/no_newton_correction}"
input_dir="./build/src/${dataset}"
pvpython="${PVPYTHON:-/usr/local/bin/pvpython}"
prefix="${output_dir}/${dataset}_no_newton_${step_ratio}"
if [[ $dataset == boussinesq_3d ]]; then
    tracking="./build/src/critical_point_tracking/critical_point_tracking"
    degenerate_input="${input_dir}/${dataset}_degenerate${step_ratio}.dat"
    model="${input_dir}/${dataset}.mfa"
    controls="${input_dir}/${dataset}_cpt.dat"
    test -s "${model}"
    test -s "${controls}"
    input_args=(-f "${model}" -a "${controls}" -j "${output_dir}/${dataset}_boundary")
else
    tracking="./build/src/critical_point_tracking/critical_point_tracking_explicit"
    degenerate_input="${input_dir}/${dataset}_degenerate_${step_ratio}.dat"
    input_args=(-f "${dataset}")
fi

test -x "${tracking}"
test -s "${degenerate_input}"
mkdir -p "${output_dir}"
# The tracker writes seed sidecars beside this file; keep them in the experiment.
degenerate_copy="${output_dir}/${dataset}_degenerate_${step_ratio}.dat"
cp "${degenerate_input}" "${degenerate_copy}"

tracking_args=("${input_args[@]}" -b "${prefix}"
    -z "${step_ratio}" -g "${step_ratio}" -x 1e-10
    -s "${degenerate_copy}" -p 4.0 -n 0
    -e "${prefix}_edge_type.csv")
printf '%q ' "${tracking}" "${tracking_args[@]}" > "${prefix}_command.txt"
printf '\n' >> "${prefix}_command.txt"
"${tracking}" "${tracking_args[@]}" 2>&1 | tee "${prefix}.log"

"${pvpython}" --no-mpi ./src/python/merge_obj_edge_type.py \
    -i "${prefix}.obj" -j "${prefix}_edge_type.csv" -o "${prefix}.vtp"
printf 'Exported: %s\n' "${prefix}.vtp"

if [[ $dataset == boussinesq_3d ]]; then
    reference="${input_dir}/${dataset}_16_remove_obstacle_regionid_colorid.vtp"
    test -s "${reference}"
    # The matching interpreter needs VTK, NumPy and SciPy. It can differ from
    # the ParaView installation used for the OBJ-to-VTP conversion above.
    postprocess_python="${POSTPROCESS_PYTHON:-pvpython}"
    post=("${postprocess_python}")
    if [[ ${postprocess_python##*/} == pvpython* ]]; then
        post+=(--no-mpi)
    fi
    "${post[@]}" ./src/python/remove_obstacle.py --data "${dataset}" \
        -i "${prefix}.vtp" -o "${prefix}_remove_obstacle.vtp" \
        2>&1 | tee "${prefix}_remove_obstacle.log"
    # Give every branch the closest reference branch's color. Keep distances
    # so geometric assignments can be inspected even when a track has drifted.
    "${post[@]}" ./src/python/seperate_branch_assign_color2.py \
        -i "${prefix}_remove_obstacle.vtp" -r "${reference}" \
        -o "${prefix}_remove_obstacle_regionid_colorid.vtp" \
        --distance-mode mean --match-distance-threshold inf --seed 0 \
        --ref-region-association cell \
        --matches-csv "${prefix}_color_matches.csv" \
        2>&1 | tee "${prefix}_color_matches.log"
    printf 'Colored result: %s\n' "${prefix}_remove_obstacle_regionid_colorid.vtp"
fi
