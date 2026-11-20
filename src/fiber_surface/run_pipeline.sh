#!/usr/bin/env bash
# Shared by the repository-root fiber-surface scripts.

fiber_setup() {
    cd "$fiber_repo"
    if [[ $(type -t spack || true) != function ]]; then
        local setup="${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"
        if [[ ! -f "$setup" ]]; then
            printf 'Initialize Spack shell support or set SPACK_ROOT before running this script.\n' >&2
            return 1
        fi
        source "$setup" || return
    fi
    # load-env.sh intentionally tolerates deactivation when no env is active.
    # Source it in a conditional so Bash's errexit does not stop at that probe.
    source ./load-env.sh || return
    if [[ ! -d "$MFA_PATH/include" || ! -d "$TBB_PATH" || ! -d "$EIGEN_PATH" || ! -d "$TORCH_PATH" ]]; then
        printf 'Spack activation did not provide the required dependency paths.\n' >&2
        return 1
    fi
    fiber_binary_dir="$fiber_repo/build/src/fiber_surface"
    fiber_results_root="${FIBER_RESULTS_ROOT:-$fiber_repo/build/fiber-surface-results}"
    if [[ "${FIBER_SKIP_BUILD:-0}" != 1 ]]; then
        cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-flto=auto" -Dmfa_thread=tbb
        cmake --build build --target fiber_find_seeds fiber_trace_endpoint_fibers \
            fiber_march_surface fiber_stitch_surface fiber_evaluate_surface -j 3
    fi
}

fiber_stage_command() {
    local run="$1" stage="$2" executable="$3"
    shift 3
    mkdir -p "$run"
    local command=("$fiber_binary_dir/$executable" --run-dir "$run" "$@")
    if [[ ! -x "${command[0]}" ]]; then
        printf 'Missing executable: %s\n' "${command[0]}" >&2
        return 1
    fi
    # Shell-escaped commands can be copied directly into Bash, including paths
    # with spaces. Each stage log is replaced when that stage is rerun.
    {
        printf '#!/usr/bin/env bash\nset -eo pipefail\nfiber_repo=%q\n' "$fiber_repo"
        printf 'source "$fiber_repo/src/fiber_surface/run_pipeline.sh"\nFIBER_SKIP_BUILD=1 fiber_setup\n'
        printf '%q ' "${command[@]}"
        printf '\n'
    } > "$run/${stage}_command.sh"
    local status=0
    "${command[@]}" 2>&1 | tee "$run/${stage}.log" || status=$?
    return "$status"
}

fiber_run_pipeline() {
    local run="$1" selection="$2"
    shift 2
    local names=(seeds fibers march stitch evaluate)
    local executables=(fiber_find_seeds fiber_trace_endpoint_fibers fiber_march_surface fiber_stitch_surface fiber_evaluate_surface)
    local index status overall=0
    for index in 0 1 2 3 4; do
        [[ "$selection" == all || "$selection" == "${names[index]}" ]] || continue
        local args=()
        if [[ "$index" == 0 ]]; then args=("$@"); fi
        status=0
        fiber_stage_command "$run" "${names[index]}" "${executables[index]}" "${args[@]}" || status=$?
        if [[ "$status" != 0 && "$status" != 2 ]]; then return "$status"; fi
        # Keep going after an incomplete extraction so evaluation can write the
        # final partial mesh and report. Status 2 remains visible to the caller.
        if [[ "$status" == 2 ]]; then overall=2; fi
    done
    printf 'Results: %s\n' "$run"
    return "$overall"
}

fiber_check_stage() {
    case "$1" in all|seeds|fibers|march|stitch|evaluate) return 0;;
        *) printf 'Unknown stage: %s (all, seeds, fibers, march, stitch, evaluate)\n' "$1" >&2; return 1;;
    esac
}
