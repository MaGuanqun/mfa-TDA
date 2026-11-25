# Shared by the Isabel outer-ring and vortex street experiment launchers.
# Run from the repository root with the build environment loaded.
run_triangle_stages() (
    set -eo pipefail
    run_dir=$1
    stage=$2
    step=$3
    shift 3
    bin=build/src/fiber_surface
    stages="$run_dir/stages"
    seed=("$bin/fiber_find_seeds" --run-dir "$stages" --step "$step"
          --fiber-seed-radius "$step" --spatial-tolerance "$step"
          --max-trace-steps 100000 "$@")
    trace=("$bin/fiber_trace_endpoint_fibers" --run-dir "$stages")
    march=("$bin/fiber_march_surface" --run-dir "$stages")
    stitch=("$bin/fiber_stitch_surface" --run-dir "$stages")
    assemble=("$bin/fiber_volume_from_patches" --run-dir "$stages" --output "$run_dir/volume")
    if [[ $stage == surface ]]; then
        assemble+=(--stage surface)
        audit=("$bin/fiber_grid_recovery" audit "$stages" "$run_dir/volume_boundary.ply"
               "$run_dir/audit" --all-face-samples)
    fi

    record() { printf '%q ' "$@"; printf '\n'; }
    {
        record "${seed[@]}"
        record "${trace[@]}"
        for side in 0 1 2; do
            record "${march[@]}" --segment "$side"
            record "${stitch[@]}" --segment "$side"
        done
        record "${assemble[@]}"
        if [[ $stage == surface ]]; then record "${audit[@]}"; fi
    } > "$run_dir/run.txt"

    run_stage() {
        # Exit 2 means the stage saved a partial result for the next stage.
        local status=0
        "$@" || status=$?
        (( status == 0 || status == 2 )) || return "$status"
    }
    {
        "${seed[@]}"
        run_stage "${trace[@]}"
        printf 'March and stitch AB, BC and CA with separate checkpoints.\n'
        labels=(AB BC CA)
        pids=()
        for side in 0 1 2; do
            (
                "${march[@]}" --segment "$side" &&
                    run_stage "${stitch[@]}" --segment "$side"
            ) 2>&1 | tee "$run_dir/${labels[side]}.log" | sed -u "s/^/[${labels[side]}] /" &
            pids+=("$!")
        done
        failed=0
        for pid in "${pids[@]}"; do
            wait "$pid" || failed=1
        done
        (( failed == 0 )) || exit 1
        run_stage "${assemble[@]}"
        if [[ $stage == surface ]]; then "${audit[@]}"; fi
    } 2>&1 | tee "$run_dir/extract.log"
)
