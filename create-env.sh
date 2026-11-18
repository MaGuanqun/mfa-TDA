#!/bin/bash

# Provision in a subshell so errors and activation do not change the caller.
# Safe to rerun: existing environments and installed packages are retained.
(
    set -euo pipefail
    project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd) || exit 1
    environment=mfa-TDA-env

    if [[ $(type -t spack || true) != function ]]; then
        spack_setup=${SPACK_SETUP_ENV:-"${SPACK_ROOT:-$HOME/spack}/share/spack/setup-env.sh"}
        if [[ ! -f "$spack_setup" ]]; then
            printf 'Initialize Spack shell support or set SPACK_SETUP_ENV to share/spack/setup-env.sh\n' >&2
            exit 1
        fi
        source "$spack_setup" || exit 1
    fi

    if ! spack location -e "$environment" >/dev/null 2>&1; then
        spack env create "$environment" "$project_root/env.yaml" || exit 1
    fi
    spack env activate "$environment" || exit 1

    # Keep project recipes in this environment's configuration, with no
    # machine-specific paths committed to env.yaml.
    recipe_repo="$project_root/spack-repo"
    registered_repos=$(spack repo list --scope "env:$environment") || exit 1
    if [[ "$registered_repos" != *"$recipe_repo"* ]]; then
        spack repo add --scope "env:$environment" "$recipe_repo" || exit 1
    fi

    # env.yaml is the single source of root specs, including when updating an
    # older environment. spack add leaves already-present specs unchanged.
    # spack python -c uses an interactive parser: keep this a single statement.
    specs=$(MFA_TDA_MANIFEST="$project_root/env.yaml" spack python -c 'import os; import spack.util.spack_yaml as yaml; print("\n".join(yaml.load(open(os.environ["MFA_TDA_MANIFEST"]))["spack"]["specs"]))') || exit 1
    if [[ -z "$specs" ]]; then
        printf 'Could not read dependency specs from %s/env.yaml\n' "$project_root" >&2
        exit 1
    fi
    while IFS= read -r spec; do
        spack add "$spec" || exit 1
    done <<< "$specs"

    # Retain the existing MFA-before-LibTorch installation order.
    spack install -j 3 'mfa~examples~tests thread=tbb' || exit 1
    MFA_PATH=$(spack location -i mfa) || exit 1
    export MFA_PATH
    spack install -j 3 || exit 1
    printf 'Environment %s is ready. Run: source ./load-env.sh\n' "$environment"
)

