#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
for component in u v; do
    dataset="vortex_street_3d_${component}"
    mkdir -p "../Exp/$dataset"
    python -u main.py --train train --dataset "$dataset" \
        --application super-spatial-temporal --num_res 5 --active sine \
        --init 64 --num_epochs 300 --omega_0 30.0 \
        --factor 1 --batch_size 16000 --lap_weight 0.0 \
        | tee "../Exp/$dataset/training.log"
done
