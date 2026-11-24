#!/usr/bin/env bash
set -uo pipefail
cd ~/Code/virtual-cell-challenge-2026
B=s3://<your-s3-bucket>
for ARM in od1 ks010 ts070 kt16000 amp070 amp085 amp120 amp145; do
  echo "=== launching $ARM at $(date -u +%FT%TZ) ==="
  VCC_CLOUD_SCORE_CLEARANCE=YES scripts/cloud_score.sh \
    --real $B/input/harness/h1-test.h5ad \
    --prologue experiments/317-operating-point/aws/prologue-$ARM.sh \
    --source-include README.md --source-include pyproject.toml --source-include uv.lock \
    --source-include scripts --source-include src \
    --source-include experiments/061-perturbmatch-hct-count/build_candidate.py \
    --source-include experiments/062-cd4-post-whiten/run.py \
    --source-include experiments/317-operating-point/build317.py \
    --pred-from-prologue "$ARM" \
    --out experiments/317-operating-point/cloud-output --run-id "s317-$ARM" \
    --no-spot --max-minutes 215 2>&1 | tail -5
done
echo "=== all launches attempted at $(date -u +%FT%TZ) ==="
