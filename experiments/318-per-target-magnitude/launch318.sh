#!/usr/bin/env bash
# 318 — one disposable worker per arm. `base` is NOT scored (317's KA-1 row is its row, and
# gate 5 asserts byte identity), so six workers cover the whole design inside the $12 cap.
set -uo pipefail
cd ~/Code/virtual-cell-challenge-2026
B=s3://<your-s3-bucket>
for ARM in "$@"; do
  echo "=== launching $ARM at $(date -u +%FT%TZ) ==="
  VCC_CLOUD_SCORE_CLEARANCE=YES scripts/cloud_score.sh \
    --real $B/input/harness/h1-test.h5ad \
    --prologue experiments/318-per-target-magnitude/aws/prologue-$ARM.sh \
    --source-include README.md --source-include pyproject.toml --source-include uv.lock \
    --source-include scripts --source-include src \
    --source-include experiments/061-perturbmatch-hct-count/build_candidate.py \
    --source-include experiments/062-cd4-post-whiten/run.py \
    --source-include experiments/318-per-target-magnitude/build318.py \
    --source-include experiments/318-per-target-magnitude/work/panel_profile.json \
    --pred-from-prologue "$ARM" \
    --out experiments/318-per-target-magnitude/cloud-output --run-id "s318-$ARM" \
    --no-spot --max-minutes 215 2>&1 | tail -4
done
echo "=== all launches attempted at $(date -u +%FT%TZ) ==="
