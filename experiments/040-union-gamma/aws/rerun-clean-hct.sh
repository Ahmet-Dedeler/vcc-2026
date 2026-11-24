#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
WORK="/work"
REPO="$WORK/repo"
LOG="/var/log/vcc-clean-hct-rerun.log"

exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  aws s3 cp "$LOG" \
    "s3://${BUCKET}/output-clean-hct/cloud.log" --only-show-errors || true
  printf '%s\n' "$status" > "$WORK/clean-hct-exit-code.txt"
  aws s3 cp "$WORK/clean-hct-exit-code.txt" \
    "s3://${BUCKET}/output-clean-hct/exit-code.txt" --only-show-errors || true
  shutdown -h now
  exit "$status"
}
trap finish EXIT

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="$WORK/uv-cache"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

test -s "$REPO/data/harness/orion-hct116.h5ad"
test -s "$REPO/data/donors/k562_gwps/delta.npy"
test -s "$REPO/data/donors/orion_hek293t/delta.npy"

# The minimal source archive omitted the PEP 621 readme referenced by
# pyproject.toml. It is metadata-only and need not contain prose for this run.
touch "$REPO/README.md"

cd "$REPO"
/usr/bin/time -v uv run python experiments/031-somatic/run.py \
  --real data/harness/orion-hct116.h5ad \
  --donors \
    k562_gwps=data/donors/k562_gwps \
    orion_hek293t=data/donors/orion_hek293t \
  --configs experiments/040-union-gamma/grid_hct_clean.json \
  --outdir experiments/040-union-gamma/out_hct_clean \
  --whiten-axis union \
  --exclude-panel challenge \
  --jobs 1 \
  2>&1 | tee "$WORK/score-clean-hct.log"

aws s3 cp experiments/040-union-gamma/out_hct_clean/results.json \
  "s3://${BUCKET}/output-clean-hct/results.json" --only-show-errors
aws s3 cp "$WORK/score-clean-hct.log" \
  "s3://${BUCKET}/output-clean-hct/score.log" --only-show-errors
for run in experiments/040-union-gamma/out_hct_clean/runs/*
do
  tag="$(basename "$run")"
  aws s3 cp "$run/run/agg_results.csv" \
    "s3://${BUCKET}/output-clean-hct/${tag}/agg_results.csv" \
    --only-show-errors
done
