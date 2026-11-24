#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
WORK="/work"
REPO="$WORK/repo"
LOG="/var/log/vcc-clean-hct-score.log"

mkdir -p "$WORK"
exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  aws s3 cp "$LOG" "s3://${BUCKET}/output-clean-hct/cloud.log" || true
  printf '%s\n' "$status" > "$WORK/exit-code.txt"
  aws s3 cp "$WORK/exit-code.txt" \
    "s3://${BUCKET}/output-clean-hct/exit-code.txt" || true
  shutdown -h now
  exit "$status"
}
trap finish EXIT

wait_for_object() {
  key="$1"
  until aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null 2>&1
  do
    echo "waiting for s3://${BUCKET}/${key}"
    sleep 30
  done
}

dnf install -y gzip tar time
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="$WORK/uv-cache"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

for key in \
  input/clean-hct/source.tar.gz \
  input/orion-hct116.h5ad \
  input/clean-hct/k562_gwps/delta.npy \
  input/clean-hct/orion_hek293t/delta.npy
do
  wait_for_object "$key"
done

mkdir -p "$REPO"
aws s3 cp "s3://${BUCKET}/input/clean-hct/source.tar.gz" "$WORK/source.tar.gz"
tar -xzf "$WORK/source.tar.gz" -C "$REPO"
mkdir -p "$REPO/data/harness" "$REPO/data/donors/k562_gwps" \
  "$REPO/data/donors/orion_hek293t"
aws s3 cp "s3://${BUCKET}/input/orion-hct116.h5ad" \
  "$REPO/data/harness/orion-hct116.h5ad"
aws s3 cp "s3://${BUCKET}/input/clean-hct/k562_gwps/" \
  "$REPO/data/donors/k562_gwps/" --recursive
aws s3 cp "s3://${BUCKET}/input/clean-hct/orion_hek293t/" \
  "$REPO/data/donors/orion_hek293t/" --recursive

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
  2>&1 | tee "$WORK/score.log"

aws s3 cp experiments/040-union-gamma/out_hct_clean/results.json \
  "s3://${BUCKET}/output-clean-hct/results.json"
aws s3 cp "$WORK/score.log" \
  "s3://${BUCKET}/output-clean-hct/score.log"
for run in experiments/040-union-gamma/out_hct_clean/runs/*
do
  tag="$(basename "$run")"
  aws s3 cp "$run/run/agg_results.csv" \
    "s3://${BUCKET}/output-clean-hct/${tag}/agg_results.csv"
done
