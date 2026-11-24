#!/bin/bash
set -Eeuo pipefail

TAG="${1:?run tag required}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
REPO="${REPO:-/work/repo}"
REAL="${REAL:-$REPO/data/harness/orion-hct116.h5ad}"
RUN_DIR="$REPO/experiments/040-union-gamma/out_hct_clean/runs/$TAG"
LOG="/work/score-${TAG}.log"
EXIT_FILE="/work/score-${TAG}.exit-code.txt"
S3_PREFIX="s3://${BUCKET}/output-clean-hct/${TAG}"

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="/work/uv-cache"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXIT_FILE"
  aws s3 cp "$LOG" "$S3_PREFIX/score.log" --only-show-errors || true
  aws s3 cp "$EXIT_FILE" "$S3_PREFIX/exit-code.txt" --only-show-errors || true
  if [[ -s "$RUN_DIR/run/agg_results.csv" ]]; then
    aws s3 cp "$RUN_DIR/run/agg_results.csv" \
      "$S3_PREFIX/agg_results.csv" --only-show-errors || true
  fi
  exit "$status"
}
trap finish EXIT

test -s "$REAL"
test -s "$RUN_DIR/pred.h5ad"
rm -rf "$RUN_DIR/run"

cd "$REPO"
/usr/bin/time -v uv run cell-eval2 run \
  -ap "$RUN_DIR/pred.h5ad" \
  -ar "$REAL" \
  -o "$RUN_DIR/run" \
  --preset vcc2026 \
  --pert-col target_gene \
  --control non-targeting \
  2>&1 | tee "$LOG"
