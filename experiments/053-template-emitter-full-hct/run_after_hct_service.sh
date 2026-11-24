#!/bin/bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
REPO="${REPO:-/work/repo}"
HCT_SERVICE="${HCT_SERVICE:-vcc-hct-g065.service}"
GAMMA_TAG="${GAMMA_TAG:-g0.65_q1.0_a1.08_spread0.25}"
S3_INPUT_URI="${S3_INPUT_URI:-s3://${BUCKET}/input/template053}"
S3_OUTPUT_URI="${S3_OUTPUT_URI:-s3://${BUCKET}/output/template-hct}"
LOG="${LOG:-/work/template-${GAMMA_TAG}.log}"
EXIT_FILE="${EXIT_FILE:-/work/template-${GAMMA_TAG}.exit-code.txt}"

export PATH="/root/.local/bin:${REPO}/.venv/bin:${PATH}"

exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXIT_FILE"
  aws s3 cp "$LOG" "${S3_OUTPUT_URI%/}/${GAMMA_TAG}/runner.log" \
    --only-show-errors || true
  aws s3 cp "$EXIT_FILE" "${S3_OUTPUT_URI%/}/${GAMMA_TAG}/exit-code.txt" \
    --only-show-errors || true
  exit "$status"
}
trap finish EXIT

mkdir -p "$REPO/experiments/053-template-emitter-full-hct" "$REPO/scripts"
aws s3 cp "$S3_INPUT_URI/" "$REPO/experiments/053-template-emitter-full-hct/" \
  --recursive --exclude "run_with_memory_guard.py" --only-show-errors
aws s3 cp "$S3_INPUT_URI/run_with_memory_guard.py" \
  "$REPO/scripts/run_with_memory_guard.py" --only-show-errors

while systemctl is-active --quiet "$HCT_SERVICE"
do
  echo "waiting for ${HCT_SERVICE} at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  sleep 30
done

test -s \
  "$REPO/experiments/040-union-gamma/out_hct_clean/runs/$GAMMA_TAG/run/agg_results.csv"

REPO="$REPO" \
GAMMA_TAG="$GAMMA_TAG" \
S3_OUTPUT_URI="$S3_OUTPUT_URI" \
bash "$REPO/experiments/053-template-emitter-full-hct/run_cloud_after_gamma.sh"
