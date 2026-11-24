#!/usr/bin/env bash
set -euo pipefail

ROOT=/work/exp063
REPO="$ROOT/repo"
INPUT="$ROOT/input"
OUT="$REPO/experiments/063-pds-mse-pareto/out"
S3_INPUT_PREFIX="${S3_INPUT_PREFIX:?set S3_INPUT_PREFIX}"
S3_OUTPUT_PREFIX="${S3_OUTPUT_PREFIX:?set S3_OUTPUT_PREFIX}"
INSTANCE_ID="${INSTANCE_ID:-unknown}"
INSTANCE_TYPE="${INSTANCE_TYPE:-r7i.2xlarge}"
REGION="${AWS_REGION:-us-east-1}"

mkdir -p "$ROOT" "$REPO" "$INPUT" "$OUT"
aws s3 cp "$S3_INPUT_PREFIX/source.tar.gz" "$ROOT/source.tar.gz"
tar -xzf "$ROOT/source.tar.gz" -C "$REPO"

for donor in k562_gwps orion_hct116 orion_hek293t; do
  mkdir -p "$INPUT/$donor"
  aws s3 cp \
    "s3://<your-s3-bucket>/input/clean-hct/$donor/" \
    "$INPUT/$donor/" \
    --recursive \
    --exclude '*' \
    --include 'delta.npy' \
    --include 'b_ctrl.npy' \
    --include 'genes.txt' \
    --include 'perts.txt' \
    --include 'meta.json'
done

aws s3 cp \
  s3://<your-s3-bucket>/input/exp058/gene_names.csv \
  "$INPUT/gene_names.csv"
aws s3 cp \
  s3://<your-s3-bucket>/input/exp058/pert_counts.csv \
  "$INPUT/pert_counts.csv"

cd "$REPO"
set +e
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 55 \
  --min-free-percent 5 \
  --poll-seconds 5 \
  -- \
  python3 experiments/063-pds-mse-pareto/run.py \
    --k562 "$INPUT/k562_gwps" \
    --hct116 "$INPUT/orion_hct116" \
    --hek293t "$INPUT/orion_hek293t" \
    --genes-csv "$INPUT/gene_names.csv" \
    --panel-csv "$INPUT/pert_counts.csv" \
    --output-dir "$OUT" \
    --max-rss-gib 55 \
    --instance-id "$INSTANCE_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --region "$REGION" \
    --s3-prefix "$S3_OUTPUT_PREFIX" \
  >"$OUT/run.log" 2>&1

status=$?
set -e
printf '%s\n' "$status" >"$OUT/exit-code.txt"
(
  cd "$REPO"
  sha256sum \
    experiments/063-pds-mse-pareto/README.md \
    experiments/063-pds-mse-pareto/run.py \
    experiments/063-pds-mse-pareto/test_run.py \
    experiments/063-pds-mse-pareto/aws/run.sh \
    experiments/063-pds-mse-pareto/out/results.json \
    experiments/063-pds-mse-pareto/out/summary.csv \
    experiments/063-pds-mse-pareto/out/per_target_metrics.csv \
    experiments/063-pds-mse-pareto/out/gap_translation.csv \
    >"$OUT/hashes.sha256"
)
aws s3 cp "$OUT/" "$S3_OUTPUT_PREFIX/" --recursive
exit "$status"
