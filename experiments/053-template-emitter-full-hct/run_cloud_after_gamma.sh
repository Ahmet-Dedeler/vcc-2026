#!/usr/bin/env bash
set -Eeuo pipefail

# Run this on the EC2 worker after experiment 040's clean HCT116 gamma sweep.
# Transformations are capped at 4 GiB process-tree RSS. cell-eval2 is always
# invoked serially and is intentionally not subject to that transformer cap.

REPO="${REPO:-/work/repo}"
REAL="${REAL:-${REPO}/data/harness/orion-hct116.h5ad}"
GAMMA_ROOT="${GAMMA_ROOT:-${REPO}/experiments/040-union-gamma/out_hct_clean}"
OUT_ROOT="${OUT_ROOT:-${REPO}/experiments/053-template-emitter-full-hct/cloud-output}"
CELL_EVAL_VERSION="${CELL_EVAL_VERSION:-0.16.0}"
SCORE_QUASI="${SCORE_QUASI:-0}"
FORCE="${FORCE:-0}"
S3_OUTPUT_URI="${S3_OUTPUT_URI:-}"
GAMMA_TAG="${GAMMA_TAG:-}"

cd "$REPO"
test -s "$REAL"
test -d "$GAMMA_ROOT/runs"
test -f "experiments/053-template-emitter-full-hct/transform.py"
test -f "scripts/run_with_memory_guard.py"
mkdir -p "$OUT_ROOT"

score_one() {
  local prediction="$1"
  local output_dir="$2"
  if [[ -s "$output_dir/agg_results.csv" && "$FORCE" != "1" ]]; then
    echo "reuse score: $output_dir/agg_results.csv"
    return
  fi
  if [[ -e "$output_dir" ]]; then
    if [[ "$FORCE" != "1" ]]; then
      echo "incomplete score directory exists: $output_dir" >&2
      exit 2
    fi
    rm -rf "$output_dir"
  fi
  /usr/bin/time -v uvx --python 3.12 --from "cell-eval2==${CELL_EVAL_VERSION}" \
    cell-eval2 run \
    -ap "$prediction" \
    -ar "$REAL" \
    -o "$output_dir" \
    --preset vcc2026 \
    --pert-col target_gene \
    --control non-targeting
}

found=0
if [[ -n "$GAMMA_TAG" ]]; then
  gamma_runs=("$GAMMA_ROOT/runs/$GAMMA_TAG")
else
  gamma_runs=("$GAMMA_ROOT"/runs/*)
fi

for gamma_run in "${gamma_runs[@]}"
do
  [[ -d "$gamma_run" ]] || continue
  quasi_prediction="$gamma_run/pred.h5ad"
  [[ -s "$quasi_prediction" ]] || continue
  found=$((found + 1))
  tag="$(basename "$gamma_run")"
  target_dir="$OUT_ROOT/$tag"
  template_prediction="$target_dir/pred-template.h5ad"
  transform_report="$target_dir/transform-report.json"
  mkdir -p "$target_dir"

  if [[ ! -s "$template_prediction" || ! -s "$transform_report" || "$FORCE" == "1" ]]; then
    overwrite_args=()
    if [[ "$FORCE" == "1" ]]; then
      overwrite_args+=(--overwrite)
    fi
    uv run python scripts/run_with_memory_guard.py \
      --max-rss-gib 4.0 \
      --min-free-percent 12 \
      --poll-seconds 1 \
      -- \
      uv run python experiments/053-template-emitter-full-hct/transform.py \
        --prediction "$quasi_prediction" \
        --real "$REAL" \
        --output "$template_prediction" \
        --report "$transform_report" \
        --max-rss-gib 3.8 \
        "${overwrite_args[@]}"
  else
    echo "reuse transform: $template_prediction"
  fi

  if [[ "$SCORE_QUASI" == "1" ]]; then
    score_one "$quasi_prediction" "$target_dir/quasi-run"
    cp "$target_dir/quasi-run/agg_results.csv" \
      "$target_dir/quasi-agg_results.csv"
  elif [[ -s "$gamma_run/run/agg_results.csv" ]]; then
    cp "$gamma_run/run/agg_results.csv" \
      "$target_dir/quasi-agg_results.csv"
  else
    echo "missing quasi score for $tag; set SCORE_QUASI=1 to generate it" >&2
    exit 2
  fi

  score_one "$template_prediction" "$target_dir/template-run"

  if [[ -n "$S3_OUTPUT_URI" ]]; then
    aws s3 cp "$transform_report" \
      "${S3_OUTPUT_URI%/}/$tag/transform-report.json" --only-show-errors
    aws s3 cp "$target_dir/quasi-agg_results.csv" \
      "${S3_OUTPUT_URI%/}/$tag/quasi-agg_results.csv" --only-show-errors
    aws s3 cp "$target_dir/template-run/agg_results.csv" \
      "${S3_OUTPUT_URI%/}/$tag/template-agg_results.csv" --only-show-errors
  fi
done

if [[ "$found" == "0" ]]; then
  echo "no gamma prediction files found under $GAMMA_ROOT/runs" >&2
  exit 2
fi

uv run python experiments/053-template-emitter-full-hct/compare_scores.py \
  --root "$OUT_ROOT" \
  --output "$OUT_ROOT/comparison.json"

if [[ -n "$S3_OUTPUT_URI" ]]; then
  aws s3 cp "$OUT_ROOT/comparison.json" \
    "${S3_OUTPUT_URI%/}/comparison.json" --only-show-errors
fi

echo "completed $found gamma configurations serially"
