#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_PREFIX="${RUN_PREFIX:?RUN_PREFIX is required}"
INPUT="s3://${BUCKET}/input/${RUN_PREFIX}"
OUTPUT="${OUTPUT_PREFIX:-s3://${BUCKET}/output/${RUN_PREFIX}}"
BASELINE_PREFIX="s3://${BUCKET}/output-clean-hct/g0.65_q1.0_a1.08_spread0.25"
WORK="${WORK:-/work/exp061}"
REPO="${WORK}/repo"
OUT="${REPO}/experiments/061-perturbmatch-hct-count/cloud-output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

mkdir -p "$WORK" "$OUT"
rm -f "$EXIT_FILE"
exec 9>"${WORK}.lock"
if ! flock -n 9; then
  echo "another experiment-061 runner already holds ${WORK}.lock" >&2
  exit 70
fi
exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  aws s3 cp "$LOG" "${OUTPUT}/cloud.log" --only-show-errors || true
  if [[ -s "$EXIT_FILE" ]]; then
    aws s3 cp "$EXIT_FILE" "${OUTPUT}/exit-code.txt" --only-show-errors || true
  fi
  if [[ -d "$OUT" ]]; then
    aws s3 sync "$OUT/" "${OUTPUT}/" --only-show-errors || true
  fi
}

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXIT_FILE"
  upload_outputs
  shutdown -h now || true
  exit "$status"
}
trap finish EXIT

dnf install -y gzip tar time
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

aws s3 cp "${INPUT}/source.tar.gz" "${WORK}/source.tar.gz" --only-show-errors
mkdir -p "$REPO"
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"
mkdir -p \
  "$REPO/data/harness" \
  "$REPO/data/challenge" \
  "$REPO/data/donors/k562_gwps" \
  "$REPO/data/donors/orion_hek293t" \
  "$REPO/data/external/perturbmatch" \
  "$REPO/experiments/040-union-gamma/out_hct_clean"

aws s3 cp "${INPUT}/orion-hct116.h5ad" \
  "$REPO/data/harness/orion-hct116.h5ad" --only-show-errors
aws s3 cp "${INPUT}/challenge/" "$REPO/data/challenge/" \
  --recursive --only-show-errors
aws s3 cp "${INPUT}/k562_gwps/" "$REPO/data/donors/k562_gwps/" \
  --recursive --only-show-errors
aws s3 cp "${INPUT}/orion_hek293t/" "$REPO/data/donors/orion_hek293t/" \
  --recursive --only-show-errors
aws s3 cp "${INPUT}/perturbmatch/" "$REPO/data/external/perturbmatch/" \
  --recursive --only-show-errors
aws s3 cp "${INPUT}/whitener_union.npz" \
  "$REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.npz" \
  --only-show-errors
aws s3 cp "${INPUT}/whitener_union.json" \
  "$REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.json" \
  --only-show-errors

cd "$REPO"
uv sync --frozen
uv run python - <<'PY' > "$OUT/runtime.json"
import importlib.metadata
import json
import platform
print(json.dumps({
    "python": platform.python_version(),
    "cell_eval2": importlib.metadata.version("cell-eval2"),
}, indent=2))
PY

rm -f "$OUT/candidate.h5ad" "$OUT/candidate.h5ad.sha256" "$OUT/build-report.json"
/usr/bin/time -v uv run python \
  experiments/061-perturbmatch-hct-count/build_candidate.py \
  --real data/harness/orion-hct116.h5ad \
  --k562 data/donors/k562_gwps \
  --hek293t data/donors/orion_hek293t \
  --perturbmatch data/external/perturbmatch \
  --panel-csv data/challenge/pert_counts.csv \
  --genes-csv data/challenge/gene_names.csv \
  --whitener experiments/040-union-gamma/out_hct_clean/whitener_union.npz \
  --whitener-metadata experiments/040-union-gamma/out_hct_clean/whitener_union.json \
  --output "$OUT/candidate.h5ad" \
  --report "$OUT/build-report.json" \
  --max-rss-gib 28 \
  2>&1 | tee "$OUT/build.log"

sha256sum "$OUT/candidate.h5ad" > "$OUT/candidate.h5ad.sha256"
aws s3 cp "$OUT/candidate.h5ad" "${OUTPUT}/candidate.h5ad" --only-show-errors
aws s3 cp "$OUT/candidate.h5ad.sha256" \
  "${OUTPUT}/candidate.h5ad.sha256" --only-show-errors
upload_outputs

echo "candidate built; waiting for the existing clean gamma-0.65 baseline score"
until aws s3api head-object \
  --bucket "$BUCKET" \
  --key "output-clean-hct/g0.65_q1.0_a1.08_spread0.25/agg_results.csv" \
  >/dev/null 2>&1
do
  date -u
  sleep 30
done
aws s3 cp "${BASELINE_PREFIX}/agg_results.csv" \
  "$OUT/baseline-agg_results.csv" --only-show-errors
if aws s3api head-object \
  --bucket "$BUCKET" \
  --key "output-clean-hct/g0.65_q1.0_a1.08_spread0.25/pred.h5ad.sha256" \
  >/dev/null 2>&1
then
  aws s3 cp "${BASELINE_PREFIX}/pred.h5ad.sha256" \
    "$OUT/baseline-pred.h5ad.sha256" --only-show-errors
fi

rm -rf "$OUT/candidate-run"
/usr/bin/time -v uvx --python 3.12 --from "cell-eval2==0.16.0" \
  cell-eval2 run \
  -ap "$OUT/candidate.h5ad" \
  -ar data/harness/orion-hct116.h5ad \
  -o "$OUT/candidate-run" \
  --preset vcc2026 \
  --pert-col target_gene \
  --control non-targeting \
  2>&1 | tee "$OUT/score.log"

uv run python experiments/061-perturbmatch-hct-count/compare_scores.py \
  --baseline "$OUT/baseline-agg_results.csv" \
  --candidate "$OUT/candidate-run/agg_results.csv" \
  --output "$OUT/results.json" \
  2>&1 | tee "$OUT/compare.log"

sha256sum \
  "$OUT/build-report.json" \
  "$OUT/baseline-agg_results.csv" \
  "$OUT/candidate-run/agg_results.csv" \
  "$OUT/results.json" \
  "$OUT/runtime.json" \
  > "$OUT/artifact-sha256.txt"
upload_outputs
