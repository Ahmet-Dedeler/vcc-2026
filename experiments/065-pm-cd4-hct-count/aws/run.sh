#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp065-pm-cd4-hct-count/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp065-pm-cd4-hct-count/${RUN_ID}}"
EXP061_REPO="${EXP061_REPO:-/work/exp061/repo}"
EXP061_EXIT_FILE="${EXP061_EXIT_FILE:-/work/exp061/exit-code.txt}"
WORK="${WORK:-/work/exp065}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

mkdir -p "$WORK" "$OUT"
rm -f "$EXIT_FILE"

# Hold the experiment-061 lock for this complete build+score, so a retry of
# experiment 061 cannot overlap. The separate global lock protects other
# memory-heavy VCC work that follows the same convention.
exec 8>/work/exp061.lock
if ! flock -n 8; then
  echo "experiment 061 still holds /work/exp061.lock" >&2
  exit 71
fi
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 72
fi

if pgrep -af \
  'cell-eval2|build_candidate.py|build_joint_candidate.py|build_submission.py' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing build or scorer:" >&2
  cat "${WORK}/preexisting-heavy-processes.txt" >&2
  exit 73
fi

exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUTPUT_PREFIX}/cloud.log" \
    --only-show-errors || true
  if [[ -s "$EXIT_FILE" ]]; then
    aws s3 cp "$EXIT_FILE" "s3://${BUCKET}/${OUTPUT_PREFIX}/exit-code.txt" \
      --only-show-errors || true
  fi
  if [[ -d "$OUT" ]]; then
    aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
      --only-show-errors || true
  fi
}

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "$status" > "$EXIT_FILE"
  (
    cd "$OUT"
    find . -type f ! -name artifact-sha256.txt -print0 |
      sort -z |
      xargs -0 sha256sum > artifact-sha256.txt
  ) || true
  upload_outputs
  exit "$status"
}
trap finish EXIT

exp061_exit="$(tr -d '[:space:]' < "$EXP061_EXIT_FILE")"
if [[ "$exp061_exit" != "0" ]]; then
  echo "experiment 061 exit code is ${exp061_exit}, not 0" >&2
  exit 74
fi

for required in \
  "$EXP061_REPO/data/harness/orion-hct116.h5ad" \
  "$EXP061_REPO/data/challenge/pert_counts.csv" \
  "$EXP061_REPO/data/challenge/gene_names.csv" \
  "$EXP061_REPO/data/donors/k562_gwps/delta.npy" \
  "$EXP061_REPO/data/donors/orion_hek293t/delta.npy" \
  "$EXP061_REPO/data/external/perturbmatch/delta.npy" \
  "$EXP061_REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.npz" \
  "$EXP061_REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.json" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/build-report.json" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/candidate.h5ad" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/candidate.h5ad.sha256" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/results.json" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/candidate-run/agg_results.csv" \
  "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/baseline-agg_results.csv"
do
  test -s "$required"
done
(
  cd "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output"
  sha256sum -c candidate.h5ad.sha256
)

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

dnf install -y gzip tar time
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

rm -rf "$REPO"
mkdir -p "$REPO" "$INPUT/cd4" "$OUT"
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"
aws s3 cp \
  "s3://${BUCKET}/input/exp060-cd4-union-k562/cd4_compact.tar.gz" \
  "${INPUT}/cd4_compact.tar.gz" --only-show-errors
tar -xzf "${INPUT}/cd4_compact.tar.gz" -C "${INPUT}/cd4"

cd "$REPO"
uv sync --frozen
uv run python experiments/065-pm-cd4-hct-count/test_candidate.py \
  2>&1 | tee "$OUT/tests.log"

uv run python - <<'PY' > "$OUT/runtime.json"
import importlib.metadata
import json
import platform
print(json.dumps({
    "python": platform.python_version(),
    "cell_eval2": importlib.metadata.version("cell-eval2"),
}, indent=2))
PY

cp "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/baseline-agg_results.csv" \
  "$OUT/baseline-agg_results.csv"
cp "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/results.json" \
  "$OUT/exp061-results.json"
cp "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/candidate-run/agg_results.csv" \
  "$OUT/exp061-candidate-agg_results.csv"

rm -f "$OUT/candidate.h5ad" "$OUT/candidate.h5ad.sha256" \
  "$OUT/build-report.json"
uv run python experiments/065-pm-cd4-hct-count/memory_guard.py \
  --max-rss-gib 32 \
  --poll-seconds 2 \
  --report "$OUT/build-memory-guard.json" \
  -- \
  /usr/bin/time -v uv run python \
    experiments/065-pm-cd4-hct-count/build_joint_candidate.py \
    --real "$EXP061_REPO/data/harness/orion-hct116.h5ad" \
    --k562 "$EXP061_REPO/data/donors/k562_gwps" \
    --hek293t "$EXP061_REPO/data/donors/orion_hek293t" \
    --perturbmatch "$EXP061_REPO/data/external/perturbmatch" \
    --cd4-root "${INPUT}/cd4/cd4_compact" \
    --panel-csv "$EXP061_REPO/data/challenge/pert_counts.csv" \
    --genes-csv "$EXP061_REPO/data/challenge/gene_names.csv" \
    --whitener "$EXP061_REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.npz" \
    --whitener-metadata "$EXP061_REPO/experiments/040-union-gamma/out_hct_clean/whitener_union.json" \
    --exp061-build-report "$EXP061_REPO/experiments/061-perturbmatch-hct-count/cloud-output/build-report.json" \
    --output "$OUT/candidate.h5ad" \
    --report "$OUT/build-report.json" \
    --max-rss-gib 28 \
  2>&1 | tee "$OUT/build.log"

sha256sum "$OUT/candidate.h5ad" > "$OUT/candidate.h5ad.sha256"
upload_outputs

# The candidate build has exited before the one and only scorer starts.
rm -rf "$OUT/candidate-run"
uv run python experiments/065-pm-cd4-hct-count/memory_guard.py \
  --max-rss-gib 55 \
  --poll-seconds 2 \
  --report "$OUT/score-memory-guard.json" \
  -- \
  /usr/bin/time -v uvx --python 3.12 --from "cell-eval2==0.16.0" \
    cell-eval2 run \
    -ap "$OUT/candidate.h5ad" \
    -ar "$EXP061_REPO/data/harness/orion-hct116.h5ad" \
    -o "$OUT/candidate-run" \
    --preset vcc2026 \
    --pert-col target_gene \
    --control non-targeting \
  2>&1 | tee "$OUT/score.log"

uv run python experiments/065-pm-cd4-hct-count/compare_scores.py \
  --baseline "$OUT/baseline-agg_results.csv" \
  --candidate "$OUT/candidate-run/agg_results.csv" \
  --output "$OUT/results.json" \
  2>&1 | tee "$OUT/compare.log"

upload_outputs
