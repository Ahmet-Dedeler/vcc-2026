#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp069-perturbmatch-aligned-whitener/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp069-perturbmatch-aligned-whitener/${RUN_ID}}"
WORK="${WORK:-/work/exp069}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

mkdir -p "$WORK" "$OUT"
rm -f "$EXIT_FILE"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi

if pgrep -af \
  '[c]ell-eval2|[b]uild_candidate.py|[b]uild_joint_candidate.py|[e]xperiments/069-perturbmatch-aligned-whitener/run.py' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC build, scorer, or exp069 run" >&2
  cat "${WORK}/preexisting-heavy-processes.txt" >&2
  exit 72
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

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

dnf install -y gzip tar time
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

rm -rf "$REPO" "$INPUT"
mkdir -p \
  "$REPO" \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/baseline-whiteners" \
  "$INPUT/exp061-gate" \
  "$OUT"

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"

for donor in k562_gwps orion_hct116 orion_hek293t
do
  mkdir -p "${INPUT}/donors/${donor}"
  aws s3 cp \
    "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "${INPUT}/donors/${donor}/" \
    --recursive --only-show-errors
done
aws s3 cp \
  "s3://${BUCKET}/input/exp061-perturbmatch-hct-count/perturbmatch/" \
  "${INPUT}/perturbmatch/" \
  --recursive --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/exp060-cd4-union-k562/cd4_compact.tar.gz" \
  "${INPUT}/cd4_compact.tar.gz" --only-show-errors
tar -xzf "${INPUT}/cd4_compact.tar.gz" -C "${INPUT}/cd4"
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/baseline-whiteners/" \
  "${INPUT}/baseline-whiteners/" \
  --recursive --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/exp061-gate/" \
  "${INPUT}/exp061-gate/" \
  --recursive --only-show-errors

cd "$REPO"
uv sync --frozen

# Defense in depth: validate the same signed-off experiment-061 artifacts
# again on the worker before tests or any eigendecomposition.
uv run python \
  experiments/069-perturbmatch-aligned-whitener/validate_exp061_gate.py \
  --results "${INPUT}/exp061-gate/results.json" \
  --baseline-agg "${INPUT}/exp061-gate/baseline-agg_results.csv" \
  --candidate-agg "${INPUT}/exp061-gate/candidate-agg_results.csv" \
  --build-report "${INPUT}/exp061-gate/build-report.json" \
  --candidate-checksum "${INPUT}/exp061-gate/candidate.h5ad.sha256" \
  --exit-code "${INPUT}/exp061-gate/exit-code.txt" \
  --candidate-object-bytes "$(
    tr -d '[:space:]' < "${INPUT}/exp061-gate/candidate-object-bytes.txt"
  )" \
  | tee "$OUT/exp061-gate.json"

PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/069-perturbmatch-aligned-whitener/memory_guard.py \
  --max-rss-gib 1.0 \
  --poll-seconds 0.25 \
  --report "$OUT/test-memory-guard.json" \
  -- \
  uv run python experiments/069-perturbmatch-aligned-whitener/test_run.py \
  2>&1 | tee "$OUT/tests.log"

TOKEN="$(curl -fsS -X PUT \
  -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600' \
  http://169.254.169.254/latest/api/token)"
INSTANCE_ID="$(curl -fsS \
  -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-id)"
MARKET="$(curl -fsS \
  -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-life-cycle || true)"

sha256sum \
  experiments/069-perturbmatch-aligned-whitener/run.py \
  experiments/069-perturbmatch-aligned-whitener/test_run.py \
  experiments/069-perturbmatch-aligned-whitener/validate_exp061_gate.py \
  experiments/069-perturbmatch-aligned-whitener/memory_guard.py \
  experiments/069-perturbmatch-aligned-whitener/aws/run.sh \
  > "$OUT/code.sha256"

rm -rf "$OUT/whiteners" "$WORK/scratch"
uv run python \
  experiments/069-perturbmatch-aligned-whitener/memory_guard.py \
  --max-rss-gib 55 \
  --poll-seconds 2 \
  --report "$OUT/memory-guard.json" \
  -- \
  /usr/bin/time -v uv run python \
    experiments/069-perturbmatch-aligned-whitener/run.py \
    --k562 "${INPUT}/donors/k562_gwps" \
    --hct116 "${INPUT}/donors/orion_hct116" \
    --hek293t "${INPUT}/donors/orion_hek293t" \
    --perturbmatch "${INPUT}/perturbmatch" \
    --cd4-root "${INPUT}/cd4/cd4_compact" \
    --panel-csv data/challenge/pert_counts.csv \
    --genes-csv data/challenge/gene_names.csv \
    --baseline-whitener-dir "${INPUT}/baseline-whiteners" \
    --candidate-whitener-dir "$OUT/whiteners" \
    --scratch "$WORK/scratch" \
    --output "$OUT/results.json" \
    --max-rss-gib 55 \
    --instance-id "$INSTANCE_ID" \
    --instance-type "r7i.2xlarge" \
    --market "${MARKET:-unknown}" \
    --s3-prefix "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
    --guarded \
  2>&1 | tee "$OUT/run.log"

upload_outputs
