#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
BASE_RUN_ID="${BASE_RUN_ID:?BASE_RUN_ID is required}"
RETRY_ID="${RETRY_ID:?RETRY_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:?SOURCE_PREFIX is required}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:?OUTPUT_PREFIX is required}"
NULL_REPS="${NULL_REPS:-1000}"
MAX_RSS_GIB="${MAX_RSS_GIB:-55}"
BASE_WORK="/work/exp070-${BASE_RUN_ID}"
BASE_INPUT="${BASE_WORK}/input"
COMPACT="${BASE_WORK}/output/viperturb_compact"
RETRY_WORK="${BASE_WORK}/evaluation-retry-${RETRY_ID}"
RETRY_REPO="${RETRY_WORK}/repo"
OUT="${RETRY_WORK}/output"
LOG="${OUT}/evaluate-only.log"
EXIT_FILE="${OUT}/exit-code.txt"

if [[ -e "$RETRY_WORK" ]]; then
  echo "refusing to overwrite existing retry directory ${RETRY_WORK}" >&2
  exit 66
fi
mkdir -p "$RETRY_REPO" "$OUT"

exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi
if pgrep -af \
  '[c]ell-eval2|[b]uild_candidate.py|[b]uild_joint_candidate.py|[e]xtract_viperturb.R|[e]xperiments/070-viperturb-k562/run.py' \
  >"${RETRY_WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC build, scorer, or exp070 run" >&2
  cat "${RETRY_WORK}/preexisting-heavy-processes.txt" >&2
  exit 72
fi

exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
    --only-show-errors || true
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

required_paths=(
  "${BASE_WORK}/repo/.venv/bin/python"
  "${BASE_INPUT}/donors/k562_gwps/delta.npy"
  "${BASE_INPUT}/donors/orion_hct116/delta.npy"
  "${BASE_INPUT}/donors/orion_hek293t/delta.npy"
  "${BASE_INPUT}/perturbmatch/delta.npy"
  "${BASE_INPUT}/cd4/cd4_compact"
  "${BASE_INPUT}/whiteners"
  "${COMPACT}/meta.json"
)
for path in "${required_paths[@]}"; do
  if [[ ! -e "$path" ]]; then
    echo "missing preserved evaluation input: ${path}" >&2
    exit 73
  fi
done

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${RETRY_WORK}/source.tar.gz" --only-show-errors
if tar -tzf "${RETRY_WORK}/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata files" >&2
  exit 77
fi
tar -xzf "${RETRY_WORK}/source.tar.gz" -C "$RETRY_REPO"

PYTHON="${BASE_WORK}/repo/.venv/bin/python"
cd "$RETRY_REPO"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

PYTHONDONTWRITEBYTECODE=1 "$PYTHON" \
  experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib 1.0 \
  --poll-seconds 0.25 \
  --report "$OUT/test-memory-guard.json" \
  -- \
  "$PYTHON" experiments/070-viperturb-k562/test_run.py \
  2>&1 | tee "$OUT/tests.log"

TOKEN="$(curl -fsS -X PUT \
  -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600' \
  http://169.254.169.254/latest/api/token)"
INSTANCE_ID="$(curl -fsS \
  -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-id)"
INSTANCE_TYPE="$(curl -fsS \
  -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-type)"
MARKET="$(curl -fsS \
  -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-life-cycle || true)"

sha256sum \
  experiments/070-viperturb-k562/run.py \
  experiments/070-viperturb-k562/test_run.py \
  experiments/070-viperturb-k562/memory_guard.py \
  experiments/070-viperturb-k562/aws/evaluate_only.sh \
  experiments/070-viperturb-k562/aws/evaluate_only_remote.sh \
  > "$OUT/code.sha256"

"$PYTHON" experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib "$MAX_RSS_GIB" \
  --poll-seconds 2 \
  --report "$OUT/evaluate-memory-guard.json" \
  -- \
  /usr/bin/time -v "$PYTHON" \
    experiments/070-viperturb-k562/run.py \
    --k562 "${BASE_INPUT}/donors/k562_gwps" \
    --hct116 "${BASE_INPUT}/donors/orion_hct116" \
    --hek293t "${BASE_INPUT}/donors/orion_hek293t" \
    --perturbmatch "${BASE_INPUT}/perturbmatch" \
    --viperturb "$COMPACT" \
    --cd4-root "${BASE_INPUT}/cd4/cd4_compact" \
    --panel-csv data/challenge/pert_counts.csv \
    --genes-csv data/challenge/gene_names.csv \
    --whitener-dir "${BASE_INPUT}/whiteners" \
    --output "$OUT/results.json" \
    --max-rss-gib "$MAX_RSS_GIB" \
    --null-reps "$NULL_REPS" \
    --minimum-modified-targets 70 \
    --minimum-non-low-targets 20 \
    --instance-id "$INSTANCE_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --market "${MARKET:-unknown}" \
    --s3-prefix "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

upload_outputs
