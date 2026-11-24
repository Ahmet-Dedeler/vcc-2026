#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
EXPECTED_SOURCE_SHA256="${EXPECTED_SOURCE_SHA256:?EXPECTED_SOURCE_SHA256 is required}"
EXPECTED_RUNNER_SHA256="${EXPECTED_RUNNER_SHA256:?EXPECTED_RUNNER_SHA256 is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp080-listwise-metric-adapter/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp080-listwise-metric-adapter/${RUN_ID}}"
WORK_ROOT="${WORK_ROOT:-/work/exp080}"
PYTHON_VERSION="3.12.12"
NUMPY_VERSION="2.5.2"
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 63
fi
WORK="${WORK_ROOT}/${RUN_ID}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

for command in aws awk curl df flock free grep pgrep sha256sum systemctl tar tee; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required command is missing: ${command}" >&2
    exit 69
  fi
done
if [[ ! -x /usr/bin/time ]]; then
  if ! command -v dnf >/dev/null 2>&1; then
    echo "/usr/bin/time is missing and dnf is unavailable" >&2
    exit 69
  fi
  dnf install -y time
fi

actual_runner_sha256="$(sha256sum "$0" | awk '{print $1}')"
if [[ "$actual_runner_sha256" != "$EXPECTED_RUNNER_SHA256" ]]; then
  echo "runner SHA-256 ${actual_runner_sha256} != ${EXPECTED_RUNNER_SHA256}" >&2
  exit 76
fi

mkdir -p "$WORK_ROOT"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi
if pgrep -af \
  '[c]ell-eval|[p]ython(3)? .*experiments/|[u]v (run|tool|x) |[t]orchrun|[a]ccelerate launch|[r]ay::|[m]ake_submission|[b]uild_submission|[r]un_with_memory_guard|[s]tage_cache.py' \
  >"${WORK_ROOT}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC heavy process" >&2
  cat "${WORK_ROOT}/preexisting-heavy-processes.txt" >&2
  exit 72
fi
active_heavy_units="$(
  systemctl list-units --type=service --state=running --no-legend |
    awk '{print $1}' |
    grep -E '^(vcc-|cell-eval)' |
    grep -v "^vcc-exp080-${RUN_ID}\\.service$" || true
)"
if [[ -n "$active_heavy_units" ]]; then
  echo "refusing to overlap an active VCC/system scorer service" >&2
  printf '%s\n' "$active_heavy_units" >&2
  exit 77
fi

if [[ "$WORK" != "${WORK_ROOT}/${RUN_ID}" || "$WORK" == "/" ]]; then
    echo "unsafe per-run work path: ${WORK}" >&2
    exit 70
fi
if [[ -e "$WORK" ]]; then
  echo "refusing to reuse existing work directory ${WORK}" >&2
  exit 79
fi
mkdir -p "$REPO" "$INPUT" "$OUT"
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
  printf '%s\n' "$status" >"$EXIT_FILE"
  (
    cd "$OUT"
    find . -type f ! -name artifact-sha256.txt -print0 |
      sort -z |
      xargs -0 sha256sum >artifact-sha256.txt
  ) || true
  upload_outputs
  exit "$status"
}
trap finish EXIT

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
visible_mem_kib="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
available_mem_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
minimum_visible_kib=$((60 * 1024 * 1024))
minimum_available_kib=$((48 * 1024 * 1024))
if [[ -z "$visible_mem_kib" || "$visible_mem_kib" -lt "$minimum_visible_kib" ]]; then
  echo "worker exposes ${visible_mem_kib:-0} KiB RAM; require nominal >=64 GiB" >&2
  exit 74
fi
if [[ -z "$available_mem_kib" || "$available_mem_kib" -lt "$minimum_available_kib" ]]; then
  echo "worker has ${available_mem_kib:-0} KiB available RAM; require >=48 GiB" >&2
  exit 78
fi
free_kib="$(df -Pk "$WORK_ROOT" | awk 'NR==2 {print $4}')"
minimum_free_kib=$((35 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  echo "worker has ${free_kib:-0} KiB free; require at least 35 GiB" >&2
  exit 75
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-8}"

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install "$PYTHON_VERSION"
uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" \
  python -c \
  "import platform,numpy; assert platform.python_version() == '${PYTHON_VERSION}'; assert numpy.__version__ == '${NUMPY_VERSION}'"

mkdir -p \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/truth-cache" \
  "$INPUT/kolf2" \
  "$INPUT/operators" \
  "$OUT/arm-d-query-cache" \
  "$OUT/evaluation"

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "$WORK/source.tar.gz" --only-show-errors
printf '%s  %s\n' "$EXPECTED_SOURCE_SHA256" "$WORK/source.tar.gz" |
  sha256sum --check -
if tar -tzf "$WORK/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata" >&2
  exit 73
fi
tar -xzf "$WORK/source.tar.gz" -C "$REPO"
cd "$REPO"
sha256sum -c experiments/080-listwise-metric-adapter/dependency-sha256.txt

for donor in k562_gwps orion_hct116 orion_hek293t; do
  mkdir -p "$INPUT/donors/$donor"
  aws s3 cp \
    "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "$INPUT/donors/$donor/" \
    --recursive \
    --exclude '*' \
    --include 'delta.npy' \
    --include 'b_ctrl.npy' \
    --include 'genes.txt' \
    --include 'perts.txt' \
    --include 'meta.json' \
    --only-show-errors
done

aws s3 cp \
  "s3://${BUCKET}/input/exp061-perturbmatch-hct-count/perturbmatch/" \
  "$INPUT/perturbmatch/" \
  --recursive \
  --exclude '*' \
  --include 'delta.npy' \
  --include 'genes.txt' \
  --include 'perts.txt' \
  --include 'meta.json' \
  --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/exp060-cd4-union-k562/cd4_compact.tar.gz" \
  "$INPUT/cd4_compact.tar.gz" --only-show-errors
tar -xzf "$INPUT/cd4_compact.tar.gz" -C "$INPUT/cd4"

for filename in \
  meta.json \
  master_genes.txt \
  train_perts.txt \
  outer_perts.txt \
  train_mask.npy \
  rpe1_raw_train.npy \
  nadig_jurkat_raw_train.npy \
  orion_hct116_raw_outer.npy \
  orion_hek293t_raw_outer.npy
do
  aws s3 cp \
    "s3://${BUCKET}/input/exp052-cache/${filename}" \
    "$INPUT/truth-cache/${filename}" \
    --only-show-errors
done
for filename in delta.npy b_ctrl.npy genes.txt perts.txt meta.json; do
  aws s3 cp \
    "s3://${BUCKET}/input/exp079-kolf-local-centering/20260823T064635Z/kolf2/${filename}" \
    "$INPUT/kolf2/${filename}" \
    --only-show-errors
done
for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_${held}.${suffix}" \
      "$INPUT/operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.npz" \
  "$INPUT/operators/held_kolf2.npz" \
  --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.json" \
  "$INPUT/operators/held_kolf2.json" \
  --only-show-errors

EXP080_REAL_DONOR_ROOT="$INPUT/donors" \
EXP080_REAL_KOLF_ROOT="$INPUT/kolf2" \
EXP080_REAL_TRUTH_CACHE="$INPUT/truth-cache" \
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.2 \
  --report "$OUT/test-memory-guard.json" \
  -- \
  uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
    experiments/080-listwise-metric-adapter/test_run.py \
  2>&1 | tee "$OUT/tests.log"

sha256sum \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/080-listwise-metric-adapter/stage_cache.py \
  experiments/080-listwise-metric-adapter/run.py \
  experiments/080-listwise-metric-adapter/test_run.py \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  experiments/080-listwise-metric-adapter/artifact_manifest.json \
  experiments/080-listwise-metric-adapter/dependency-sha256.txt \
  experiments/080-listwise-metric-adapter/aws/run.sh \
  >"$OUT/code.sha256"

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 48 \
  --poll-seconds 2 \
  --report "$OUT/stage-memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
    experiments/080-listwise-metric-adapter/stage_cache.py \
    --donor-root "$INPUT/donors" \
    --perturbmatch "$INPUT/perturbmatch" \
    --cd4-root "$INPUT/cd4/cd4_compact" \
    --truth-cache "$INPUT/truth-cache" \
    --operator-dir "$INPUT/operators" \
    --exp064-results \
      experiments/064-pm-cd4-composition/cloud-output/results.json \
    --artifact-manifest \
      experiments/080-listwise-metric-adapter/artifact_manifest.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output "$OUT/arm-d-query-cache" \
    --max-rss-gib 48 \
    --guarded \
  2>&1 | tee "$OUT/stage.log"

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 16 \
  --poll-seconds 1 \
  --report "$OUT/evaluate-memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python "$PYTHON_VERSION" --with "numpy==${NUMPY_VERSION}" python \
    experiments/080-listwise-metric-adapter/run.py \
    --donor-root "$INPUT/donors" \
    --truth-cache "$INPUT/truth-cache" \
    --kolf-root "$INPUT/kolf2" \
    --query-cache "$OUT/arm-d-query-cache" \
    --artifact-manifest \
      experiments/080-listwise-metric-adapter/artifact_manifest.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output-dir "$OUT/evaluation" \
    --max-rss-gib 16 \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${INSTANCE_TYPE}",
  "market": "${MARKET:-unknown}",
  "source_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "minimum_instance_memory_gib": 64,
  "minimum_available_memory_gib": 48,
  "source_archive_sha256": "${EXPECTED_SOURCE_SHA256}",
  "runner_sha256": "${EXPECTED_RUNNER_SHA256}",
  "stage_guard_gib": 48,
  "evaluation_guard_gib": 16,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

upload_outputs
