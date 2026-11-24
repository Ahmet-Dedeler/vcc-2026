#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp073-postwhiten-csls-distill/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp073-postwhiten-csls-distill/${RUN_ID}}"
WORK_ROOT="${WORK_ROOT:-/work/exp073}"
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

for command in aws df flock pgrep curl tar gzip sha256sum awk grep tee; do
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

mkdir -p "$WORK_ROOT"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi

if [[ "$WORK" != "${WORK_ROOT}/${RUN_ID}" || "$WORK" == "/" ]]; then
  echo "unsafe per-run work path: ${WORK}" >&2
  exit 70
fi
rm -rf "$WORK"
mkdir -p "$WORK" "$OUT"

if pgrep -af \
  '[c]ell-eval2|[e]xperiments/073-postwhiten-csls-distill/run.py' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing scorer or experiment-073 run" >&2
  cat "${WORK}/preexisting-heavy-processes.txt" >&2
  exit 72
fi

exec > >(tee -a "$LOG") 2>&1

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
minimum_visible_kib=$((60 * 1024 * 1024))
if [[ -z "$visible_mem_kib" || "$visible_mem_kib" -lt "$minimum_visible_kib" ]]; then
  echo "worker exposes ${visible_mem_kib:-0} KiB RAM; require a nominal >=64 GiB instance" >&2
  exit 74
fi
free_kib="$(df -Pk "$WORK_ROOT" | awk 'NR==2 {print $4}')"
minimum_free_kib=$((30 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  echo "worker has ${free_kib:-0} KiB free; require at least 30 GiB" >&2
  exit 75
fi

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

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

mkdir -p \
  "$REPO" \
  "$INPUT/donors" \
  "$INPUT/truth-cache" \
  "$INPUT/operators" \
  "$OUT/query-cache"

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
if tar -tzf "${WORK}/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata files" >&2
  exit 73
fi
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"

for donor in k562_gwps orion_hct116 orion_hek293t; do
  mkdir -p "${INPUT}/donors/${donor}"
  aws s3 cp \
    "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "${INPUT}/donors/${donor}/" \
    --recursive \
    --exclude '*' \
    --include 'delta.npy' \
    --include 'genes.txt' \
    --include 'perts.txt' \
    --only-show-errors
done

for filename in \
  meta.json \
  master_genes.txt \
  train_perts.txt \
  outer_perts.txt \
  train_mask.npy \
  rpe1_raw_train.npy \
  nadig_jurkat_raw_train.npy \
  k562_gwps_raw_outer.npy \
  orion_hct116_raw_outer.npy \
  orion_hek293t_raw_outer.npy
do
  aws s3 cp \
    "s3://${BUCKET}/input/exp052-cache/${filename}" \
    "${INPUT}/truth-cache/${filename}" \
    --only-show-errors
done

for held in k562_gwps orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_${held}.${suffix}" \
      "${INPUT}/operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done

cd "$REPO"
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 --with numpy python \
  experiments/073-postwhiten-csls-distill/test_run.py \
  2>&1 | tee "$OUT/tests.log"

sha256sum \
  experiments/073-postwhiten-csls-distill/run.py \
  experiments/073-postwhiten-csls-distill/test_run.py \
  experiments/073-postwhiten-csls-distill/memory_guard.py \
  experiments/073-postwhiten-csls-distill/artifact_manifest.json \
  experiments/073-postwhiten-csls-distill/aws/run.sh \
  >"$OUT/code.sha256"

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 --with numpy python \
  experiments/073-postwhiten-csls-distill/memory_guard.py \
  --max-rss-gib 48 \
  --poll-seconds 2 \
  --report "$OUT/memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy python \
    experiments/073-postwhiten-csls-distill/run.py \
    --donor-root "$INPUT/donors" \
    --truth-cache "$INPUT/truth-cache" \
    --operator-dir "$INPUT/operators" \
    --artifact-manifest \
      experiments/073-postwhiten-csls-distill/artifact_manifest.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --query-cache "$OUT/query-cache" \
    --output "$OUT/results.json" \
    --per-target-output "$OUT/per_target_pds.csv" \
    --max-rss-gib 48 \
    --guarded \
  2>&1 | tee "$OUT/run.log"

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${INSTANCE_TYPE}",
  "market": "${MARKET:-unknown}",
  "source_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "submission": false
}
EOF

upload_outputs
