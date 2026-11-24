#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
EXPECTED_SOURCE_SHA256="${EXPECTED_SOURCE_SHA256:?EXPECTED_SOURCE_SHA256 is required}"
EXPECTED_RUNNER_SHA256="${EXPECTED_RUNNER_SHA256:?EXPECTED_RUNNER_SHA256 is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp082-biological-target-kernel/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp082-biological-target-kernel/${RUN_ID}}"
WORK_ROOT="${WORK_ROOT:-/work/exp082}"
if [[ ! "$RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 63
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "required command is missing: python3" >&2
  exit 69
fi
CANONICAL_WORK_ROOT="$(python3 -c \
  'import os, sys; print(os.path.realpath(sys.argv[1]))' \
  "$WORK_ROOT")"
if [[ "$WORK_ROOT" != "$CANONICAL_WORK_ROOT" || "$CANONICAL_WORK_ROOT" == "/" ]]; then
  echo "WORK_ROOT must be an absolute canonical non-root path: ${WORK_ROOT}" >&2
  exit 70
fi
WORK="$(python3 -c \
  'import os, sys; print(os.path.realpath(sys.argv[1]))' \
  "${CANONICAL_WORK_ROOT}/${RUN_ID}")"
if [[ "$WORK" != "${CANONICAL_WORK_ROOT}/${RUN_ID}" ]]; then
  echo "non-canonical per-run work path: ${WORK}" >&2
  exit 70
fi
case "$WORK" in
  "${CANONICAL_WORK_ROOT}/"*) ;;
  *)
    echo "per-run work path escapes WORK_ROOT: ${WORK}" >&2
    exit 70
    ;;
esac
if [[ -e "$WORK" || -L "$WORK" ]]; then
  echo "refusing existing local run directory: ${WORK}" >&2
  exit 70
fi
WORK_ROOT="$CANONICAL_WORK_ROOT"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

for command in aws awk curl df find flock free grep pgrep python3 sha256sum systemctl tar tee; do
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

output_key_count="$(aws s3api list-objects-v2 \
  --bucket "$BUCKET" \
  --prefix "${OUTPUT_PREFIX%/}/" \
  --max-keys 1 \
  --query 'KeyCount' \
  --output text)"
if [[ ! "$output_key_count" =~ ^[0-9]+$ || "$output_key_count" -ne 0 ]]; then
  echo "refusing reused or unverifiable output prefix: s3://${BUCKET}/${OUTPUT_PREFIX}/" >&2
  exit 78
fi

mkdir -p "$WORK_ROOT"
if [[ "$(python3 -c \
  'import os, sys; print(os.path.realpath(sys.argv[1]))' \
  "$WORK_ROOT")" != "$CANONICAL_WORK_ROOT" ]]; then
  echo "WORK_ROOT changed during validation: ${WORK_ROOT}" >&2
  exit 70
fi
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
    grep -v "^vcc-exp082-${RUN_ID}\\.service$" || true
)"
if [[ -n "$active_heavy_units" ]]; then
  echo "refusing to overlap an active VCC/system scorer service" >&2
  printf '%s\n' "$active_heavy_units" >&2
  exit 77
fi

if ! mkdir "$WORK"; then
  echo "refusing existing or uncreatable local run directory: ${WORK}" >&2
  exit 70
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
  echo "worker has ${available_mem_kib:-0} KiB MemAvailable; require >=48 GiB" >&2
  exit 76
fi
free_kib="$(df -Pk "$WORK_ROOT" | awk 'NR==2 {print $4}')"
minimum_free_kib=$((40 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  echo "worker has ${free_kib:-0} KiB free; require at least 40 GiB" >&2
  exit 75
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-16}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-16}"

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

mkdir -p \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/truth-cache" \
  "$INPUT/kolf2" \
  "$INPUT/operators" \
  "$INPUT/embedding" \
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
sha256sum -c experiments/082-biological-target-kernel/dependency-sha256.txt

for donor in k562_gwps orion_hct116 orion_hek293t rpe1 nadig_jurkat; do
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
  "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_kolf2.npz" \
  "$INPUT/operators/held_kolf2.npz" \
  --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_kolf2.json" \
  "$INPUT/operators/held_kolf2.json" \
  --only-show-errors

# Download exactly one Hugging Face file. No dataset snapshot or expression
# artifact is requested.
EMBEDDING="$INPUT/embedding/ESM2_pert_features.parquet"
curl -fL --retry 5 --retry-all-errors \
  "https://huggingface.co/datasets/altoslabs/primeflow-vcc-datasets/resolve/d773a3476d15e5d32d385b13cebf87e62580097e/ESM2_pert_features.parquet" \
  -o "$EMBEDDING"
printf '%s  %s\n' \
  'f05fe94e3a592bfa0973b3816b956af8c0a99387484d4c68ee89e0eac0b21db3' \
  "$EMBEDDING" | sha256sum -c -
if find "$INPUT" -type f \( -iname '*.h5ad' -o -iname '*.h5ad.gz' \) |
  grep -q .
then
  echo "forbidden expression bundle was downloaded" >&2
  exit 77
fi

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy --with pyarrow python \
  experiments/082-biological-target-kernel/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.2 \
  --report "$OUT/test-memory-guard.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy --with pyarrow python \
    experiments/082-biological-target-kernel/test_run.py \
  2>&1 | tee "$OUT/tests.log"

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy --with pyarrow python \
  experiments/082-biological-target-kernel/memory_guard.py \
  --max-rss-gib 48 \
  --poll-seconds 2 \
  --report "$OUT/stage-memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy --with pyarrow python \
    experiments/082-biological-target-kernel/stage_cache.py \
    --donor-root "$INPUT/donors" \
    --perturbmatch "$INPUT/perturbmatch" \
    --cd4-root "$INPUT/cd4/cd4_compact" \
    --truth-cache "$INPUT/truth-cache" \
    --operator-dir "$INPUT/operators" \
    --exp064-results \
      experiments/064-pm-cd4-composition/cloud-output/results.json \
    --artifact-manifest \
      experiments/082-biological-target-kernel/artifact_manifest.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output "$OUT/arm-d-query-cache" \
    --exp064-atol 1e-6 \
    --max-rss-gib 48 \
    --guarded \
  2>&1 | tee "$OUT/stage.log"

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy --with pyarrow python \
  experiments/082-biological-target-kernel/memory_guard.py \
  --max-rss-gib 24 \
  --poll-seconds 1 \
  --report "$OUT/evaluate-memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy --with pyarrow python \
    experiments/082-biological-target-kernel/run.py \
    --truth-cache "$INPUT/truth-cache" \
    --donor-root "$INPUT/donors" \
    --kolf-root "$INPUT/kolf2" \
    --query-cache "$OUT/arm-d-query-cache" \
    --embedding-parquet "$EMBEDDING" \
    --artifact-manifest \
      experiments/082-biological-target-kernel/artifact_manifest.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output-dir "$OUT/evaluation" \
    --max-rss-gib 24 \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

sha256sum \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/082-biological-target-kernel/stage_cache.py \
  experiments/082-biological-target-kernel/run.py \
  experiments/082-biological-target-kernel/test_run.py \
  experiments/082-biological-target-kernel/memory_guard.py \
  experiments/082-biological-target-kernel/artifact_manifest.json \
  experiments/082-biological-target-kernel/dependency-sha256.txt \
  experiments/082-biological-target-kernel/aws/launch.sh \
  experiments/082-biological-target-kernel/aws/run.sh \
  >"$OUT/code.sha256"

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
  "available_memory_kib_at_worker_start": ${available_mem_kib},
  "source_archive_sha256": "${EXPECTED_SOURCE_SHA256}",
  "runner_sha256": "${EXPECTED_RUNNER_SHA256}",
  "stage_guard_gib": 48,
  "evaluation_guard_gib": 24,
  "downloaded_only_esm2_parquet": true,
  "expression_bundle_downloaded": false,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

upload_outputs
