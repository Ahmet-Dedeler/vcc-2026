#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
EXPECTED_SOURCE_SHA256="${EXPECTED_SOURCE_SHA256:?EXPECTED_SOURCE_SHA256 is required}"
EXPECTED_RUNNER_SHA256="${EXPECTED_RUNNER_SHA256:?EXPECTED_RUNNER_SHA256 is required}"
EXP085_ONE_TIME_CLEARANCE="${EXP085_ONE_TIME_CLEARANCE:?clearance is required}"
EXP085_DURABLE_MARKER_URI="${EXP085_DURABLE_MARKER_URI:?durable marker URI is required}"
EXP085_DURABLE_FREEZE_URI="${EXP085_DURABLE_FREEZE_URI:?durable freeze URI is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp085-kolf-pm-gate/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp085-kolf-pm-gate/${RUN_ID}}"
WORK_ROOT="${WORK_ROOT:-/work/exp085}"
if [[ "$EXP085_ONE_TIME_CLEARANCE" != "YES" ]]; then
  echo "invalid one-time clearance" >&2
  exit 64
fi
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 63
fi
WORK="${WORK_ROOT}/${RUN_ID}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${REPO}/experiments/085-kolf-pm-gate/cloud-output/${RUN_ID}"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

for command in aws awk curl df find flock free grep pgrep sha256sum systemctl tar tee; do
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
    grep -v "^vcc-exp085-${RUN_ID}\\.service$" || true
)"
if [[ -n "$active_heavy_units" ]]; then
  echo "refusing to overlap an active VCC/system scorer service" >&2
  printf '%s\n' "$active_heavy_units" >&2
  exit 77
fi
if [[ -e "$WORK" ]]; then
  echo "refusing reused work path ${WORK}" >&2
  exit 78
fi
remote_output_count="$(aws s3api list-objects-v2 \
  --bucket "$BUCKET" \
  --prefix "${OUTPUT_PREFIX}/" \
  --max-keys 1 \
  --query KeyCount \
  --output text)"
if [[ "$remote_output_count" != "0" ]]; then
  echo "refusing reused output prefix s3://${BUCKET}/${OUTPUT_PREFIX}/" >&2
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
      xargs -0 sha256sum >"${WORK}/artifact-sha256.tmp"
    mv "${WORK}/artifact-sha256.tmp" artifact-sha256.txt
  ) || true
  upload_outputs
  exit "$status"
}
trap finish EXIT

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
minimum_free_kib=$((30 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  echo "worker has ${free_kib:-0} KiB free; require at least 30 GiB" >&2
  exit 75
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-16}"

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

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
sha256sum -c experiments/085-kolf-pm-gate/dependency-sha256.txt

mkdir -p \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/held-operators" \
  "$INPUT/production-operator" \
  "$INPUT/kolf2"

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

for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/${SOURCE_PREFIX}/held-operators/held_${held}.${suffix}" \
      "$INPUT/held-operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done
for suffix in npz json; do
  aws s3 cp \
    "s3://${BUCKET}/input/production-challenge/whitener_union.${suffix}" \
    "$INPUT/production-operator/whitener_union.${suffix}" \
    --only-show-errors
done

# Stage the compact KOLF artifact without hashing, importing, or parsing it.
# The Python runner authenticates and opens it only after the pre-KOLF freeze.
for filename in \
  b_ctrl.npy \
  delta.npy \
  gene_map.csv \
  genes.txt \
  meta.json \
  perts.txt \
  provenance.json \
  source_rows_int8.npy \
  target_row_map.csv
do
  aws s3 cp \
    "s3://${BUCKET}/${SOURCE_PREFIX}/kolf2/${filename}" \
    "$INPUT/kolf2/${filename}" \
    --only-show-errors
done

PYTHON_ENV=(
  uv run --no-project --python 3.12
  --with numpy==2.5.2
  --with scipy==1.18.0
  --with anndata==0.13.2
  --with h5py==3.16.0
  --with pandas==3.0.5
)

PYTHONDONTWRITEBYTECODE=1 \
"${PYTHON_ENV[@]}" python experiments/085-kolf-pm-gate/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.05 \
  --report "$OUT/test-memory.json" \
  -- \
  "${PYTHON_ENV[@]}" python experiments/085-kolf-pm-gate/test_run.py \
  2>&1 | tee "$OUT/tests.log"

sha256sum \
  src/vcc2026/donor_transfer.py \
  src/vcc2026/multi_donor.py \
  src/vcc2026/whiten.py \
  scripts/build_submission.py \
  scripts/make_submission.py \
  scripts/make_submission_pm_cd4.py \
  experiments/061-perturbmatch-hct-count/build_candidate.py \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv \
  experiments/066-leakage-safe-target-gate/RESULTS.md \
  experiments/066-leakage-safe-target-gate/metrics.json \
  experiments/085-kolf-pm-gate/run.py \
  experiments/085-kolf-pm-gate/test_run.py \
  experiments/085-kolf-pm-gate/memory_guard.py \
  experiments/085-kolf-pm-gate/preregistered_protocol.json \
  experiments/085-kolf-pm-gate/artifact_manifest.json \
  experiments/085-kolf-pm-gate/dependency-sha256.txt \
  experiments/085-kolf-pm-gate/aws/run.sh \
  >"$OUT/code.sha256"

PYTHONDONTWRITEBYTECODE=1 \
EXP085_ONE_TIME_CLEARANCE=YES \
EXP085_DURABLE_MARKER_URI="$EXP085_DURABLE_MARKER_URI" \
EXP085_DURABLE_FREEZE_URI="$EXP085_DURABLE_FREEZE_URI" \
"${PYTHON_ENV[@]}" python experiments/085-kolf-pm-gate/memory_guard.py \
  --max-rss-gib 40 \
  --poll-seconds 1 \
  --report "$OUT/evaluate-memory.json" \
  -- \
  /usr/bin/time -v \
  "${PYTHON_ENV[@]}" python experiments/085-kolf-pm-gate/run.py \
    --k562 "$INPUT/donors/k562_gwps" \
    --hct116 "$INPUT/donors/orion_hct116" \
    --hek293t "$INPUT/donors/orion_hek293t" \
    --perturbmatch "$INPUT/perturbmatch" \
    --cd4-root "$INPUT/cd4/cd4_compact" \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --held-whitener-dir "$INPUT/held-operators" \
    --production-whitener "$INPUT/production-operator/whitener_union.npz" \
    --production-whitener-metadata "$INPUT/production-operator/whitener_union.json" \
    --exp064-results experiments/064-pm-cd4-composition/cloud-output/results.json \
    --exp064-per-target experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv \
    --kolf "$INPUT/kolf2" \
    --artifact-manifest experiments/085-kolf-pm-gate/artifact_manifest.json \
    --output-dir "$OUT/result" \
    --max-rss-gib 40 \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

if [[ -f "$OUT/result/results.json" ]]; then
  cp "$OUT/result/results.json" "$OUT/results.json"
fi
if [[ -f "$OUT/result/per_target_pds.csv" ]]; then
  cp "$OUT/result/per_target_pds.csv" "$OUT/per_target_pds.csv"
fi
if [[ -f "$OUT/result/KOLF_ACCESS_CONSUMED.json" ]]; then
  cp "$OUT/result/KOLF_ACCESS_CONSUMED.json" "$OUT/KOLF_ACCESS_CONSUMED.json"
fi

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

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${INSTANCE_TYPE}",
  "market": "${MARKET:-unknown}",
  "source_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "worker_was_created_or_started": false,
  "minimum_available_memory_gib": 48,
  "systemd_memory_max_gib": 44,
  "systemd_swap_max_gib": 0,
  "one_time_kolf_evaluation": true,
  "count_generation": false,
  "submission": false
}
EOF

upload_outputs
