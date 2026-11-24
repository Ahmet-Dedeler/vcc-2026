#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp079-kolf-local-centering/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp079-kolf-local-centering/${RUN_ID}}"
WORK="${WORK:-/work/exp079/${RUN_ID}}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${REPO}/experiments/079-kolf-local-centering/cloud-output/${RUN_ID}"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

mkdir -p /work "$WORK" "$OUT"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi

if pgrep -af \
  '[c]ell-eval2|[e]xperiments/07[0-9].*/run.py|[m]ake_submission|[f]it.*whiten' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC experiment or scorer" >&2
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

dnf install -y gzip tar time util-linux
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

mkdir -p "$REPO" "$INPUT/donors" "$INPUT/kolf2"
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
if tar -tzf "${WORK}/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata files" >&2
  exit 73
fi
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"
cd "$REPO"
sha256sum -c experiments/079-kolf-local-centering/dependency-sha256.txt

for donor in k562_gwps orion_hct116 orion_hek293t; do
  mkdir -p "${INPUT}/donors/${donor}"
  aws s3 cp \
    "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "${INPUT}/donors/${donor}/" \
    --recursive \
    --exclude '*' \
    --include 'delta.npy' \
    --include 'b_ctrl.npy' \
    --include 'genes.txt' \
    --include 'perts.txt' \
    --include 'meta.json' \
    --only-show-errors
done

for filename in delta.npy b_ctrl.npy genes.txt perts.txt meta.json; do
  aws s3 cp \
    "s3://${BUCKET}/${SOURCE_PREFIX}/kolf2/${filename}" \
    "${INPUT}/kolf2/${filename}" \
    --only-show-errors
done

aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.npz" \
  "${INPUT}/whitener_union.npz" \
  --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.json" \
  "${INPUT}/whitener_union.json" \
  --only-show-errors

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy==2.5.2 \
python experiments/079-kolf-local-centering/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.25 \
  --report "$OUT/test-memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy==2.5.2 \
  python experiments/079-kolf-local-centering/test_run.py \
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
  src/vcc2026/donor_transfer.py \
  src/vcc2026/multi_donor.py \
  src/vcc2026/whiten.py \
  experiments/075-postwhiten-local-centering/run.py \
  experiments/075-postwhiten-local-centering/cloud-output/20260823T062845Z/results.json \
  experiments/079-kolf-local-centering/run.py \
  experiments/079-kolf-local-centering/test_run.py \
  experiments/079-kolf-local-centering/memory_guard.py \
  experiments/079-kolf-local-centering/preregistered_protocol.json \
  experiments/079-kolf-local-centering/dependency-sha256.txt \
  experiments/079-kolf-local-centering/aws/run.sh \
  >"$OUT/code.sha256"

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy==2.5.2 \
python experiments/079-kolf-local-centering/memory_guard.py \
  --max-rss-gib 32 \
  --poll-seconds 1 \
  --report "$OUT/evaluate-memory.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy==2.5.2 \
  python experiments/079-kolf-local-centering/run.py \
    --donor-root "$INPUT/donors" \
    --kolf "$INPUT/kolf2" \
    --whitener "$INPUT/whitener_union.npz" \
    --whitener-metadata "$INPUT/whitener_union.json" \
    --exp075-results \
      experiments/075-postwhiten-local-centering/cloud-output/20260823T062845Z/results.json \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output-dir "$OUT/result" \
    --max-rss-gib 32 \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

cp "$OUT/result/results.json" "$OUT/results.json"
cp "$OUT/result/per_target_pds.csv" "$OUT/per_target_pds.csv"
rm -rf "$OUT/result"

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${INSTANCE_TYPE}",
  "market": "${MARKET:-unknown}",
  "source_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "worker_was_created_or_started": false,
  "fits_eigensystem": false,
  "cell_eval2": false,
  "count_emission": false,
  "submission": false
}
EOF

upload_outputs
