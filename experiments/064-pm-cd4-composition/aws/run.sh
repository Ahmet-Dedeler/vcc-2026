#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
PREFIX="exp064-pm-cd4-composition"
WORK="/work/exp064"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUTPUT="${WORK}/output"
LOG="${WORK}/cloud.log"

mkdir -p "${WORK}" "${OUTPUT}"
exec > >(tee -a "${LOG}") 2>&1

finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "${status}" > "${OUTPUT}/exit-code.txt"
  cp "${LOG}" "${OUTPUT}/cloud.log" || true
  (
    cd "${OUTPUT}"
    find . -type f ! -name manifest.sha256 -print0 |
      sort -z |
      xargs -0 sha256sum > manifest.sha256
  ) || true
  aws s3 cp "${OUTPUT}/" \
    "s3://${BUCKET}/output/${PREFIX}/" \
    --recursive --only-show-errors || true
  exit "${status}"
}
trap finish EXIT

export PATH="/root/.local/bin:${PATH}"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

dnf install -y gzip tar time
if ! command -v uv >/dev/null 2>&1
then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

mkdir -p \
  "${REPO}" \
  "${INPUT}/donors" \
  "${INPUT}/perturbmatch" \
  "${INPUT}/cd4" \
  "${INPUT}/whiteners"
aws s3 cp \
  "s3://${BUCKET}/input/${PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
tar -xzf "${WORK}/source.tar.gz" -C "${REPO}"

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
for held in orion_hct116 orion_hek293t
do
  aws s3 cp \
    "s3://${BUCKET}/input/${PREFIX}/whiteners/held_${held}.npz" \
    "${INPUT}/whiteners/held_${held}.npz" --only-show-errors
  aws s3 cp \
    "s3://${BUCKET}/input/${PREFIX}/whiteners/held_${held}.json" \
    "${INPUT}/whiteners/held_${held}.json" --only-show-errors
done

cd "${REPO}"
uv run --python 3.12 --with numpy python \
  experiments/064-pm-cd4-composition/test_run.py

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
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/test_run.py \
  experiments/064-pm-cd4-composition/memory_guard.py \
  /root/exp064-run.sh \
  > "${OUTPUT}/code.sha256"

uv run --python 3.12 --with numpy python \
  experiments/064-pm-cd4-composition/memory_guard.py \
  --max-rss-gib 55 \
  --poll-seconds 2 \
  --report "${OUTPUT}/memory-guard.json" \
  -- \
  /usr/bin/time -v \
  uv run --python 3.12 --with numpy python \
    experiments/064-pm-cd4-composition/run.py \
    --k562 "${INPUT}/donors/k562_gwps" \
    --hct116 "${INPUT}/donors/orion_hct116" \
    --hek293t "${INPUT}/donors/orion_hek293t" \
    --perturbmatch "${INPUT}/perturbmatch" \
    --cd4-root "${INPUT}/cd4/cd4_compact" \
    --panel-csv data/challenge/pert_counts.csv \
    --genes-csv data/challenge/gene_names.csv \
    --whitener-dir "${INPUT}/whiteners" \
    --output "${OUTPUT}/results.json" \
    --max-rss-gib 55 \
    --instance-id "${INSTANCE_ID}" \
    --instance-type "r7i.2xlarge" \
    --market "${MARKET:-unknown}" \
    --s3-prefix "s3://${BUCKET}/output/${PREFIX}/" \
    --guarded
