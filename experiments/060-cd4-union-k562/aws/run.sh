#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
PREFIX="exp060-cd4-union-k562"
WORK="/work/exp060"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUTPUT="${REPO}/experiments/060-cd4-union-k562"
LOG="${WORK}/cloud.log"

mkdir -p "${WORK}"
exec > >(tee -a "${LOG}") 2>&1

finish() {
  status=$?
  trap - EXIT
  mkdir -p "${OUTPUT}"
  printf '%s\n' "${status}" > "${OUTPUT}/exit-code.txt"
  aws s3 cp "${OUTPUT}/" \
    "s3://${BUCKET}/output/${PREFIX}/" \
    --recursive --only-show-errors || true
  aws s3 cp "${LOG}" \
    "s3://${BUCKET}/output/${PREFIX}/cloud.log" \
    --only-show-errors || true
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

mkdir -p "${REPO}" "${INPUT}/donors" "${INPUT}/cd4"
aws s3 cp \
  "s3://${BUCKET}/input/${PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
tar -xzf "${WORK}/source.tar.gz" -C "${REPO}"
aws s3 cp \
  "s3://${BUCKET}/input/${PREFIX}/root-README.md" \
  "${REPO}/README.md" --only-show-errors

for donor in orion_hct116 orion_hek293t k562_gwps
do
  mkdir -p "${INPUT}/donors/${donor}"
  aws s3 cp \
    "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "${INPUT}/donors/${donor}/" \
    --recursive --only-show-errors
done

aws s3 cp \
  "s3://${BUCKET}/input/${PREFIX}/cd4_compact.tar.gz" \
  "${INPUT}/cd4_compact.tar.gz" --only-show-errors
tar -xzf "${INPUT}/cd4_compact.tar.gz" -C "${INPUT}/cd4"

cd "${REPO}"
uv run python experiments/060-cd4-union-k562/test_run.py

INSTANCE_ID="$(curl -fsS -X PUT \
  -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600' \
  http://169.254.169.254/latest/api/token | \
  xargs -I{} curl -fsS \
  -H 'X-aws-ec2-metadata-token: {}' \
  http://169.254.169.254/latest/meta-data/instance-id)"

/usr/bin/time -v uv run python experiments/060-cd4-union-k562/run.py \
  --hct116 "${INPUT}/donors/orion_hct116" \
  --hek293t "${INPUT}/donors/orion_hek293t" \
  --k562 "${INPUT}/donors/k562_gwps" \
  --cd4-root "${INPUT}/cd4/cd4_compact" \
  --panel-csv data/challenge/pert_counts.csv \
  --genes-csv data/challenge/gene_names.csv \
  --output "${OUTPUT}/results.json" \
  --whitener "${OUTPUT}/whitener_union.npz" \
  --scratch "${WORK}/scratch" \
  --max-rss-gib 55 \
  --instance-id "${INSTANCE_ID}" \
  --instance-type "r7i.2xlarge" \
  --s3-prefix "s3://${BUCKET}/output/${PREFIX}/"
