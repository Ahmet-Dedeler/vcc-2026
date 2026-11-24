#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
WORK="/work"
LOG="/var/log/vcc-cloud-score.log"

mkdir -p "$WORK"
exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  aws s3 cp "$LOG" "s3://${BUCKET}/output/cloud.log" || true
  printf '%s\n' "$status" > "$WORK/exit-code.txt"
  aws s3 cp "$WORK/exit-code.txt" "s3://${BUCKET}/output/exit-code.txt" || true
  shutdown -h now
  exit "$status"
}
trap finish EXIT

dnf install -y gzip tar time
curl -LsSf https://astral.sh/uv/install.sh | sh
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="$WORK/uv-cache"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

wait_for_object() {
  key="$1"
  until aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null 2>&1
  do
    echo "waiting for s3://${BUCKET}/${key}"
    sleep 30
  done
}

wait_for_object "input/orion-hct116.h5ad"
aws s3 cp "s3://${BUCKET}/input/orion-hct116.h5ad" "$WORK/real.h5ad"

for tag in \
  g0.45_q1.0_a1.08_spread0.25 \
  g0.55_q1.0_a1.08_spread0.25 \
  g0.75_q1.0_a1.08_spread0.25 \
  g0.85_q1.0_a1.08_spread0.25
do
  mkdir -p "$WORK/$tag"
  wait_for_object "input/preds/${tag}/pred.h5ad"
  aws s3 cp "s3://${BUCKET}/input/preds/${tag}/pred.h5ad" "$WORK/$tag/pred.h5ad"
  /usr/bin/time -v uvx --python 3.12 --from cell-eval2==0.16.0 \
    cell-eval2 run \
    -ap "$WORK/$tag/pred.h5ad" \
    -ar "$WORK/real.h5ad" \
    -o "$WORK/$tag/run" \
    --preset vcc2026 \
    --pert-col target_gene \
    --control non-targeting \
    2>&1 | tee "$WORK/$tag/score.log"
  aws s3 cp "$WORK/$tag/run" "s3://${BUCKET}/output/${tag}/run/" --recursive
  aws s3 cp "$WORK/$tag/score.log" "s3://${BUCKET}/output/${tag}/score.log"
  rm -rf "$WORK/$tag"
done
