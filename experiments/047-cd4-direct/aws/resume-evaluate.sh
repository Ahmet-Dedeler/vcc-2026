#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
WORK="/work/cd4"
REPO="$WORK/repo"
COMPACT="$WORK/cd4_compact"
LOG="/var/log/vcc-cd4-resume.log"

exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  aws s3 cp "$LOG" \
    "s3://${BUCKET}/output-cd4/resume.log" --only-show-errors || true
  printf '%s\n' "$status" > "$WORK/resume-exit-code.txt"
  aws s3 cp "$WORK/resume-exit-code.txt" \
    "s3://${BUCKET}/output-cd4/exit-code.txt" --only-show-errors || true
  shutdown -h now
  exit "$status"
}
trap finish EXIT

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="/work/uv-cache"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

test -s "$COMPACT/manifest.json"
for donor in k562_gwps orion_hct116 orion_hek293t
do
  test -s "$WORK/donors/$donor/delta.npy"
done

aws s3 cp "s3://${BUCKET}/input/cd4/evaluate_lodo.py" \
  "$REPO/experiments/047-cd4-direct/evaluate_lodo.py" --only-show-errors

cd "$REPO"
/usr/bin/time -v uv run python experiments/047-cd4-direct/evaluate_lodo.py \
  --cd4-root "$COMPACT" \
  --k562 "$WORK/donors/k562_gwps" \
  --hct116 "$WORK/donors/orion_hct116" \
  --hek293t "$WORK/donors/orion_hek293t" \
  --max-rss-gib 4 \
  --output "$WORK/results.json"

tar -czf "$WORK/cd4_compact.tar.gz" -C "$WORK" cd4_compact
aws s3 cp "$WORK/cd4_compact.tar.gz" \
  "s3://${BUCKET}/output-cd4/cd4_compact.tar.gz" --only-show-errors
aws s3 cp "$WORK/results.json" \
  "s3://${BUCKET}/output-cd4/results.json" --only-show-errors
aws s3 cp "$WORK/results.csv" \
  "s3://${BUCKET}/output-cd4/results.csv" --only-show-errors
aws s3 cp "$COMPACT/manifest.json" \
  "s3://${BUCKET}/output-cd4/manifest.json" --only-show-errors
