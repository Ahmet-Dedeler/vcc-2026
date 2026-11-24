#!/bin/bash
set -Eeuo pipefail

BUCKET="<your-s3-bucket>"
WORK="/work/cd4"
REPO="$WORK/repo"
SOURCE="$WORK/GWCD4i.DE_stats.h5ad"
COMPACT="$WORK/cd4_compact"
LOG="/var/log/vcc-cd4-direct.log"

mkdir -p "$WORK"
exec > >(tee -a "$LOG") 2>&1

finish() {
  status=$?
  trap - EXIT
  aws s3 cp "$LOG" "s3://${BUCKET}/output-cd4/cloud.log" || true
  printf '%s\n' "$status" > "$WORK/exit-code.txt"
  aws s3 cp "$WORK/exit-code.txt" \
    "s3://${BUCKET}/output-cd4/exit-code.txt" || true
  shutdown -h now
  exit "$status"
}
trap finish EXIT

dnf install -y gzip tar time
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="/work/uv-cache"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

wait_for_object() {
  key="$1"
  until aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null 2>&1
  do
    echo "waiting for s3://${BUCKET}/${key}"
    sleep 30
  done
}

# The previous HCT sweep can leave several large prediction files. Its compact
# scores have already been uploaded before this script is started.
rm -rf /work/repo /work/source.tar.gz /work/score.log /work/exit-code.txt
rm -rf /work/g0.*
rm -rf "$WORK"
mkdir -p "$REPO" "$WORK/donors"

wait_for_object "input/cd4/source.tar.gz"
aws s3 cp "s3://${BUCKET}/input/cd4/source.tar.gz" "$WORK/source.tar.gz"
tar -xzf "$WORK/source.tar.gz" -C "$REPO"
# The compact archive intentionally carries only task-relevant files, while
# pyproject.toml references a root README as package metadata.
touch "$REPO/README.md"

aws s3 cp --no-sign-request \
  "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad" \
  "$SOURCE"

cd "$REPO"
/usr/bin/time -v uv run python experiments/047-cd4-direct/extract_cd4.py \
  --input "$SOURCE" \
  --output "$COMPACT" \
  --max-rss-gib 4

# The 16.8 GB source is deliberately ephemeral. Keep only compact panel slices.
rm -f "$SOURCE"

for donor in k562_gwps orion_hct116 orion_hek293t
do
  wait_for_object "input/clean-hct/${donor}/delta.npy"
  mkdir -p "$WORK/donors/$donor"
  aws s3 cp "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "$WORK/donors/$donor/" --recursive --only-show-errors
done

/usr/bin/time -v uv run python experiments/047-cd4-direct/evaluate_lodo.py \
  --cd4-root "$COMPACT" \
  --k562 "$WORK/donors/k562_gwps" \
  --hct116 "$WORK/donors/orion_hct116" \
  --hek293t "$WORK/donors/orion_hek293t" \
  --max-rss-gib 4 \
  --output "$WORK/results.json"

tar -czf "$WORK/cd4_compact.tar.gz" -C "$WORK" cd4_compact
aws s3 cp "$WORK/cd4_compact.tar.gz" \
  "s3://${BUCKET}/output-cd4/cd4_compact.tar.gz"
aws s3 cp "$WORK/results.json" \
  "s3://${BUCKET}/output-cd4/results.json"
aws s3 cp "$WORK/results.csv" \
  "s3://${BUCKET}/output-cd4/results.csv"
aws s3 cp "$COMPACT/manifest.json" \
  "s3://${BUCKET}/output-cd4/manifest.json"
