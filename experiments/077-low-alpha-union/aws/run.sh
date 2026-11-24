#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp077-low-alpha-union/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp077-low-alpha-union/${RUN_ID}}"
WORK="${WORK:-/work/exp077-${RUN_ID}}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"
REAL_SHA256="d20a1e9da59debca8db7d83b955d12c3e1384b530ae6154191d95039b1655863"
WHITENER_SHA256="a46264ebb6cb4b877032a6e4c68e76069eb9112c6f1bb7f15e1133fae04a8a77"
WHITENER_JSON_SHA256="cf5c1de31d03875b15f6fd11d3f628aadfccf04e3346493c831c3a6f441b4793"
BASELINE_SHA256="f1011cd12786d2196949daa6ce3bf328d52c0afc9d46966e8d3e882f6200cccf"

mkdir -p "$WORK"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi
if pgrep -af \
  '[c]ell-eval2|[e]xperiments/031-somatic/run.py|[e]xperiments/077-low-alpha-union' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap another heavy VCC process" >&2
  exit 72
fi

rm -rf "$REPO" "$INPUT" "$OUT" "$LOG" "$EXIT_FILE"
mkdir -p "$REPO" "$INPUT/donors" "$OUT/runs"
exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUTPUT_PREFIX}/cloud.log" \
    --only-show-errors || true
  if [[ -s "$EXIT_FILE" ]]; then
    aws s3 cp "$EXIT_FILE" \
      "s3://${BUCKET}/${OUTPUT_PREFIX}/exit-code.txt" --only-show-errors || true
  fi
  if [[ -d "$OUT" ]]; then
    aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
      --exclude '*.h5ad' --exclude 'whitener_union.npz' --only-show-errors || true
  fi
}
finish() {
  status=$?
  trap - EXIT
  printf '%s\n' "$status" >"$EXIT_FILE"
  upload_outputs
  exit "$status"
}
trap finish EXIT

mem_kib="$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
if (( mem_kib < 120000000 )); then
  echo "worker has less than 120000000 KiB RAM" >&2
  exit 73
fi
avail_kib="$(df -Pk /work | awk 'NR==2 {print $4}')"
if (( avail_kib < 45 * 1024 * 1024 )); then
  echo "worker has less than 45 GiB free disk" >&2
  exit 74
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

dnf install -y gzip tar time
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
if tar -tzf "${WORK}/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata files" >&2
  exit 75
fi
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/orion-hct116.h5ad" \
  "$INPUT/orion-hct116.h5ad" --only-show-errors
for donor in k562_gwps orion_hek293t; do
  mkdir -p "$INPUT/donors/$donor"
  aws s3 cp "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "$INPUT/donors/$donor/" --recursive --only-show-errors
done
for filename in whitener_union.npz whitener_union.json; do
  aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/${filename}" \
    "$OUT/${filename}" --only-show-errors
done
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/baseline-agg_results.csv" \
  "$INPUT/baseline-agg_results.csv" --only-show-errors
printf '%s  %s\n' "$REAL_SHA256" "$INPUT/orion-hct116.h5ad" |
  sha256sum --check -
printf '%s  %s\n' "$WHITENER_SHA256" "$OUT/whitener_union.npz" |
  sha256sum --check -
printf '%s  %s\n' "$WHITENER_JSON_SHA256" "$OUT/whitener_union.json" |
  sha256sum --check -
printf '%s  %s\n' "$BASELINE_SHA256" "$INPUT/baseline-agg_results.csv" |
  sha256sum --check -

cd "$REPO"
uv sync --frozen

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 110 \
  --min-free-percent 0 \
  --poll-seconds 2 \
  -- \
  /usr/bin/time -v uv run python experiments/031-somatic/run.py \
    --real "$INPUT/orion-hct116.h5ad" \
    --donors \
      "k562_gwps=$INPUT/donors/k562_gwps" \
      "orion_hek293t=$INPUT/donors/orion_hek293t" \
    --configs experiments/077-low-alpha-union/grid.json \
    --outdir "$OUT" \
    --whiten-axis union \
    --exclude-panel challenge \
    --jobs 2 \
  2>&1 | tee "$OUT/run.log"

uv run python experiments/077-low-alpha-union/summarize.py \
  --baseline "$INPUT/baseline-agg_results.csv" \
  --results "$OUT/results.json" \
  --output "$OUT/summary.json" \
  2>&1 | tee "$OUT/summary.log"

(
  cd "$OUT"
  find . -type f \
    ! -name '*.h5ad' \
    ! -name 'whitener_union.npz' \
    ! -name artifact-sha256.txt \
    -print0 |
    sort -z |
    xargs -0 sha256sum >artifact-sha256.txt
)
upload_outputs
