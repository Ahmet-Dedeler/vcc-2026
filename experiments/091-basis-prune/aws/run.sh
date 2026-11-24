#!/usr/bin/env bash
# Experiment 088 worker runner. Parameterised by env so one script serves all
# three tracks and both harnesses.
#
#   HARNESS   h1-test | orion-hct116      (which held-out truth to score against)
#   GRID      grid-ab.json | grid-c.json | grid-smoke.json
#   DONORS    space-separated name=path pairs, relative to $INPUT/donors
#   JOBS      concurrent cell-eval2 workers
#
# Adapted from experiments/077-low-alpha-union/aws/run.sh, which is the
# provenance for the staging layout and the safety guards.
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
HARNESS="${HARNESS:?HARNESS is required}"
GRID="${GRID:?GRID is required}"
DONOR_NAMES="${DONOR_NAMES:?DONOR_NAMES is required}"
JOBS="${JOBS:-3}"
BASIS_PRUNE="${BASIS_PRUNE:-none}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp091-basis-prune/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp091-basis-prune/${RUN_ID}}"
WORK="${WORK:-/work/exp091-${RUN_ID}}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

# Pinned truth hashes. h1-test is uploaded by launch.sh; orion-hct116 is the
# same object experiment 077 verified.
H1_TEST_SHA256="06c6baea8c62ef3d38cacdbf58e104908eb93f0ef3b25aca56b86d28240f9692"
ORION_HCT116_SHA256="d20a1e9da59debca8db7d83b955d12c3e1384b530ae6154191d95039b1655863"

mkdir -p "$WORK"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi
if pgrep -af '[c]ell-eval2|[e]xperiments/031-somatic/run.py' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap another heavy VCC process" >&2
  exit 72
fi

rm -rf "$REPO" "$INPUT" "$OUT" "$LOG" "$EXIT_FILE"
mkdir -p "$REPO" "$INPUT/donors" "$OUT/runs"
exec > >(tee -a "$LOG") 2>&1

upload_outputs() {
  aws s3 cp "$LOG" "s3://${BUCKET}/${OUTPUT_PREFIX}/cloud.log" --only-show-errors || true
  if [[ -s "$EXIT_FILE" ]]; then
    aws s3 cp "$EXIT_FILE" "s3://${BUCKET}/${OUTPUT_PREFIX}/exit-code.txt" \
      --only-show-errors || true
  fi
  if [[ -d "$OUT" ]]; then
    # Predictions are 2.3 GB each and reproducible from the pinned grid; the
    # fitted whitener is 1.3 GB and re-derivable. Neither goes to S3.
    aws s3 sync "$OUT/" "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
      --exclude '*.h5ad' --exclude 'whitener_*.npz' --only-show-errors || true
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
# Guard at ~85% of physical RAM. Experiment 071 died three times against a
# hardcoded ceiling that did not track the instance it was running on.
guard_gib="$(( mem_kib / 1024 / 1024 * 85 / 100 ))"
avail_kib="$(df -Pk /work | awk 'NR==2 {print $4}')"
if (( avail_kib < 60 * 1024 * 1024 )); then
  echo "worker has less than 60 GiB free disk" >&2
  exit 74
fi
echo "worker: $((mem_kib/1024/1024)) GiB RAM, guard ${guard_gib} GiB, jobs ${JOBS}"

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

dnf install -y gzip tar time >/dev/null
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" "${WORK}/source.tar.gz" \
  --only-show-errors
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"

case "$HARNESS" in
  h1-test)
    REAL_SHA256="$H1_TEST_SHA256"
    aws s3 cp "s3://${BUCKET}/input/harness/h1-test.h5ad" \
      "$INPUT/${HARNESS}.h5ad" --only-show-errors ;;
  orion-hct116)
    REAL_SHA256="$ORION_HCT116_SHA256"
    aws s3 cp "s3://${BUCKET}/input/harness/orion-hct116.h5ad" \
      "$INPUT/${HARNESS}.h5ad" --only-show-errors ;;
  *) echo "unknown HARNESS ${HARNESS}" >&2; exit 75 ;;
esac
printf '%s  %s\n' "$REAL_SHA256" "$INPUT/${HARNESS}.h5ad" | sha256sum --check -

donor_args=()
for donor in $DONOR_NAMES; do
  mkdir -p "$INPUT/donors/$donor"
  aws s3 cp "s3://${BUCKET}/input/clean-hct/${donor}/" "$INPUT/donors/$donor/" \
    --recursive --only-show-errors
  donor_args+=("${donor}=$INPUT/donors/$donor")
done

# A whitener fitted for this exact (axis, donor set, excluded panel) may already
# exist from an earlier wave. run.py validates the cache against its own
# metadata and refits on any mismatch, so a stale object cannot be used silently.
aws s3 cp "s3://${BUCKET}/input/exp091-whitener/${HARNESS}/whitener_union.npz" \
  "$OUT/whitener_union.npz" --only-show-errors 2>/dev/null || true
aws s3 cp "s3://${BUCKET}/input/exp091-whitener/${HARNESS}/whitener_union.json" \
  "$OUT/whitener_union.json" --only-show-errors 2>/dev/null || true

cd "$REPO"
uv sync --frozen

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib "$guard_gib" \
  --min-free-percent 0 \
  --poll-seconds 5 \
  -- \
  /usr/bin/time -v uv run python experiments/031-somatic/run.py \
    --real "$INPUT/${HARNESS}.h5ad" \
    --donors "${donor_args[@]}" \
    --configs "experiments/091-basis-prune/${GRID}" \
    --outdir "$OUT" \
    --whiten-axis union \
    --exclude-panel challenge \
    --basis-prune "$BASIS_PRUNE" \
    --jobs "$JOBS" \
  2>&1 | tee "$OUT/run.log"

# The worker role can write output/ but not input/, so publishing the fitted
# whitener from here is not possible. It refits in about 4 minutes per worker,
# which is cheaper than plumbing a wider IAM policy for it.

(
  cd "$OUT"
  find . -type f ! -name '*.h5ad' ! -name 'whitener_*.npz' \
    ! -name artifact-sha256.txt -print0 |
    sort -z | xargs -0 sha256sum >artifact-sha256.txt
)
echo "done: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
