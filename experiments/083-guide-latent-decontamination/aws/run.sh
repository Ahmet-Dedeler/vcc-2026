#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp083-guide-latent-decontamination/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp083-guide-latent-decontamination/${RUN_ID}}"
SOURCE_ARCHIVE_SHA256="${SOURCE_ARCHIVE_SHA256:?SOURCE_ARCHIVE_SHA256 is required}"
WORK_ROOT="${WORK_ROOT:-/work/exp083}"
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 63
fi
WORK="${WORK_ROOT}/${RUN_ID}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

for command in aws awk curl df flock grep pgrep sha256sum tar tee; do
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

mkdir -p "$WORK_ROOT"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi
if pgrep -af \
  '[c]ell-eval2|[e]xperiments/0(62|64|65|69|70|71|73|75|77|79|80|83)-.*\\.py' \
  >"${WORK_ROOT}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC heavy process" >&2
  cat "${WORK_ROOT}/preexisting-heavy-processes.txt" >&2
  exit 72
fi

if [[ "$WORK" != "${WORK_ROOT}/${RUN_ID}" || "$WORK" == "/" ]]; then
  echo "unsafe per-run path: ${WORK}" >&2
  exit 70
fi
rm -rf "$WORK"
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

visible_mem_kib="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
minimum_visible_kib=$((60 * 1024 * 1024))
if [[ -z "$visible_mem_kib" || "$visible_mem_kib" -lt "$minimum_visible_kib" ]]; then
  echo "worker exposes ${visible_mem_kib:-0} KiB RAM; require nominal >=64 GiB" >&2
  exit 74
fi
available_mem_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
minimum_available_kib=$((55 * 1024 * 1024))
if [[ -z "$available_mem_kib" || "$available_mem_kib" -lt "$minimum_available_kib" ]]; then
  echo "worker exposes ${available_mem_kib:-0} KiB available RAM; require >=55 GiB" >&2
  exit 76
fi
free_kib="$(df -Pk "$WORK_ROOT" | awk 'NR==2 {print $4}')"
minimum_free_kib=$((80 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  echo "worker has ${free_kib:-0} KiB free; require at least 80 GiB" >&2
  exit 75
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-16}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-16}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-16}"
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "$WORK/source.tar.gz" --only-show-errors
printf '%s  %s\n' "$SOURCE_ARCHIVE_SHA256" "$WORK/source.tar.gz" |
  sha256sum -c -
if tar -tzf "$WORK/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata" >&2
  exit 73
fi
tar -xzf "$WORK/source.tar.gz" -C "$REPO"
cd "$REPO"
sha256sum -c experiments/083-guide-latent-decontamination/dependency-sha256.txt

mkdir -p \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/operators" \
  "$INPUT/kolf" \
  "$OUT/preflight"
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/ReplogleWeissman2022_K562_gwps.h5ad" \
  "$INPUT/K562_gwps.h5ad" --only-show-errors
curl -fL \
  "https://raw.githubusercontent.com/AustinHartman/perturb_seed_reproduce/198fe3d84263002e9cfa17e1a2178d4e8dd0a669/data/replogle_k562_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz" \
  -o "$INPUT/replogle-seed-candidates.csv.gz"

# The same structural gate is repeated remotely before any expression block is
# read. With the currently audited source this exits 78 and the trap uploads
# the evidence; extraction and fitting remain unreachable.
PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.2 \
  --report "$OUT/preflight/memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/inspect_source.py \
    --h5ad "$INPUT/K562_gwps.h5ad" \
    --panel-csv data/challenge/pert_counts.csv \
    --perturbmatch-perts data/external/perturbmatch/perts.txt \
    --seed-table "$INPUT/replogle-seed-candidates.csv.gz" \
    --hash-h5ad \
    --output "$OUT/preflight/source-audit.json" \
    --require-promotion-coverage

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
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/kolf2/" "$INPUT/kolf/" \
  --recursive --only-show-errors
aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/operators/" "$INPUT/operators/" \
  --recursive --only-show-errors

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 python \
  experiments/083-guide-latent-decontamination/verify_inputs.py \
  --manifest experiments/083-guide-latent-decontamination/artifact_manifest.json \
  --h5ad "$INPUT/K562_gwps.h5ad" \
  --seed-table "$INPUT/replogle-seed-candidates.csv.gz" \
  --genes-csv data/challenge/gene_names.csv \
  --panel-csv data/challenge/pert_counts.csv \
  --donor-root "$INPUT/donors" \
  --perturbmatch "$INPUT/perturbmatch" \
  --kolf "$INPUT/kolf" \
  --cd4-root "$INPUT/cd4/cd4_compact" \
  --operator-dir "$INPUT/operators" \
  2>&1 | tee "$OUT/input-verification.log"

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.2 \
  --report "$OUT/test-memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/test_model.py \
  2>&1 | tee "$OUT/tests.log"

mkdir -p "$OUT/guide-cache" "$OUT/feature-cache" "$OUT/frozen" "$OUT/evaluation"
PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 24 \
  --poll-seconds 1 \
  --report "$OUT/extract-memory.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/extract_guide_effects.py \
    --h5ad "$INPUT/K562_gwps.h5ad" \
    --compact-k562 "$INPUT/donors/k562_gwps" \
    --panel-csv data/challenge/pert_counts.csv \
    --output "$OUT/guide-cache" \
    --guarded \
  2>&1 | tee "$OUT/extract.log"

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/build_seed_features.py \
  --guide-ids "$OUT/guide-cache/guide_ids.txt" \
  --source "$INPUT/replogle-seed-candidates.csv.gz" \
  --output "$OUT/feature-cache" \
  2>&1 | tee "$OUT/features.log"

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 48 \
  --poll-seconds 1 \
  --report "$OUT/fit-memory.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/fit.py \
    --guide-cache "$OUT/guide-cache" \
    --feature-cache "$OUT/feature-cache" \
    --panel-csv data/challenge/pert_counts.csv \
    --perturbmatch-perts "$INPUT/perturbmatch/perts.txt" \
    --source-audit "$OUT/preflight/source-audit.json" \
    --output "$OUT/frozen" \
    --guarded \
  2>&1 | tee "$OUT/fit.log"

if [[ ! -f "$OUT/frozen/panel-evaluation-authorized.json" ]]; then
  echo "fit phase did not authorize held-panel evaluation" >&2
  exit 79
fi

PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 48 \
  --poll-seconds 1 \
  --report "$OUT/evaluate-memory.json" \
  -- \
  /usr/bin/time -v \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/evaluate.py \
    --frozen-dir "$OUT/frozen" \
    --source-audit "$OUT/preflight/source-audit.json" \
    --guide-cache "$OUT/guide-cache" \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --k562 "$INPUT/donors/k562_gwps" \
    --hct116 "$INPUT/donors/orion_hct116" \
    --hek293t "$INPUT/donors/orion_hek293t" \
    --perturbmatch "$INPUT/perturbmatch" \
    --kolf "$INPUT/kolf" \
    --cd4-root "$INPUT/cd4/cd4_compact" \
    --operator-dir "$INPUT/operators" \
    --output "$OUT/evaluation" \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "source_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "minimum_instance_memory_gib": 64,
  "minimum_free_disk_gib": 80,
  "extract_guard_gib": 24,
  "fit_guard_gib": 48,
  "evaluation_guard_gib": 48,
  "experiment_080_artifacts_used": false,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

upload_outputs
