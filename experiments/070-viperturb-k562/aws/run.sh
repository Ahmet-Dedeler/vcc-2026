#!/usr/bin/env bash
set -Eeuo pipefail

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
SOURCE_PREFIX="${SOURCE_PREFIX:-input/exp070-viperturb-k562/${RUN_ID}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-output/exp070-viperturb-k562/${RUN_ID}}"
WORK="${WORK:-/work/exp070-${RUN_ID}}"
REPO="${WORK}/repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"
MIN_FREE_GIB="${MIN_FREE_GIB:-24}"
MAX_RSS_GIB="${MAX_RSS_GIB:-55}"
NULL_REPS="${NULL_REPS:-1000}"
RDS_URL="https://zenodo.org/api/records/18460279/files/genome_wide_filtered.rds/content"
RDS_BYTES="3611371127"
RDS_MD5="ac52c41233fc7a5a5c3810e0a075007f"
MANIFEST_URL="https://zenodo.org/api/records/18460279/files/genome_wide_manifest.txt/content"
MANIFEST_MD5="c095c8f2d16136753c22ef5d9b30a22b"

mkdir -p "$WORK" "$OUT"
rm -f "$EXIT_FILE"
exec 9>/work/vcc-heavy.lock
if ! flock -n 9; then
  echo "another heavy VCC process holds /work/vcc-heavy.lock" >&2
  exit 71
fi

if pgrep -af \
  '[c]ell-eval2|[b]uild_candidate.py|[b]uild_joint_candidate.py|[e]xtract_viperturb.R|[e]xperiments/070-viperturb-k562/run.py' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  echo "refusing to overlap an existing VCC build, scorer, or exp070 run" >&2
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
  printf '%s\n' "$status" > "$EXIT_FILE"
  (
    cd "$OUT"
    find . -type f ! -name artifact-sha256.txt -print0 |
      sort -z |
      xargs -0 sha256sum > artifact-sha256.txt
  ) || true
  upload_outputs
  exit "$status"
}
trap finish EXIT

require_free_disk() {
  path="$1"
  minimum_gib="$2"
  available_kib="$(df -Pk "$path" | awk 'NR==2 {print $4}')"
  required_kib="$((minimum_gib * 1024 * 1024))"
  if (( available_kib < required_kib )); then
    echo "insufficient disk at ${path}: $((available_kib / 1024 / 1024)) GiB free, ${minimum_gib} GiB required" >&2
    exit 73
  fi
  echo "disk preflight: $((available_kib / 1024 / 1024)) GiB free at ${path}"
}

install_host_tools() {
  if command -v dnf >/dev/null 2>&1; then
    # Amazon Linux 2023 ships curl-minimal. Never let bootstrap replace it.
    command -v curl >/dev/null 2>&1 || {
      echo "Amazon Linux worker lacks its preinstalled curl command" >&2
      exit 74
    }
    retained_curl_minimal=false
    if rpm -q curl-minimal >/dev/null 2>&1; then
      retained_curl_minimal=true
    fi
    while IFS=: read -r command package; do
      if [[ "$command" == /* ]]; then
        [[ -x "$command" ]] && continue
      else
        command -v "$command" >/dev/null 2>&1 && continue
      fi
      dnf install -y \
        --setopt=install_weak_deps=False \
        --exclude=curl \
        --exclude=curl-minimal \
        --exclude=libcurl \
        --exclude=libcurl-minimal \
        "$package"
    done <<'EOF'
bzip2:bzip2
gzip:gzip
tar:tar
flock:util-linux
/usr/bin/time:time
EOF
    if [[ "$retained_curl_minimal" == true ]] &&
      ! rpm -q curl-minimal >/dev/null 2>&1
    then
      echo "host bootstrap replaced curl-minimal; refusing to continue" >&2
      exit 74
    fi
  elif command -v apt-get >/dev/null 2>&1; then
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
      bzip2 ca-certificates curl gzip tar time util-linux
  else
    for command in bzip2 curl tar gzip; do
      command -v "$command" >/dev/null 2>&1 || {
        echo "missing required host command: ${command}" >&2
        exit 74
      }
    done
  fi
  command -v curl >/dev/null 2>&1 || {
    echo "missing required host command: curl" >&2
    exit 74
  }
}

install_micromamba() {
  mkdir -p "${WORK}/bin"
  if [[ ! -x "${WORK}/bin/micromamba" ]]; then
    machine="$(uname -m)"
    case "$machine" in
      x86_64) platform="linux-64" ;;
      aarch64|arm64) platform="linux-aarch64" ;;
      *)
        echo "unsupported micromamba architecture: ${machine}" >&2
        exit 75
        ;;
    esac
    curl -fsSL --retry 5 \
      "https://micro.mamba.pm/api/micromamba/${platform}/latest" |
      tar -xj -C "${WORK}/bin" --strip-components=1 bin/micromamba
  fi
}

export PATH="${WORK}/bin:/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export MAMBA_ROOT_PREFIX="${WORK}/micromamba"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

install_host_tools
require_free_disk "$WORK" "$MIN_FREE_GIB"
install_micromamba

rm -rf "$REPO" "$INPUT" "$OUT"
mkdir -p \
  "$REPO" \
  "$INPUT/donors" \
  "$INPUT/perturbmatch" \
  "$INPUT/cd4" \
  "$INPUT/whiteners" \
  "$INPUT/viperturb-source" \
  "$OUT"

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "${WORK}/source.tar.gz" --only-show-errors
if tar -tzf "${WORK}/source.tar.gz" |
  grep -Eq '(^|/)\._|(^|/)\.DS_Store$'
then
  echo "source archive contains forbidden macOS metadata files" >&2
  exit 77
fi
tar -xzf "${WORK}/source.tar.gz" -C "$REPO"

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
aws s3 cp \
  "s3://${BUCKET}/${SOURCE_PREFIX}/whiteners/" \
  "${INPUT}/whiteners/" \
  --recursive --only-show-errors

cd "$REPO"
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="/root/.local/bin:$PATH"
fi
uv python install 3.12
uv sync --frozen

PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib 1.0 \
  --poll-seconds 0.25 \
  --report "$OUT/test-memory-guard.json" \
  -- \
  uv run python experiments/070-viperturb-k562/test_run.py \
  2>&1 | tee "$OUT/tests.log"

require_free_disk "$WORK" "$MIN_FREE_GIB"
curl -fsSL --retry 8 --retry-all-errors --continue-at - \
  "$MANIFEST_URL" \
  --output "${INPUT}/viperturb-source/genome_wide_manifest.txt"
printf '%s  %s\n' \
  "$MANIFEST_MD5" \
  "${INPUT}/viperturb-source/genome_wide_manifest.txt" |
  md5sum --check -

uv run python experiments/070-viperturb-k562/audit_manifest.py \
  --manifest "${INPUT}/viperturb-source/genome_wide_manifest.txt" \
  --panel-csv data/challenge/pert_counts.csv \
  --expected-filtered-coverage 124 \
  --low-cell-threshold 30 \
  --k562-perts "${INPUT}/donors/k562_gwps/perts.txt" \
  --perturbmatch-perts "${INPUT}/perturbmatch/perts.txt" \
  --hct116-perts "${INPUT}/donors/orion_hct116/perts.txt" \
  --hek293t-perts "${INPUT}/donors/orion_hek293t/perts.txt" \
  --minimum-modified-targets 70 \
  --output "$OUT/manifest-audit.json"

R_ENV="${WORK}/r-env"
if [[ ! -x "${R_ENV}/bin/Rscript" ]]; then
  micromamba create -y --strict-channel-priority \
    -p "$R_ENV" \
    -f experiments/070-viperturb-k562/environment.yml
fi
micromamba run -p "$R_ENV" Rscript -e \
  'stopifnot(as.character(packageVersion("Seurat")) == "5.3.1"); sessionInfo()' \
  > "$OUT/r-environment.txt"
micromamba list -p "$R_ENV" --json > "$OUT/r-environment-packages.json"
uv run python experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib 1.0 \
  --poll-seconds 0.25 \
  --report "$OUT/r-schema-test-memory-guard.json" \
  -- \
  micromamba run -p "$R_ENV" Rscript \
    experiments/070-viperturb-k562/test_extract_schema.R \
  2>&1 | tee "$OUT/r-schema-tests.log"

require_free_disk "$WORK" "$MIN_FREE_GIB"
curl -fL --retry 8 --retry-all-errors --continue-at - \
  "$RDS_URL" \
  --output "${INPUT}/viperturb-source/genome_wide_filtered.rds.part"
actual_bytes="$(stat -c '%s' "${INPUT}/viperturb-source/genome_wide_filtered.rds.part")"
if [[ "$actual_bytes" != "$RDS_BYTES" ]]; then
  echo "RDS byte count ${actual_bytes} != ${RDS_BYTES}" >&2
  exit 76
fi
printf '%s  %s\n' \
  "$RDS_MD5" \
  "${INPUT}/viperturb-source/genome_wide_filtered.rds.part" |
  md5sum --check -
mv \
  "${INPUT}/viperturb-source/genome_wide_filtered.rds.part" \
  "${INPUT}/viperturb-source/genome_wide_filtered.rds"

uv run python experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib "$MAX_RSS_GIB" \
  --poll-seconds 2 \
  --report "$OUT/extract-memory-guard.json" \
  -- \
  /usr/bin/time -v micromamba run -p "$R_ENV" Rscript \
    experiments/070-viperturb-k562/extract_viperturb.R \
    --input "${INPUT}/viperturb-source/genome_wide_filtered.rds" \
    --manifest "${INPUT}/viperturb-source/genome_wide_manifest.txt" \
    --genes-csv data/challenge/gene_names.csv \
    --panel-csv data/challenge/pert_counts.csv \
    --output "$OUT/viperturb_compact" \
    --verified-source-md5 "$RDS_MD5" \
    --verified-manifest-md5 "$MANIFEST_MD5" \
    --low-cell-threshold 30 \
    --group-chunk-size 256 \
  2>&1 | tee "$OUT/extract.log"

rm -f "${INPUT}/viperturb-source/genome_wide_filtered.rds"

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
  experiments/070-viperturb-k562/run.py \
  experiments/070-viperturb-k562/test_run.py \
  experiments/070-viperturb-k562/audit_manifest.py \
  experiments/070-viperturb-k562/memory_guard.py \
  experiments/070-viperturb-k562/extract_viperturb.R \
  experiments/070-viperturb-k562/test_extract_schema.R \
  experiments/070-viperturb-k562/environment.yml \
  experiments/070-viperturb-k562/aws/run.sh \
  > "$OUT/code.sha256"

uv run python experiments/070-viperturb-k562/memory_guard.py \
  --max-rss-gib "$MAX_RSS_GIB" \
  --poll-seconds 2 \
  --report "$OUT/evaluate-memory-guard.json" \
  -- \
  /usr/bin/time -v uv run python \
    experiments/070-viperturb-k562/run.py \
    --k562 "${INPUT}/donors/k562_gwps" \
    --hct116 "${INPUT}/donors/orion_hct116" \
    --hek293t "${INPUT}/donors/orion_hek293t" \
    --perturbmatch "${INPUT}/perturbmatch" \
    --viperturb "$OUT/viperturb_compact" \
    --cd4-root "${INPUT}/cd4/cd4_compact" \
    --panel-csv data/challenge/pert_counts.csv \
    --genes-csv data/challenge/gene_names.csv \
    --whitener-dir "${INPUT}/whiteners" \
    --output "$OUT/results.json" \
    --max-rss-gib "$MAX_RSS_GIB" \
    --null-reps "$NULL_REPS" \
    --minimum-modified-targets 70 \
    --minimum-non-low-targets 20 \
    --instance-id "$INSTANCE_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --market "${MARKET:-unknown}" \
    --s3-prefix "s3://${BUCKET}/${OUTPUT_PREFIX}/" \
    --guarded \
  2>&1 | tee "$OUT/evaluate.log"

upload_outputs
