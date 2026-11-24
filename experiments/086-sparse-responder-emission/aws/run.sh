#!/usr/bin/env bash
set -Eeuo pipefail

fail() {
  echo "$1" >&2
  exit "${2:-64}"
}

require_match() {
  local label="$1"
  local value="$2"
  local pattern="$3"
  [[ "$value" =~ $pattern ]] || fail "unsafe ${label}: ${value}" 67
}

BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:?RUN_ID is required}"
EXP077_RUN_ID="${EXP077_RUN_ID:?EXP077_RUN_ID is required}"
EXPECTED_SOURCE_SHA256="${EXPECTED_SOURCE_SHA256:?EXPECTED_SOURCE_SHA256 is required}"
EXPECTED_RUNNER_SHA256="${EXPECTED_RUNNER_SHA256:?EXPECTED_RUNNER_SHA256 is required}"
EXPECTED_RESERVATION_SHA256="${EXPECTED_RESERVATION_SHA256:?EXPECTED_RESERVATION_SHA256 is required}"
SOURCE_PREFIX_RAW="${SOURCE_PREFIX:-}"
OUTPUT_PREFIX_RAW="${OUTPUT_PREFIX:-}"
WORK_ROOT="${WORK_ROOT:-/work/exp086}"

# Validate every supplied identifier before using it in a path, S3 key,
# process match, or systemd unit.
require_match "run id" "$RUN_ID" \
  '^[0-9]{8}T[0-9]{6}Z(-[a-z0-9][a-z0-9-]{0,31})?$'
require_match "exp077 run id" "$EXP077_RUN_ID" \
  '^[0-9]{8}T[0-9]{6}Z(-[a-z0-9][a-z0-9-]{0,31})?$'
require_match "bucket" "$BUCKET" '^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$'
require_match "source SHA-256" "$EXPECTED_SOURCE_SHA256" '^[0-9a-f]{64}$'
require_match "runner SHA-256" "$EXPECTED_RUNNER_SHA256" '^[0-9a-f]{64}$'
require_match "reservation SHA-256" "$EXPECTED_RESERVATION_SHA256" \
  '^[0-9a-f]{64}$'
[[ "$BUCKET" != *..* && "$BUCKET" != *.-* && "$BUCKET" != *-.* ]] ||
  fail "unsafe bucket punctuation" 67
[[ "$BUCKET" == "<your-s3-bucket>" ]] ||
  fail "bucket differs from frozen provenance" 67
[[ "$EXP077_RUN_ID" == "20260823T060938Z" ]] ||
  fail "exp077 run differs from frozen provenance" 67

SOURCE_PREFIX="input/exp086-sparse-responder-emission/${RUN_ID}"
OUTPUT_PREFIX="output/exp086-sparse-responder-emission/${RUN_ID}"
if [[ -n "$SOURCE_PREFIX_RAW" && "$SOURCE_PREFIX_RAW" != "$SOURCE_PREFIX" ]]; then
  fail "source prefix is not canonical for RUN_ID" 67
fi
if [[ -n "$OUTPUT_PREFIX_RAW" && "$OUTPUT_PREFIX_RAW" != "$OUTPUT_PREFIX" ]]; then
  fail "output prefix is not canonical for RUN_ID" 67
fi
EXP077_INPUT_PREFIX="input/exp077-low-alpha-union/${EXP077_RUN_ID}"
EXP077_OUTPUT_PREFIX="output/exp077-low-alpha-union/${EXP077_RUN_ID}"

python3 - "$WORK_ROOT" <<'PY'
import os
import sys
from pathlib import Path

raw = sys.argv[1]
path = Path(raw)
if not path.is_absolute() or os.path.normpath(raw) != raw:
    raise SystemExit(f"WORK_ROOT is not an absolute normalized path: {raw!r}")
probe = Path("/")
for part in path.parts[1:]:
    probe /= part
    if probe.is_symlink():
        raise SystemExit(f"WORK_ROOT contains a symlink component: {probe}")
if str(path.resolve(strict=False)) != raw:
    raise SystemExit(f"WORK_ROOT is not canonical: {raw!r}")
PY

WORK="${WORK_ROOT}/${RUN_ID}"
LOCK_FILE="${WORK_ROOT%/*}/vcc-heavy.lock"
REPO="${WORK}/repo"
EXP077_REPO="${WORK}/exp077-repo"
INPUT="${WORK}/input"
OUT="${WORK}/output"
LOG="${WORK}/cloud.log"
EXIT_FILE="${WORK}/exit-code.txt"

if [[ -e "$WORK" || -L "$WORK" ]]; then
  fail "refusing existing or symlink work directory: ${WORK}" 79
fi

for command in aws awk curl df flock free grep pgrep python3 sha256sum systemctl tar tee; do
  if ! command -v "$command" >/dev/null 2>&1; then
    fail "required command is missing: ${command}" 69
  fi
done

actual_runner_sha256="$(sha256sum "$0" | awk '{print $1}')"
if [[ "$actual_runner_sha256" != "$EXPECTED_RUNNER_SHA256" ]]; then
  fail "runner SHA-256 ${actual_runner_sha256} != ${EXPECTED_RUNNER_SHA256}" 76
fi

mkdir -p "$WORK_ROOT"
[[ ! -L "$WORK_ROOT" ]] || fail "WORK_ROOT is a symlink" 79
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
  fail "another heavy VCC process holds ${LOCK_FILE}" 71
fi
if [[ -e "$WORK" || -L "$WORK" ]]; then
  fail "refusing existing or symlink work directory: ${WORK}" 79
fi
mkdir -m 700 "$WORK"
[[ ! -L "$WORK" && "$(realpath "$WORK")" == "$WORK" ]] ||
  fail "created work directory is noncanonical or a symlink" 79
mkdir -p "$REPO" "$EXP077_REPO" "$INPUT/exp077-output" "$INPUT/donors" "$OUT"

if pgrep -af \
  '[c]ell-eval|[p]ython(3)? .*experiments/|[u]v (run|tool|x) |[t]orchrun|[a]ccelerate launch|[r]ay::|[m]ake_submission|[b]uild_submission|[r]un_with_memory_guard' \
  >"${WORK}/preexisting-heavy-processes.txt"
then
  fail "refusing to overlap an existing VCC-heavy process" 72
fi

active_heavy_units="$(
  systemctl list-units --type=service --state=running --no-legend |
    awk '{print $1}' |
    grep -E '^(vcc-|cell-eval)' |
    grep -v "^vcc-exp086-${RUN_ID}\\.service$" || true
)"
if [[ -n "$active_heavy_units" ]]; then
  printf '%s\n' "$active_heavy_units" >&2
  fail "refusing to overlap an active VCC/system scorer service" 77
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
      --exclude '*.h5ad' --only-show-errors || true
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

# Both prefixes were conditionally reserved by launch.sh. Verify ownership
# before downloading source or writing any result.
for prefix in "$SOURCE_PREFIX" "$OUTPUT_PREFIX"; do
  marker="$INPUT/$(basename "$prefix")-reservation.json"
  aws s3 cp "s3://${BUCKET}/${prefix}/launch-reservation.json" \
    "$marker" --only-show-errors
  printf '%s  %s\n' "$EXPECTED_RESERVATION_SHA256" "$marker" |
    sha256sum --check -
done

visible_mem_kib="$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)"
available_mem_kib="$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)"
minimum_visible_kib=$((120 * 1024 * 1024))
minimum_available_kib=$((100 * 1024 * 1024))
if [[ -z "$visible_mem_kib" || "$visible_mem_kib" -lt "$minimum_visible_kib" ]]; then
  fail "worker exposes ${visible_mem_kib:-0} KiB RAM; require >=120 GiB" 74
fi
if [[ -z "$available_mem_kib" || "$available_mem_kib" -lt "$minimum_available_kib" ]]; then
  fail "worker has ${available_mem_kib:-0} KiB available; require >=100 GiB" 78
fi
free_kib="$(df -Pk "$WORK_ROOT" | awk 'NR==2 {print $4}')"
minimum_free_kib=$((45 * 1024 * 1024))
if [[ -z "$free_kib" || "$free_kib" -lt "$minimum_free_kib" ]]; then
  fail "worker has ${free_kib:-0} KiB free; require >=45 GiB" 75
fi

export PATH="/root/.local/bin:$PATH"
export UV_CACHE_DIR="${WORK}/uv-cache"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-8}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-8}"

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi
uv python install 3.12

aws s3 cp "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  "$WORK/source.tar.gz" --only-show-errors
printf '%s  %s\n' "$EXPECTED_SOURCE_SHA256" "$WORK/source.tar.gz" |
  sha256sum --check -
if tar -tzf "$WORK/source.tar.gz" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
  fail "source archive contains forbidden macOS metadata" 73
fi
tar -xzf "$WORK/source.tar.gz" -C "$REPO"
cd "$REPO"
sha256sum -c experiments/086-sparse-responder-emission/dependency-sha256.txt

# Authenticate exp077 using AWS-authoritative launch and terminal SSM records.
# The current pinned policy intentionally fails closed until a dedicated
# terminal attestation command id has been recorded after exp077 completes.
trusted="experiments/086-sparse-responder-emission/exp077-trusted-provenance.json"
terminal_command_id="$(
  python3 - "$trusted" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1]))
terminal = record.get("terminal_attestation")
print("" if not isinstance(terminal, dict) else terminal.get("ssm_command_id", ""))
PY
)"
if [[ -z "$terminal_command_id" ]]; then
  python3 experiments/086-sparse-responder-emission/gate.py \
    --exp077-output "$INPUT/exp077-output" \
    --run-id "$EXP077_RUN_ID" \
    --trusted-provenance "$trusted" \
    --launch-metadata-snapshot \
      experiments/086-sparse-responder-emission/exp077-launch-metadata.snapshot.json \
    --output "$OUT/exp077-gate.json"
  fail "exp077 provenance unexpectedly passed without terminal evidence" 80
fi
require_match "terminal SSM command id" "$terminal_command_id" \
  '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
launch_command_id="$(
  python3 - "$trusted" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1]))["launch"]["ssm_command_id"])
PY
)"
require_match "launch SSM command id" "$launch_command_id" \
  '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'

for filename in exit-code.txt summary.json results.json artifact-sha256.txt; do
  aws s3 cp \
    "s3://${BUCKET}/${EXP077_OUTPUT_PREFIX}/${filename}" \
    "$INPUT/exp077-output/${filename}" \
    --only-show-errors
done
aws ssm list-commands \
  --region us-east-1 \
  --command-id "$launch_command_id" \
  --output json >"$INPUT/exp077-launch-command.json"
aws ssm list-commands \
  --region us-east-1 \
  --command-id "$terminal_command_id" \
  --output json >"$INPUT/exp077-terminal-command.json"
aws ssm get-command-invocation \
  --region us-east-1 \
  --command-id "$terminal_command_id" \
  --instance-id i-0e4070204e83188fb \
  --output json >"$INPUT/exp077-terminal-invocation.json"
python3 experiments/086-sparse-responder-emission/gate.py \
  --exp077-output "$INPUT/exp077-output" \
  --run-id "$EXP077_RUN_ID" \
  --trusted-provenance "$trusted" \
  --launch-metadata-snapshot \
    experiments/086-sparse-responder-emission/exp077-launch-metadata.snapshot.json \
  --launch-command-evidence "$INPUT/exp077-launch-command.json" \
  --terminal-command-evidence "$INPUT/exp077-terminal-command.json" \
  --terminal-invocation-evidence "$INPUT/exp077-terminal-invocation.json" \
  --output "$OUT/exp077-gate.json"

aws s3 cp "s3://${BUCKET}/${EXP077_INPUT_PREFIX}/source.tar.gz" \
  "$INPUT/exp077-source.tar.gz" --only-show-errors
printf '%s  %s\n' \
  "e1ec0cdb2b5647224431ce334931d9210d1d038494a92660cbff63bcecee3292" \
  "$INPUT/exp077-source.tar.gz" |
  sha256sum --check -
tar -xzf "$INPUT/exp077-source.tar.gz" -C "$EXP077_REPO"

aws s3 cp "s3://${BUCKET}/${EXP077_INPUT_PREFIX}/orion-hct116.h5ad" \
  "$INPUT/orion-hct116.h5ad" --only-show-errors
aws s3 cp "s3://${BUCKET}/${EXP077_INPUT_PREFIX}/whitener_union.npz" \
  "$INPUT/whitener_union.npz" --only-show-errors
aws s3 cp "s3://${BUCKET}/${EXP077_INPUT_PREFIX}/whitener_union.json" \
  "$INPUT/whitener_union.json" --only-show-errors
printf '%s  %s\n' \
  "d20a1e9da59debca8db7d83b955d12c3e1384b530ae6154191d95039b1655863" \
  "$INPUT/orion-hct116.h5ad" |
  sha256sum --check -
printf '%s  %s\n' \
  "a46264ebb6cb4b877032a6e4c68e76069eb9112c6f1bb7f15e1133fae04a8a77" \
  "$INPUT/whitener_union.npz" |
  sha256sum --check -
printf '%s  %s\n' \
  "cf5c1de31d03875b15f6fd11d3f628aadfccf04e3346493c831c3a6f441b4793" \
  "$INPUT/whitener_union.json" |
  sha256sum --check -

for donor in k562_gwps orion_hek293t; do
  mkdir -p "$INPUT/donors/$donor"
  aws s3 cp "s3://${BUCKET}/input/clean-hct/${donor}/" \
    "$INPUT/donors/$donor/" --recursive --only-show-errors
done

uv sync --frozen
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  python experiments/086-sparse-responder-emission/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report "$OUT/test-memory.json" \
  -- \
  uv run python \
    experiments/086-sparse-responder-emission/test_emitter.py \
  2>&1 | tee "$OUT/tests.log"

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 110 \
  --min-free-percent 0 \
  --poll-seconds 2 \
  -- \
  uv run python \
    experiments/086-sparse-responder-emission/run.py \
    --exp077-output "$INPUT/exp077-output" \
    --exp077-run-id "$EXP077_RUN_ID" \
    --exp077-source-archive "$INPUT/exp077-source.tar.gz" \
    --exp077-repo "$EXP077_REPO" \
    --exp077-launch-command-evidence "$INPUT/exp077-launch-command.json" \
    --exp077-terminal-command-evidence "$INPUT/exp077-terminal-command.json" \
    --exp077-terminal-invocation-evidence \
      "$INPUT/exp077-terminal-invocation.json" \
    --real "$INPUT/orion-hct116.h5ad" \
    --donors \
      "k562_gwps=$INPUT/donors/k562_gwps" \
      "orion_hek293t=$INPUT/donors/orion_hek293t" \
    --whitener "$INPUT/whitener_union.npz" \
    --whitener-metadata "$INPUT/whitener_union.json" \
    --panel-csv data/challenge/pert_counts.csv \
    --output "$OUT/evaluation" \
    --fractions 0.35 0.50 0.65 \
    --maximum-scorers 2 \
    --seed 8602026 \
    --guarded \
  2>&1 | tee "$OUT/run.log"

cat >"$OUT/infrastructure.json" <<EOF
{
  "provider": "AWS EC2",
  "exp077_run_id": "${EXP077_RUN_ID}",
  "source_archive_sha256": "${EXPECTED_SOURCE_SHA256}",
  "runner_sha256": "${EXPECTED_RUNNER_SHA256}",
  "reservation_sha256": "${EXPECTED_RESERVATION_SHA256}",
  "memory_guard_gib": 110,
  "maximum_parallel_scorers": 1,
  "submission": false
}
EOF

upload_outputs
