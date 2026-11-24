#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES for exp086" >&2
  exit 64
fi

fail() {
  echo "$1" >&2
  exit "${2:-64}"
}

require_match() {
  local label="$1"
  local value="$2"
  local pattern="$3"
  if [[ ! "$value" =~ $pattern ]]; then
    printf 'unsafe %s: %q\n' "$label" "$value" >&2
    exit 67
  fi
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/086-sparse-responder-emission"
TRUSTED_PROVENANCE="${HERE}/exp077-trusted-provenance.json"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an existing >=128 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
EXP077_RUN_ID="${EXP077_RUN_ID:-20260823T060938Z}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"

# Validate every externally supplied identifier before constructing paths,
# S3 keys, systemd unit names, JSON, or SSM shell text from it.
require_match "run id" "$RUN_ID" \
  '^[0-9]{8}T[0-9]{6}Z(-[a-z0-9][a-z0-9-]{0,31})?$'
require_match "exp077 run id" "$EXP077_RUN_ID" \
  '^[0-9]{8}T[0-9]{6}Z(-[a-z0-9][a-z0-9-]{0,31})?$'
require_match "instance id" "$INSTANCE_ID" '^i-[0-9a-f]{8,17}$'
require_match "AWS region" "$REGION" '^[a-z]{2}(-gov)?-[a-z]+-[0-9]$'
require_match "bucket" "$BUCKET" '^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$'
[[ "$BUCKET" != *..* && "$BUCKET" != *.-* && "$BUCKET" != *-.* ]] ||
  fail "unsafe bucket punctuation" 67
[[ "$REGION" == "us-east-1" ]] || fail "region differs from frozen provenance" 67
[[ "$BUCKET" == "<your-s3-bucket>" ]] ||
  fail "bucket differs from frozen provenance" 67
[[ "$EXP077_RUN_ID" == "20260823T060938Z" ]] ||
  fail "exp077 run differs from frozen provenance" 67

SOURCE_PREFIX="input/exp086-sparse-responder-emission/${RUN_ID}"
OUTPUT_PREFIX="output/exp086-sparse-responder-emission/${RUN_ID}"
EXP077_OUTPUT_PREFIX="output/exp077-low-alpha-union/${EXP077_RUN_ID}"
METADATA="${HERE}/launch-metadata-${RUN_ID}.json"
[[ ! -e "$METADATA" && ! -L "$METADATA" ]] ||
  fail "local launch metadata already exists: ${METADATA}" 79

PREFLIGHT="$(mktemp -d "${TMPDIR:-/tmp}/exp086-preflight.${RUN_ID}.XXXXXX")"
ARCHIVE="${PREFLIGHT}/source.${RUN_ID}.tar.gz"
RESERVATION="${PREFLIGHT}/launch-reservation.json"
cleanup() {
  find "$PREFLIGHT" -depth -type f -delete 2>/dev/null || true
  find "$PREFLIGHT" -depth -type d -empty -delete 2>/dev/null || true
}
trap cleanup EXIT

attestor_sha256="$(
  python3 - "$TRUSTED_PROVENANCE" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1]))["required_terminal_attestation"]["attestor_sha256"])
PY
)"
require_match "attestor SHA-256" "$attestor_sha256" '^[0-9a-f]{64}$'
printf '%s  %s\n' "$attestor_sha256" "$HERE/aws/attest_exp077.py" |
  shasum -a 256 -c -

terminal_command_id="$(
  python3 - "$TRUSTED_PROVENANCE" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1]))
terminal = record.get("terminal_attestation")
print("" if not isinstance(terminal, dict) else terminal.get("ssm_command_id", ""))
PY
)"
if [[ -z "$terminal_command_id" ]]; then
  # Let gate.py provide the single authoritative fail-closed explanation.
  python3 "$HERE/gate.py" \
    --exp077-output "$PREFLIGHT" \
    --run-id "$EXP077_RUN_ID" \
    --trusted-provenance "$TRUSTED_PROVENANCE" \
    --launch-metadata-snapshot \
      "$HERE/exp077-launch-metadata.snapshot.json" \
    --output "$PREFLIGHT/gate.json"
  fail "exp077 provenance unexpectedly passed without terminal evidence" 80
fi
require_match "terminal SSM command id" "$terminal_command_id" \
  '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
launch_command_id="$(
  python3 - "$TRUSTED_PROVENANCE" <<'PY'
import json
import sys

print(json.load(open(sys.argv[1]))["launch"]["ssm_command_id"])
PY
)"
require_match "launch SSM command id" "$launch_command_id" \
  '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'

# Fetch terminal result files and AWS-authoritative command records. The gate
# rejects any mismatch before packaging or dispatch.
for filename in exit-code.txt summary.json results.json artifact-sha256.txt; do
  aws s3 cp \
    "s3://${BUCKET}/${EXP077_OUTPUT_PREFIX}/${filename}" \
    "$PREFLIGHT/${filename}" \
    --only-show-errors
done
aws ssm list-commands \
  --region "$REGION" \
  --command-id "$launch_command_id" \
  --output json >"$PREFLIGHT/launch-command.json"
aws ssm list-commands \
  --region "$REGION" \
  --command-id "$terminal_command_id" \
  --output json >"$PREFLIGHT/terminal-command.json"
aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$terminal_command_id" \
  --instance-id "$INSTANCE_ID" \
  --output json >"$PREFLIGHT/terminal-invocation.json"
python3 "$HERE/gate.py" \
  --exp077-output "$PREFLIGHT" \
  --run-id "$EXP077_RUN_ID" \
  --trusted-provenance "$TRUSTED_PROVENANCE" \
  --launch-metadata-snapshot "$HERE/exp077-launch-metadata.snapshot.json" \
  --launch-command-evidence "$PREFLIGHT/launch-command.json" \
  --terminal-command-evidence "$PREFLIGHT/terminal-command.json" \
  --terminal-invocation-evidence "$PREFLIGHT/terminal-invocation.json" \
  --output "$PREFLIGHT/gate.json"

prefix_is_empty() {
  local prefix="$1"
  local count
  count="$(
    aws s3api list-objects-v2 \
      --bucket "$BUCKET" \
      --prefix "${prefix}/" \
      --max-keys 1 \
      --query KeyCount \
      --output text
  )"
  [[ "$count" == "0" ]]
}
prefix_is_empty "$SOURCE_PREFIX" ||
  fail "S3 input prefix already exists: s3://${BUCKET}/${SOURCE_PREFIX}/" 74
prefix_is_empty "$OUTPUT_PREFIX" ||
  fail "S3 output prefix already exists: s3://${BUCKET}/${OUTPUT_PREFIX}/" 74

read -r state instance_type memory_mib < <(
  aws ec2 describe-instances \
    --region "$REGION" \
    --instance-ids "$INSTANCE_ID" \
    --query 'Reservations[0].Instances[0].[State.Name,InstanceType]' \
    --output text |
  while read -r instance_state type; do
    require_match "instance type" "$type" '^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$'
    memory="$(aws ec2 describe-instance-types \
      --region "$REGION" \
      --instance-types "$type" \
      --query 'InstanceTypes[0].MemoryInfo.SizeInMiB' \
      --output text)"
    printf '%s %s %s\n' "$instance_state" "$type" "$memory"
  done
)
[[ "$state" == "running" ]] ||
  fail "worker ${INSTANCE_ID} is ${state}; refusing to start or create it" 65
[[ "$memory_mib" =~ ^[0-9]+$ ]] || fail "invalid worker memory value" 66
(( memory_mib >= 120000 )) ||
  fail "worker ${instance_type} has ${memory_mib} MiB; require >=120000" 66

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/086-sparse-responder-emission/launch-metadata-*.json' \
  -czf "$ARCHIVE" \
  pyproject.toml \
  uv.lock \
  data/challenge/pert_counts.csv \
  experiments/043-template-emitter/template_emitter.py \
  experiments/086-sparse-responder-emission \
  scripts/run_with_memory_guard.py

source_sha256="$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')"
runner_sha256="$(shasum -a 256 "$HERE/aws/run.sh" | awk '{print $1}')"
require_match "source SHA-256" "$source_sha256" '^[0-9a-f]{64}$'
require_match "runner SHA-256" "$runner_sha256" '^[0-9a-f]{64}$'
printf '{"run_id":"%s","source_sha256":"%s","runner_sha256":"%s"}\n' \
  "$RUN_ID" "$source_sha256" "$runner_sha256" >"$RESERVATION"
reservation_sha256="$(shasum -a 256 "$RESERVATION" | awk '{print $1}')"

# Conditional writes make the empty-prefix check race-safe. A failed partial
# launch deliberately leaves its content-addressed reservation behind.
aws s3api put-object \
  --bucket "$BUCKET" \
  --key "${SOURCE_PREFIX}/launch-reservation.json" \
  --body "$RESERVATION" \
  --if-none-match '*' >/dev/null
aws s3api put-object \
  --bucket "$BUCKET" \
  --key "${OUTPUT_PREFIX}/launch-reservation.json" \
  --body "$RESERVATION" \
  --if-none-match '*' >/dev/null
aws s3api put-object \
  --bucket "$BUCKET" \
  --key "${SOURCE_PREFIX}/source.tar.gz" \
  --body "$ARCHIVE" \
  --if-none-match '*' >/dev/null
aws s3api put-object \
  --bucket "$BUCKET" \
  --key "${SOURCE_PREFIX}/run.sh" \
  --body "$HERE/aws/run.sh" \
  --if-none-match '*' >/dev/null

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp086-run-${RUN_ID}.sh --only-show-errors\",
    \"chmod +x /root/exp086-run-${RUN_ID}.sh\",
    \"systemd-run --unit=vcc-exp086-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=EXP077_RUN_ID=${EXP077_RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=EXPECTED_SOURCE_SHA256=${source_sha256} --setenv=EXPECTED_RUNNER_SHA256=${runner_sha256} --setenv=EXPECTED_RESERVATION_SHA256=${reservation_sha256} /bin/bash /root/exp086-run-${RUN_ID}.sh\"
  ]" \
  --comment "Launch gated exp086 sparse-responder HCT study; no submission" \
  --query 'Command.CommandId' \
  --output text)"
require_match "SSM command id" "$command_id" \
  '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'

set -o noclobber
cat >"$METADATA" <<EOF
{
  "run_id": "${RUN_ID}",
  "exp077_run_id": "${EXP077_RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${instance_type}",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "source_archive_sha256": "${source_sha256}",
  "runner_sha256": "${runner_sha256}",
  "reservation_sha256": "${reservation_sha256}",
  "ssm_command_id": "${command_id}",
  "submission": false
}
EOF
set +o noclobber

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\nmetadata=%s\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX" "$METADATA"
