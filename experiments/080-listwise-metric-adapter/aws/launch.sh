#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after explicit approval" >&2
  exit 64
fi

for command in aws python3 sha256sum tar grep; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/080-listwise-metric-adapter"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an idle running >=64 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 69
fi
SOURCE_PREFIX="input/exp080-listwise-metric-adapter/${RUN_ID}"
OUTPUT_PREFIX="output/exp080-listwise-metric-adapter/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

for prefix in "$SOURCE_PREFIX" "$OUTPUT_PREFIX"; do
  existing_key="$(aws s3api list-objects-v2 \
    --region "$REGION" \
    --bucket "$BUCKET" \
    --prefix "${prefix}/" \
    --max-keys 1 \
    --query 'Contents[0].Key' \
    --output text)"
  if [[ "$existing_key" != "None" && -n "$existing_key" ]]; then
    echo "refusing reused RUN_ID ${RUN_ID}: s3://${BUCKET}/${prefix}/ exists" >&2
    exit 74
  fi
done

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to create or start it" >&2
  exit 65
fi
instance_type="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].InstanceType' \
  --output text)"
memory_mib="$(aws ec2 describe-instance-types \
  --region "$REGION" \
  --instance-types "$instance_type" \
  --query 'InstanceTypes[0].MemoryInfo.SizeInMiB' \
  --output text)"
if [[ ! "$memory_mib" =~ ^[0-9]+$ || "$memory_mib" -lt 65536 ]]; then
  echo "worker ${INSTANCE_ID} (${instance_type}) has ${memory_mib} MiB; require >=65536" >&2
  exit 66
fi
ping_status="$(aws ssm describe-instance-information \
  --region "$REGION" \
  --filters "Key=InstanceIds,Values=${INSTANCE_ID}" \
  --query 'InstanceInformationList[0].PingStatus' \
  --output text)"
if [[ "$ping_status" != "Online" ]]; then
  echo "worker ${INSTANCE_ID} SSM status is ${ping_status}; require Online" >&2
  exit 68
fi

preflight_parameters="$(python3 - <<'PY'
import json

commands = [
    "set -Eeuo pipefail",
    "available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)",
    "test -n \"$available_kib\" && test \"$available_kib\" -ge $((48 * 1024 * 1024))",
    "mkdir -p /work",
    "flock -n /work/vcc-heavy.lock -c true",
    (
        "if pgrep -af '[c]ell-eval|[p]ython(3)? .*experiments/|"
        "[u]v (run|tool|x) |[t]orchrun|[a]ccelerate launch|[r]ay::|"
        "[m]ake_submission|[b]uild_submission|[r]un_with_memory_guard|"
        "[s]tage_cache.py'; then exit 72; fi"
    ),
    (
        "if systemctl list-units --type=service --state=running --no-legend "
        "| awk '{print $1}' | grep -E '^(vcc-|cell-eval)'; then exit 73; fi"
    ),
    "printf 'MemAvailableKiB=%s\\n' \"$available_kib\"",
]
print(json.dumps({"commands": commands}))
PY
)"
preflight_command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters "$preflight_parameters" \
  --comment "Read-only exp080 memory and process-collision preflight" \
  --query 'Command.CommandId' \
  --output text)"
if ! aws ssm wait command-executed \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID"
then
  aws ssm get-command-invocation \
    --region "$REGION" \
    --command-id "$preflight_command_id" \
    --instance-id "$INSTANCE_ID" \
    --output json >&2 || true
  echo "worker preflight failed; refusing dispatch" >&2
  exit 70
fi
preflight_status="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" \
  --query Status \
  --output text)"
preflight_output="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" \
  --query StandardOutputContent \
  --output text)"
if [[ "$preflight_status" != "Success" ]]; then
  echo "worker preflight status is ${preflight_status}; refusing dispatch" >&2
  exit 71
fi
available_memory_kib="${preflight_output##*=}"
available_memory_kib="${available_memory_kib//$'\n'/}"
if [[ ! "$available_memory_kib" =~ ^[0-9]+$ ]]; then
  echo "could not parse worker available memory: ${preflight_output}" >&2
  exit 72
fi

cd "$ROOT"
sha256sum -c "$HERE/dependency-sha256.txt"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/080-listwise-metric-adapter/aws/source.tar.gz' \
  --exclude='experiments/080-listwise-metric-adapter/local-test-memory.json' \
  -czf "$ARCHIVE" \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/080-listwise-metric-adapter

if tar -tzf "$ARCHIVE" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
  echo "refusing upload: source archive contains macOS metadata" >&2
  exit 67
fi
source_archive_sha256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"
runner_sha256="$(sha256sum "$HERE/aws/run.sh" | awk '{print $1}')"
aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" \
  --only-show-errors

# The worker role is input-only. Copy the exact frozen fold operators into the
# run's own input prefix before dispatch.
for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/output/exp062-cd4-post-whiten/whiteners/held_${held}.${suffix}" \
      "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp080-run.sh --only-show-errors\",
    \"printf '${runner_sha256}  /root/exp080-run.sh\\\\n' | sha256sum --check -\",
    \"chmod +x /root/exp080-run.sh\",
    \"systemd-run --unit=vcc-exp080-${RUN_ID} --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=EXPECTED_SOURCE_SHA256=${source_archive_sha256} --setenv=EXPECTED_RUNNER_SHA256=${runner_sha256} /bin/bash /root/exp080-run.sh\"
  ]" \
  --comment "Launch exp080 leakage-safe listwise adapter; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat >"${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${instance_type}",
  "instance_memory_mib": ${memory_mib},
  "worker_preflight_command_id": "${preflight_command_id}",
  "worker_available_memory_kib_at_dispatch": ${available_memory_kib},
  "minimum_available_memory_gib": 48,
  "source_archive_sha256": "${source_archive_sha256}",
  "runner_sha256": "${runner_sha256}",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "creates_or_starts_instance": false,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
