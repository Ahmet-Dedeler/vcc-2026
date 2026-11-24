#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after explicit approval" >&2
  exit 64
fi

for command in aws awk grep python3 sha256sum tar uv; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/085-kolf-pm-gate"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an idle running >=64 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 69
fi
SOURCE_PREFIX="input/exp085-kolf-pm-gate/${RUN_ID}"
OUTPUT_PREFIX="output/exp085-kolf-pm-gate/${RUN_ID}"
GLOBAL_MARKER_KEY="output/exp085-kolf-pm-gate/KOLF_ACCESS_CONSUMED.json"
ARCHIVE="${HERE}/aws/source.tar.gz"

global_marker_count="$(aws s3api list-objects-v2 \
  --bucket "$BUCKET" \
  --prefix "$GLOBAL_MARKER_KEY" \
  --max-keys 1 \
  --query KeyCount \
  --output text)"
if [[ "$global_marker_count" != "0" ]]; then
  echo "experiment 085 global KOLF access marker exists; one-time evaluation is consumed" >&2
  exit 78
fi
prior_output_keys="$(aws s3api list-objects-v2 \
  --bucket "$BUCKET" \
  --prefix "output/exp085-kolf-pm-gate/" \
  --query 'Contents[].Key' \
  --output text)"
if tr '\t' '\n' <<<"$prior_output_keys" |
  grep -q '/KOLF_ACCESS_CONSUMED.json$'
then
  echo "experiment 085 already has a KOLF access marker; one-time evaluation is consumed" >&2
  exit 78
fi
for prefix in "$SOURCE_PREFIX" "$OUTPUT_PREFIX"; do
  key_count="$(aws s3api list-objects-v2 \
    --bucket "$BUCKET" \
    --prefix "${prefix}/" \
    --max-keys 1 \
    --query KeyCount \
    --output text)"
  if [[ "$key_count" != "0" ]]; then
    echo "refusing reused S3 prefix s3://${BUCKET}/${prefix}/" >&2
    exit 79
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

preflight_parameters="$(python3 - "$RUN_ID" <<'PY'
import json
import sys

run_id = sys.argv[1]
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
        "| awk '{print $1}' | grep -E '^(vcc-|cell-eval)' "
        f"| grep -v '^vcc-exp085-{run_id}\\\\.service$'; then exit 73; fi"
    ),
    f"test ! -e /work/exp085/{run_id}",
    "printf 'MemAvailableKiB=%s\\n' \"$available_kib\"",
]
print(json.dumps({"commands": commands}))
PY
)"
preflight_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters "$preflight_parameters" \
  --comment "Read-only exp085 memory, lock, and collision preflight" \
  --query 'Command.CommandId' \
  --output text)"
if ! aws ssm wait command-executed \
  --region "$REGION" \
  --command-id "$preflight_id" \
  --instance-id "$INSTANCE_ID"
then
  aws ssm get-command-invocation \
    --region "$REGION" \
    --command-id "$preflight_id" \
    --instance-id "$INSTANCE_ID" \
    --output json >&2 || true
  echo "worker preflight failed; refusing dispatch" >&2
  exit 70
fi
preflight_status="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_id" \
  --instance-id "$INSTANCE_ID" \
  --query Status \
  --output text)"
preflight_output="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_id" \
  --instance-id "$INSTANCE_ID" \
  --query StandardOutputContent \
  --output text)"
if [[ "$preflight_status" != "Success" ]]; then
  echo "worker preflight status is ${preflight_status}; refusing dispatch" >&2
  exit 71
fi
available_kib="${preflight_output##*=}"
available_kib="${available_kib//$'\n'/}"
if [[ ! "$available_kib" =~ ^[0-9]+$ ]]; then
  echo "could not parse worker available memory: ${preflight_output}" >&2
  exit 72
fi

cd "$ROOT"
sha256sum -c "$HERE/dependency-sha256.txt"
python3 "$HERE/memory_guard.py" \
  --max-rss-gib 1 \
  --poll-seconds 0.05 \
  --report "$HERE/local-launch-test-memory.json" \
  -- \
  uv run --frozen python "$HERE/test_run.py"

export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/085-kolf-pm-gate/aws/source.tar.gz' \
  --exclude='experiments/085-kolf-pm-gate/cloud-output/*' \
  --exclude='experiments/085-kolf-pm-gate/launch-metadata.json' \
  --exclude='experiments/085-kolf-pm-gate/*memory.json' \
  -czf "$ARCHIVE" \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  src/vcc2026/__init__.py \
  src/vcc2026/donor_transfer.py \
  src/vcc2026/multi_donor.py \
  src/vcc2026/whiten.py \
  scripts/build_submission.py \
  scripts/make_submission.py \
  scripts/make_submission_pm_cd4.py \
  experiments/061-perturbmatch-hct-count/build_candidate.py \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/code.sha256 \
  experiments/064-pm-cd4-composition/cloud-output/manifest.sha256 \
  experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/066-leakage-safe-target-gate/RESULTS.md \
  experiments/066-leakage-safe-target-gate/metrics.json \
  experiments/085-kolf-pm-gate

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
for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/output/exp062-cd4-post-whiten/whiteners/held_${held}.${suffix}" \
      "s3://${BUCKET}/${SOURCE_PREFIX}/held-operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done
for filename in \
  b_ctrl.npy \
  delta.npy \
  gene_map.csv \
  genes.txt \
  meta.json \
  perts.txt \
  provenance.json \
  source_rows_int8.npy \
  target_row_map.csv
do
  aws s3 cp \
    "data/external/kolf2/${filename}" \
    "s3://${BUCKET}/${SOURCE_PREFIX}/kolf2/${filename}" \
    --only-show-errors
done

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp085-run.sh --only-show-errors\",
    \"printf '${runner_sha256}  /root/exp085-run.sh\\\\n' | sha256sum --check -\",
    \"chmod +x /root/exp085-run.sh\",
    \"systemd-run --unit=vcc-exp085-${RUN_ID} --property=MemoryMax=44G --property=MemorySwapMax=0 --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=EXPECTED_SOURCE_SHA256=${source_archive_sha256} --setenv=EXPECTED_RUNNER_SHA256=${runner_sha256} --setenv=EXP085_ONE_TIME_CLEARANCE=YES --setenv=EXP085_DURABLE_FREEZE_URI=s3://${BUCKET}/${OUTPUT_PREFIX}/pre_kolf/ --setenv=EXP085_DURABLE_MARKER_URI=s3://${BUCKET}/${GLOBAL_MARKER_KEY} /bin/bash /root/exp085-run.sh\"
  ]" \
  --comment "Launch frozen exp085 one-time KOLF PM gate; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat >"${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${instance_type}",
  "instance_memory_mib": ${memory_mib},
  "worker_preflight_command_id": "${preflight_id}",
  "worker_available_memory_kib_at_dispatch": ${available_kib},
  "minimum_available_memory_gib": 48,
  "source_archive_sha256": "${source_archive_sha256}",
  "runner_sha256": "${runner_sha256}",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "durable_pre_kolf_freeze_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/pre_kolf/",
  "global_kolf_access_marker": "s3://${BUCKET}/${GLOBAL_MARKER_KEY}",
  "ssm_command_id": "${command_id}",
  "systemd_memory_max_gib": 44,
  "systemd_swap_max_gib": 0,
  "creates_or_starts_instance": false,
  "one_time_kolf_evaluation": true,
  "count_generation": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
