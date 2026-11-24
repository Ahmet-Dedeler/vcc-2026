#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after explicit approval" >&2
  exit 64
fi

for command in aws awk grep python3 sha256sum tar; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/087-string-target-kernel"
HELD_OPERATOR_DIR="${ROOT}/experiments/062-cd4-post-whiten/cloud-output/whiteners"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an idle running >=64 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
if [[ ! "$RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 69
fi
SOURCE_PREFIX="input/exp087-string-target-kernel/${RUN_ID}"
OUTPUT_PREFIX="output/exp087-string-target-kernel/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

assert_unused_s3_prefix() {
  local prefix="$1"
  local key_count
  key_count="$(aws s3api list-objects-v2 \
    --region "$REGION" \
    --bucket "$BUCKET" \
    --prefix "${prefix%/}/" \
    --max-keys 1 \
    --query 'KeyCount' \
    --output text)"
  if [[ ! "$key_count" =~ ^[0-9]+$ ]]; then
    echo "could not verify S3 prefix is unused: s3://${BUCKET}/${prefix}/" >&2
    exit 78
  fi
  if [[ "$key_count" -ne 0 ]]; then
    echo "refusing reused S3 prefix: s3://${BUCKET}/${prefix}/" >&2
    exit 79
  fi
}

assert_unused_s3_prefix "$SOURCE_PREFIX"
assert_unused_s3_prefix "$OUTPUT_PREFIX"

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
preflight_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters "$preflight_parameters" \
  --comment "Read-only exp087 memory and process-collision preflight" \
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
  --query 'StandardOutputContent' \
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
(
  cd "$HELD_OPERATOR_DIR"
  printf '%s  %s\n' \
    "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc" \
    "held_orion_hct116.npz" \
    "c4e00ad1fa5fbba7b9b0f93e4dae91b76da73638ef9431316923bb6ca4b585bc" \
    "held_orion_hct116.json" \
    "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1" \
    "held_orion_hek293t.npz" \
    "cb22406a7a529380c4348bf636f2532c29a43ab0fc4b00d6f6d98fb5e02360c6" \
    "held_orion_hek293t.json" |
    sha256sum --check -
)
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/087-string-target-kernel/aws/source.tar.gz' \
  --exclude='experiments/087-string-target-kernel/local-test-memory.json' \
  -czf "$ARCHIVE" \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  data/donors/orion_hct116/genes.txt \
  data/donors/orion_hct116/perts.txt \
  data/donors/orion_hek293t/genes.txt \
  data/donors/orion_hek293t/perts.txt \
  data/donors/rpe1/genes.txt \
  data/donors/rpe1/perts.txt \
  data/donors/nadig_jurkat/genes.txt \
  data/donors/nadig_jurkat/perts.txt \
  data/external/slim/README.md \
  data/external/slim/gene_string_embeddings.v0.3.h5 \
  experiments/044-csls-distill/cache/meta.json \
  experiments/044-csls-distill/cache/master_genes.txt \
  experiments/044-csls-distill/cache/train_perts.txt \
  experiments/044-csls-distill/cache/outer_perts.txt \
  experiments/044-csls-distill/cache/train_mask.npy \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/087-string-target-kernel

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
      "${HELD_OPERATOR_DIR}/held_${held}.${suffix}" \
      "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_${held}.${suffix}" \
      --only-show-errors
  done
done
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.npz" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_kolf2.npz" \
  --only-show-errors
aws s3 cp \
  "s3://${BUCKET}/input/production-challenge/whitener_union.json" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/operators/held_kolf2.json" \
  --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp087-run.sh --only-show-errors\",
    \"printf '${runner_sha256}  /root/exp087-run.sh\\\\n' | sha256sum --check -\",
    \"chmod +x /root/exp087-run.sh\",
    \"systemd-run --unit=vcc-exp087-${RUN_ID} --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=EXPECTED_SOURCE_SHA256=${source_archive_sha256} --setenv=EXPECTED_RUNNER_SHA256=${runner_sha256} /bin/bash /root/exp087-run.sh\"
  ]" \
  --comment "Launch exp087 leakage-safe STRING residual; no submission" \
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
  "ssm_command_id": "${command_id}",
  "creates_or_starts_instance": false,
  "bundles_only_pinned_string_embedding": true,
  "downloads_expression_bundle": false,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
