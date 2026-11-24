#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES only after Ahmet approves exp073" >&2
  exit 64
fi

for command in aws tar grep; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/073-postwhiten-csls-distill"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:-i-0e05556f09944111f}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp073-postwhiten-csls-distill/${RUN_ID}"
OUTPUT_PREFIX="output/exp073-postwhiten-csls-distill/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create an instance" >&2
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
  echo "worker ${INSTANCE_ID} (${instance_type}) has ${memory_mib} MiB; require >=65536 MiB" >&2
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

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/073-postwhiten-csls-distill/aws/source.tar.gz' \
  --exclude='experiments/073-postwhiten-csls-distill/aws/source.test.tar.gz' \
  --exclude='experiments/073-postwhiten-csls-distill/cloud-output/*' \
  --exclude='experiments/073-postwhiten-csls-distill/query-cache/*' \
  --exclude='experiments/073-postwhiten-csls-distill/launch-metadata.json' \
  -czf "$ARCHIVE" \
  README.md \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  experiments/073-postwhiten-csls-distill

if tar -tzf "$ARCHIVE" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
  echo "refusing upload: source archive contains macOS metadata files" >&2
  exit 67
fi

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" \
  --only-show-errors
for held in k562_gwps orion_hct116 orion_hek293t; do
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
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp073-run.sh --only-show-errors\",
    \"chmod +x /root/exp073-run.sh\",
    \"systemd-run --unit=vcc-exp073-${RUN_ID} --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp073-run.sh\"
  ]" \
  --comment "Launch exp073 frozen-post-whitening CSLS distillation; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat >"${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "instance_type": "${instance_type}",
  "instance_memory_mib": ${memory_mib},
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "reuses_exp062_operators": true,
  "refits_whitening_or_full_covariance_eigensystem": false,
  "uses_randomized_low_rank_distillation_svd": true,
  "cell_eval2": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
