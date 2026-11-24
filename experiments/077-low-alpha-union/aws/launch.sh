#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES for exp077" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/077-low-alpha-union"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to a running >=128 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp077-low-alpha-union/${RUN_ID}"
OUTPUT_PREFIX="output/exp077-low-alpha-union/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

read -r state memory_mib < <(
  aws ec2 describe-instances \
    --region "$REGION" \
    --instance-ids "$INSTANCE_ID" \
    --query 'Reservations[0].Instances[0].[State.Name,InstanceType]' \
    --output text |
  while read -r instance_state instance_type; do
    memory="$(aws ec2 describe-instance-types \
      --region "$REGION" \
      --instance-types "$instance_type" \
      --query 'InstanceTypes[0].MemoryInfo.SizeInMiB' \
      --output text)"
    printf '%s %s\n' "$instance_state" "$memory"
  done
)
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create it" >&2
  exit 65
fi
if (( memory_mib < 120000 )); then
  echo "worker has only ${memory_mib} MiB RAM; need at least 120000 MiB" >&2
  exit 66
fi

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/077-low-alpha-union/aws/source.tar.gz' \
  -czf "$ARCHIVE" \
  README.md \
  pyproject.toml \
  uv.lock \
  src \
  data/challenge/pert_counts.csv \
  experiments/031-somatic/run.py \
  experiments/077-low-alpha-union \
  scripts/run_with_memory_guard.py

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" --only-show-errors
aws s3 cp data/harness/orion-hct116.h5ad \
  "s3://${BUCKET}/${SOURCE_PREFIX}/orion-hct116.h5ad" --only-show-errors
for filename in whitener_union.npz whitener_union.json; do
  aws s3 cp \
    "s3://${BUCKET}/output-clean-hct/${filename}" \
    "s3://${BUCKET}/${SOURCE_PREFIX}/${filename}" \
    --only-show-errors
done
aws s3 cp \
  "s3://${BUCKET}/output-clean-hct/g0.65_q1.0_a1.08_spread0.25/agg_results.csv" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/baseline-agg_results.csv" \
  --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp077-run.sh --only-show-errors\",
    \"chmod +x /root/exp077-run.sh\",
    \"systemd-run --unit=vcc-exp077-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp077-run.sh\"
  ]" \
  --comment "Launch low-alpha union-whitened HCT count test; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat >"${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "cell_eval2_jobs": 2,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"

