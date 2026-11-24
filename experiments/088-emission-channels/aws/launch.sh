#!/usr/bin/env bash
# Stage experiment 088 onto one already-running worker and start it via SSM.
#
#   INSTANCE_ID=i-... HARNESS=h1-test GRID=grid-ab.json \
#   DONOR_NAMES="k562_gwps orion_hct116 orion_hek293t" JOBS=3 ./launch.sh
#
# Refuses to create or start instances: provisioning is a separate, deliberate
# act. Same contract as experiments/077-low-alpha-union/aws/launch.sh.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/088-emission-channels"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to a running worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
HARNESS="${HARNESS:?set HARNESS}"
GRID="${GRID:?set GRID}"
DONOR_NAMES="${DONOR_NAMES:?set DONOR_NAMES}"
JOBS="${JOBS:-3}"
SOURCE_PREFIX="input/exp088-emission-channels/${RUN_ID}"
OUTPUT_PREFIX="output/exp088-emission-channels/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

state="$(aws ec2 describe-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[].Instances[].State.Name' --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create it" >&2
  exit 65
fi

if [[ ! -f "${HERE}/${GRID}" ]]; then
  echo "grid ${HERE}/${GRID} does not exist" >&2
  exit 66
fi

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' --exclude='*.pyc' \
  --exclude='experiments/088-emission-channels/aws/source.tar.gz' \
  -czf "$ARCHIVE" \
  README.md pyproject.toml uv.lock src \
  data/challenge/pert_counts.csv \
  experiments/031-somatic/run.py \
  experiments/088-emission-channels \
  scripts/run_with_memory_guard.py

aws s3 cp "$ARCHIVE" "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" --only-show-errors
aws s3 cp "${HERE}/aws/run.sh" "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --comment "exp088 ${HARNESS} ${GRID} ${RUN_ID}" \
  --timeout-seconds 3600 \
  --parameters "commands=[\
\"set -Eeuo pipefail\",\
\"mkdir -p /root\",\
\"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp088-run.sh --only-show-errors\",\
\"chmod +x /root/exp088-run.sh\",\
\"systemd-run --unit=vcc-exp088-${RUN_ID} --collect \
--setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} \
--setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} \
--setenv=HARNESS=${HARNESS} --setenv=GRID=${GRID} \
--setenv=DONOR_NAMES='${DONOR_NAMES}' --setenv=JOBS=${JOBS} \
/bin/bash /root/exp088-run.sh\"\
]" \
  --query 'Command.CommandId' --output text)"

cat >"${HERE}/launch-metadata-${RUN_ID}.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "harness": "${HARNESS}",
  "grid": "${GRID}",
  "donors": "${DONOR_NAMES}",
  "jobs": ${JOBS},
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "submission": false
}
EOF
echo "launched ${RUN_ID} on ${INSTANCE_ID} (${HARNESS}, ${GRID}, jobs=${JOBS})"
echo "  ssm command ${command_id}"
echo "  outputs s3://${BUCKET}/${OUTPUT_PREFIX}/"
