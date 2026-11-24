#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after approval" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/075-postwhiten-local-centering"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an idle running >=64 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp075-postwhiten-local-centering/${RUN_ID}"
OUTPUT_PREFIX="output/exp075-postwhiten-local-centering/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create it" >&2
  exit 65
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
  --exclude='experiments/075-postwhiten-local-centering/aws/source.tar.gz' \
  -czf "$ARCHIVE" \
  README.md \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  experiments/073-postwhiten-csls-distill/run.py \
  experiments/073-postwhiten-csls-distill/memory_guard.py \
  experiments/073-postwhiten-csls-distill/artifact_manifest.json \
  experiments/075-postwhiten-local-centering

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" \
  --only-show-errors

# The worker role is deliberately input-only. Stage the exact operators that
# experiment 073 already copied and hash-verified instead of asking the worker
# to read the shared output prefix.
for held in k562_gwps orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/input/exp073-postwhiten-csls-distill/20260823T060938Z/operators/held_${held}.${suffix}" \
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
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp075-run.sh --only-show-errors\",
    \"chmod +x /root/exp075-run.sh\",
    \"systemd-run --unit=vcc-exp075-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp075-run.sh\"
  ]" \
  --comment "Launch exp075 post-whitening reciprocal local centering; no submission" \
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
  "reuses_exp062_operators": true,
  "reuses_exp073_staging": true,
  "fits_eigensystem": false,
  "cell_eval2": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
