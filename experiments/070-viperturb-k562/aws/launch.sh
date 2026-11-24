#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES only after Ahmet clears exp070" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/070-viperturb-k562"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?INSTANCE_ID must name an already-running worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
EXP062_OUTPUT_PREFIX="${EXP062_OUTPUT_PREFIX:-output/exp062-cd4-post-whiten}"
NULL_REPS="${NULL_REPS:-1000}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp070-viperturb-k562/${RUN_ID}"
OUTPUT_PREFIX="output/exp070-viperturb-k562/${RUN_ID}"
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

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/070-viperturb-k562/aws/source.tar.gz' \
  --exclude='experiments/070-viperturb-k562/cloud-output/*' \
  -czf "$ARCHIVE" \
  README.md \
  pyproject.toml \
  uv.lock \
  src \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/070-viperturb-k562

aws s3 cp "$ARCHIVE" "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/${EXP062_OUTPUT_PREFIX}/whiteners/held_${held}.${suffix}" \
      "s3://${BUCKET}/${SOURCE_PREFIX}/whiteners/held_${held}.${suffix}" \
      --only-show-errors
  done
done
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp070-run.sh --only-show-errors\",
    \"chmod +x /root/exp070-run.sh\",
    \"systemd-run --unit=vcc-exp070-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=NULL_REPS=${NULL_REPS} /bin/bash /root/exp070-run.sh\"
  ]" \
  --comment "Launch VIPerturb extraction and frozen effect-level test; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat > "${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "zenodo_record": "18460279",
  "source_md5": "ac52c41233fc7a5a5c3810e0a075007f",
  "cell_eval2": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
