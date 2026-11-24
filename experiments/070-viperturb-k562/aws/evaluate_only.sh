#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing evaluation retry: set I_HAVE_EXPLICIT_CLEARANCE=YES" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/070-viperturb-k562"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?INSTANCE_ID must name the existing exp070 worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
BASE_RUN_ID="${BASE_RUN_ID:-20260823T063300Z}"
RETRY_ID="${RETRY_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
NULL_REPS="${NULL_REPS:-1000}"
SOURCE_PREFIX="input/exp070-viperturb-k562/${BASE_RUN_ID}/evaluation-retry/${RETRY_ID}"
OUTPUT_PREFIX="output/exp070-viperturb-k562/${BASE_RUN_ID}/evaluation-retry/${RETRY_ID}"
ARCHIVE="$(mktemp "${TMPDIR:-/tmp}/exp070-eval-source.XXXXXX")"

cleanup() {
  rm -f "$ARCHIVE"
}
trap cleanup EXIT

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

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/evaluate_only_remote.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/evaluate_only_remote.sh" \
  --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/evaluate_only_remote.sh /root/exp070-evaluate-only-${RETRY_ID}.sh --only-show-errors\",
    \"chmod +x /root/exp070-evaluate-only-${RETRY_ID}.sh\",
    \"systemd-run --unit=vcc-exp070-eval-${RETRY_ID} --setenv=BUCKET=${BUCKET} --setenv=BASE_RUN_ID=${BASE_RUN_ID} --setenv=RETRY_ID=${RETRY_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=NULL_REPS=${NULL_REPS} /bin/bash /root/exp070-evaluate-only-${RETRY_ID}.sh\"
  ]" \
  --comment "Evaluation-only exp070 retry using preserved compact; no RDS read or submission" \
  --query 'Command.CommandId' \
  --output text)"

printf 'base_run_id=%s\nretry_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$BASE_RUN_ID" "$RETRY_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
