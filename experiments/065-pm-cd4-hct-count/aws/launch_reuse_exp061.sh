#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES only after Ahmet clears exp065" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/065-pm-cd4-hct-count"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${EXP061_INSTANCE_ID:-i-0343d04439564b709}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
EXP061_OUTPUT_PREFIX="${EXP061_OUTPUT_PREFIX:-output/exp061-perturbmatch-hct-count/attempt2}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp065-pm-cd4-hct-count/${RUN_ID}"
OUTPUT_PREFIX="output/exp065-pm-cd4-hct-count/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

for key in \
  "${EXP061_OUTPUT_PREFIX}/results.json" \
  "${EXP061_OUTPUT_PREFIX}/candidate-run/agg_results.csv" \
  "${EXP061_OUTPUT_PREFIX}/exit-code.txt"
do
  aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null
done
exp061_exit="$(aws s3 cp \
  "s3://${BUCKET}/${EXP061_OUTPUT_PREFIX}/exit-code.txt" - \
  --only-show-errors | tr -d '[:space:]')"
if [[ "$exp061_exit" != "0" ]]; then
  echo "experiment 061 is not successful: exit code ${exp061_exit}" >&2
  exit 65
fi

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
case "$state" in
  stopped)
    aws ec2 start-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
      >/dev/null
    aws ec2 wait instance-status-ok --region "$REGION" \
      --instance-ids "$INSTANCE_ID"
    ;;
  running)
    ;;
  *)
    echo "exp061 worker ${INSTANCE_ID} is ${state}; refusing to create a worker" >&2
    exit 66
    ;;
esac

cd "$ROOT"
tar \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  -czf "$ARCHIVE" \
  README.md \
  pyproject.toml \
  uv.lock \
  src \
  experiments/061-perturbmatch-hct-count/build_candidate.py \
  experiments/062-cd4-post-whiten/run.py \
  experiments/065-pm-cd4-hct-count/build_joint_candidate.py \
  experiments/065-pm-cd4-hct-count/compare_scores.py \
  experiments/065-pm-cd4-hct-count/memory_guard.py \
  experiments/065-pm-cd4-hct-count/test_candidate.py \
  experiments/065-pm-cd4-hct-count/aws/launch_reuse_exp061.sh \
  experiments/065-pm-cd4-hct-count/aws/run.sh

aws s3 cp "$ARCHIVE" "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/run.sh" "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" \
  --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp065-run.sh --only-show-errors\",
    \"chmod +x /root/exp065-run.sh\",
    \"systemd-run --unit=vcc-exp065-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp065-run.sh\"
  ]" \
  --comment "Launch frozen exp065 on completed exp061 worker; no submission" \
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
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
