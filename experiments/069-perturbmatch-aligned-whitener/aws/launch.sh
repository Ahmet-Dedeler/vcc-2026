#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES only after Ahmet clears exp069" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/069-perturbmatch-aligned-whitener"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:-${EXP061_INSTANCE_ID:-i-0343d04439564b709}}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
EXP061_OUTPUT_PREFIX="${EXP061_OUTPUT_PREFIX:-output/exp061-perturbmatch-hct-count/attempt2}"
EXP062_OUTPUT_PREFIX="${EXP062_OUTPUT_PREFIX:-output/exp062-cd4-post-whiten}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp069-perturbmatch-aligned-whitener/${RUN_ID}"
OUTPUT_PREFIX="output/exp069-perturbmatch-aligned-whitener/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"
GATE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/exp069-gate.XXXXXX")"
trap 'rm -rf "$GATE_DIR"' EXIT

GATE_FILES=(
  "results.json|${EXP061_OUTPUT_PREFIX}/results.json"
  "baseline-agg_results.csv|${EXP061_OUTPUT_PREFIX}/baseline-agg_results.csv"
  "candidate-agg_results.csv|${EXP061_OUTPUT_PREFIX}/candidate-run/agg_results.csv"
  "build-report.json|${EXP061_OUTPUT_PREFIX}/build-report.json"
  "candidate.h5ad.sha256|${EXP061_OUTPUT_PREFIX}/candidate.h5ad.sha256"
  "exit-code.txt|${EXP061_OUTPUT_PREFIX}/exit-code.txt"
)
for mapping in "${GATE_FILES[@]}"; do
  local_name="${mapping%%|*}"
  key="${mapping#*|}"
  aws s3api head-object --bucket "$BUCKET" --key "$key" >/dev/null
  aws s3 cp "s3://${BUCKET}/${key}" "${GATE_DIR}/${local_name}" \
    --only-show-errors
done
candidate_object_bytes="$(aws s3api head-object \
  --bucket "$BUCKET" \
  --key "${EXP061_OUTPUT_PREFIX}/candidate.h5ad" \
  --query ContentLength \
  --output text)"
printf '%s\n' "$candidate_object_bytes" \
  > "${GATE_DIR}/candidate-object-bytes.txt"

if ! uv run python "$HERE/validate_exp061_gate.py" \
  --results "${GATE_DIR}/results.json" \
  --baseline-agg "${GATE_DIR}/baseline-agg_results.csv" \
  --candidate-agg "${GATE_DIR}/candidate-agg_results.csv" \
  --build-report "${GATE_DIR}/build-report.json" \
  --candidate-checksum "${GATE_DIR}/candidate.h5ad.sha256" \
  --exit-code "${GATE_DIR}/exit-code.txt" \
  --candidate-object-bytes "$candidate_object_bytes"
then
  echo "experiment 061 prerequisite failed; refusing experiment 069" >&2
  exit 65
fi

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create an instance" >&2
  exit 66
fi

cd "$ROOT"
tar \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/069-perturbmatch-aligned-whitener/aws/source.tar.gz' \
  -czf "$ARCHIVE" \
  README.md \
  pyproject.toml \
  uv.lock \
  src \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/069-perturbmatch-aligned-whitener

aws s3 cp "$ARCHIVE" "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
for mapping in "${GATE_FILES[@]}"; do
  local_name="${mapping%%|*}"
  aws s3 cp "${GATE_DIR}/${local_name}" \
    "s3://${BUCKET}/${SOURCE_PREFIX}/exp061-gate/${local_name}" \
    --only-show-errors
done
aws s3 cp "${GATE_DIR}/candidate-object-bytes.txt" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/exp061-gate/candidate-object-bytes.txt" \
  --only-show-errors

for held in orion_hct116 orion_hek293t; do
  for suffix in npz json; do
    aws s3 cp \
      "s3://${BUCKET}/${EXP062_OUTPUT_PREFIX}/whiteners/held_${held}.${suffix}" \
      "s3://${BUCKET}/${SOURCE_PREFIX}/baseline-whiteners/held_${held}.${suffix}" \
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
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp069-run.sh --only-show-errors\",
    \"chmod +x /root/exp069-run.sh\",
    \"systemd-run --unit=vcc-exp069-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp069-run.sh\"
  ]" \
  --comment "Launch gated exp069 effect-level matched-whitener test; no submission" \
  --query 'Command.CommandId' \
  --output text)"

cat > "${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "exp061_output_prefix": "s3://${BUCKET}/${EXP061_OUTPUT_PREFIX}/",
  "exp061_count_level_avg6_gate": "nonnegative",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "ssm_command_id": "${command_id}",
  "cell_eval2": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
