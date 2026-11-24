#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after explicit approval" >&2
  exit 64
fi

for command in aws tar grep sha256sum python3 uv; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/079-kolf-local-centering"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:-i-0e05556f09944111f}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp079-kolf-local-centering/${RUN_ID}"
OUTPUT_PREFIX="output/exp079-kolf-local-centering/${RUN_ID}"
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
  exit 67
fi

cd "$ROOT"
sha256sum -c "$HERE/dependency-sha256.txt"
for filename in delta.npy b_ctrl.npy genes.txt perts.txt meta.json; do
  expected="$(python3 - "$filename" <<'PY'
import ast
import sys
from pathlib import Path

tree = ast.parse(Path("experiments/079-kolf-local-centering/run.py").read_text())
for node in tree.body:
    if isinstance(node, ast.Assign):
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == "EXPECTED_KOLF_HASHES":
                values = ast.literal_eval(node.value)
                print(values[sys.argv[1]])
                raise SystemExit
raise SystemExit("EXPECTED_KOLF_HASHES not found")
PY
)"
  printf '%s  %s\n' "$expected" "data/external/kolf2/${filename}" |
    sha256sum -c -
done

PYTHONDONTWRITEBYTECODE=1 \
python3 "$HERE/memory_guard.py" \
  --max-rss-gib 1 \
  --poll-seconds 0.25 \
  --report "$HERE/local-launch-test-memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy==2.5.2 \
  python "$HERE/test_run.py"

cd "$ROOT"
export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/079-kolf-local-centering/aws/source.tar.gz' \
  --exclude='experiments/079-kolf-local-centering/cloud-output/*' \
  --exclude='experiments/079-kolf-local-centering/launch-metadata.json' \
  --exclude='experiments/079-kolf-local-centering/*memory.json' \
  -czf "$ARCHIVE" \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  src/vcc2026/__init__.py \
  src/vcc2026/donor_transfer.py \
  src/vcc2026/multi_donor.py \
  src/vcc2026/whiten.py \
  experiments/075-postwhiten-local-centering/run.py \
  experiments/075-postwhiten-local-centering/cloud-output/20260823T062845Z/results.json \
  experiments/079-kolf-local-centering

if tar -tzf "$ARCHIVE" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
  echo "refusing upload: source archive contains macOS metadata files" >&2
  exit 68
fi

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" \
  --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" \
  --only-show-errors
for filename in delta.npy b_ctrl.npy genes.txt perts.txt meta.json; do
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
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp079-run.sh --only-show-errors\",
    \"chmod +x /root/exp079-run.sh\",
    \"systemd-run --unit=vcc-exp079-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} /bin/bash /root/exp079-run.sh\"
  ]" \
  --comment "Launch frozen exp079 KOLF truth-only confirmation; no submission" \
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
  "worker_was_created_or_started": false,
  "retunes_on_kolf": false,
  "fits_eigensystem": false,
  "cell_eval2": false,
  "count_emission": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
