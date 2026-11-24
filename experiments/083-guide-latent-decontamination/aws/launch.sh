#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES after explicit approval" >&2
  exit 64
fi

for command in aws grep python3 sha256sum tar uv; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required local command is missing: ${command}" >&2
    exit 63
  fi
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/083-guide-latent-decontamination"
H5AD="${H5AD:-${ROOT}/data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad}"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:?set INSTANCE_ID to an idle running >=64 GiB worker}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
SOURCE_PREFIX="input/exp083-guide-latent-decontamination/${RUN_ID}"
OUTPUT_PREFIX="output/exp083-guide-latent-decontamination/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"

if [[ ! "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "unsafe RUN_ID: ${RUN_ID}" >&2
  exit 69
fi

cd "$ROOT"
sha256sum -c "$HERE/dependency-sha256.txt"

# This gate intentionally precedes all AWS state inspection, archive creation,
# H5AD upload, and remote dispatch. The current source exits 78 here (8 < 100).
PYTHONDONTWRITEBYTECODE=1 \
python3 "$HERE/memory_guard.py" \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report "$HERE/local-preflight-memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    "$HERE/inspect_source.py" \
    --h5ad "$H5AD" \
    --panel-csv "$ROOT/data/challenge/pert_counts.csv" \
    --allow-network \
    --hash-h5ad \
    --output "$HERE/source-audit-launch.json" \
    --require-promotion-coverage

PYTHONDONTWRITEBYTECODE=1 \
python3 "$HERE/memory_guard.py" \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report "$HERE/local-test-memory.json" \
  -- \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    "$HERE/test_model.py"

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

preflight_command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"awk '/^MemAvailable:/ {print \\\\$2}' /proc/meminfo\",
    \"if pgrep -af '[c]ell-eval2|[e]xperiments/0[0-9][0-9]-.*\\\\.py|[e]xtract_guide_effects.py|[f]it.py|[e]valuate.py'; then exit 42; fi\"
  ]" \
  --comment "Read-only exp083 worker availability preflight" \
  --query 'Command.CommandId' \
  --output text)"
aws ssm wait command-executed \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" || true
preflight_status="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" \
  --query 'Status' \
  --output text)"
preflight_code="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" \
  --query 'ResponseCode' \
  --output text)"
preflight_output="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$preflight_command_id" \
  --instance-id "$INSTANCE_ID" \
  --query 'StandardOutputContent' \
  --output text)"
if [[ "$preflight_status" != "Success" || "$preflight_code" != "0" ]]; then
  echo "worker availability preflight failed (status=${preflight_status}, code=${preflight_code})" >&2
  printf '%s\n' "$preflight_output" >&2
  exit 70
fi
available_kib="$(printf '%s\n' "$preflight_output" | head -n 1)"
minimum_available_kib=$((55 * 1024 * 1024))
if [[ ! "$available_kib" =~ ^[0-9]+$ || "$available_kib" -lt "$minimum_available_kib" ]]; then
  echo "worker has ${available_kib:-0} KiB available; require >=55 GiB" >&2
  exit 71
fi

export COPYFILE_DISABLE=1
export COPY_EXTENDED_ATTRIBUTES_DISABLE=1
tar \
  --exclude='._*' \
  --exclude='.DS_Store' \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/083-guide-latent-decontamination/aws/source.tar.gz' \
  --exclude='experiments/083-guide-latent-decontamination/source-audit-launch.json' \
  --exclude='experiments/083-guide-latent-decontamination/local-*-memory.json' \
  -czf "$ARCHIVE" \
  data/challenge/gene_names.csv \
  data/challenge/pert_counts.csv \
  data/external/perturbmatch/perts.txt \
  experiments/062-cd4-post-whiten/run.py \
  experiments/064-pm-cd4-composition/run.py \
  experiments/064-pm-cd4-composition/cloud-output/results.json \
  experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv \
  experiments/068-seed-offtarget-ablation/run.py \
  experiments/083-guide-latent-decontamination

if tar -tzf "$ARCHIVE" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
  echo "refusing upload: source archive contains macOS metadata" >&2
  exit 67
fi
archive_sha256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"
runner_sha256="$(sha256sum "$HERE/aws/run.sh" | awk '{print $1}')"

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" --only-show-errors
aws s3 cp "$H5AD" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/ReplogleWeissman2022_K562_gwps.h5ad" \
  --only-show-errors
for filename in delta.npy b_ctrl.npy genes.txt perts.txt meta.json; do
  aws s3 cp \
    "$ROOT/data/external/kolf2/${filename}" \
    "s3://${BUCKET}/${SOURCE_PREFIX}/kolf2/${filename}" \
    --only-show-errors
done
for held in orion_hct116 orion_hek293t; do
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
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp083-run.sh --only-show-errors\",
    \"printf '%s  %s\\\\n' '${runner_sha256}' '/root/exp083-run.sh' | sha256sum -c -\",
    \"chmod +x /root/exp083-run.sh\",
    \"systemd-run --unit=vcc-exp083-${RUN_ID} --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=SOURCE_ARCHIVE_SHA256=${archive_sha256} /bin/bash /root/exp083-run.sh\"
  ]" \
  --comment "Launch exp083 guide latent decontamination; no submission" \
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
  "worker_available_memory_kib_at_dispatch": ${available_kib},
  "source_archive_sha256": "${archive_sha256}",
  "remote_runner_sha256": "${runner_sha256}",
  "creates_or_starts_instance": false,
  "preflight_required_before_upload": true,
  "cell_eval2": false,
  "count_generation": false,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
