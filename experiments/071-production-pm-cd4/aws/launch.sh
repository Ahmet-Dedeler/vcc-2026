#!/usr/bin/env bash
set -Eeuo pipefail

if [[ "${I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES" ]]; then
  echo "refusing launch: set I_HAVE_EXPLICIT_CLEARANCE=YES for exploratory packaging" >&2
  exit 64
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HERE="${ROOT}/experiments/071-production-pm-cd4"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_ID="${INSTANCE_ID:-i-0343d04439564b709}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
EXP065_RUN_ID="${EXP065_RUN_ID:?set the exact terminal experiment-065 run ID}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
EXP065_PREFIX="output/exp065-pm-cd4-hct-count/${EXP065_RUN_ID}"
SOURCE_PREFIX="input/exp071-production-pm-cd4/${RUN_ID}"
OUTPUT_PREFIX="output/exp071-production-pm-cd4/${RUN_ID}"
ARCHIVE="${HERE}/aws/source.tar.gz"
GATE_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/exp071-gate.XXXXXX")"
AUTH_STAGE="${GATE_ROOT}/authorization"
EXP065_STAGE="${AUTH_STAGE}/exp065"
HEAD_STAGE="${GATE_ROOT}/heads"

cleanup() {
  rm -rf "$GATE_ROOT"
}
trap cleanup EXIT

if [[ ! "$EXP065_RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "invalid EXP065_RUN_ID: ${EXP065_RUN_ID}" >&2
  exit 65
fi
if [[ ! "$RUN_ID" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "invalid RUN_ID: ${RUN_ID}" >&2
  exit 66
fi

state="$(aws ec2 describe-instances \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].State.Name' \
  --output text)"
if [[ "$state" != "running" ]]; then
  echo "worker ${INSTANCE_ID} is ${state}; refusing to start or create it" >&2
  exit 67
fi

exp065_status="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[\"systemctl is-active vcc-exp065-${EXP065_RUN_ID}.service || true\"]" \
  --comment "Check exact experiment 065 run is no longer active" \
  --query 'Command.CommandId' \
  --output text)"
sleep 3
active="$(aws ssm get-command-invocation \
  --region "$REGION" \
  --command-id "$exp065_status" \
  --instance-id "$INSTANCE_ID" \
  --query 'StandardOutputContent' \
  --output text | tr -d '[:space:]')"
if [[ "$active" == "active" || "$active" == "activating" ]]; then
  echo "experiment 065 run ${EXP065_RUN_ID} is still ${active}" >&2
  exit 68
fi

required_artifacts=(
  "exit-code.txt"
  "artifact-sha256.txt"
  "results.json"
  "build-report.json"
  "baseline-agg_results.csv"
  "candidate-run/agg_results.csv"
  "candidate.h5ad.sha256"
  "score-memory-guard.json"
  "runtime.json"
)
mkdir -p "$EXP065_STAGE" "$HEAD_STAGE"
for relative in "${required_artifacts[@]}"; do
  mkdir -p "$(dirname "${EXP065_STAGE}/${relative}")"
  mkdir -p "$(dirname "${HEAD_STAGE}/${relative}.json")"
  aws s3 cp \
    "s3://${BUCKET}/${EXP065_PREFIX}/${relative}" \
    "${EXP065_STAGE}/${relative}" \
    --only-show-errors
  aws s3api head-object \
    --bucket "$BUCKET" \
    --key "${EXP065_PREFIX}/${relative}" \
    --output json > "${HEAD_STAGE}/${relative}.json"
done

python3 - \
  "$EXP065_STAGE" \
  "$HEAD_STAGE" \
  "$AUTH_STAGE/terminal-s3-manifest.json" \
  "$BUCKET" \
  "$EXP065_PREFIX" \
  "$EXP065_RUN_ID" \
  "${required_artifacts[@]}" <<'PY'
import hashlib
import json
import pathlib
import sys

artifact_root = pathlib.Path(sys.argv[1])
head_root = pathlib.Path(sys.argv[2])
output = pathlib.Path(sys.argv[3])
bucket, prefix, run_id = sys.argv[4:7]
relative_paths = sys.argv[7:]

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()

artifacts = []
for relative in relative_paths:
    path = artifact_root / relative
    head = json.loads((head_root / f"{relative}.json").read_text())
    if int(head["ContentLength"]) != path.stat().st_size:
        raise SystemExit(f"{relative}: S3/local size mismatch")
    artifacts.append({
        "path": relative,
        "key": f"{prefix}/{relative}",
        "size": path.stat().st_size,
        "sha256": sha256(path),
        "etag": head.get("ETag"),
        "last_modified": head.get("LastModified"),
        "expiration": head.get("Expiration"),
        "version_id": head.get("VersionId"),
    })

manifest = {
    "schema": 1,
    "experiment": "065-pm-cd4-hct-count",
    "run_id": run_id,
    "s3": {"bucket": bucket, "prefix": prefix},
    "artifacts": artifacts,
}
output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY

python3 "$HERE/validate_authorization.py" \
  --artifact-root "$EXP065_STAGE" \
  --terminal-manifest "$AUTH_STAGE/terminal-s3-manifest.json" \
  --exp065-run-id "$EXP065_RUN_ID" \
  --expected-prefix "$EXP065_PREFIX" \
  --output "$AUTH_STAGE/authorization.json"

cd "$ROOT"
rm -f "$ARCHIVE"
COPYFILE_DISABLE=1 tar \
  --exclude='*/__pycache__/*' \
  --exclude='*.pyc' \
  --exclude='experiments/071-production-pm-cd4/aws/source.tar.gz' \
  -czf "$ARCHIVE" \
  README.md \
  pyproject.toml \
  uv.lock \
  src \
  scripts/build_submission.py \
  scripts/make_submission.py \
  scripts/make_submission_pm_cd4.py \
  scripts/run_with_memory_guard.py \
  scripts/test_make_submission_pm_cd4.py \
  experiments/061-perturbmatch-hct-count/build_candidate.py \
  experiments/062-cd4-post-whiten/run.py \
  experiments/071-production-pm-cd4

python3 - \
  "$ROOT" \
  "$ARCHIVE" \
  "$HERE/aws/run.sh" \
  "$AUTH_STAGE" \
  "$GATE_ROOT/source-manifest.json" \
  "$RUN_ID" \
  "$SOURCE_PREFIX" \
  "$EXP065_RUN_ID" \
  "$EXP065_PREFIX" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
archive = pathlib.Path(sys.argv[2])
runner = pathlib.Path(sys.argv[3])
authorization = pathlib.Path(sys.argv[4])
output = pathlib.Path(sys.argv[5])
run_id, source_prefix, exp065_run_id, exp065_prefix = sys.argv[6:10]

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()

files = [
    {
        "path": "source.tar.gz",
        "role": "source_archive",
        "sha256": sha256(archive),
        "size": archive.stat().st_size,
    },
    {
        "path": "run.sh",
        "role": "runner",
        "sha256": sha256(runner),
        "size": runner.stat().st_size,
    },
]
for path in sorted(authorization.rglob("*")):
    if path.is_file():
        relative = pathlib.PurePosixPath("authorization") / path.relative_to(
            authorization
        )
        files.append({
            "path": relative.as_posix(),
            "role": "authorization_evidence",
            "sha256": sha256(path),
            "size": path.stat().st_size,
        })

members = []
for relative in (
    "experiments/071-production-pm-cd4/validate_authorization.py",
    "scripts/make_submission_pm_cd4.py",
    "scripts/test_make_submission_pm_cd4.py",
    "scripts/run_with_memory_guard.py",
):
    path = root / relative
    members.append({
        "path": relative,
        "sha256": sha256(path),
        "size": path.stat().st_size,
    })

manifest = {
    "schema": 1,
    "experiment": "071-production-pm-cd4",
    "run_id": run_id,
    "source_prefix": source_prefix,
    "exp065_run_id": exp065_run_id,
    "exp065_prefix": exp065_prefix,
    "challenge_submission_authorized": False,
    "files": files,
    "archive_members": members,
}
output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
PY

SOURCE_ARCHIVE_SHA256="$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')"
RUNNER_SHA256="$(shasum -a 256 "$HERE/aws/run.sh" | awk '{print $1}')"
SOURCE_MANIFEST_SHA256="$(
  shasum -a 256 "$GATE_ROOT/source-manifest.json" | awk '{print $1}'
)"
AUTHORIZATION_SHA256="$(
  shasum -a 256 "$AUTH_STAGE/authorization.json" | awk '{print $1}'
)"

aws s3 cp "$ARCHIVE" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source.tar.gz" --only-show-errors
aws s3 cp "$HERE/aws/run.sh" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/run.sh" --only-show-errors
aws s3 cp "$GATE_ROOT/source-manifest.json" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/source-manifest.json" --only-show-errors
aws s3 cp "$AUTH_STAGE/" \
  "s3://${BUCKET}/${SOURCE_PREFIX}/authorization/" \
  --recursive --only-show-errors

command_id="$(aws ssm send-command \
  --region "$REGION" \
  --instance-ids "$INSTANCE_ID" \
  --document-name AWS-RunShellScript \
  --parameters commands="[
    \"set -Eeuo pipefail\",
    \"aws s3 cp s3://${BUCKET}/${SOURCE_PREFIX}/run.sh /root/exp071-run.sh --only-show-errors\",
    \"echo '${RUNNER_SHA256}  /root/exp071-run.sh' | sha256sum --check -\",
    \"chmod +x /root/exp071-run.sh\",
    \"systemd-run --unit=vcc-exp071-${RUN_ID} --setenv=BUCKET=${BUCKET} --setenv=RUN_ID=${RUN_ID} --setenv=SOURCE_PREFIX=${SOURCE_PREFIX} --setenv=OUTPUT_PREFIX=${OUTPUT_PREFIX} --setenv=EXP065_RUN_ID=${EXP065_RUN_ID} --setenv=EXP065_PREFIX=${EXP065_PREFIX} --setenv=SOURCE_ARCHIVE_SHA256=${SOURCE_ARCHIVE_SHA256} --setenv=RUNNER_SHA256=${RUNNER_SHA256} --setenv=SOURCE_MANIFEST_SHA256=${SOURCE_MANIFEST_SHA256} --setenv=AUTHORIZATION_SHA256=${AUTHORIZATION_SHA256} /bin/bash /root/exp071-run.sh\"
  ]" \
  --comment "Build exploratory PM+CD4 package only; never submit" \
  --query 'Command.CommandId' \
  --output text)"

cat > "${HERE}/launch-metadata.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "instance_id": "${INSTANCE_ID}",
  "exp065_run_id": "${EXP065_RUN_ID}",
  "exp065_s3_prefix": "s3://${BUCKET}/${EXP065_PREFIX}/",
  "source_s3_prefix": "s3://${BUCKET}/${SOURCE_PREFIX}/",
  "output_s3_prefix": "s3://${BUCKET}/${OUTPUT_PREFIX}/",
  "source_archive_sha256": "${SOURCE_ARCHIVE_SHA256}",
  "runner_sha256": "${RUNNER_SHA256}",
  "source_manifest_sha256": "${SOURCE_MANIFEST_SHA256}",
  "authorization_sha256": "${AUTHORIZATION_SHA256}",
  "ssm_command_id": "${command_id}",
  "authorization_scope": "exploratory_packaging_only",
  "challenge_submission_authorized": false,
  "vcc_prep": true,
  "submission": false
}
EOF

printf 'run_id=%s\ncommand_id=%s\noutput=s3://%s/%s/\n' \
  "$RUN_ID" "$command_id" "$BUCKET" "$OUTPUT_PREFIX"
