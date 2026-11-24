#!/usr/bin/env bash
# cloud_score.sh — score `.h5ad` candidates through the REAL `cell-eval2` on a disposable
# EC2 worker, then terminate it.
#
# WHY THIS EXISTS
# ---------------
# `AGENTS.md` §6.4's memory-and-validation note has said since 2026-08-22: *"never run more
# than one `cell-eval2` process locally. A scorer can consume about 25 GB RSS … Use
# `scripts/run_with_memory_guard.py`, **or score on an isolated cloud worker.**"* No worker was
# ever built, so in practice the full six-member scorer was almost never run and members were
# priced piecemeal through proxies. §8.4 (now docs/FINDINGS.md §8.4) measured what that cost: `cd4f034` was priced at
# +0.00437 and delivered +0.0130557, with **60% of the realised gain coming from the four
# members carried at zero**. This is the tool that removes the excuse.
#
# WHAT IT DOES
#   1. uploads the reference and every candidate to S3 (or reuses S3 URIs already there);
#   2. launches ONE right-sized instance — spot by default, on-demand fallback;
#   3. runs the full pipeline: `baseline` -> `prep-real-bundle` -> `run` -> `score`,
#      `--preset vcc2026`, on every candidate against one shared bundle;
#   4. pulls back `agg_results.csv` / `score.csv` with sha256s;
#   5. TERMINATES the instance and proves it terminated.
#
# SELF-TERMINATION, three independent guards (see `cloud_score_worker.sh`):
#   the instance is launched `instance-initiated-shutdown-behavior=terminate`; user-data arms
#   `shutdown -h +$MAX_MINUTES` before any work starts; and the worker's EXIT trap powers off.
#   A crashed coordinator therefore cannot leave a running instance behind. `scripts/aws_reap.sh`
#   is the outer sweep for anything that still escapes.
#
# USAGE
#   VCC_CLOUD_SCORE_CLEARANCE=YES scripts/cloud_score.sh \
#       --real data/harness/h1-test.h5ad \
#       --pred a15=experiments/xxx/a15.h5ad \
#       --pred rc_hybrid=s3://bucket/key.h5ad \
#       --out experiments/221-cloud-scorer/cloud-output --wait
#
#   # with a build step on the worker (heavy staging stays off the laptop):
#   ... --prologue experiments/221-cloud-scorer/aws/build_h1_candidates.sh \
#       --source-include src --source-include scripts \
#       --payload experiments/170-gamma-mse/work/pre.npz \
#       --pred-from-prologue a15 --pred-from-prologue cd4f034
#
#   scripts/cloud_score.sh --collect <run-id> --out <dir>     # resume a collection
set -Eeuo pipefail

REGION="${AWS_REGION:-us-east-1}"
BUCKET="${BUCKET:?set BUCKET to your S3 bucket}"
INSTANCE_TYPE="${INSTANCE_TYPE:-r7i.2xlarge}"
IAM_PROFILE="${IAM_PROFILE:-vcc-scoring-20260822}"
CE_VERSION="${CE_VERSION:-0.16.0}"
PRESET="vcc2026"
PERT_COL="target_gene"
CONTROL="non-targeting"
ROOT_GB=200
MAX_MINUTES=240
USE_SPOT=1
WAIT=0
DRY_RUN=0
RUN_ID=""
COLLECT=""
OUT_DIR=""
REAL=""
BUNDLE_URI=""
BUNDLE_WAIT=0
BASELINE_PRED=""
BUNDLE_ONLY=0
PROLOGUE=""
declare -a PRED_SPECS=() PROLOGUE_PREDS=() PAYLOADS=() SOURCE_INCLUDES=()

usage() { sed -n '2,46p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --real) REAL="$2"; shift 2 ;;
    --bundle) BUNDLE_URI="$2"; shift 2 ;;
    --baseline-pred) BASELINE_PRED="$2"; shift 2 ;;
    --bundle-only) BUNDLE_ONLY=1; shift ;;
    --bundle-wait-minutes) BUNDLE_WAIT="$2"; shift 2 ;;
    --pred) PRED_SPECS+=("$2"); shift 2 ;;
    --pred-from-prologue) PROLOGUE_PREDS+=("$2"); shift 2 ;;
    --prologue) PROLOGUE="$2"; shift 2 ;;
    --payload) PAYLOADS+=("$2"); shift 2 ;;
    --source-include) SOURCE_INCLUDES+=("$2"); shift 2 ;;
    --out) OUT_DIR="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --collect) COLLECT="$2"; shift 2 ;;
    --instance-type) INSTANCE_TYPE="$2"; shift 2 ;;
    --region) REGION="$2"; shift 2 ;;
    --root-gb) ROOT_GB="$2"; shift 2 ;;
    --max-minutes) MAX_MINUTES="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    --pert-col) PERT_COL="$2"; shift 2 ;;
    --control) CONTROL="$2"; shift 2 ;;
    --no-spot) USE_SPOT=0; shift ;;
    --wait) WAIT=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 64 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for c in aws python3 jq shasum tar; do
  command -v "$c" >/dev/null 2>&1 || { echo "required command missing: $c" >&2; exit 63; }
done
sha256_of() { shasum -a 256 "$1" | awk '{print $1}'; }

[[ -n "$OUT_DIR" ]] || { echo "--out is required" >&2; exit 64; }
mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------- collect helper
collect() {
  local run_id="$1" out_prefix="output/cloud-score/${run_id}"
  echo "[collect] s3://${BUCKET}/${out_prefix}/ -> ${OUT_DIR}/${run_id}/"
  mkdir -p "${OUT_DIR}/${run_id}"
  aws s3 sync "s3://${BUCKET}/${out_prefix}/" "${OUT_DIR}/${run_id}/" \
    --region "$REGION" --exclude '*.h5ad' --only-show-errors
  echo "[collect] status: $(cat "${OUT_DIR}/${run_id}/STATUS" 2>/dev/null || echo 'no STATUS')"
}

if [[ -n "$COLLECT" ]]; then collect "$COLLECT"; exit 0; fi

# ---------------------------------------------------------------------------- clearance gate
if [[ "${VCC_CLOUD_SCORE_CLEARANCE:-}" != "YES" && "$DRY_RUN" -eq 0 ]]; then
  echo "refusing launch: set VCC_CLOUD_SCORE_CLEARANCE=YES (this creates a billable instance)" >&2
  exit 65
fi
[[ -n "$REAL" ]] || { echo "--real is required" >&2; exit 64; }
[[ ${#PRED_SPECS[@]} -gt 0 || ${#PROLOGUE_PREDS[@]} -gt 0 || "$BUNDLE_ONLY" -eq 1 ]] || {
  echo "at least one --pred or --pred-from-prologue is required (or --bundle-only)" >&2; exit 64; }

RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
[[ "$RUN_ID" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "unsafe run id: $RUN_ID" >&2; exit 69; }
PREFIX="input/cloud-score/${RUN_ID}"
OUT_PREFIX="output/cloud-score/${RUN_ID}"
NAME_TAG="vcc-cloud-score-${RUN_ID}"

# A reused prefix would silently mix two runs' artifacts. 082's launcher refused that; so do we.
for p in "$PREFIX" "$OUT_PREFIX"; do
  n="$(aws s3api list-objects-v2 --region "$REGION" --bucket "$BUCKET" \
        --prefix "${p}/" --max-keys 1 --query 'KeyCount' --output text)"
  [[ "$n" == "0" ]] || { echo "refusing reused S3 prefix: s3://${BUCKET}/${p}/" >&2; exit 79; }
done

# ---------------------------------------------------------------------------- upload inputs
put() {  # local-or-s3 -> s3 uri the WORKER can read
  local src="$1" key="$2"
  # The instance role may only GetObject under `input/*`. An `s3://` source already there is
  # passed through untouched (that is how a re-used harness avoids a seventh upload); one
  # anywhere else — a previous run's parked `output/.../built-preds/` — is server-side copied
  # into this run's input prefix, which costs no egress and no laptop bandwidth.
  if [[ "$src" == "s3://${BUCKET}/input/"* ]]; then echo "$src"; return; fi
  if [[ "$src" == s3://* ]]; then
    echo "[stage] $src -> s3://${BUCKET}/${key} (server-side copy)" >&2
    [[ "$DRY_RUN" -eq 1 ]] || aws s3 cp "$src" "s3://${BUCKET}/${key}" --region "$REGION" --only-show-errors >&2
    echo "s3://${BUCKET}/${key}"; return
  fi
  [[ -f "$src" ]] || { echo "no such file: $src" >&2; exit 66; }
  echo "[upload] $src -> s3://${BUCKET}/${key} ($(du -h "$src" | cut -f1))" >&2
  [[ "$DRY_RUN" -eq 1 ]] || aws s3 cp "$src" "s3://${BUCKET}/${key}" --region "$REGION" --only-show-errors >&2
  echo "s3://${BUCKET}/${key}"
}

REAL_URI="$(put "$REAL" "${PREFIX}/real.h5ad")"

PRED_JSON="[]"
for spec in "${PRED_SPECS[@]+"${PRED_SPECS[@]}"}"; do
  name="${spec%%=*}"; path="${spec#*=}"
  [[ "$name" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "unsafe pred name: $name" >&2; exit 69; }
  uri="$(put "$path" "${PREFIX}/preds/${name}.h5ad")"
  PRED_JSON="$(jq -c --arg n "$name" --arg u "$uri" '. + [{name:$n, uri:$u}]' <<<"$PRED_JSON")"
done
for name in "${PROLOGUE_PREDS[@]+"${PROLOGUE_PREDS[@]}"}"; do
  [[ "$name" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "unsafe pred name: $name" >&2; exit 69; }
  PRED_JSON="$(jq -c --arg n "$name" '. + [{name:$n, uri:"PROLOGUE"}]' <<<"$PRED_JSON")"
done

BASELINE_PRED_URI=""
if [[ -n "$BASELINE_PRED" ]]; then
  BASELINE_PRED_URI="$(put "$BASELINE_PRED" "${PREFIX}/baseline-pred.h5ad")"
fi

PROLOGUE_NAME=""
if [[ -n "$PROLOGUE" ]]; then
  [[ -f "$PROLOGUE" ]] || { echo "no such prologue: $PROLOGUE" >&2; exit 66; }
  PROLOGUE_NAME="$(basename "$PROLOGUE")"
  [[ ${#SOURCE_INCLUDES[@]} -gt 0 ]] || SOURCE_INCLUDES=(src scripts pyproject.toml uv.lock)
  ARCHIVE="${OUT_DIR}/source-${RUN_ID}.tar.gz"
  COPYFILE_DISABLE=1 COPY_EXTENDED_ATTRIBUTES_DISABLE=1 tar \
    --exclude='._*' --exclude='.DS_Store' --exclude='__pycache__' --exclude='*.pyc' \
    -czf "$ARCHIVE" -C "$ROOT" "${SOURCE_INCLUDES[@]}"
  if tar -tzf "$ARCHIVE" | grep -Eq '(^|/)\._|(^|/)\.DS_Store$'; then
    echo "refusing upload: source archive carries macOS metadata" >&2; exit 67
  fi
  put "$ARCHIVE" "${PREFIX}/source.tar.gz" >/dev/null
  put "$PROLOGUE" "${PREFIX}/prologue.sh" >/dev/null
  for pay in "${PAYLOADS[@]+"${PAYLOADS[@]}"}"; do
    rel="${pay#"$ROOT"/}"; rel="${rel#./}"
    put "$pay" "${PREFIX}/payload/${rel}" >/dev/null
  done
fi

jq -n \
  --arg real "$REAL_URI" --arg preset "$PRESET" --arg pc "$PERT_COL" --arg ctl "$CONTROL" \
  --arg cev "$CE_VERSION" --arg pro "$PROLOGUE_NAME" --arg bun "$BUNDLE_URI" \
  --argjson bw "${BUNDLE_WAIT:-0}" --arg bp "$BASELINE_PRED_URI" \
  --argjson preds "$PRED_JSON" \
  '{real:$real, preds:$preds, preset:$preset, pert_col:$pc, control:$ctl,
    cell_eval2_version:$cev, prologue:$pro, bundle:$bun,
    bundle_wait_minutes:$bw, baseline_pred:$bp}' >"${OUT_DIR}/job-${RUN_ID}.json"
cat "${OUT_DIR}/job-${RUN_ID}.json"
put "${OUT_DIR}/job-${RUN_ID}.json" "${PREFIX}/job.json" >/dev/null

WORKER="${ROOT}/scripts/cloud_score_worker.sh"
WORKER_SHA="$(sha256_of "$WORKER")"
put "$WORKER" "${PREFIX}/worker.sh" >/dev/null

# ---------------------------------------------------------------------------- launch
AMI="$(aws ssm get-parameter --region "$REGION" \
  --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
  --query 'Parameter.Value' --output text)"
SPOT_PRICE="$(aws ec2 describe-spot-price-history --region "$REGION" \
  --instance-types "$INSTANCE_TYPE" --product-descriptions 'Linux/UNIX' \
  --start-time "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  --query 'SpotPriceHistory[0].SpotPrice' --output text 2>/dev/null | head -1 || echo unknown)"

USER_DATA="$(cat <<EOF
#!/bin/bash
# GUARD (b): arm the dead-man switch BEFORE anything can hang. Combined with
# instance-initiated-shutdown-behavior=terminate, this instance cannot outlive ${MAX_MINUTES} min.
shutdown -h +${MAX_MINUTES}
set -Eeuo pipefail
export AWS_REGION=${REGION}
mkdir -p /work
aws s3 cp s3://${BUCKET}/${PREFIX}/worker.sh /root/worker.sh --only-show-errors
printf '${WORKER_SHA}  /root/worker.sh\n' | sha256sum --check -
chmod +x /root/worker.sh
systemd-run --unit=vcc-cloud-score --collect \
  --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} \
  --setenv=PREFIX=${PREFIX} --setenv=OUT=${OUT_PREFIX} \
  /bin/bash /root/worker.sh
EOF
)"

cat >"${OUT_DIR}/launch-${RUN_ID}.json" <<EOF
{
  "run_id": "${RUN_ID}",
  "region": "${REGION}",
  "bucket": "${BUCKET}",
  "instance_type": "${INSTANCE_TYPE}",
  "ami": "${AMI}",
  "spot": ${USE_SPOT},
  "spot_price_at_launch": "${SPOT_PRICE}",
  "root_gb": ${ROOT_GB},
  "max_minutes": ${MAX_MINUTES},
  "worker_sha256": "${WORKER_SHA}",
  "input_prefix": "s3://${BUCKET}/${PREFIX}/",
  "output_prefix": "s3://${BUCKET}/${OUT_PREFIX}/",
  "shutdown_behavior": "terminate",
  "submission": false,
  "vcc_submit": false
}
EOF

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "[dry-run] would launch ${INSTANCE_TYPE} (${AMI}) in ${REGION}; spot=${USE_SPOT}"
  cat "${OUT_DIR}/launch-${RUN_ID}.json"; exit 0
fi

run_instances() {  # $1 = "spot" | "ondemand"
  local -a market=()
  [[ "$1" == "spot" ]] && market=(--instance-market-options \
    'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}')
  aws ec2 run-instances --region "$REGION" \
    --image-id "$AMI" --instance-type "$INSTANCE_TYPE" --count 1 \
    --iam-instance-profile "Name=${IAM_PROFILE}" \
    --instance-initiated-shutdown-behavior terminate \
    --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":${ROOT_GB},\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]" \
    --metadata-options "HttpTokens=required,HttpEndpoint=enabled" \
    --tag-specifications \
      "ResourceType=instance,Tags=[{Key=Name,Value=${NAME_TAG}},{Key=project,Value=vcc2026},{Key=purpose,Value=cell-eval2-scoring},{Key=reap-after-minutes,Value=${MAX_MINUTES}}]" \
    --user-data "$USER_DATA" \
    "${market[@]+"${market[@]}"}" \
    --query 'Instances[0].InstanceId' --output text
}

MARKET="spot"
if [[ "$USE_SPOT" -eq 1 ]]; then
  INSTANCE_ID="$(run_instances spot 2>"${OUT_DIR}/spot-error-${RUN_ID}.txt" || true)"
else
  INSTANCE_ID=""
fi
if [[ -z "${INSTANCE_ID:-}" || "$INSTANCE_ID" == "None" ]]; then
  [[ "$USE_SPOT" -eq 1 ]] && echo "[launch] spot unavailable, falling back to on-demand:" \
    && cat "${OUT_DIR}/spot-error-${RUN_ID}.txt" >&2
  MARKET="ondemand"
  INSTANCE_ID="$(run_instances ondemand)"
fi
echo "[launch] ${INSTANCE_ID} (${INSTANCE_TYPE}, ${MARKET}) in ${REGION}; name ${NAME_TAG}"
python3 - "${OUT_DIR}/launch-${RUN_ID}.json" "$INSTANCE_ID" "$MARKET" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["instance_id"], d["market"] = sys.argv[2], sys.argv[3]
json.dump(d, open(p, "w"), indent=1)
PY

echo "run_id=${RUN_ID}"
echo "instance_id=${INSTANCE_ID}"
echo "output=s3://${BUCKET}/${OUT_PREFIX}/"

[[ "$WAIT" -eq 1 ]] || { echo "not waiting; collect later with --collect ${RUN_ID}"; exit 0; }

# ---------------------------------------------------------------------------- wait + collect
deadline=$(( $(date -u +%s) + MAX_MINUTES * 60 + 900 ))
status=""
while :; do
  if aws s3 cp "s3://${BUCKET}/${OUT_PREFIX}/STATUS" - --region "$REGION" >"${OUT_DIR}/.status" 2>/dev/null; then
    status="$(tr -d '\n' <"${OUT_DIR}/.status")"
    [[ -n "$status" ]] && { echo "[wait] worker reported: ${status}"; break; }
  fi
  state="$(aws ec2 describe-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
      --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null || echo gone)"
  if [[ "$state" == "terminated" || "$state" == "shutting-down" || "$state" == "gone" ]]; then
    echo "[wait] instance is ${state} with no STATUS marker — collecting what exists"; break
  fi
  [[ "$(date -u +%s)" -lt "$deadline" ]] || { echo "[wait] deadline passed"; break; }
  sleep 60
done

collect "$RUN_ID"

echo "[terminate] ${INSTANCE_ID}"
aws ec2 terminate-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
  --query 'TerminatingInstances[].{id:InstanceId,to:CurrentState.Name}' --output text || true
aws ec2 wait instance-terminated --region "$REGION" --instance-ids "$INSTANCE_ID" || true
aws ec2 describe-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
  --query 'Reservations[0].Instances[0].{id:InstanceId,state:State.Name}' --output json \
  | tee "${OUT_DIR}/terminated-${RUN_ID}.json"
[[ "$status" == DONE* ]] || { echo "worker did not report DONE" >&2; exit 75; }
