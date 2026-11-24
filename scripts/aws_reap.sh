#!/usr/bin/env bash
# aws_reap.sh — find and stop idle EC2 workers this project left running.
#
# WHY THIS EXISTS
# ---------------
# Experiment 061 finished on 2026-08-22 and its worker (`i-0af19b9a537e08d4b`, `r7i.2xlarge`,
# us-west-2) stayed RUNNING for 13 days at 0.27% average CPU. Nothing noticed, because nothing
# was looking. Two independent guards now exist and this is the outer one:
#
#   1. every launcher registers a boot-time `shutdown -h +N` and sets
#      `instance-initiated-shutdown-behavior=terminate`, so a worker dies on its own even if
#      the coordinator crashes  (`scripts/cloud_score.sh`);
#   2. this script sweeps EVERY enabled region for anything that outlived its guard.
#
# WHAT IT DOES
#   - lists every RUNNING instance in every enabled region, with age, type and `Name` tag;
#   - pulls 24 h of CloudWatch `CPUUtilization` (1 h periods) for each;
#   - marks REAP when age > --max-age-hours AND the instance looks idle on BOTH readings:
#     max CPU% under --cpu-max, and max BUSY CORES (cpu% x vCPUs / 100) under --busy-cores.
#
# WHY TWO READINGS. CloudWatch reports CPU as a PERCENTAGE OF THE WHOLE INSTANCE, and
# `cell-eval2` scores single-threaded: measured on this project's own workers, a fully busy
# scorer reads 13.7% on an 8-vCPU r7i.2xlarge — and would read 3.4% on a 32-vCPU one, UNDER a
# 5% threshold. A percentage-only rule therefore reaps a working scorer as soon as the
# instance is big enough. The cores reading is instance-size independent: 0.152 busy cores for
# the 13-day orphan this script exists for, 1.10 for a busy 8-vCPU scorer, 1.06 for the same
# job on 32 vCPUs. Both must say idle before anything is stopped.
#   - DRY RUN by default. `--apply` actually stops them. It never terminates: an EBS-backed
#     stop is reversible and the root volume survives, so a false positive costs nothing.
#
# An instance younger than the age threshold is never reaped, whatever its CPU — a worker that
# has just booted and is still pulling 60 GB from S3 reads as 0% CPU.
#
# USAGE
#   scripts/aws_reap.sh                          # dry run, all regions, 24 h / 5 %
#   scripts/aws_reap.sh --apply                  # actually stop the flagged ones
#   scripts/aws_reap.sh --max-age-hours 6 --cpu-max 2 --apply
#   scripts/aws_reap.sh --regions us-east-1,us-west-2 --json work/reap.json
#   scripts/aws_reap.sh --s3 <your-s3-bucket>   # also report bucket state
set -Eeuo pipefail

APPLY=0
MAX_AGE_HOURS=24
CPU_MAX=5.0
BUSY_CORES=0.25
REGIONS=""
JSON_OUT=""
S3_BUCKET=""

usage() { sed -n '2,40p' "$0"; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply) APPLY=1; shift ;;
    --max-age-hours) MAX_AGE_HOURS="$2"; shift 2 ;;
    --cpu-max) CPU_MAX="$2"; shift 2 ;;
    --busy-cores) BUSY_CORES="$2"; shift 2 ;;
    --regions) REGIONS="$2"; shift 2 ;;
    --json) JSON_OUT="$2"; shift 2 ;;
    --s3) S3_BUCKET="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 64 ;;
  esac
done

for c in aws python3; do
  command -v "$c" >/dev/null 2>&1 || { echo "required command missing: $c" >&2; exit 63; }
done

if [[ -z "$REGIONS" ]]; then
  REGIONS="$(aws ec2 describe-regions --query 'Regions[].RegionName' --output text | tr '\t' ',')"
fi

NOW_UTC="$(python3 -c 'import datetime;print(datetime.datetime.now(datetime.timezone.utc).isoformat())')"
START_UTC="$(python3 -c "
import datetime
print((datetime.datetime.now(datetime.timezone.utc)-datetime.timedelta(hours=24)).strftime('%Y-%m-%dT%H:%M:%SZ'))")"
END_UTC="$(python3 -c "
import datetime
print(datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'))")"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
: >"$TMP/rows.tsv"

printf '%-14s  %-21s  %-13s  %7s  %7s  %6s  %-30s  %s\n' \
  REGION INSTANCE TYPE AGE_H CPU_PCT CORES NAME VERDICT

IFS=',' read -r -a REGION_LIST <<<"$REGIONS"
for region in "${REGION_LIST[@]}"; do
  [[ -n "$region" ]] || continue
  # `|| true`: an opted-in-but-unreachable region must not abort the sweep.
  insts="$(aws ec2 describe-instances --region "$region" \
      --filters Name=instance-state-name,Values=running \
      --query 'Reservations[].Instances[].[InstanceId,InstanceType,LaunchTime,(Tags[?Key==`Name`].Value|[0]),CpuOptions.CoreCount,CpuOptions.ThreadsPerCore]' \
      --output text 2>/dev/null || true)"
  [[ -n "$insts" ]] || continue
  while IFS=$'\t' read -r iid itype launch name cores threads; do
    [[ -n "${iid:-}" ]] || continue
    [[ "$name" == "None" ]] && name="(untagged)"
    vcpus=1
    [[ "$cores" =~ ^[0-9]+$ && "$threads" =~ ^[0-9]+$ ]] && vcpus=$(( cores * threads ))
    age_h="$(python3 -c "
import datetime,sys
launch=datetime.datetime.fromisoformat(sys.argv[1].replace('Z','+00:00'))
now=datetime.datetime.fromisoformat(sys.argv[2])
print(f'{(now-launch).total_seconds()/3600:.1f}')" "$launch" "$NOW_UTC")"
    cpu="$(aws cloudwatch get-metric-statistics --region "$region" \
        --namespace AWS/EC2 --metric-name CPUUtilization \
        --dimensions "Name=InstanceId,Value=${iid}" \
        --start-time "$START_UTC" --end-time "$END_UTC" \
        --period 3600 --statistics Maximum \
        --query 'max(Datapoints[].Maximum)' --output text 2>/dev/null || echo None)"
    [[ -z "$cpu" || "$cpu" == "None" || "$cpu" == "null" ]] && cpu="nan"
    verdict="$(python3 -c "
import sys
age, cpu, max_age = float(sys.argv[1]), sys.argv[2], float(sys.argv[3])
cpu_max, vcpus, busy_cores = float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
if age <= max_age:
    print('keep (young)')
elif cpu == 'nan':
    print('keep (no CPU data)')
elif float(cpu) < cpu_max and float(cpu) * vcpus / 100.0 < busy_cores:
    print('REAP')
else:
    print('keep (busy)')" "$age_h" "$cpu" "$MAX_AGE_HOURS" "$CPU_MAX" "$vcpus" "$BUSY_CORES")"
    busy="$(python3 -c "
import sys
print('nan' if sys.argv[1]=='nan' else f'{float(sys.argv[1])*int(sys.argv[2])/100.0:.3f}')" "$cpu" "$vcpus")"
    printf '%-14s  %-21s  %-13s  %7s  %7s  %6s  %-30s  %s\n' \
      "$region" "$iid" "$itype" "$age_h" "${cpu:0:7}" "$busy" "${name:0:30}" "$verdict"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$region" "$iid" "$itype" "$age_h" "$cpu" "$name" "$verdict" "$busy" >>"$TMP/rows.tsv"
  done <<<"$insts"
done

n_total="$(wc -l <"$TMP/rows.tsv" | tr -d ' ')"
n_reap="$(grep -c $'\tREAP$' "$TMP/rows.tsv" || true)"
echo
echo "running instances: ${n_total}   flagged for reaping: ${n_reap}   " \
     "(threshold: age > ${MAX_AGE_HOURS} h AND 24 h max CPU < ${CPU_MAX} % AND < ${BUSY_CORES} busy cores)"

if [[ "$n_reap" -gt 0 ]]; then
  if [[ "$APPLY" -eq 1 ]]; then
    while IFS=$'\t' read -r region iid _ _ _ _ verdict _; do
      [[ "$verdict" == "REAP" ]] || continue
      echo "stopping ${iid} in ${region}"
      aws ec2 stop-instances --region "$region" --instance-ids "$iid" \
        --query 'StoppingInstances[].{id:InstanceId,from:PreviousState.Name,to:CurrentState.Name}' \
        --output text
    done <"$TMP/rows.tsv"
  else
    echo "DRY RUN — re-run with --apply to stop them."
  fi
fi

if [[ -n "$S3_BUCKET" ]]; then
  echo
  echo "=== s3://${S3_BUCKET} ==="
  aws s3 ls "s3://${S3_BUCKET}/" --recursive --summarize | tail -2
fi

if [[ -n "$JSON_OUT" ]]; then
  python3 - "$TMP/rows.tsv" "$JSON_OUT" "$NOW_UTC" "$MAX_AGE_HOURS" "$CPU_MAX" <<'PY'
import json, sys
rows = []
with open(sys.argv[1]) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) != 8:
            continue
        rows.append(dict(region=f[0], instance_id=f[1], instance_type=f[2],
                         age_hours=float(f[3]),
                         cpu_max_24h=None if f[4] == "nan" else float(f[4]),
                         name=f[5], verdict=f[6],
                         busy_cores_24h=None if f[7] == "nan" else float(f[7])))
with open(sys.argv[2], "w") as fh:
    json.dump({"swept_utc": sys.argv[3], "max_age_hours": float(sys.argv[4]),
               "cpu_max_pct": float(sys.argv[5]),
               "running": len(rows),
               "flagged": sum(r["verdict"] == "REAP" for r in rows),
               "instances": rows}, fh, indent=1)
PY
  echo "wrote ${JSON_OUT}"
fi
