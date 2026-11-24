#!/usr/bin/env bash
# labcloud.sh — the lab's EC2 runner (experiments 8NN): one instance per batch, code tarball + S3 data mirror + a job list.
#   stage <relpath>... | launch --name N --jobs FILE [--type T] [--minutes M] [--parallel P] [--spot] [--root-gb G]
#     [--sync RELPATH]... [--collect RELPATH]... [--no-venv] [--code-extra RELPATH]... | wait <run> | collect <run> | status | reap
# Guards (charter + Amendment 1): Project=lab tags (instance + volume), terminate-on-shutdown, user-data `shutdown -h +M`,
# worker EXIT trap powers off; caps $50/launch, $150/UTC day, $800 total; no Mac stage > 1 GB in a reserved window
# (window800.sh). Every launch/termination goes to cloud/runs.jsonl and docs/LAB-LOG.md. Full notes: DESIGN.md section 11.
# Runs from a frozen per-content copy (cloud/.frozen/), so editing this file never disturbs a live `wait` (bash reads
# scripts incrementally; AGENTS 6.2).
if [[ -z "${LABCLOUD_FROZEN:-}" ]]; then
  _src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
  _f="$(dirname "$_src")/.frozen/labcloud-$(shasum -a 256 "$_src" | cut -c1-12).sh"
  mkdir -p "$(dirname "$_f")"; [[ -f "$_f" ]] || cp "$_src" "$_f"
  LABCLOUD_FROZEN=1 LABCLOUD_ROOT="$(cd "$(dirname "$_src")/../../.." && pwd)" exec bash "$_f" "$@"
fi
set -Eeuo pipefail
REGION="${AWS_REGION:-us-east-1}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
IAM_PROFILE="${IAM_PROFILE:-vcc-scoring-20260822}"
ROOT="${LABCLOUD_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
HERE="$ROOT/experiments/800-checker/cloud"
RUNS="$HERE/runs.jsonl"
LABLOG="$ROOT/docs/LAB-LOG.md"
BUDGET="${LAB_BUDGET_USD:-800}"
LAUNCH_CAP="${LAB_LAUNCH_CAP_USD:-50}"
DAY_CAP="${LAB_DAY_CAP_USD:-150}"
mkdir -p "$HERE/runs"
utc() { date -u +%FT%TZ; }
die() { echo "labcloud: $*" >&2; exit 2; }

price_of() {  # on-demand USD/h for an instance type, from the Pricing API (fallback table)
  local t="$1" p
  p="$(aws pricing get-products --region us-east-1 --service-code AmazonEC2 \
      --filters Type=TERM_MATCH,Field=instanceType,Value="$t" Type=TERM_MATCH,Field=location,Value="US East (N. Virginia)" \
      Type=TERM_MATCH,Field=operatingSystem,Value=Linux Type=TERM_MATCH,Field=tenancy,Value=Shared \
      Type=TERM_MATCH,Field=preInstalledSw,Value=NA Type=TERM_MATCH,Field=capacitystatus,Value=Used --max-results 1 --output json 2>/dev/null \
      | python3 -c 'import json,sys
d=json.load(sys.stdin)
for p in d.get("PriceList",[]):
    p=json.loads(p)
    for t in p["terms"]["OnDemand"].values():
        for x in t["priceDimensions"].values(): print(x["pricePerUnit"]["USD"]); raise SystemExit' 2>/dev/null || true)"
  [[ -n "$p" ]] && { echo "$p"; return; }
  case "$t" in r7i.large) echo 0.1323;; r7i.xlarge) echo 0.2646;; r7i.2xlarge) echo 0.5292;; r7i.4xlarge) echo 1.0584;;
    r7i.8xlarge) echo 2.1168;; r7i.12xlarge) echo 3.1752;; r7i.16xlarge) echo 4.2336;; c7i.4xlarge) echo 0.714;;
    c7i.8xlarge) echo 1.428;; m7i.4xlarge) echo 0.8064;; *) echo 5.0;; esac
}

spent() {  # USD logged so far (actual where known, else the launch cap)
  [[ -f "$RUNS" ]] || { echo 0; return; }
  python3 - "$RUNS" <<'PY'
import json, sys
tot = {}
for ln in open(sys.argv[1]):
    if not ln.strip(): continue
    r = json.loads(ln)
    rid = r["run_id"]
    if r.get("event") == "launch": tot.setdefault(rid, r["cap_usd"])
    if r.get("event") == "terminated": tot[rid] = r["actual_usd"]
print(round(sum(tot.values()), 2))
PY
}

spent_today() {  # USD of launches started this UTC day (actual where known, else the cap)
  [[ -f "$RUNS" ]] || { echo 0; return; }
  python3 - "$RUNS" "$(date -u +%Y-%m-%d)" <<'PY'
import json, sys
tot, day = {}, sys.argv[2]
for ln in open(sys.argv[1]):
    if not ln.strip(): continue
    r = json.loads(ln)
    if r.get("event") == "launch" and r["utc"].startswith(day): tot[r["run_id"]] = r["cap_usd"]
    if r.get("event") == "terminated" and r["run_id"] in tot: tot[r["run_id"]] = r["actual_usd"]
print(round(sum(tot.values()), 2))
PY
}

lablog_cost() {  # append one row to docs/LAB-LOG.md's cost table (the table is the file's last section)
  printf '| %s | %s | %s | %s | %s | %s |\n' "$@" >> "$LABLOG"
}

cmd_stage() {
  [[ $# -gt 0 ]] || die "stage: give repo-relative paths"
  # reserved windows (window800.sh): no transfer > 1 GB from the Mac
  if ! bash "$ROOT/experiments/800-checker/window800.sh" >/dev/null; then
    local kb=0; for rel in "$@"; do kb=$(( kb + $(du -sk "$ROOT/${rel%/}" | cut -f1) )); done
    (( kb <= 1048576 )) || die "reserved window ($(bash "$ROOT/experiments/800-checker/window800.sh")): ${kb} KB > 1 GB from the Mac refused"
  fi
  for rel in "$@"; do
    rel="${rel%/}"; local p="$ROOT/$rel"
    [[ -e "$p" ]] || die "no such path: $rel"
    if [[ -d "$p" ]]; then
      aws s3 sync "$p/" "s3://${BUCKET}/input/lab/repo/${rel}/" --region "$REGION" --only-show-errors --exclude '*.pyc' --exclude '__pycache__/*' --exclude '.DS_Store'
    else
      aws s3 cp "$p" "s3://${BUCKET}/input/lab/repo/${rel}" --region "$REGION" --only-show-errors
    fi
    printf '%s\t%s\t%s\n' "$(utc)" "$rel" "$(du -sk "$p" | cut -f1)" >> "$HERE/staged.tsv"
    echo "staged $rel ($(du -sh "$p" | cut -f1))"
  done
}

cmd_launch() {
  local NAME="" JOBS="" TYPE="r7i.4xlarge" MIN=120 PAR=4 SPOT=0 ROOTGB=200 VENV=true
  local -a SYNC=() COLLECT=() EXTRA=()
  while [[ $# -gt 0 ]]; do case "$1" in
    --name) NAME="$2"; shift 2;; --jobs) JOBS="$2"; shift 2;; --type) TYPE="$2"; shift 2;; --minutes) MIN="$2"; shift 2;;
    --parallel) PAR="$2"; shift 2;; --spot) SPOT=1; shift;; --root-gb) ROOTGB="$2"; shift 2;; --no-venv) VENV=false; shift;;
    --sync) SYNC+=("${2%/}"); shift 2;; --collect) COLLECT+=("$2"); shift 2;; --code-extra) EXTRA+=("$2"); shift 2;;
    *) die "launch: unknown argument $1";; esac; done
  [[ "$NAME" =~ ^[a-z0-9-]+$ ]] || die "--name must be [a-z0-9-]+"
  [[ -f "$JOBS" ]] || die "--jobs file missing"
  [[ "$SPOT" -eq 0 || "$MIN" -le 30 ]] || die "spot only for jobs capped at <= 30 min (charter)"
  local RID; RID="$(date -u +%Y%m%dT%H%M%SZ)-${NAME}"
  local PREFIX="input/lab/runs/${RID}" OUT="output/lab/${RID}" D="$HERE/runs/${RID}"
  mkdir -p "$D"
  local PRICE CAP SP
  PRICE="$(price_of "$TYPE")"
  CAP="$(python3 -c "print(round(${PRICE}*(${MIN}+10)/60 + ${ROOTGB}*0.08/730*(${MIN}+10)/60, 2))")"
  SP="$(spent)"
  python3 -c "import sys; sys.exit(0 if ${SP}+${CAP} <= ${BUDGET} else 1)" || die "budget: spent ${SP} + cap ${CAP} > ${BUDGET} USD"
  python3 -c "import sys; sys.exit(0 if ${CAP} <= ${LAUNCH_CAP} else 1)" || die "per-launch cap: ${CAP} > ${LAUNCH_CAP} USD (shorter --minutes or a smaller --type)"
  local SD; SD="$(spent_today)"
  python3 -c "import sys; sys.exit(0 if ${SD}+${CAP} <= ${DAY_CAP} else 1)" || die "daily cap: today ${SD} + ${CAP} > ${DAY_CAP} USD"
  # every staged sync path must exist in the mirror (a missing input would fail 10 minutes into the run)
  for rel in "${SYNC[@]+"${SYNC[@]}"}"; do
    aws s3 ls "s3://${BUCKET}/input/lab/repo/${rel}" --region "$REGION" >/dev/null 2>&1 || die "not staged: ${rel} (run: labcloud.sh stage ${rel})"
  done
  # code: src, scripts, the lock, every experiment's *.py / *.sh (small), plus --code-extra paths
  ( cd "$ROOT" && { printf '%s\n' pyproject.toml uv.lock .python-version README.md; find src scripts -type f ! -name '*.pyc' ! -path '*/__pycache__/*' -size -5M;
      find experiments -type f \( -name '*.py' -o -name '*.sh' \) -size -300k ! -path '*/__pycache__/*' ! -path '*/.venv/*' \
        ! -path '*/site-packages/*' ! -path '*/node_modules/*' ! -path 'experiments/800-checker/cloud/runs/*' \
        ! -name '.frozen-*' ! -name 'pin-*.py' ! -name '*.builder-before.py';     # venvs and builder byte copies out (was 141 MB)
      find data/challenge -maxdepth 1 -type f \( -name '*.csv' -o -name '*.json' \);
      for e in "${EXTRA[@]+"${EXTRA[@]}"}"; do find "$e" -type f; done; } | sort -u > "$D/code.list"
    COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs -czf "$D/code.tar.gz" -T "$D/code.list" 2>/dev/null )
  python3 - "$D/job.json" "$NAME" "$PAR" "$VENV" "${#SYNC[@]}" "${SYNC[@]+"${SYNC[@]}"}" "${COLLECT[@]+"${COLLECT[@]}"}" <<'PY'
import json, sys
out, name, par, venv, ns = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] == "true", int(sys.argv[5])
rest = sys.argv[6:]
json.dump(dict(name=name, parallel=par, venv=venv, sync=rest[:ns], collect=rest[ns:]), open(out, "w"), indent=1)
PY
  cp "$JOBS" "$D/jobs.txt"
  cp "$HERE/labworker.sh" "$D/labworker.sh"
  local WSHA; WSHA="$(shasum -a 256 "$D/labworker.sh" | cut -d' ' -f1)"
  for f in job.json jobs.txt code.tar.gz labworker.sh; do
    aws s3 cp "$D/$f" "s3://${BUCKET}/${PREFIX}/${f}" --region "$REGION" --only-show-errors
  done
  rm -f "$D/code.tar.gz"          # it lives on S3 now (s3://…/${PREFIX}/code.tar.gz); keep only its list locally
  local AMI; AMI="$(aws ssm get-parameter --region "$REGION" --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 --query 'Parameter.Value' --output text)"
  local DEADLINE; DEADLINE="$(date -u -v+"${MIN}"M +%FT%TZ 2>/dev/null || date -u -d "+${MIN} min" +%FT%TZ)"
  local UD
  UD="$(cat <<EOF
#!/bin/bash
shutdown -h +${MIN}
export AWS_REGION=${REGION}
aws s3 cp s3://${BUCKET}/${PREFIX}/labworker.sh /root/labworker.sh --only-show-errors
printf '${WSHA}  /root/labworker.sh\n' | sha256sum --check - || shutdown -h now
systemd-run --unit=vcc-lab --collect --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} --setenv=PREFIX=${PREFIX} --setenv=OUT=${OUT} /bin/bash /root/labworker.sh
EOF
)"
  local -a MARKET=()
  [[ "$SPOT" -eq 1 ]] && MARKET=(--instance-market-options 'MarketType=spot,SpotOptions={SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}')
  local IID
  IID="$(aws ec2 run-instances --region "$REGION" --image-id "$AMI" --instance-type "$TYPE" --count 1 \
    --iam-instance-profile "Name=${IAM_PROFILE}" --instance-initiated-shutdown-behavior terminate \
    --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":${ROOTGB},\"VolumeType\":\"gp3\",\"Throughput\":500,\"Iops\":6000,\"DeleteOnTermination\":true}}]" \
    --metadata-options "HttpTokens=required,HttpEndpoint=enabled" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=vcc-lab-${RID}},{Key=Project,Value=lab},{Key=purpose,Value=${NAME}},{Key=deadline,Value=${DEADLINE}},{Key=reap-after-minutes,Value=${MIN}}]" "ResourceType=volume,Tags=[{Key=Project,Value=lab},{Key=Name,Value=vcc-lab-${RID}}]" \
    --user-data "$UD" "${MARKET[@]+"${MARKET[@]}"}" --query 'Instances[0].InstanceId' --output text)"
  python3 - "$RUNS" <<PY
import json
json.dump(dict(event="launch", run_id="${RID}", utc="$(utc)", instance_id="${IID}", type="${TYPE}", market="$([[ $SPOT -eq 1 ]] && echo spot || echo ondemand)",
               price_usd_h=float("${PRICE}"), minutes_cap=${MIN}, cap_usd=float("${CAP}"), deadline="${DEADLINE}", ami="${AMI}",
               worker_sha256="${WSHA}", out="s3://${BUCKET}/${OUT}/"), open("${RUNS}", "a")); open("${RUNS}", "a").write("\n")
PY
  lablog_cost "$(utc)" "$RID" "$TYPE $( [[ $SPOT -eq 1 ]] && echo spot || echo on-demand) \$${PRICE}/h" "cap ${MIN} min" "launch, cap \$${CAP}" "$IID"
  echo "run_id=${RID}"; echo "instance_id=${IID}"; echo "cap=\$${CAP} (spent before: \$${SP} of \$${BUDGET})"; echo "out=s3://${BUCKET}/${OUT}/"
}

cmd_collect() {
  local RID="${1:?run id}" D="$HERE/runs/$1"
  mkdir -p "$D/out"
  aws s3 sync "s3://${BUCKET}/output/lab/${RID}/" "$D/out/" --region "$REGION" --only-show-errors
  echo "status: $(cat "$D/out/STATUS" 2>/dev/null || echo 'no STATUS yet')"
  [[ -f "$D/out/jobs/summary.tsv" ]] && column -t -s $'\t' "$D/out/jobs/summary.tsv" | cut -c1-200 || true
}

cmd_wait() {
  local RID="${1:?run id}"
  local IID; IID="$(python3 -c "
import json
for l in open('$RUNS'):
    r=json.loads(l)
    if r['run_id']=='$RID' and r['event']=='launch': print(r['instance_id'])")"
  [[ -n "$IID" ]] || die "unknown run $RID"
  while :; do
    if aws s3 ls "s3://${BUCKET}/output/lab/${RID}/STATUS" --region "$REGION" >/dev/null 2>&1; then break; fi
    local st; st="$(aws ec2 describe-instances --region "$REGION" --instance-ids "$IID" --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null || echo gone)"
    [[ "$st" == terminated || "$st" == shutting-down || "$st" == gone ]] && { echo "instance ${st} without STATUS"; break; }
    sleep 60
  done
  cmd_collect "$RID"
  aws ec2 terminate-instances --region "$REGION" --instance-ids "$IID" --output text >/dev/null 2>&1 || true
  aws ec2 wait instance-terminated --region "$REGION" --instance-ids "$IID" 2>/dev/null || true
  local INFO; INFO="$(aws ec2 describe-instances --region "$REGION" --instance-ids "$IID" --query 'Reservations[0].Instances[0].[State.Name,LaunchTime,StateTransitionReason]' --output text)"
  python3 - "$RUNS" "$RID" "$IID" "$INFO" <<'PY'
import datetime as dt, json, re, sys
runs, rid, iid, info = sys.argv[1:5]
state, launch, reason = (info.split("\t") + ["", ""])[:3]
L = [json.loads(l) for l in open(runs) if l.strip()]
rec = [r for r in L if r["run_id"] == rid and r["event"] == "launch"][0]
t0 = dt.datetime.fromisoformat(launch.replace("Z", "+00:00"))
m = re.search(r"\((\d{4}-\d\d-\d\d \d\d:\d\d:\d\d) GMT\)", reason)
t1 = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S").replace(tzinfo=dt.timezone.utc) if m else dt.datetime.now(dt.timezone.utc)
h = (t1 - t0).total_seconds() / 3600
usd = round(rec["price_usd_h"] * h + 200 * 0.08 / 730 * h, 3)
out = dict(event="terminated", run_id=rid, instance_id=iid, state=state, hours=round(h, 3), actual_usd=usd,
           utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
open(runs, "a").write(json.dumps(out) + "\n")
print(json.dumps(out))
PY
  local LAST; LAST="$(tail -1 "$RUNS")"
  lablog_cost "$(utc)" "$RID" "$(python3 -c "import json;r=json.loads('''$LAST''');print(r['state'])")" "$(python3 -c "import json;r=json.loads('''$LAST''');print(str(r['hours'])+' h')")" "actual \$$(python3 -c "import json;print(json.loads('''$LAST''')['actual_usd'])")" "$IID"
  echo "spent so far: \$$(spent) of \$${BUDGET}"
}

cmd_status() {
  aws ec2 describe-instances --region "$REGION" --filters Name=tag:Project,Values=lab Name=instance-state-name,Values=pending,running,stopping,shutting-down \
    --query 'Reservations[].Instances[].[InstanceId,InstanceType,State.Name,LaunchTime,Tags[?Key==`deadline`].Value|[0],Tags[?Key==`Name`].Value|[0]]' --output text
  echo "spent (logged): \$$(spent) of \$${BUDGET}"
}

cmd_reap() {
  local now; now="$(date -u +%s)"
  aws ec2 describe-instances --region "$REGION" --filters Name=tag:Project,Values=lab Name=instance-state-name,Values=pending,running \
    --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`deadline`].Value|[0]]' --output text | while read -r iid dl; do
      [[ -n "$iid" ]] || continue
      local t; t="$(python3 -c "import datetime as d;print(int(d.datetime.fromisoformat('$dl'.replace('Z','+00:00')).timestamp()))")"
      if (( now > t + 900 )); then aws ec2 terminate-instances --region "$REGION" --instance-ids "$iid" --output text; echo "reaped $iid (deadline $dl)"; fi
    done
}

sub="${1:-}"; shift || true
case "$sub" in
  stage) cmd_stage "$@";; launch) cmd_launch "$@";; wait) cmd_wait "$@";; collect) cmd_collect "$@";;
  status) cmd_status;; reap) cmd_reap;; *) sed -n '2,20p' "$0"; exit 2;;
esac
#############################################################################################################################################################
exit $?
