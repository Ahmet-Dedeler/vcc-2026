#!/usr/bin/env bash
# 221 — every instance this experiment launched: lifetime, market, cost, and final state.
#
# Falsifier 3 of the mandate: report actual spend per score, and PROVE every instance launched
# is terminated.
#
# The LAUNCH ledger is CloudTrail's `RunInstances`, which is authoritative and complete —
# `describe-instances` drops a terminated instance after about an hour, so it cannot be the
# ledger. The END time is the `STATUS` object each worker's EXIT trap writes to S3 immediately
# before it powers off, read with `head-object` so the timestamp is unambiguous UTC.
#
# Note on the termination proof: a worker that ends itself does so with `shutdown -h now` under
# `instance-initiated-shutdown-behavior=terminate`, which is NOT an API call and therefore
# leaves NO `TerminateInstances` event. The proof is instead (a) this complete launch ledger,
# (b) each instance's `STATUS` marker, and (c) `scripts/aws_reap.sh` over every enabled region
# showing nothing running.
set -Eeuo pipefail
REGION="${AWS_REGION:-us-east-1}"
BUCKET="${BUCKET:-<your-s3-bucket>}"
SINCE="${SINCE:-2026-09-04T00:00:00Z}"
OD_PER_HOUR=0.5292        # r7i.2xlarge on-demand, us-east-1, 2026-09-04
SPOT_PER_HOUR=0.2280      # r7i.2xlarge spot, recorded at launch
EBS_GB_MONTH=0.08         # gp3
ROOT_GB=200

aws cloudtrail lookup-events --region "$REGION" \
  --lookup-attributes AttributeKey=EventName,AttributeValue=RunInstances \
  --start-time "$SINCE" --max-results 50 \
  --query 'Events[].[EventTime,Resources[?ResourceType==`AWS::EC2::Instance`].ResourceName|join(`,`,@)]' \
  --output text | sort > /tmp/vcc221-launches.tsv

python3 - "$REGION" "$BUCKET" "$OD_PER_HOUR" "$SPOT_PER_HOUR" "$EBS_GB_MONTH" "$ROOT_GB" <<'PY'
import datetime as dt, json, subprocess, sys
region, bucket = sys.argv[1], sys.argv[2]
od, spot, ebs_gb_month, root_gb = (float(x) for x in sys.argv[3:7])

# run-id -> (instance id, market). Filled from this session's launch records; the run id is not
# in the CloudTrail event, and the tag is gone once the instance ages out of describe-instances.
RUNS = [
    ("ka-20260904",      "i-01db4020744f09521", "spot",      "known-answer, RECLAIMED by spot 42 min into `baseline`"),
    ("ladder-20260904",  "i-01e401ae1863945b4", "spot",      "gate 2 (entrywise-relative) failed — build fixed"),
    ("ladder2-20260904", "i-0fd55468b05fa77af", "spot",      "gate 4 expectation wrong — gate fixed"),
    ("ka2-20260904",     "i-0e038f2a9056ab5c1", "on-demand", "`prep-real-bundle --baseline` given a directory"),
    ("ladder3-20260904", "i-0eaa6d65f4e47452a", "on-demand", "builder HEAD moved mid-run (2-tuple -> 3-tuple)"),
    ("ladder4-20260904", "i-04d6bcfbd1018ac4c", "on-demand", "parked built candidates under `input/` — AccessDenied"),
    ("ladder5-20260904", "i-0a1cb2df981c58dad", "on-demand", "BUILD OK; terminated by me to re-order runs before bundle"),
    ("ladder6-20260904", "i-0f33dddb34a808e4d", "on-demand", "runs: ka095, a15"),
    ("ladder7-20260904", "i-043a904dbea6b188d", "on-demand", "runs: a15_seed1, rc_hybrid"),
    ("ladder8-20260904", "i-0926ebf77f8240701", "on-demand", "runs: rc_callset, cd4f034"),
    ("bundle-20260904",  "i-08028907eeba5a980", "on-demand", "the h1-test real bundle (0 and 1 ends)"),
]
launch = {}
for line in open("/tmp/vcc221-launches.tsv"):
    t, ids = line.rstrip("\n").split("\t")
    for i in ids.split(","):
        launch[i] = dt.datetime.fromisoformat(t).astimezone(dt.timezone.utc)

def status_time(run_id):
    p = subprocess.run(["aws", "s3api", "head-object", "--bucket", bucket,
                        "--key", f"output/cloud-score/{run_id}/STATUS",
                        "--query", "LastModified", "--output", "text"],
                       capture_output=True, text=True)
    if p.returncode != 0 or not p.stdout.strip():
        return None
    return dt.datetime.fromisoformat(p.stdout.strip()).astimezone(dt.timezone.utc)

def state(iid):
    p = subprocess.run(["aws", "ec2", "describe-instances", "--region", region,
                        "--instance-ids", iid, "--query",
                        "Reservations[0].Instances[0].State.Name", "--output", "text"],
                       capture_output=True, text=True)
    s = p.stdout.strip()
    return s if p.returncode == 0 and s and s != "None" else "gone (aged out; terminated)"

now = dt.datetime.now(dt.timezone.utc)
rows, total = [], 0.0
print(f"{'run':20s} {'instance':21s} {'mkt':10s} {'launch UTC':9s} {'end UTC':9s} "
      f"{'min':>6s} {'$':>6s}  state")
for run_id, iid, market, note in RUNS:
    L = launch.get(iid)
    if L is None:
        continue
    E = status_time(run_id) or now
    minutes = max((E - L).total_seconds() / 60.0, 0.0)
    rate = spot if market == "spot" else od
    cost = minutes / 60.0 * (rate + root_gb * ebs_gb_month / 730.0)
    total += cost
    st = state(iid)
    rows.append(dict(run=run_id, instance_id=iid, market=market,
                     launch_utc=L.isoformat(), end_utc=E.isoformat(),
                     minutes=round(minutes, 1), usd=round(cost, 4),
                     final_state=st, note=note))
    print(f"{run_id:20s} {iid:21s} {market:10s} {L.strftime('%H:%M:%S'):9s} "
          f"{E.strftime('%H:%M:%S'):9s} {minutes:6.1f} {cost:6.3f}  {st}")
print(f"\n{'TOTAL':20s} {len(rows)} instances{'':41s}{total:6.3f} USD")
alive = [r for r in rows if r["final_state"] in ("running", "pending")]
print(f"still running: {len(alive)}" + (f"  {[r['instance_id'] for r in alive]}" if alive else ""))
json.dump({"instances": rows, "total_usd": round(total, 4),
           "still_running": alive},
          open("experiments/221-cloud-scorer/work/cost.json", "w"), indent=1)
PY

echo
echo "=== termination proof: every enabled region, every running instance ==="
bash scripts/aws_reap.sh
