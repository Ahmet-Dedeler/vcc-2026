"""317 — what the sweep cost, from EC2's own launch/termination times."""
import json, subprocess, datetime as dt

PRICE_HR = 0.5292 + 200 * 0.08 / 730.0   # r7i.2xlarge on-demand us-east-1 + 200 GB gp3
NAMES = ["vcc-cloud-score-ka1-317"] + [f"vcc-cloud-score-s317-{a}" for a in
         ("od1", "od8", "ks010", "ts070", "kt16000",
          "amp070", "amp085", "amp120", "amp145", "amp200")]
out = subprocess.run(
    ["aws", "ec2", "describe-instances", "--region", "us-east-1",
     "--filters", "Name=tag:Name,Values=vcc-cloud-score-*317*",
     "--query", "Reservations[].Instances[].[Tags[?Key=='Name']|[0].Value,InstanceId,"
                "LaunchTime,State.Name,StateTransitionReason]", "--output", "json"],
    capture_output=True, text=True).stdout
rows = json.loads(out) if out.strip() else []
now = dt.datetime.now(dt.timezone.utc)
tot = 0.0
print(f"{'name':<32} {'state':<12} {'launch UTC':<21} {'min':>7} {'$':>7}")
for name, iid, launch, state, reason in sorted(rows):
    t0 = dt.datetime.fromisoformat(launch.replace("Z", "+00:00"))
    t1 = now
    if "(" in (reason or ""):
        try:
            t1 = dt.datetime.strptime(reason.split("(")[1].rstrip(")"),
                                      "%Y-%m-%d %H:%M:%S %Z").replace(
                                          tzinfo=dt.timezone.utc)
        except Exception:
            pass
    mins = (t1 - t0).total_seconds() / 60
    cost = mins / 60 * PRICE_HR
    tot += cost
    print(f"{name:<32} {state:<12} {launch[:19]:<21} {mins:7.1f} {cost:7.3f}")
print(f"{'TOTAL (11 instances)':<32} {'':<12} {'':<21} {'':>7} {tot:7.2f}")
print(f"per six-member score: ${tot/max(len(rows),1):.2f}   (221 measured $1.48)")
