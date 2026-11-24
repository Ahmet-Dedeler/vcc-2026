"""206 — submission-time burst structure. If several accounts are one operator running a script,
their submission timestamps cluster within seconds. Over 640 live rows + 17 archived snapshots,
count how often any two DIFFERENT teams submitted within N seconds of each other, and compare
against the rate a Poisson null predicts.

Read-only.

Known-answer gate: 186 §3 recorded `Andrew Liu` (AL-2574) and `Hongjian Zhou` (pm-r14d) at the
same minute, 2026-09-02T17:47. That pair must come out of the detector.
"""

import json
import glob
import os
import sys
from collections import defaultdict
from datetime import datetime

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")

paths = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
paths.append(os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json"))
paths.append(os.path.join(HERE, "board/leaderboard_live.json"))

seen = {}
for p in paths:
    with open(p) as fh:
        d = json.load(fh)
    for e in (d["entries"] if isinstance(d, dict) else d):
        if e.get("score_avg") is None:
            continue
        key = (e["team_name"], e.get("submission_count"))
        if key not in seen:
            seen[key] = {
                "team": e["team_name"], "n": e.get("submission_count"),
                "date": e["submission_date"], "model": e.get("model_name"),
                "avg": e["score_avg"], "org": e.get("organization"),
                "mse": e.get("expr_mse_unbiased_capped_norm"),
                "pds": e.get("pds_cosine"),
            }

rows = sorted(seen.values(), key=lambda r: r["date"])
for r in rows:
    r["t"] = datetime.fromisoformat(r["date"].replace("Z", "+00:00")).timestamp()
print(f"distinct (team, submission_count) observations: {len(rows)}")

WIN = 60.0
pairs = []
for i, a in enumerate(rows):
    for b in rows[i + 1:]:
        if b["t"] - a["t"] > WIN:
            break
        if b["team"] != a["team"]:
            pairs.append((b["t"] - a["t"], a, b))

# Poisson null: rows are spread over the observed span; expected same-window cross-team pairs.
span = rows[-1]["t"] - rows[0]["t"]
n = len(rows)
exp_pairs = n * (n - 1) / 2 * (2 * WIN / span)
print(f"cross-team pairs within {WIN:.0f}s: {len(pairs)}   Poisson-null expectation: {exp_pairs:.1f}"
      f"   ratio {len(pairs)/exp_pairs:.2f}x")

gate = [p for p in pairs
        if {p[1]["team"], p[2]["team"]} == {"Andrew Liu", "Hongjian Zhou"}
        and p[1]["date"].startswith("2026-09-02T17:47")]
if not gate:
    print("GATE FAIL: 186 §3's Andrew Liu / Hongjian Zhou 09-02T17:47 pair not detected",
          file=sys.stderr)
    sys.exit(1)
print(f"GATE PASS: 186 §3's pair detected at Δt={gate[0][0]:.0f}s\n")

# ---- which team-pairs co-fire repeatedly, and how tightly
co = defaultdict(list)
for dt, a, b in pairs:
    co[tuple(sorted((a["team"], b["team"])))].append((dt, a, b))

print(f"=== TEAM PAIRS THAT CO-FIRE WITHIN {WIN:.0f}s, 2+ TIMES (sorted by count) ===")
for k, v in sorted(co.items(), key=lambda kv: (-len(kv[1]), min(d for d, _, _ in kv[1]))):
    if len(v) < 2:
        continue
    dts = sorted(d for d, _, _ in v)
    print(f"  {k[0]:<22} ~ {k[1]:<22} n={len(v):>2}  Δt={[f'{d:.0f}s' for d in dts]}")

print(f"\n=== EVERY CO-FIRE INVOLVING THE CLUSTER, Δt <= {WIN:.0f}s ===")
CL = {"Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng"}
hits = [p for p in pairs if p[1]["team"] in CL or p[2]["team"] in CL]
for dt, a, b in sorted(hits, key=lambda p: p[1]["date"]):
    print(f"  {a['date'][:19]}  Δt={dt:>5.0f}s  "
          f"{a['team']:<16}#{a['n']:<3}{str(a['model'])[:16]:<18}ovr {a['avg']:.5f}"
          f"   ||  {b['team']:<16}#{b['n']:<3}{str(b['model'])[:16]:<18}ovr {b['avg']:.5f}")

# ---- how tight are the cluster's own Δt values vs everyone else's
cl_dt = [dt for dt, a, b in pairs if a["team"] in CL and b["team"] in CL]
other = [dt for dt, a, b in pairs if not (a["team"] in CL and b["team"] in CL)]
print(f"\nwithin-cluster co-fires: n={len(cl_dt)}  Δt={sorted(f'{d:.0f}' for d in cl_dt)}")
print(f"all other co-fires:      n={len(other)}  median Δt="
      f"{sorted(other)[len(other)//2]:.0f}s")
