"""539: the leaders' board trajectories from our own leaderboard snapshots (data/leaderboard/*.json).

For each team, every distinct entry (submission_date) seen in any snapshot, with its scaled members; then per team the
decomposition of its climb since a start date into the six members (overall-equivalent = scaled delta / 6), and the
member that moved least (a frozen channel). Reads only files on disk; no network.
usage: uv run python experiments/539-rival-trajectories/traj539.py [--teams "Illumina AI,cqawesome,Testh,nailong king"] [--since 2026-09-15]
"""
import argparse
import glob
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M = ["pds", "mse", "nmae", "fid", "reach", "jac"]
ap = argparse.ArgumentParser()
ap.add_argument("--teams", default="Illumina AI,cqawesome,Testh,nailong king,Ahmet Dedeler")
ap.add_argument("--since", default="2026-09-15")
ap.add_argument("--out", default=str(ROOT / "experiments/539-rival-trajectories/traj539.json"))
a = ap.parse_args()
teams = [t.strip() for t in a.teams.split(",")]
seen = {}
snaps = sorted(glob.glob(str(ROOT / "data/leaderboard/*.json")))
for f in snaps:
    try:
        d = json.load(open(f))
    except Exception:
        continue
    rows = d["entries"] if isinstance(d, dict) and "entries" in d else d
    for r in rows if isinstance(rows, list) else []:
        if not isinstance(r, dict) or r.get("team_name") not in teams or r.get("score_avg") is None:
            continue
        k = (r["team_name"], r.get("submission_date"))
        if k not in seen:
            seen[k] = {"team": r["team_name"], "date": r.get("submission_date"), "first_seen": os.path.basename(f)[:16],
                       "n": r.get("submission_count"), "avg": r["score_avg"], **{m: r.get(f"score_{m}") for m in M}}
out = {"snapshots": len(snaps), "first_snapshot": os.path.basename(snaps[0])[:16] if snaps else None,
       "last_snapshot": os.path.basename(snaps[-1])[:16] if snaps else None, "teams": {}}
for t in teams:
    rows = sorted([v for v in seen.values() if v["team"] == t], key=lambda v: v["date"] or "")
    if not rows:
        continue
    base = next((v for v in rows if (v["date"] or "") >= a.since), rows[0])
    last = rows[-1]
    dec = {m: round(((last[m] or 0) - (base[m] or 0)) / 6, 5) for m in M}
    out["teams"][t] = {"entries_seen": len(rows), "from": {"date": base["date"], "avg": base["avg"]},
                       "to": {"date": last["date"], "avg": last["avg"], "n": last["n"]},
                       "climb": round(last["avg"] - base["avg"], 5), "by_member_overall_equiv": dec,
                       "steps": [{"date": v["date"], "avg": round(v["avg"], 5), **{m: round(v[m] or 0, 4) for m in M}}
                                 for v in rows]}
    ups = [round(rows[i]["avg"] - rows[i - 1]["avg"], 5) for i in range(1, len(rows)) if (rows[i]["date"] or "") >= a.since]
    out["teams"][t]["step_deltas_since"] = ups
    print(f"{t:14s} {base['date'][:10]} {base['avg']:.4f} -> {last['date'][:10]} {last['avg']:.4f}  climb {last['avg'] - base['avg']:+.4f}"
          f"  over {len(ups)} entries (median step {sorted(ups)[len(ups) // 2] if ups else 0:+.4f}, max {max(ups) if ups else 0:+.4f})")
    print("   by member (overall-equiv): " + "  ".join(f"{m} {dec[m]:+.4f}" for m in M))
json.dump(out, open(a.out, "w"), indent=1)
print("wrote", a.out)
