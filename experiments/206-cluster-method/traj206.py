"""206 — reconstruct the joining trajectories member by member, and mine every model_name and
description on the board with a diff against 186's snapshots.

Read-only. Writes only into experiments/206-cluster-method/work/.

Known-answer gate: 186 §1.3 recorded jw_13288 at overall 0.1810 with model `cedar_90` on the
2026-09-03 00:30 snapshot, and Andrew Liu at 0.2083 / `AL-2574`. Both must come back out of the
reconstructed trajectories.
"""

import json
import glob
import os
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")
WORK = os.path.join(HERE, "work")

MEMBERS = [
    ("pds", "pds_cosine"),
    ("mse", "expr_mse_unbiased_capped_norm"),
    ("nmae", "de_wilcoxon_lfc_nmae"),
    ("fid", "de_wilcoxon_direction_fidelity_yield_raw"),
    ("reach", "de_wilcoxon_direction_reach_raw"),
    ("jac", "de_wilcoxon_sig_jaccard"),
]
SCALED = [("pds", "score_pds"), ("mse", "score_mse"), ("nmae", "score_nmae"),
          ("fid", "score_fid"), ("reach", "score_reach"), ("jac", "score_jac")]

NUM, DEN = 1.023611, 2.074495
cos_floor = lambda m: (NUM - m) / DEN


def row(e):
    d = {"team": e.get("team_name"), "org": e.get("organization"),
         "members": tuple(e.get("team_members") or []),
         "model": e.get("model_name"), "desc": e.get("description"),
         "date": e.get("submission_date"), "n": e.get("submission_count"),
         "avg": e.get("score_avg")}
    for k, f in MEMBERS:
        d[k] = e.get(f)
    for k, f in SCALED:
        d["s_" + k] = e.get(f)
    return d


paths = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
paths.append(os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json"))
paths.append(os.path.join(HERE, "board/leaderboard_live.json"))

traj = defaultdict(dict)
snapshot_names = {}
for p in paths:
    with open(p) as fh:
        d = json.load(fh)
    ents = d["entries"] if isinstance(d, dict) else d
    tag = os.path.basename(p)
    snapshot_names[p] = tag
    for e in ents:
        r = row(e)
        if r["avg"] is None:
            continue
        key = (r["team"], r["n"])
        if key not in traj[r["team"]]:
            traj[r["team"]][key] = dict(r, first_seen=tag)

TRAJ = {t: sorted(v.values(), key=lambda x: (x["n"] or 0, x["date"] or "")) for t, v in traj.items()}

# ---- gate
g = TRAJ.get("jw_13288", [])
hit = [r for r in g if r["model"] == "cedar_90"]
if not hit or abs(hit[0]["avg"] - 0.1810) > 5e-4:
    print("GATE FAIL: jw_13288 cedar_90 not reproduced", file=sys.stderr)
    sys.exit(1)
al = [r for r in TRAJ.get("Andrew Liu", []) if r["model"] == "AL-2574"]
if not al or abs(al[0]["avg"] - 0.2083) > 5e-4:
    print("GATE FAIL: Andrew Liu AL-2574 not reproduced", file=sys.stderr)
    sys.exit(1)
print("GATE PASS: 186 §1.3 rows reproduced from the reconstructed trajectories\n")

HDR = (f"{'#':>3} {'date':<20}{'model':<20}{'ovr':>9}{'pds':>9}{'mse':>9}{'nmae':>9}"
       f"{'fid':>8}{'reach':>8}{'jac':>8}{'cosflr':>8}")


def show(team):
    rows = TRAJ.get(team, [])
    print(f"===== {team} ({rows[-1]['org']!r}, members={rows[-1]['members']}) — "
          f"{len(rows)} of {rows[-1]['n']} submissions observed =====")
    print(HDR)
    prev = None
    for r in rows:
        cf = cos_floor(r["mse"])
        print(f"{r['n']:>3} {r['date'][:19]:<20}{str(r['model'])[:19]:<20}{r['avg']:>9.5f}"
              f"{r['pds']:>9.4f}{r['mse']:>9.4f}{r['nmae']:>9.4f}{r['fid']:>8.4f}"
              f"{r['reach']:>8.4f}{r['jac']:>8.4f}{cf:>8.5f}")
        if prev is not None:
            d = {k: r[k] - prev[k] for k, _ in MEMBERS}
            ds = {k: (r['s_' + k] - prev['s_' + k]) for k, _ in SCALED}
            print(f"    Δ raw          {'':<20}{r['avg']-prev['avg']:>+9.5f}{d['pds']:>+9.4f}"
                  f"{d['mse']:>+9.4f}{d['nmae']:>+9.4f}{d['fid']:>+8.4f}{d['reach']:>+8.4f}"
                  f"{d['jac']:>+8.4f}{cf-cos_floor(prev['mse']):>+8.5f}")
            print(f"    Δ scaled(/6)   {'':<20}{'':>9}{ds['pds']/6:>+9.5f}{ds['mse']/6:>+9.5f}"
                  f"{ds['nmae']/6:>+9.5f}{ds['fid']/6:>+8.5f}{ds['reach']/6:>+8.5f}"
                  f"{ds['jac']/6:>+8.5f}")
        prev = r
    print()


for t in ["MT Zeng", "jw_13288", "Hongjian Zhou", "Andrew Liu", "cqawesome", "Ahmet Dedeler"]:
    show(t)

# ---------------------------------------------------------------- description / model-name mining
live = [row(e) for e in json.load(open(os.path.join(HERE, "board/leaderboard_live.json")))["entries"]]
prev_board = {r["team"]: r for r in
              [row(e) for e in json.load(open(os.path.join(
                  ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json")))["entries"]]}
aug28 = {r["team"]: r for r in
         [row(e) for e in json.load(open(os.path.join(
             ROOT, "data/leaderboard/lb_20260828T1325.json")))["entries"]]}

with_desc = [r for r in live if r["desc"]]
print(f"live rows with a non-empty description: {len(with_desc)} of {len(live)}")

new_desc = [r for r in with_desc
            if r["team"] not in prev_board or prev_board[r["team"]]["desc"] != r["desc"]]
print(f"  ...of which NEW or CHANGED since the 09-03 00:30 snapshot: {len(new_desc)}\n")

json.dump({"traj": {t: TRAJ[t] for t in
                    ["MT Zeng", "jw_13288", "Hongjian Zhou", "Andrew Liu", "cqawesome"]},
           "new_desc": new_desc},
          open(os.path.join(WORK, "traj206.json"), "w"), indent=1, default=str)

print("===== NEW / CHANGED DESCRIPTIONS SINCE 09-03 00:30 =====")
for r in sorted(new_desc, key=lambda x: -x["avg"]):
    print(f"\n--- {r['team']} | {r['model']} | ovr {r['avg']:.5f} | org={r['org']!r} | "
          f"mse {r['mse']:.4f} cosflr {cos_floor(r['mse']):.5f}")
    print("    " + str(r["desc"]).replace("\n", "\n    ")[:1600])
