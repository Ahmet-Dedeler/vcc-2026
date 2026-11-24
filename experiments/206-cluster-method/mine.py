"""206 — mine every model_name/description on the live board, rebuild the shared-prediction
cluster, reconstruct the two joining trajectories member by member, and invert 197's cosine floor
for the whole field.

Read-only. Writes only into experiments/206-cluster-method/work/.

Known-answer gates (abort on mismatch):
  KA-1  197's floor identity on OUR OWN displayed row: raw 0.8738045 -> floor 0.072240,
        against 197's measured cosine 0.079249 (tightness ratio 1.0970).
  KA-2  the mandate's stated cluster floor: raw 0.7515 -> 0.13117.
  KA-3  186 §1.4's Sep-03 00:30 row for Hongjian Zhou (pm-r14d) must be reproduced exactly out of
        the archived snapshot.
  KA-4  159's law + a15 board<->surface anchor must return k_eff = 3.09 for our own cosine.
"""

import json
import glob
import math
import os
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")
WORK = os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)

# ---------------------------------------------------------------- 197's inversion
NUM = 1.023611
DEN = 2.074495
A15_ANCHOR = 2.078          # board cosine / surface cosine
LAW_C, LAW_B = 66.5, 1920.0  # 1/cos_surface^2 = C + B/k


def cos_floor(raw_mse: float) -> float:
    """197 §2.1. One-sided; vacuous at or above raw 1.023611."""
    return (NUM - raw_mse) / DEN


def k_eff(board_cos: float, anchor: float = A15_ANCHOR) -> float:
    s = board_cos / anchor
    d = 1.0 / (s * s) - LAW_C
    return LAW_B / d if d > 0 else float("inf")


MEMBERS = [
    ("pds", "pds_cosine"),
    ("mse", "expr_mse_unbiased_capped_norm"),
    ("nmae", "de_wilcoxon_lfc_nmae"),
    ("fid", "de_wilcoxon_direction_fidelity_yield_raw"),
    ("reach", "de_wilcoxon_direction_reach_raw"),
    ("jac", "de_wilcoxon_sig_jaccard"),
]
SCALED = [
    ("pds", "score_pds"), ("mse", "score_mse"), ("nmae", "score_nmae"),
    ("fid", "score_fid"), ("reach", "score_reach"), ("jac", "score_jac"),
]


def row(e):
    d = {
        "team": e.get("team_name"),
        "team_id": e.get("team_id"),
        "org": e.get("organization"),
        "members": e.get("team_members") or [],
        "model": e.get("model_name"),
        "desc": e.get("description"),
        "date": e.get("submission_date"),
        "n": e.get("submission_count"),
        "avg": e.get("score_avg"),
        "status": e.get("status"),
    }
    for k, f in MEMBERS:
        d[k] = e.get(f)
    for k, f in SCALED:
        d["s_" + k] = e.get(f)
    return d


def load(path):
    with open(path) as fh:
        d = json.load(fh)
    ents = d["entries"] if isinstance(d, dict) else d
    return [row(e) for e in ents]


# ---------------------------------------------------------------- load every snapshot
snap_paths = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
snap_paths += [os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json")]
LIVE = os.path.join(HERE, "board/leaderboard_live.json")
snap_paths.append(LIVE)

snaps = {}
for p in snap_paths:
    try:
        snaps[p] = load(p)
    except Exception as exc:  # pragma: no cover
        print(f"skip {p}: {exc}", file=sys.stderr)

live = snaps[LIVE]
print(f"snapshots loaded: {len(snaps)}   live rows: {len(live)}")

# ---------------------------------------------------------------- gates
fails = []

# 197 §2 table row a15: raw 0.87381 -> floor 0.07221, true 0.07925, tightness (floor/true) 0.911
f_ours = cos_floor(0.8738045054832962)
if abs(f_ours - 0.07221) > 5e-6:
    fails.append(f"KA-1 floor(ours)={f_ours:.6f} != 197's 0.07221")
tight = f_ours / 0.079249
if abs(tight - 0.911) > 1e-3:
    fails.append(f"KA-1 tightness={tight:.4f} != 197's 0.911")
ratio = 1.0 / tight

f_cl = cos_floor(0.7515)
if abs(f_cl - 0.13117) > 5e-5:
    fails.append(f"KA-2 floor(0.7515)={f_cl:.5f} != 0.13117")

k_ours = k_eff(0.079249)
if abs(k_ours - 3.09) > 0.02:
    fails.append(f"KA-4 k_eff(ours)={k_ours:.3f} != 3.09")

s0903 = snaps[os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json")]
hz = [r for r in s0903 if r["team"] == "Hongjian Zhou"]
if not hz or hz[0]["model"] != "pm-r14d":
    fails.append(f"KA-3 Hongjian Zhou row on 09-03 00:30 is {hz[0]['model'] if hz else None}, want pm-r14d")
else:
    want = dict(pds=0.8325, mse=0.760, nmae=0.9291, fid=0.5084, reach=0.2454, jac=0.0224)
    for k, v in want.items():
        if abs(hz[0][k] - v) > 6e-4:
            fails.append(f"KA-3 {k}={hz[0][k]:.4f} != {v}")

if fails:
    print("KNOWN-ANSWER GATES FAILED:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print(f"GATES PASS  floor(ours)={f_ours:.6f} tightness={ratio:.4f} floor(cluster)={f_cl:.5f} "
      f"k_eff(ours)={k_ours:.3f}")

# ---------------------------------------------------------------- trajectories
# The API exposes only each team's LATEST row, so a trajectory is the union over snapshots keyed by
# (team, submission_count).
traj = defaultdict(dict)
for p, rows in snaps.items():
    tag = os.path.basename(p)
    for r in rows:
        if r["avg"] is None:
            continue
        key = (r["team"], r["n"])
        if key not in traj[r["team"]]:
            traj[r["team"]][key] = dict(r, seen=[tag])
        else:
            traj[r["team"]][key]["seen"].append(tag)

TRAJ = {t: sorted(v.values(), key=lambda x: (x["n"], x["date"])) for t, v in traj.items()}
print(f"teams with any observed row: {len(TRAJ)}")

# ---------------------------------------------------------------- the cluster
CLUSTER_SEED = ["Hongjian Zhou", "jw_13288", "Andrew Liu", "MT Zeng"]
seed = {r["team"]: r for r in live if r["team"] in CLUSTER_SEED}
ctr = {k: sum(seed[t][k] for t in seed) / len(seed) for k, _ in MEMBERS}

TOL = dict(pds=0.0025, mse=0.020, nmae=0.0035, fid=0.0025, reach=0.020, jac=0.0012)


def in_cluster(r):
    return all(r[k] is not None and abs(r[k] - ctr[k]) <= TOL[k] for k, _ in MEMBERS)


cluster_live = sorted([r for r in live if in_cluster(r)], key=lambda r: -r["avg"])

# and: has any team EVER been in the cluster (over all snapshots)?
ever = defaultdict(list)
for t, rows in TRAJ.items():
    for r in rows:
        if in_cluster(r):
            ever[t].append(r)

out = {
    "generated": "2026-09-03",
    "live_rows": len(live),
    "gates": {"floor_ours": f_ours, "tightness": ratio, "floor_cluster": f_cl, "k_eff_ours": k_ours},
    "cluster_centre": ctr,
    "cluster_live": cluster_live,
    "ever_in_cluster": {t: len(v) for t, v in ever.items()},
}

# ---------------------------------------------------------------- field-wide inversion
field = []
for r in live:
    if r["mse"] is None:
        continue
    fl = cos_floor(r["mse"])
    field.append({
        "team": r["team"], "model": r["model"], "avg": r["avg"], "mse": r["mse"],
        "pds": r["pds"], "reach": r["reach"], "nmae": r["nmae"], "jac": r["jac"], "fid": r["fid"],
        "org": r["org"], "n": r["n"], "date": r["date"],
        "cos_floor": fl, "informative": fl > 0, "k_eff": k_eff(fl) if fl > 0 else None,
    })
field.sort(key=lambda x: -x["cos_floor"])
out["field_top_cosine"] = field[:40]
out["n_informative"] = sum(1 for f in field if f["informative"])

with open(os.path.join(WORK, "mine.json"), "w") as fh:
    json.dump(out, fh, indent=1, default=str)

# ---------------------------------------------------------------- print
print("\n=== LIVE CLUSTER (auto-detected, tol on all six) ===")
print(f"{'team':<22}{'model':<16}{'ovr':>9}{'pds':>9}{'mse':>8}{'nmae':>8}{'fid':>8}"
      f"{'reach':>8}{'jac':>8}{'n':>4}  {'cosfloor':>8} {'k':>5}  org")
for r in cluster_live:
    fl = cos_floor(r["mse"])
    print(f"{r['team']:<22}{str(r['model'])[:15]:<16}{r['avg']:>9.5f}{r['pds']:>9.4f}"
          f"{r['mse']:>8.4f}{r['nmae']:>8.4f}{r['fid']:>8.4f}{r['reach']:>8.4f}{r['jac']:>8.4f}"
          f"{r['n']:>4}  {fl:>8.5f} {k_eff(fl):>5.1f}  {str(r['org'])[:28]}")

us = [r for r in live if r["team"] == "Ahmet Dedeler"][0]
fl = cos_floor(us["mse"])
print(f"{'US':<22}{str(us['model'])[:15]:<16}{us['avg']:>9.5f}{us['pds']:>9.4f}{us['mse']:>8.4f}"
      f"{us['nmae']:>8.4f}{us['fid']:>8.4f}{us['reach']:>8.4f}{us['jac']:>8.4f}{us['n']:>4}  "
      f"{fl:>8.5f} {k_eff(fl):>5.1f}")

print("\n=== TOP 20 BY INVERTED COSINE FLOOR (whole board) ===")
print(f"{'team':<24}{'model':<18}{'ovr':>9}{'mse':>8}{'pds':>8}{'cosfloor':>10}{'k_eff':>7}")
for f in field[:20]:
    print(f"{f['team']:<24}{str(f['model'])[:17]:<18}{f['avg']:>9.5f}{f['mse']:>8.4f}"
          f"{f['pds']:>8.4f}{f['cos_floor']:>10.5f}{f['k_eff']:>7.1f}")
print(f"\ninformative rows (raw mse < 1.023611): {out['n_informative']} of {len(field)}")

print("\n=== TEAMS EVER MATCHING THE CLUSTER TUPLE (any snapshot) ===")
for t, v in sorted(ever.items(), key=lambda kv: -len(kv[1])):
    ms = ",".join(str(r["model"]) for r in v)
    print(f"  {t:<24} {len(v)} row(s): {ms}")
