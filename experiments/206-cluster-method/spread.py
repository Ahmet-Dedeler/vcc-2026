"""206 — is the cluster's number SPREADING as a method, or is the board simply shrinking its
amplitude and revealing directions it already held?

197's floor is a function of raw `mse` alone. A team that is over-amplified has an UNINFORMATIVE
floor no matter how good its direction. So "more teams at a high floor" is ambiguous between
"more teams got a better direction" and "more teams found their amplitude optimum". The
discriminator is `pds`, which is scale-invariant: an amplitude-only arm moves the floor a lot and
`pds` by ~0.

Read-only. Known-answer gate: 186 §1.2 counted 8 teams with raw `mse` < 1.0 on the Aug-28
snapshot and 44 on the 09-03 00:30 snapshot. Both must be reproduced.
"""

import json
import glob
import os
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")
NUM, DEN = 1.023611, 2.074495
cos_floor = lambda m: (NUM - m) / DEN

paths = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
paths.append(os.path.join(HERE, "board/leaderboard_live.json"))

snaps = []
for p in paths:
    with open(p) as fh:
        d = json.load(fh)
    ents = [e for e in (d["entries"] if isinstance(d, dict) else d) if e.get("score_avg") is not None]
    snaps.append((os.path.basename(p), ents))

a28 = [e for n, e in snaps if n == "lb_20260828T1325.json"][0]
n_a28 = sum(1 for e in a28 if e["expr_mse_unbiased_capped_norm"] < 1.0)
s0903 = json.load(open(os.path.join(
    ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json")))["entries"]
n_0903 = sum(1 for e in s0903
             if e.get("score_avg") is not None and e["expr_mse_unbiased_capped_norm"] < 1.0)
if n_a28 != 8 or n_0903 != 44:
    print(f"GATE FAIL: mse<1.0 counts {n_a28} (want 8) / {n_0903} (want 44)", file=sys.stderr)
    sys.exit(1)
print(f"GATE PASS: 186 §1.2's raw mse<1.0 counts reproduced — Aug 28 = {n_a28}, "
      f"Sep 3 00:30 = {n_0903}\n")

print("=== THE COSINE-FLOOR DISTRIBUTION OVER TIME ===")
print(f"{'snapshot':<26}{'teams':>7}{'flr>0':>7}{'>0.05':>7}{'>0.09':>7}{'>0.11':>7}"
      f"{'best flr':>10}{'pds>0.83':>10}")
for n, ents in snaps:
    fl = [cos_floor(e["expr_mse_unbiased_capped_norm"]) for e in ents]
    print(f"{n:<26}{len(ents):>7}{sum(1 for x in fl if x > 0):>7}"
          f"{sum(1 for x in fl if x > 0.05):>7}{sum(1 for x in fl if x > 0.09):>7}"
          f"{sum(1 for x in fl if x > 0.11):>7}{max(fl):>10.5f}"
          f"{sum(1 for e in ents if e['pds_cosine'] > 0.83):>10}")

# ---- the discriminator, over every consecutive same-team pair we can observe
seen = {}
for p in paths + [os.path.join(ROOT, "experiments/186-intel-0903/board/leaderboard_2026-09-03.json")]:
    with open(p) as fh:
        d = json.load(fh)
    for e in (d["entries"] if isinstance(d, dict) else d):
        if e.get("score_avg") is None:
            continue
        seen.setdefault((e["team_name"], e.get("submission_count")), e)
by = defaultdict(list)
for e in seen.values():
    by[e["team_name"]].append(e)
for t in by:
    by[t].sort(key=lambda e: (e.get("submission_count") or 0, e["submission_date"]))

print("\n=== EVERY OBSERVED ARM THAT RAISED THE COSINE FLOOR BY >= 0.03 ===")
print("    Δpds is the discriminator: ~0 means the direction did not change, only the amplitude")
print(f"{'team':<22}{'arm':<30}{'date':<12}{'Δfloor':>9}{'Δpds':>9}{'Δreach':>9}{'Δnmae':>9}"
      f"{'Δjac':>8}{'Δovr':>9}")
rowsout = []
for t, es in by.items():
    for a, b in zip(es, es[1:]):
        df = cos_floor(b["expr_mse_unbiased_capped_norm"]) - cos_floor(a["expr_mse_unbiased_capped_norm"])
        if df < 0.03 or cos_floor(b["expr_mse_unbiased_capped_norm"]) < 0.05:
            continue
        rowsout.append((df, t, a, b))
for df, t, a, b in sorted(rowsout, key=lambda r: -cos_floor(r[3]["expr_mse_unbiased_capped_norm"])):
    print(f"{t:<22}{(str(a.get('model_name'))[:13] + ' -> ' + str(b.get('model_name'))[:13]):<30}"
          f"{b['submission_date'][5:10]:<12}{df:>+9.4f}{b['pds_cosine']-a['pds_cosine']:>+9.4f}"
          f"{b['de_wilcoxon_direction_reach_raw']-a['de_wilcoxon_direction_reach_raw']:>+9.4f}"
          f"{b['de_wilcoxon_lfc_nmae']-a['de_wilcoxon_lfc_nmae']:>+9.4f}"
          f"{b['de_wilcoxon_sig_jaccard']-a['de_wilcoxon_sig_jaccard']:>+8.4f}"
          f"{b['score_avg']-a['score_avg']:>+9.5f}")

dp = [abs(b["pds_cosine"] - a["pds_cosine"]) for _, _, a, b in rowsout]
print(f"\n  n = {len(rowsout)} such arms.  median |Δpds| = {sorted(dp)[len(dp)//2]:.4f}"
      f"   max |Δpds| = {max(dp):.4f}")
print("  For reference, moving cosine 0.079 -> 0.131 on 159's curve is a DIRECTION change of the")
print("  kind that lifted OUR OWN pds by 0.09 over the same cosine span (197 §3.1).")
