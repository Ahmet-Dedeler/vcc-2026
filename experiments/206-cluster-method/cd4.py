"""206 — the DATA hypothesis, priced from the board.

Several teams publish, in their own descriptions, that they added a CD4+ T-cell / "Marson" /
Jurkat-matched perturbation source. Context A is Jurkat-like (AlphaBio and peng yu both say so in
public). If the cluster's cosine edge is a CONTEXT-MATCHED donor rather than more exchangeable
donors, these teams' own board arms price it.

Read-only. Known-answer gate: GISL's live row must be `GISL v15 ...` at overall 0.17233, and its
description must contain the string "Marson".
"""

import json
import glob
import os
import re
import sys
from collections import defaultdict

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments/206-cluster-method")

MEM = [("pds", "pds_cosine"), ("mse", "expr_mse_unbiased_capped_norm"),
       ("nmae", "de_wilcoxon_lfc_nmae"), ("fid", "de_wilcoxon_direction_fidelity_yield_raw"),
       ("reach", "de_wilcoxon_direction_reach_raw"), ("jac", "de_wilcoxon_sig_jaccard")]

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
        seen.setdefault((e["team_name"], e.get("submission_count")), e)
rows = sorted(seen.values(), key=lambda e: (e["team_name"], e.get("submission_count") or 0))

by_team = defaultdict(list)
for e in rows:
    by_team[e["team_name"]].append(e)

live = {e["team_name"]: e for e in
        json.load(open(os.path.join(HERE, "board/leaderboard_live.json")))["entries"]}
g = live.get("GISL")
if not g or abs(g["score_avg"] - 0.17233) > 5e-5 or "Marson" not in (g.get("description") or ""):
    print(f"GATE FAIL: GISL live row {g and g['score_avg']}", file=sys.stderr)
    sys.exit(1)
print("GATE PASS: GISL v15 @ 0.17233 with 'Marson' in the description\n")

# ---- 1. every description mentioning a T-cell / CD4 / Marson / Jurkat-matched source
PAT = r"marson|cd4|\bT[- ]cell|jurkat|GSE314342|primary T"
print("=== EVERY DESCRIPTION EVER POSTED THAT NAMES A T-CELL / CD4 SOURCE ===")
for t, es in sorted(by_team.items()):
    for e in es:
        d = str(e.get("description") or "")
        if d and re.search(PAT, d, re.I):
            print(f"\n--- {t} | #{e['submission_count']} | {e.get('model_name')} | "
                  f"{e['submission_date'][:19]} | ovr {e['score_avg']:.5f} "
                  f"| mse {e['expr_mse_unbiased_capped_norm']:.4f}")
            print("    " + d.replace("\n", "\n    ")[:1500])

# ---- 2. GISL's full observed trajectory, member by member (they price Marson in public)
def show(team):
    es = by_team.get(team, [])
    if not es:
        return
    print(f"\n===== {team} — {len(es)} of {es[-1]['submission_count']} submissions observed =====")
    print(f"{'#':>3} {'date':<12}{'model':<34}{'ovr':>9}{'pds':>8}{'mse':>8}{'nmae':>8}"
          f"{'fid':>8}{'reach':>8}{'jac':>8}")
    prev = None
    for e in es:
        v = {k: e[f] for k, f in MEM}
        print(f"{e['submission_count']:>3} {e['submission_date'][5:10]:<12}"
              f"{str(e.get('model_name'))[:33]:<34}{e['score_avg']:>9.5f}{v['pds']:>8.4f}"
              f"{v['mse']:>8.3f}{v['nmae']:>8.4f}{v['fid']:>8.4f}{v['reach']:>8.4f}{v['jac']:>8.4f}")
        if prev:
            pv = {k: prev[f] for k, f in MEM}
            print(f"    {'Δ':<45}{e['score_avg']-prev['score_avg']:>+9.5f}"
                  f"{v['pds']-pv['pds']:>+8.4f}{v['mse']-pv['mse']:>+8.3f}"
                  f"{v['nmae']-pv['nmae']:>+8.4f}{v['fid']-pv['fid']:>+8.4f}"
                  f"{v['reach']-pv['reach']:>+8.4f}{v['jac']-pv['jac']:>+8.4f}")
        prev = e


for t in ["GISL", "jhtao", "Yusif Mukhtarov", "AlphaBio", "kriagr", "AgingLab_BI", "Aginglab.com"]:
    show(t)

# ---- 3. 197's exchange rate applied to the published Marson gain
print("\n=== 197's EXCHANGE RATE APPLIED TO PUBLISHED CD4/MARSON GAINS ===")
RATE = 0.00365  # overall per +0.01 cosine
for label, overall in [("GISL's own retrodiction of Marson (local)", 0.0120),
                       ("GISL's stated BOARD gain from Marson", 0.01325),
                       ("GISL's stated board loss, KOLF at full weight", -0.00245),
                       ("GISL's stated board gain, six-panel ensemble", 0.00045)]:
    print(f"  {label:<48} {overall:>+9.5f} overall  ->  "
          f"{overall / RATE / 100:>+8.4f} cosine (if it were all direction)")
print(f"\n  our cosine 0.079249 + 0.0363 = {0.079249 + 0.01325/RATE/100:.4f}"
      f"   (cluster floor 0.13117)")
