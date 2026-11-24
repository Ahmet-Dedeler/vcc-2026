#!/usr/bin/env python3
"""608 (board forensics for [21dd4c], 2026-10-02): every step (consecutive distinct entries of one team) by the CURRENT
top 30, with 610's columns (raw member deltas + overall-equivalent per member), since a cut time; DE-only flag as in
destep610 (raw pds and mse unchanged < 2e-4). Input: hist610.py's all_entries.json (data/leaderboard/*.json + lb_now).
usage: python3 steps608.py <all_entries.json> <lb_now.json> <since ISO, e.g. 2026-10-02T10:43>"""
import collections
import json
import sys

K = dict(pds="pds_cosine", mse="expr_mse_unbiased_capped_norm", nmae="de_wilcoxon_lfc_nmae",
         fid="de_wilcoxon_direction_fidelity_yield_raw", reach="de_wilcoxon_direction_reach_raw", jac="de_wilcoxon_sig_jaccard")
SC = dict(pds="score_pds", mse="score_mse", nmae="score_nmae", fid="score_fid", reach="score_reach", jac="score_jac")
es = json.load(open(sys.argv[1]))
now = json.load(open(sys.argv[2]))["entries"]
since = sys.argv[3]
top = sorted((e for e in now if e.get("score_avg") is not None), key=lambda e: -e["score_avg"])[:30]
rank = {e["team_name"]: i + 1 for i, e in enumerate(top)}
by = collections.defaultdict(list)
for e in es:
    if e.get("score_avg") is not None and e.get(K["fid"]) is not None and e["team_name"] in rank:
        by[e["team_name"]].append(e)
print(f"top-30 cut {top[-1]['score_avg']:.4f}; steps arriving at or after {since}Z")
n = 0
for t in sorted(by, key=lambda x: rank[x]):
    L = sorted(by[t], key=lambda e: e["submission_date"])
    for a, b in zip(L, L[1:]):
        if b["submission_date"][:16] < since:
            continue
        n += 1
        de_only = abs(a[K["pds"]] - b[K["pds"]]) < 2e-4 and abs(a[K["mse"]] - b[K["mse"]]) < 2e-4
        d = {m: b[K[m]] - a[K[m]] for m in K}
        ds = {m: (b[SC[m]] - a[SC[m]]) / 6 for m in SC}
        print(f"#{rank[t]:>2} {t[:18]:18s} {a['submission_date'][:16]} -> {b['submission_date'][:16]}  overall {a['score_avg']:.4f} -> "
              f"{b['score_avg']:.4f} ({b['score_avg'] - a['score_avg']:+.4f}) {'DE-ONLY' if de_only else 'mixed  '} | raw "
              + " ".join(f"{m} {d[m]:+.4f}" for m in K) + " | ov-eq " + " ".join(f"{m[0:2]}{ds[m]:+.4f}" for m in SC)
              + f" | {b.get('model_name') or ''!s:.28} | {(b.get('description') or '').replace(chr(10), ' ')[:160]}")
print(f"{n} steps")
