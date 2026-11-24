"""351 follow-up: competitive-zone improvements, unclamped-mse stratum, leader steps, jac max."""
import json, csv, os, math
WORK = "~/Code/virtual-cell-challenge-2026/experiments/351-board-forensics-and-2025-winners/work"
MEMBERS = ["score_pds", "score_mse", "score_nmae", "score_fid", "score_reach", "score_jac"]

import csv as csvm
rows = list(csvm.DictReader(open(os.path.join(WORK, "entries.csv"))))
print(f"entries rows: {len(rows)}")
# max jac ever
rows_j = [r for r in rows if r["score_jac"] not in ("", None)]
rows_j.sort(key=lambda r: -float(r["score_jac"]))
print("top-10 jac EVER (all entries):")
for r in rows_j[:10]:
    print(f"  jac={float(r['score_jac']):+.4f} avg={float(r['score_avg']):+.4f} {r['team_name'][:28]:28s} {r['model_name'][:26]:26s} {r['submission_date'][:16]}")
print("our 319:", [f"{m}={r[m]}" for r in rows if r["id"] == "xuDelEik7zGWMXdNu0R4" for m in ["score_avg", "score_jac"]])

# mse distribution: how many at clamp 0?
mse0 = sum(1 for r in rows if r["score_mse"] not in ("", None) and float(r["score_mse"]) == 0.0)
print(f"\nentries with mse==0.0: {mse0}/{len(rows)} ({mse0/len(rows):.2f})")

# competitive-zone improvements: rebuild pairs from entries.csv
from datetime import datetime
def pd(s):
    try: return datetime.fromisoformat(s.replace("Z", "+00:00"))
    except Exception: return None
teams = {}
for r in rows:
    if r["score_avg"] in ("", None): continue
    teams.setdefault(r["team_name"], []).append(r)
for t in teams:
    teams[t].sort(key=lambda e: (pd(e["submission_date"]) or datetime.min, e["id"]))

def f(x):
    return None if x in ("", None) else float(x)

zones = {"prev>=0.10": 0.10, "prev>=0.15": 0.15, "prev>=0.20": 0.20}
for zname, z in zones.items():
    imps = []
    for t, lst in teams.items():
        for a, b in zip(lst, lst[1:]):
            if a["id"] == b["id"]: continue
            pa, ca = f(a["score_avg"]), f(b["score_avg"])
            if pa is None or ca is None: continue
            if pa < z: continue
            d = ca - pa
            if d >= 0.005:
                imps.append((d, t, a, b))
    print(f"\n{zname}: {len(imps)} improvements>=0.005")
    # member carry: mean d_member/6 over these
    import collections
    carry = {m: sum((f(b[m]) or 0) - (f(a[m]) or 0) for _, _, a, b in imps) / max(1, len(imps)) / 6 for m in MEMBERS}
    for m in MEMBERS:
        print(f"   mean carry {m}: {carry[m]:+.6f}")
    imps.sort(reverse=True)
    print("   top-5 steps:")
    for d, t, a, b in imps[:5]:
        dm = {m: (f(b[m]) or 0) - (f(a[m]) or 0) for m in MEMBERS}
        s = " ".join(f"{m[-3:]}={dm[m]:+.4f}" for m in MEMBERS)
        print(f"    d={d:+.4f} {t[:26]:26s} {a['submission_date'][:10]}->{b['submission_date'][:10]} {a['model_name'][:18]:18s}->{b['model_name'][:18]:18s} {s}")

# unclamped stratum: pairs where BOTH prev and cur mse>0.001 (strictly positive)
both = []
for t, lst in teams.items():
    for a, b in zip(lst, lst[1:]):
        if a["id"] == b["id"]: continue
        ma, mb = f(a["score_mse"]), f(b["score_mse"])
        na, nb = f(a["score_nmae"]), f(b["score_nmae"])
        if None in (ma, mb, na, nb): continue
        if ma > 0.001 and mb > 0.001:
            both.append((nb - na, mb - ma))
print(f"\nunclamped-mse pairs (both>0.001): {len(both)}")
if both:
    import math
    xs = [p[0] for p in both]; ys = [p[1] for p in both]
    mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
    sxx = sum((x-mx)**2 for x in xs); syy = sum((y-my)**2 for y in ys); sxy = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    print(f"  corr(d_nmae,d_mse) unclamped = {sxy/math.sqrt(sxx*syy):+.4f}" if sxx and syy else "  degenerate")
    up = [p for p in both if p[0] > 0]
    print(f"  d_nmae>0: {len(up)}, with d_mse>=0: {sum(1 for p in up if p[1]>=0)} ({sum(1 for p in up if p[1]>=0)/max(1,len(up)):.2f})")
    # top30 teams only
    top30teams = set(r["team_name"] for r in list(csvm.DictReader(open(os.path.join(WORK, 'top30.csv')))))
    bothT = []
    for t, lst in teams.items():
        if t not in top30teams: continue
        for a, b in zip(lst, lst[1:]):
            if a["id"] == b["id"]: continue
            ma, mb = f(a["score_mse"]), f(b["score_mse"])
            na, nb = f(a["score_nmae"]), f(b["score_nmae"])
            if None in (ma, mb, na, nb): continue
            if ma > 0.001 and mb > 0.001:
                bothT.append((nb-na, mb-ma))
    print(f"  top30 unclamped pairs: {len(bothT)}")
    if bothT:
        xs = [p[0] for p in bothT]; ys = [p[1] for p in bothT]
        mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
        sxx = sum((x-mx)**2 for x in xs); syy = sum((y-my)**2 for y in ys); sxy = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
        print(f"  corr top30 unclamped = {sxy/math.sqrt(sxx*syy):+.4f}" if sxx and syy else "  degenerate")
        upT = [p for p in bothT if p[0] > 0]
        print(f"  top30 d_nmae>0: {len(upT)}, d_mse>=0: {sum(1 for p in upT if p[1]>=0)}")

# leader consecutive steps with member contributions
print("\nLEADER consecutive steps:")
lr = sorted([r for r in rows if r["team_name"] == "ViRtUaL cElL builders"],
            key=lambda e: (pd(e["submission_date"]) or datetime.min, e["id"]))
for a, b in zip(lr, lr[1:]):
    d = f(b["score_avg"]) - f(a["score_avg"])
    dm = {m: (f(b[m]) or 0) - (f(a[m]) or 0) for m in MEMBERS}
    s = " ".join(f"{m.replace('score_','')[:4]}={dm[m]/6:+.4f}" for m in MEMBERS)
    print(f"  {a['submission_date'][:16]} ({a['score_avg']}) -> {b['submission_date'][:16]} ({b['score_avg']}) d={d:+.4f} | {s} | {a['model_name'][:20]}->{b['model_name'][:20]}")
