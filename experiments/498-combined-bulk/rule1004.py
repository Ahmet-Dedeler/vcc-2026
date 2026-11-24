"""The 10-04 rule, FIXED 2026-09-27 (before any 09-28 entry exists). It is applied by chain1004.sh at 10-03 20:00Z; its
inputs are board entries that do not exist yet: 453dr (09-28 s1), the 479 pair (10-02) and 495 (10-03).

  P  = 479's pair pick (read1002's rule: both TRANSFERS -> larger dpds, tie within 0.0005 raw -> S_P2; one -> it; none -> NONE)
  C  = 495's bin vs 446bj2 (TRANSFERS: overall >= +0.0012 and dmse <= -0.006; FAILS: overall < +0.0004 or dmse > -0.002;
       else PARTIAL; UNKNOWN if 495 was not submitted)
  O  = the October path's status from 453dr vs 446b (453's bins): BAD if the bulk channel (pds+mse scaled, overall units)
       is below -0.0060, the DE channel (fid+reach+jac+nmae) is outside +-0.009, or overall < -0.02; else OK; UNKNOWN
       if 453dr was not submitted.

  s1 = 498<P> (cap-6 dump + the picked per-target form, re-levelled; vs 495 = the per-target increment on the cap-6
       base, vs 446bj2 = the total) if P != NONE and C != FAILS.
  s2 = 497dr (the October path + percentile cap + S_P2; vs 460dr) if O != BAD and C != FAILS and P == S_P2.
  Empty slots take, in order, the first READY of: 498<other form> (if C != FAILS and that form's own bin was at least
  PARTIAL), 460dr, 461dr, 469dr, 458e. READY (package, gates, novelty) is checked by the chain, not here.

usage: rule1004.py <453dr id|none> <479p2g id|none> <479sp2 id|none> <495 id|none>
prints one line: "<s1 candidates space-separated>|<s2 candidates>|P=..,C=..,O=.."  (a "-" means none)
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
BJ2, B446 = "YvVru3DAF89Jiyr34olR", "UDgjOKL9gI1udAexq6ue"


def _find(x):
    if isinstance(x, dict):
        if "score_avg" in x:
            return x
        for v in x.values():
            y = _find(v)
            if y:
                return y
    if isinstance(x, list):
        for v in x:
            y = _find(v)
            if y:
                return y


def load(i):
    p = ROOT / f"data/board-entries/{i}.json"
    return _find(json.loads(p.read_text())) if i != "none" and p.exists() else None


def d(a, b):
    return {m: a[RAW[m]] - b[RAW[m]] for m in RAW}, {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}, a["score_avg"] - b["score_avg"]


def bin479(e, ref):
    a = load(e)
    if a is None:
        return None
    dr, _, dov = d(a, ref)
    if dov >= 0.0008 and dr["pds"] >= 0.0015:
        return ("TRANSFERS", dr["pds"])
    if 0.0003 <= dov < 0.0008 or 0.0005 <= dr["pds"] < 0.0015:
        return ("PARTIAL", dr["pds"])
    return ("FAILS", dr["pds"])


e453, ep, es, e495 = (sys.argv[1:5] + ["none"] * 4)[:4]
bj2, b446 = load(BJ2), load(B446)
bp, bs = bin479(ep, bj2), bin479(es, bj2)
tp, ts = bp is not None and bp[0] == "TRANSFERS", bs is not None and bs[0] == "TRANSFERS"
if tp and ts:
    P = "S_P2" if abs(bp[1] - bs[1]) < 0.0005 else ("P2g" if bp[1] > bs[1] else "S_P2")
elif tp or ts:
    P = "P2g" if tp else "S_P2"
else:
    P = "NONE"
a495 = load(e495)
if a495 is None:
    C = "UNKNOWN"
else:
    dr, _, dov = d(a495, bj2)
    C = "TRANSFERS" if (dov >= 0.0012 and dr["mse"] <= -0.006) else "FAILS" if (dov < 0.0004 or dr["mse"] > -0.002) else "PARTIAL"
a453 = load(e453)
if a453 is None:
    O = "UNKNOWN"
else:
    dr, ds, dov = d(a453, b446)
    bulk = (ds["pds"] + ds["mse"]) / 6
    de = (ds["fid"] + ds["reach"] + ds["jac"] + ds["nmae"]) / 6
    O = "BAD" if (bulk < -0.0060 or abs(de) > 0.009 or dov < -0.02) else "OK"
tag = {"P2g": "498p2g", "S_P2": "498sp2"}
s1, s2 = [], []
if P != "NONE" and C != "FAILS":
    s1.append(tag[P])
if O != "BAD" and C != "FAILS" and P == "S_P2":
    s2.append("497dr")
fall = []
if C != "FAILS":
    other = {"P2g": ("S_P2", bs), "S_P2": ("P2g", bp)}.get(P)
    if P == "NONE":
        # no pick: a combined form still reads the cap-6 x per-target stack if its own 479 bin was at least PARTIAL
        for nm, bb in (("S_P2", bs), ("P2g", bp)):
            if bb is not None and bb[0] in ("TRANSFERS", "PARTIAL"):
                fall.append(tag[nm])
    elif other[1] is not None and other[1][0] in ("TRANSFERS", "PARTIAL"):
        fall.append(tag[other[0]])
fall += ["460dr"]   # 2026-09-30 06:38Z amendment (before any 10-01..10-04 entry): 461dr/469dr/458e carry the Arc programme,
                     # which read HARM on 09-30 (458d -0.0085, fid -0.0055 raw); they no longer fill a slot
c1 = s1 + [f for f in fall if f not in s1 + s2]
c2 = s2 + [f for f in fall if f not in s1 + s2]
print(f"{' '.join(c1) or '-'}|{' '.join(c2) or '-'}|P={P},C={C},O={O}")
