"""The 10-01 rule as REVISED 2026-09-26 (docs/NEXT-SLOT-DECISION.md), before any 09-27 entry exists.
Inputs: 09-27 s1 (446bj2) -> key_ok; 09-30 s1 (458d / 458) -> P; 09-30 s2 (458g) -> G. Prints "<candidates> | <detail>".
  key_ok : 446bj2 - 446b overall >= +0.001 and pds >= 0 (433's PREREGISTRATION-0926)
  P      : 458d - 446b by fid raw: positive >= +0.0015 (positive_reach if also reach <= -0.005 raw and overall <= 0),
           partial +0.0005 .. +0.0015, else null (458's bins + addendum 11:31Z); none if no 458d/458 entry
  G      : 458g - 458d (468's bins): bug if pds/mse differ > 1e-6 raw or |jac| > 0.0010; pos if overall >= +0.0005;
           null if |overall| < 0.0005; harm if <= -0.0005; none if no 458g entry (or s1 was not 458d)
  carried: 458g - 446b fid raw >= +0.0015 and overall >= 0 (468's table, row 3)
  candidates, first two that are ready (AMENDED 2026-09-27 07:59Z, before any 09-30 entry: 487 §2 / 468's addendum,
  October keeps le1 unless the gate reads HARM, so a NULL gate read now selects the gated October path):
    P in {positive, positive_reach, partial} and G in {pos, null}  -> 469dr, (458e unless P == positive), 460dr (if key_ok)
    P == positive and G in {harm, none, bug}                       -> 461dr, 458e, 460dr (if key_ok)
    P in {positive_reach, partial} and G in {harm, none, bug}      -> 458e, 460dr (if key_ok)
    P in {null, none} and G == pos and carried                     -> 469dr, 460dr (if key_ok)
    otherwise                                                      -> 460dr (if key_ok)
  AMENDED 2026-09-29 10:56Z (before any 09-30 entry exists): when P is null or none and 458d/458 did not read HARM
  (fid - 446b > -0.0010), 446bj2zmz (527: 446bj2z + the Arc-sourced uncalled-gene hook, 517 arm HK-N) is appended as
  the LAST candidate: 517 RESULTS §6 makes it arm 1 in the null branch, where the list above leaves slot 2 empty.
usage: uv run python rule1001b.py <id27|none> <id30s1|none> <tag30s1> <id30s2|none> <tag30s2>   (BOARD_DIR overrides data/board-entries)
"""
import json
import os
import sys
from pathlib import Path

BOARD = Path(os.environ.get("BOARD_DIR", "data/board-entries"))
FID, REACH, PDS, MSE, JAC = ("de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw", "pds_cosine",
                             "expr_mse_unbiased_capped_norm", "de_wilcoxon_sig_jaccard")


def find(x):
    if isinstance(x, dict):
        if "score_avg" in x:
            return x
        for v in x.values():
            y = find(v)
            if y:
                return y
    if isinstance(x, list):
        for v in x:
            y = find(v)
            if y:
                return y


def load(i):
    p = BOARD / f"{i}.json"
    return find(json.loads(p.read_text())) if i != "none" and p.exists() else None


id27, id1, tag1, id2, tag2 = sys.argv[1:6]
b = load("UDgjOKL9gI1udAexq6ue")
k, a, g = load(id27), load(id1), load(id2)
key_ok = bool(k) and (k["score_avg"] - b["score_avg"]) >= 0.001 and (k[PDS] - b[PDS]) >= 0
P, det = "none", ["no programme entry"]
if a and tag1 in ("458d", "458"):
    dfid, dreach, dov = a[FID] - b[FID], a[REACH] - b[REACH], a["score_avg"] - b["score_avg"]
    P = ("positive_reach" if (dreach <= -0.005 and dov <= 0) else "positive") if dfid >= 0.0015 else \
        "partial" if dfid >= 0.0005 else "null"
    det = [f"{tag1}-446b dfid {dfid:+.5f} dreach {dreach:+.5f} dov {dov:+.7f}"]
G, carried = "none", False
if g and tag2 == "458g" and a and tag1 == "458d":
    dov = g["score_avg"] - a["score_avg"]
    bug = abs(g[PDS] - a[PDS]) > 1e-6 or abs(g[MSE] - a[MSE]) > 1e-6 or abs(g[JAC] - a[JAC]) > 0.0010
    G = "bug" if bug else "pos" if dov >= 0.0005 else "harm" if dov <= -0.0005 else "null"
    carried = (g[FID] - b[FID]) >= 0.0015 and (g["score_avg"] - b["score_avg"]) >= 0
    det.append(f"458g-458d dov {dov:+.7f} dfid {g[FID] - a[FID]:+.5f} dreach {g[REACH] - a[REACH]:+.5f}"
               f" (bug {bug}); 458g-446b dfid {g[FID] - b[FID]:+.5f} dov {g['score_avg'] - b['score_avg']:+.7f}")
if P in ("positive", "positive_reach", "partial") and G in ("pos", "null"):
    cands = ["469dr"] + ([] if P == "positive" else ["458e"])
elif P == "positive":
    cands = ["461dr", "458e"]
elif P in ("positive_reach", "partial"):
    cands = ["458e"]
elif G == "pos" and carried:
    cands = ["469dr"]
else:
    cands = []
if key_ok:
    cands.append("460dr")
# 2026-09-29 amendment: the Arc-sourced hook arm fills the null branch (not after a HARM 458d read)
harm = bool(a) and tag1 in ("458d", "458") and (a[FID] - b[FID]) <= -0.0010
if P in ("null", "none") and not harm:
    cands.append("446bj2zmz")
print(" ".join(cands) if cands else "-", f"| key_ok {key_ok} P {P} G {G} carried {carried} :: " + " ; ".join(det))
