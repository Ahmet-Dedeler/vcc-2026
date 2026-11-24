"""Read 458g against its preregistered bins (written 2026-09-26 19:56Z, before any 09-27…09-30 entry exists; bins in
experiments/468-prior-gate/PREREGISTRATION.md, including the addenda: the three-way split and the combined October table).

  uv run python experiments/468-prior-gate/read468.py <458g entry> <458d entry> [<446b entry>]   (446b defaults to 446b's id)

Prints 458g − 458d (primary), its falsifiers and mechanism, the exact three-way split (458g − 446b = the programme on the
design-strong two thirds; 458d − 458g = on the weak third; valid because compare458g.py found every 458g block equal to
458d's or 446b's), 458d's own P by 458's bins, and the October decision from the combined table.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
B446 = "UDgjOKL9gI1udAexq6ue"


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
    return _find(json.loads((ROOT / f"data/board-entries/{i}.json").read_text()))


def delta(a, b, tag):
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {tag}: overall {a['score_avg']:.7f} vs {b['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    return d, dov


g, dd, b = load(sys.argv[1]), load(sys.argv[2]), load(sys.argv[3] if len(sys.argv) > 3 else B446)
d, dov = delta(g, dd, "458g vs 458d (the gate: the weak third's prior, Arc -> CLEAN4)")
bad = [m for m in ("pds", "mse") if abs(d[m]) > 1e-6]
jac_bad = abs(d["jac"]) > 0.0010
print("  FALSIFIERS:", "OK (pds/mse identical to 1e-6, jac within +-0.0010)" if not bad and not jac_bad
      else f"FIRED {bad} {'jac' if jac_bad else ''}: a bug; diagnose before reading anything else")
G = "bug" if (bad or jac_bad) else "pos" if dov >= 0.0005 else "harm" if dov <= -0.0005 else "null"
print("  PRIMARY:", {"pos": ("TRANSFERS (>= +0.0015): October sets DEF_DAY_PRIOR_GATE=q33 with any programme source" if dov >= 0.0015
                             else "WEAK POSITIVE (+0.0005 .. +0.0015): keep the gate, flagged unconfirmed"),
                     "null": "NULL (|d| < 0.0005): keep le1 with the programme (amended 09-27, 487 section 2; only HARM closes the gate)",
                     "harm": "HARM (<= -0.0005): close the gate", "bug": "BUG"}[G])
print("  MECHANISM:", "the H1 mechanism transfers (reach up, fid down)" if d["reach"] >= 0.003 and d["fid"] <= 0 else
      "the gate only removes signal (reach and fid both down)" if d["reach"] <= 0 and d["fid"] < 0 else
      f"mixed (d reach {d['reach']:+.5f}, d fid {d['fid']:+.5f})")
s, sov = delta(g, b, "SPLIT strong: 458g vs 446b = the programme on the design-strong two thirds")
w, wov = delta(dd, g, "SPLIT weak: 458d vs 458g = the programme on the design-weak third")
p, pov = delta(dd, b, "458d vs 446b (458's bins; fid primary)")
P = ("positive_reach" if (p["reach"] <= -0.005 and pov <= 0) else "positive") if p["fid"] >= 0.0015 else \
    "partial" if p["fid"] >= 0.0005 else "harm" if p["fid"] <= -0.0010 else "null"
carried = s["fid"] >= 0.0015 and sov >= 0
print(f"\n  458d P = {P};  G = {G};  458g - 446b carried (fid >= +0.0015 raw, overall >= 0): {carried}")
if G == "bug":
    dec = "none yet: a falsifier fired, diagnose first"
elif P in ("positive", "positive_reach", "partial"):
    # 468 addendum 2026-09-27 07:58Z (487 section 2): a NULL gate read keeps le1; only HARM closes the gate.
    dec = "programme source + gate (le1)" if G in ("pos", "null") else "programme source, no gate"
else:
    dec = "programme source + gate (the gated construction carries; the ungated one does not)" if carried else "no programme"
print(f"  OCTOBER (468's combined table): {dec}")
print("  10-01 (rule1001b.py) reads the same inputs at run time; this reader is for the record.")
