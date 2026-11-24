"""Read the 09-26 entries against their preregistered bins: 446b vs 433k (446's PREREGISTRATION.md) and the slot-2 arm
(446bja vs 446b: 450/PREREGISTRATION-446bja.md; or the fallback key rider vs its reference: 433/PREREGISTRATION-0926.md).

usage: uv run python experiments/446-q175-level/read0926.py <446b_entry> [<slot2_entry> <slot2_tag>]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
REF_433K = "T76JWD3tOajk5IYvzwUW"


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


def show(tag, a, b):
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {tag}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    return d, dov


b446 = load(sys.argv[1])
k = load(REF_433K)
d, dov = show("446b vs 433k", b446, k)
p, ms = d["pds"], d["mse"]
print("  BIN pds:", "TRANSFERS (>= +0.002)" if p >= 0.002 else "partial (+0.0005..+0.002)" if p >= 0.0005 else
      "flat (|d| < 0.0005)" if p > -0.001 else "OVERSHOOTS (<= -0.001): q stays 1.5")
print("  BIN mse:", "TRANSFERS (<= -0.005)" if ms <= -0.005 else "partial (-0.005..-0.001)" if ms <= -0.001 else
      "flat (-0.001..0)" if ms <= 0 else "OVERSHOOT or mis-level (> 0): re-read ntclevel's rho")
leak = [m for m, lim in (("nmae", 0.001), ("fid", 0.001), ("jac", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
print("  BIN DE members:", "pinned" if not leak else f"OUTSIDE ballast noise {leak}: stop and diagnose")
print("  BIN overall:", "q 1.75 becomes the default (A/B/C and October)" if (dov >= 0.001 and p >= 0) else
      "flat: keep q 1.5" if abs(dov) < 0.001 else "close the q curve at 1.5" if dov <= -0.001 else "partial/ambiguous")
if len(sys.argv) > 3:
    s2, tag = load(sys.argv[2]), sys.argv[3]
    if tag == "446bja":
        d2, dov2 = show("446bja vs 446b", s2, b446)
        f1 = abs(d2["pds"]) <= 1e-6 and abs(d2["mse"]) <= 1e-6
        f2 = abs(d2["jac"]) <= 0.001
        print(f"  FALSIFIERS: pds/mse pinned {f1}; jac within 0.001 {f2}" + ("" if (f1 and f2) else "  -> READING VOID"))
        f = d2["fid"]
        print("  BIN fid (primary):", "TRANSFERS (>= +0.0015)" if f >= 0.0015 else "partial (+0.0005..+0.0015)" if f >= 0.0005
              else "no transfer (|d| < 0.0005)" if f > -0.001 else "HARM (<= -0.0010)")
        n = d2["nmae"]
        print("  BIN nmae:", "confirms (<= -0.0015)" if n <= -0.0015 else "neutral" if n < 0.0015 else "contradicts fid (>= +0.0015)")
        print(f"  reach (reported, not binned): {d2['reach']:+.6f}")
        print("  BIN overall:", "446bja becomes the incumbent (>= +0.0010)" if dov2 >= 0.001 else
              "flat" if dov2 > -0.001 else "harm (<= -0.0010)")
    else:
        ref = b446 if tag == "446bj2" else k
        show(f"{tag} vs {'446b' if tag == '446bj2' else '433k'}", s2, ref)
