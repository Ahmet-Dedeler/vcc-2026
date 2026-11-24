"""Read 10-02's per-target sharpening pair against 479's preregistered bins (PREREGISTRATION.md, b692c2c5), written
2026-09-27 before any 479 arm exists.
usage: uv run python experiments/479-pertarget-sharpen/read1002.py <P2g entry> <S_P2 entry> [<446bj2 entry>]
Either entry may be "none" (the arm did not run). Prints each arm vs 446bj2, its bin, the ballast and leak checks, and
the pair's October pick (larger dpds; tie within 0.0005 raw -> S_P2; one transfer -> it; none -> closed).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
BJ2 = "YvVru3DAF89Jiyr34olR"


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


def arm(tag, e, b):
    if e == "none":
        print(f"\n== {tag}: did not run"); return None
    a = load(e)
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {tag} {e}: overall {a['score_avg']:.7f}  d {dov:+.7f} vs 446bj2")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    if dov >= 0.0008 and d["pds"] >= 0.0015:
        b_ = "TRANSFERS"
    elif 0.0003 <= dov < 0.0008 or 0.0005 <= d["pds"] < 0.0015:
        b_ = "PARTIAL"
    else:
        b_ = "FAILS"
    print(f"  BIN: {b_}  (transfers: overall >= +0.0008 and dpds >= +0.0015; partial: overall +0.0003..+0.0008 or dpds +0.0005..+0.0015)")
    n = d["nmae"]
    print("  ballast check: dnmae", f"{n:+.5f}", "expected +0.0006..+0.0010" if n <= 0.0015 else "ABOVE +0.0015: the low tercile's ballast is under-priced (record it)")
    leak = [m for m, lim in (("fid", 0.001), ("jac", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
    print("  leak check:", "clean" if not leak else f"DIAGNOSE {leak}")
    return b_, d["pds"]


b = load(sys.argv[3] if len(sys.argv) > 3 else BJ2)
p = arm("P2g (tercile:2.5,2.0,1.75)", sys.argv[1], b)
s = arm("S_P2 (stercile:2.5,2.0,1.75)", sys.argv[2], b)
tp = p is not None and p[0] == "TRANSFERS"; ts = s is not None and s[0] == "TRANSFERS"
if tp and ts:
    pick = "S_P2 (tie within 0.0005 raw)" if abs(p[1] - s[1]) < 0.0005 else ("P2g" if p[1] > s[1] else "S_P2")
elif tp or ts:
    pick = "P2g" if tp else "S_P2"
else:
    pick = "NONE: per-target sharpening is closed on the board" + (" (a PARTIAL arm is recorded, not shipped)" if any(x and x[0] == "PARTIAL" for x in (p, s)) else "")
print(f"\nOCTOBER PICK (479's pair rule): {pick}")
print("If P2g is picked, October needs 493h's carried-only cut points before DEF_DAY_SHARPEN_PT=tercile:... is used.")
