"""Read 10-03's two bulk arms against their preregistered bins, written 2026-09-27 before any 495 / 471c entry exists.
- 495 (LODO key at cap 6) vs 446bj2 (YvVru3DAF89Jiyr34olR): 495/PREREGISTRATION.md (706ad5f0).
- 471c (h1_2025_full union swap, re-levelled) vs 446b (UDgjOKL9gI1udAexq6ue): 471/PREREGISTRATION.md + addenda.
usage: uv run python experiments/495-lodo-cap6/read1003.py <495 entry|none> <471c entry|none>
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
    return _find(json.loads((ROOT / f"data/board-entries/{i}.json").read_text()))


def delta(tag, e, ref, refname):
    a, b = load(e), load(ref)
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {tag} {e}: overall {a['score_avg']:.7f}  d {dov:+.7f} vs {refname}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    return d, dov


if sys.argv[1] != "none":
    d, dov = delta("495 (cap 6)", sys.argv[1], BJ2, "446bj2")
    leak = abs(d["nmae"]) > 0.002 or abs(d["fid"]) > 0.001 or abs(d["jac"]) > 0.001 or abs(d["reach"]) > 0.01
    print("  DE leak check:", "DIAGNOSE before reading" if leak else "pinned")
    if dov >= 0.0012 and d["mse"] <= -0.006:
        print("  BIN: TRANSFERS -> October: DEF_DAY_LODO_HI_PCT=96.64697522769474 (the percentile form, 543; == cap 6 on A/B/C)")
    elif dov < 0.0004 or d["mse"] > -0.002:
        print("  BIN: FAILS -> keep cap 3 (no DEF_DAY_LODO_HI_PCT)")
    else:
        print("  BIN: PARTIAL (+0.0004 .. +0.0012) -> keep cap 3 unless a second arm confirms (498's combined arm is one)")
    print("  prediction: overall +0.0019 .. +0.0029, mse -0.0140 raw")
if len(sys.argv) > 2 and sys.argv[2] != "none":
    d, dov = delta("471c (h1_2025_full union)", sys.argv[2], B446, "446b")
    leak = d["nmae"] > 0.0015 or abs(d["fid"]) > 0.001 or abs(d["jac"]) > 0.001 or abs(d["reach"]) > 0.01
    print("  DE leak check:", "DIAGNOSE before reading" if leak else "pinned")
    print("  BIN:", "ADOPT (h1_2025_full in the October stack union; needs its own price on the LODO stack)" if dov >= 0.0002 and d["mse"] <= 0
          else "HARM: keep the default union (noise reduction in a union member does not pay)" if dov <= -0.0002
          else "FLAT: keep the default union")
    print("  prediction: mse -0.0005 .. -0.0030 raw, pds -0.0005 .. +0.0015, overall +0.0001 .. +0.0006")
