"""Read the 09-25 entries against PREREGISTRATION-0925.md's bins (433h vs 433a; 433k vs 433h) and 446's gate.

usage: uv run python experiments/433-coarse-mde-arm/read0925.py <433h_entry> [<433k_entry>]
Prints raw / scaled / overall deltas per member (board_delta arithmetic) and the bin each preregistered reading falls in.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import json  # noqa: E402

RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
SC = {m: f"score_{m}" for m in RAW}


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

REF_433A = "FCbBcvGxF1HOFPxPISUC"


def deltas(a, b):
    return {m: (a[RAW[m]] - b[RAW[m]], a[SC[m]] - b[SC[m]]) for m in RAW}, a["score_avg"] - b["score_avg"]


def show(tag, d, dov, a):
    print(f"\n== {tag}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m, (dr, ds) in d.items():
        print(f"  {m:5s} raw d {dr:+.6f}  scaled d {ds:+.6f}  -> overall {ds / 6:+.6f}")


def bin_433h(d, dov):
    nm, ms = d["nmae"][0], d["mse"][0]
    out = []
    out.append("nmae: " + ("REPAIR MATTERED (<= -0.002)" if nm <= -0.002 else "neutral (-0.002..+0.001)" if nm <= 0.001
                           else "STOP AND DIAGNOSE (> +0.001)"))
    out.append("mse: " + ("RE-LEVEL TRANSFERS (<= -0.006)" if ms <= -0.006 else "partial (-0.006..-0.002)" if ms <= -0.002
                          else "NO TRANSFER (> -0.002): re-read the cliff model"))
    out.append(f"overall: {'inside' if 0.001 <= dov <= 0.005 else 'OUTSIDE'} the predicted +0.001..+0.005")
    return out


def bin_433k(d, dov):
    p, ms = d["pds"][0], d["mse"][0]
    out = []
    out.append("pds: " + ("DOSE STEP TRANSFERS (>= +0.004)" if p >= 0.004 else "partial (+0.001..+0.004)" if p >= 0.001
                          else "flat (< +0.001)"))
    out.append("mse: " + ("TRANSFERS (<= -0.010)" if ms <= -0.010 else "partial (-0.010..0)" if ms <= 0 else "OVERSHOOTS A's cliff (> 0)"))
    leak = [m for m, lim in (("nmae", 0.001), ("reach", 0.01), ("fid", 0.001), ("jac", 0.001)) if abs(d[m][0]) > lim]
    out.append("DE members: " + ("pinned" if not leak else f"OUTSIDE ballast noise: {leak}"))
    out.append(f"overall: {'inside' if 0.002 <= dov <= 0.006 else 'OUTSIDE'} the predicted +0.002..+0.006")
    go446 = p >= 0.001 or ms <= -0.004
    out.append(f"446b gate (pds >= +0.001 or mse <= -0.004 vs 433h): {'GO' if go446 else 'NO-GO (the q curve closes at 1.5)'}")
    return out


h = load(sys.argv[1])
a = load(REF_433A)
d, dov = deltas(h, a)
show("433h vs 433a", d, dov, h)
for line in bin_433h(d, dov):
    print("  BIN", line)
if len(sys.argv) > 2:
    k = load(sys.argv[2])
    d2, dov2 = deltas(k, h)
    show("433k vs 433h", d2, dov2, k)
    for line in bin_433k(d2, dov2):
        print("  BIN", line)
    d3, dov3 = deltas(k, a)
    print(f"\n  433k vs 433a overall d {dov3:+.7f}")
