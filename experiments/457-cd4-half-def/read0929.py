"""Read the 09-29 entries against their preregistered bins (reader written 2026-09-28 14:10Z, before any 09-29 entry
exists; the bins themselves are older and unchanged).

- cd4 : 457dr (the October path with the stack at CD4 x1.5) vs 453dr (q3VuTImh4FRGKXH2waCp).
        Bins: experiments/457-cd4-half-def/PREREGISTRATION.md (2026-09-26 02:45Z).
- dld1: 446bd2 (446b + DLD-1 as the bulk's second union residual, g 0.2) vs 446b (UDgjOKL9gI1udAexq6ue).
        Bins: experiments/445-dld1-bulk-union/PREREGISTRATION.md (2026-09-24 21:40Z).

usage: uv run python experiments/457-cd4-half-def/read0929.py cd4 <457dr entry> [dld1 <446bd2 entry>]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
BASE = {"cd4": "q3VuTImh4FRGKXH2waCp", "dld1": "UDgjOKL9gI1udAexq6ue"}


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


args = sys.argv[1:]
while args:
    kind, e, args = args[0], args[1], args[2:]
    a, b = load(e), load(BASE[kind])
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {kind} {e} vs {BASE[kind]}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    if kind == "cd4":
        leak = [m for m, lim in (("nmae", 0.0015), ("fid", 0.001), ("jac", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
        print("  DE leak check:", "pinned" if not leak else f"OUTSIDE {leak}: a leak into the carriers -- diagnose")
        if dov >= 0.0040:
            print("  BIN: >= +0.0040: 438 UNDER-reads CD4; re-price x1.25 as an October hedge before Oct 22")
        elif dov >= 0.0010:
            print("  BIN: >= +0.0010: the CD4 half transfers as 438 prices it; October keeps CD4 x1 (438's panel reading stands)")
        elif dov > -0.0010:
            print("  BIN: |d| < 0.0010: flat on the board; October keeps CD4 x1, with more reason")
        else:
            print("  BIN: <= -0.0010: x1.5 loses even on A/B/C; October keeps CD4 x1; re-examine 438's CD4 residual rows")
        print("  predicted +0.0017 (k 1) .. +0.0022 (k 1.3), band +0.0005 .. +0.0035")
    elif kind == "dld1":
        leak = [m for m, lim in (("nmae", 0.001), ("fid", 0.001), ("jac", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
        print("  DE leak check:", "pinned" if not leak else f"OUTSIDE {leak}: a bulk-to-carrier leak -- STOP")
        if dov >= 0.0005 and d["mse"] <= -0.001:
            print("  BIN: DLD-1 TRANSFERS (overall >= +0.0005 with mse <= -0.001 raw): it joins the bulk and the October stack")
        elif abs(dov) < 0.0005:
            print("  BIN: flat (|d| < 0.0005): keep DLD-1 out (simpler October recipe)")
        elif dov <= -0.0005:
            print("  BIN: <= -0.0005: close DLD-1")
        else:
            print("  BIN: overall >= +0.0005 but mse above -0.001 raw: not the preregistered transfer signature -- read by hand")
        print("  predicted overall +0.0005 .. +0.0015; mse -0.001 .. -0.004 raw; pds -0.0005 .. +0.002 raw")
    else:
        raise SystemExit(f"unknown kind {kind}")
