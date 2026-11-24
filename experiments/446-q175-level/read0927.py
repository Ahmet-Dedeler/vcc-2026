"""Read the 09-27 riders against their preregistered bins (written 2026-09-26 00:50Z, before any 09-27 entry exists).

- key   : the LODO-key rider (446bj2 / 446bjaj2 / 433kj) vs its base. Bins: 433-coarse-mde-arm/PREREGISTRATION-0926.md.
- floor : the zero-repair floor -2 rider (446br2 / 446bjar2 / 433kr2) vs its base. Bins: 449-repair-floor/PREREGISTRATION.md.

usage: uv run python experiments/446-q175-level/read0927.py <kind> <entry> <base_entry> [<kind> <entry> <base_entry>]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}


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


args = sys.argv[1:]
while args:
    kind, e, base, args = args[0], args[1], args[2], args[3:]
    d, dov = show(f"{kind} rider {e} vs base {base}", load(e), load(base))
    if kind == "key":
        leak = [m for m, lim in (("nmae", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
        print("  BIN DE members:", "pinned" if not leak else f"OUTSIDE noise {leak} (ballast K may have moved; diagnose)")
        print("  BIN overall:", "LODO key TRANSFERS: it replaces 410's key in the October stack (needs its day-of fit)"
              if dov >= 0.001 and d["pds"] >= 0 else "flat: keep 410's key" if abs(dov) < 0.001
              else "close the LODO key" if dov <= -0.001 else "ambiguous (overall >= +0.001 but pds < 0): keep 410's key")
        print(f"  prediction: pds +0.001 .. +0.004 raw (got {d['pds']:+.4f}); mse -0.004 .. -0.012 (got {d['mse']:+.4f})")
    elif kind == "floor":
        n = d["nmae"]
        print("  BIN nmae:", "the shallower floor PAYS: ship -2 in the October stack" if n <= -0.002 else
              "flat: keep -5" if abs(n) < 0.001 else "deep gated truths: keep -5 (consider a cis-only deep floor)"
              if n > 0.001 else "between bins (-0.002 < d <= -0.001): partial; keep -5 unless overall >= +0.0006")
        calls = [m for m, lim in (("fid", 0.001), ("jac", 0.001), ("reach", 0.01)) if abs(d[m]) > lim]
        print("  BIN calls:", "unmoved" if not calls else f"MOVED {calls}: diagnose")
        bulk = [m for m in ("pds", "mse") if abs(d[m]) > 0.0005]
        print("  BIN bulk:", "pinned (|d| <= 0.0005)" if not bulk else f"MOVED {bulk}: the repair reached the bulk; diagnose")
        print(f"  prediction: overall +0.0006 .. +0.0026 (got {dov:+.4f}); nmae -0.002 .. -0.008 (got {n:+.4f})")
    else:
        raise SystemExit(f"unknown kind {kind}")
