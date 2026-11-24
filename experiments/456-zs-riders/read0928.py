"""Read the 09-28 entries against their preregistered bins (written 2026-09-26 02:10Z, before any 09-28 entry exists).

- dress : 453 (the October path, blind, on A/B/C) vs 446b. Bins: experiments/453-abc-dress-rehearsal/PREREGISTRATION.md
          (bulk channel = pds + mse in overall units; DE channel = fid + reach + jac + nmae in overall units).
- zs    : the zs rider vs its base. Bins: experiments/433-coarse-mde-arm/PREREGISTRATION-0926.md.

usage: uv run python experiments/456-zs-riders/read0928.py <kind> <entry> <base_entry> [<kind> <entry> <base_entry>]
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


args = sys.argv[1:]
while args:
    kind, e, base, args = args[0], args[1], args[2], args[3:]
    a, b = load(e), load(base)
    d = {m: a[RAW[m]] - b[RAW[m]] for m in RAW}
    ds = {m: a[f"score_{m}"] - b[f"score_{m}"] for m in RAW}
    dov = a["score_avg"] - b["score_avg"]
    print(f"\n== {kind} {e} vs {base}: overall {a['score_avg']:.7f}  d {dov:+.7f}")
    for m in RAW:
        print(f"  {m:5s} raw d {d[m]:+.6f}  scaled d {ds[m]:+.6f}  -> overall {ds[m] / 6:+.6f}")
    if kind == "dress":
        bulk = (ds["pds"] + ds["mse"]) / 6
        de = (ds["fid"] + ds["reach"] + ds["jac"] + ds["nmae"]) / 6
        print(f"  bulk channel (pds+mse) {bulk:+.5f}; DE channel {de:+.5f}; predicted overall -0.0024 (band -0.0045..-0.0005)")
        print("  BIN bulk:", "VALIDATED (-0.0045..-0.0005): the October path is priced right end to end" if -0.0045 <= bulk <= -0.0005
              else "MIS-PRICED (< -0.0060): attribute by local rebuilds (A/B/C operator, pinned key, 446b's factors)" if bulk < -0.0060
              else "BETTER than priced (> +0.0005): re-read 438's CD4 half" if bulk > 0.0005 else "between bins (read with the attribution)")
        print("  BIN DE:", "pinned (|d| <= 0.009)" if abs(de) <= 0.009 else "OUTSIDE +-0.009: diagnose the DEF carriers before October")
        print("  HARD:", "overall below -0.02: STOP the October plan until explained" if dov < -0.02 else "none")
    elif kind == "zs":
        n = d["nmae"]
        print("  BIN nmae:", "zs TRANSFERS: it rides the October stack" if n <= -0.001 else
              "below readability (-0.001 < d <= 0): keep only if <= 0" if n <= 0 else "worse (d > 0): close zs")
        other = [m for m in ("pds", "mse", "fid", "jac") if abs(d[m]) > (1e-6 if m in ("pds", "mse") else 0.001)]
        print("  other members:", "pinned" if not other else f"MOVED {other}: diagnose (zs is a pure-nmae pass)")
    else:
        raise SystemExit(f"unknown kind {kind}")
