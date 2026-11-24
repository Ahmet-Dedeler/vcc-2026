#!/usr/bin/env python3
"""803 sanity table: per pair (arms/pairs.json), surface and regime, the mean raw member delta over seeds (paired by seed)
and its seed SD, from the collected members.npz files only (light; no board units, no bins: check800 does that).
Also pds and mse raw deltas from the scorer's own per-target values (ce_pds mean; mse = sum num / sum den).
    uv run python experiments/803-de-engine/pairs803.py <out dir> [--md work/pairs803.md] [--json work/pairs803.json]
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REGS = ("f1", "f0.61", "f0.5", "f0.35", "cmp")
MEM = ("nmae", "fid", "jac", "reach")
CE = ("pds", "nmae", "fid", "jac", "reach")   # the scorer's own f 1 members (ce_*), beside 521's


def means(z):
    out = {f"{m}_{r}": float(np.nanmean(z[f"{m}_{r}"])) for m in MEM for r in REGS}
    out["pds"] = float(np.nanmean(z["ce_pds"]))
    for m in CE[1:]:
        out[f"ce_{m}"] = float(np.nanmean(z[f"ce_{m}"]))
    out["mse"] = float(np.nansum(z["mse_num"]) / np.nansum(z["mse_den"]))
    out["K"] = float(np.mean(z["K"]))
    out["ncall"] = float(np.mean(z["ncall"]))
    return out


def main():
    root = Path(sys.argv[1])
    md_path = Path(sys.argv[sys.argv.index("--md") + 1]) if "--md" in sys.argv else HERE / "work/pairs803.md"
    js_path = Path(sys.argv[sys.argv.index("--json") + 1]) if "--json" in sys.argv else HERE / "work/pairs803.json"
    pf = Path(sys.argv[sys.argv.index("--pairs") + 1]) if "--pairs" in sys.argv else HERE / "arms/pairs.json"
    pairs = json.loads(pf.read_text())["pairs"]
    res, md = {}, ["# 803 pair deltas (raw, arm − base, mean over seeds paired by seed; ± = seed SD)\n"]
    for sf in ("h1test", "h1tv"):
        for p in pairs:
            if sf not in p.get("surfaces", ["h1test", "h1tv"]):
                continue
            ds = []
            for s in (0, 1, 2):
                fa, fb = root / sf / p["arm"] / f"s{s}/members.npz", root / sf / p["base"] / f"s{s}/members.npz"
                if fa.exists() and fb.exists():
                    a, b = means(np.load(fa)), means(np.load(fb))
                    ds.append({k: a[k] - b[k] for k in a})
            if not ds:
                continue
            keys = list(ds[0])
            rec = {k: dict(mean=float(np.mean([d[k] for d in ds])), sd=float(np.std([d[k] for d in ds], ddof=1)) if len(ds) > 1 else None,
                           per_seed=[d[k] for d in ds]) for k in keys}
            res[f"{sf}|{p['name']}"] = dict(pair=p, seeds=len(ds), deltas=rec)
    js_path.write_text(json.dumps(res, indent=1))
    for sf in ("h1test", "h1tv"):
        md += [f"\n## {sf}\n", "| pair (arm − base) | seeds | regime | Δnmae | Δfid | Δjac | Δreach | board |", "|---|---|---|---|---|---|---|---|"]
        for k, v in res.items():
            if not k.startswith(sf + "|"):
                continue
            p, d = v["pair"], v["deltas"]
            for r in REGS:
                cells = []
                for m in MEM:
                    x = d[f"{m}_{r}"]
                    cells.append(f"{x['mean']:+.5f}" + (f" ± {x['sd']:.5f}" if x["sd"] is not None else ""))
                bd = ", ".join(f"{m} {val:+.5f}" for m, val in (p.get("board") or {}).items()) if r == "f1" else ""
                md.append(f"| {p['name']} ({p['arm']} − {p['base']}) | {v['seeds']} | {r} | " + " | ".join(cells) + f" | {bd} |")
            sd = lambda k: f" ± {d[k]['sd']:.2e}" if d[k]["sd"] is not None else ""
            bd = p.get("board") or {}
            md.append(f"| ↳ bulk / K | | | Δpds {d['pds']['mean']:+.2e}{sd('pds')} | Δmse {d['mse']['mean']:+.2e}{sd('mse')} | ΔK {d['K']['mean']:+.3f} | "
                      f"Δcalls {d['ncall']['mean']:+.1f} | " + (f"pds {bd['pds']:+.6f}, mse {bd['mse']:+.6f}" if "pds" in bd else "") + " |")
    md_path.write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
