#!/usr/bin/env python3
"""553 — the rehearsal's between-draw spread vs the access's within-surface band. Reads the access553.json of N
R7b-shaped rehearsal draws (h1r7, h1r7s2, ...: the same shape, independent controls + cells) and, per arm (vs `cpm`) and
per switch-on-recipe pair: the mean effect over draws, the between-draw sd (ddof 1), the mean within-surface band
(cells only: `band_sd`; cells + controls: `band_sd_ctrl` when the draw ran csplit) and the ratio between-draw sd / band
-- the factor by which a single surface's band understates how far a reading moves when the controls and cells are drawn
again. With N draws the sd has N-1 degrees of freedom (a chi-square 90 % interval is printed for the ratio).

usage: uv run python experiments/553-r7b-harness/draws553.py h1r7 h1r7s2 h1r7s3 ... [--out work/draws553.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent


def load(name: str) -> dict:
    return json.loads((HERE / "work" / f"access553_{name}" / "access553.json").read_text())


def rows_of(rec: dict) -> dict:
    out = {}
    for arm, row in rec["arms"].items():
        if arm != "cpm" and row["effect_scaled"]:
            out[f"{arm} vs cpm"] = row["effect_scaled"]["overall"]
    for name, pr in rec.get("pairs", {}).items():
        out[name] = pr["effect_scaled"]["overall"]
    return out


def main():
    argv = sys.argv[1:]
    args = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or argv[i - 1] not in ("--out", "--mean-out"))]
    out_path = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else HERE / "work" / "draws553.json"
    draws = {n: rows_of(load(n)) for n in args}
    keys = [k for k in draws[args[0]] if all(k in d for d in draws.values())]
    res = {}
    print(f"| effect | draws | mean | between-draw sd | mean band (cells) | ratio | mean band (cells + controls) | ratio | "
          f"ratio 90 % interval (cells + controls) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for k in keys:
        v = np.array([draws[n][k]["mean"] for n in args])
        b = np.array([draws[n][k]["band_sd"] for n in args])
        bc = [draws[n][k].get("band_sd_ctrl") for n in args]
        bc = np.array([np.nan if x is None else x for x in bc], dtype=float)
        n = len(v)
        sd = float(v.std(ddof=1)) if n > 1 else float("nan")
        # 90 % interval of the sd from its chi-square(n-1) sampling law
        lo, hi = (float(sd * np.sqrt((n - 1) / chi2.ppf(.95, n - 1))), float(sd * np.sqrt((n - 1) / chi2.ppf(.05, n - 1)))) if n > 1 else (np.nan, np.nan)
        mb, mbc = float(np.mean(b)), float(np.nanmean(bc)) if np.isfinite(bc).any() else float("nan")
        res[k] = dict(values=v.tolist(), mean=float(v.mean()), between_sd=sd, between_sd_90=[lo, hi], band_mean=mb,
                      band_ctrl_mean=mbc, ratio=sd / mb, ratio_ctrl=sd / mbc if np.isfinite(mbc) else None,
                      ratio_ctrl_90=[lo / mbc, hi / mbc] if np.isfinite(mbc) else None, draws=args)
        rc = "" if not np.isfinite(mbc) else f"{sd / mbc:.2f}"
        rci = "" if not np.isfinite(mbc) else f"[{lo / mbc:.2f}, {hi / mbc:.2f}]"
        print(f"| {k} | {n} | {v.mean():+.5f} | {sd:.5f} | {mb:.5f} | {sd / mb:.2f} | "
              f"{'' if not np.isfinite(mbc) else f'{mbc:.5f}'} | {rc} | {rci} |")
    out_path.write_text(json.dumps(res, indent=1))
    if "--mean-out" in sys.argv:
        # the rehearsal's mean raw member effects in access553.json's layout, so `access553.py --retro <dir>` can use the
        # shape-matched rehearsal (not H1-full) as the reference for R7b's board-retrodiction ratios
        recs = [load(n) for n in args]
        mean_rec = {"arms": {}, "pairs": {}, "draws": args}
        for arm in recs[0]["arms"]:
            if arm != "cpm" and all(arm in r["arms"] and r["arms"][arm]["effect"] for r in recs):
                mean_rec["arms"][arm] = {"effect": {m: {"mean": float(np.mean([r["arms"][arm]["effect"][m]["mean"] for r in recs]))}
                                                    for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}}
        for k in recs[0].get("pairs", {}):
            if all(k in r.get("pairs", {}) for r in recs):
                mean_rec["pairs"][k] = {"effect": {m: {"mean": float(np.mean([r["pairs"][k]["effect"][m]["mean"] for r in recs]))}
                                                   for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}}
        d = HERE / sys.argv[sys.argv.index("--mean-out") + 1]
        d.mkdir(parents=True, exist_ok=True)
        (d / "access553.json").write_text(json.dumps(mean_rec, indent=1))
        print(f"rehearsal mean -> {d / 'access553.json'}")


if __name__ == "__main__":
    main()
