#!/usr/bin/env python3
"""553 — the control re-draw band of a surface, from its CONTROL HALVES (`prep553.py chalf <base> <k>`, each half's
prediction side rebuilt and emitted again, each half accessed against the base's target cells with the same half of the
controls). Per split k the two halves are independent draws of half the controls (prediction AND truth side), so the
difference of an effect between them, d_k = E(kA) - E(kB), has variance 2 x the half-control variance = 4 x the
full-control variance (variance ~ 1 / n_controls): the control component of the base's band is sigma_c = rms(d_k) / 2.
The base's own band (emission seeds + split-half re-draws of the target cells, `band_sd`) does not contain it; the total
is band_total = sqrt(band_sd^2 + sigma_c^2). Each half carries one emission seed, so sigma_c also holds some seed noise
(corrected with the base's single-seed sd of the same effect; the raw value is kept). If a draws553.json is given (the rehearsal's between-draw spread), the comparison is printed.

usage: uv run python experiments/553-r7b-harness/band553.py <base> --splits 0,1,2 [--draws work/draws553.json]
writes work/band553_<base>.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def effects(rec: dict) -> dict:
    """overall-equivalent and raw-member effects of every arm vs cpm and every pair, keyed like draws553."""
    out = {}
    for arm, row in rec["arms"].items():
        if arm != "cpm" and row["effect_scaled"]:
            out[f"{arm} vs cpm"] = row
    for name, pr in rec.get("pairs", {}).items():
        out[name] = pr
    return out


def main():
    base = sys.argv[1]
    splits = [int(x) for x in sys.argv[sys.argv.index("--splits") + 1].split(",")]
    dp = Path(sys.argv[sys.argv.index("--draws") + 1]) if "--draws" in sys.argv else None
    draws = json.loads((dp if dp.is_absolute() else HERE / dp).read_text()) if dp is not None else {}
    tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else ""
    W = HERE / "work"
    b = effects(json.loads((W / f"access553_{base}{tag}" / "access553.json").read_text()))
    halves = {(k, h): effects(json.loads((W / f"access553_{base}c{k}{h}{tag}" / "access553.json").read_text()))
              for k in splits for h in "AB"}
    res = {}
    print(f"| effect | {base} effect | band_sd (seeds + cells) | control re-draw sigma_c ({len(splits)} splits; seed-corrected / raw) | "
          f"band_total | effect / band_total | rehearsal between-draw sd | between / band_total |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for k_, row in b.items():
        miss = [f"{k}{h}" for k in splits for h in "AB" if k_ not in halves[(k, h)]]
        if miss:   # 555: recorded, never dropped silently: no band_total, so this effect cannot carry a band-based rule
            res[k_] = dict(effect=row["effect_scaled"]["overall"]["mean"], band_sd=row["effect_scaled"]["overall"].get("band_sd"),
                           band_total=None, no_band=f"not on the control halves {miss}")
            print(f"| {k_} | {row['effect_scaled']['overall']['mean']:+.5f} | NO BAND: not on the halves {miss} | | | | | |")
            continue
        d = np.array([halves[(k, "A")][k_]["effect_scaled"]["overall"]["mean"]
                      - halves[(k, "B")][k_]["effect_scaled"]["overall"]["mean"] for k in splits])
        dm = {m: [halves[(k, "A")][k_]["effect"][m]["mean"] - halves[(k, "B")][k_]["effect"][m]["mean"] for k in splits]
              for m in ("nmae", "fid", "reach", "jac")}
        o = row["effect_scaled"]["overall"]
        band = o.get("band_sd")
        # each half carries ONE emission seed, so E[d^2] = 4 sigma_c^2 + 2 seed_sd^2; the base's single-seed sd of the same
        # paired effect estimates the second term (sigma_c_raw keeps it: an upper bound)
        sig_c_raw = float(np.sqrt(np.mean(d ** 2)) / 2)
        ssd = o.get("seed_sd") or 0.0
        sig_c = float(np.sqrt(max(np.mean(d ** 2) - 2 * ssd ** 2, 0.0)) / 2)
        tot = float(np.sqrt(band ** 2 + sig_c ** 2)) if band is not None else None
        dr = draws.get(k_, {})
        res[k_] = dict(effect=o["mean"], band_sd=band, d_halves=d.tolist(), d_members=dm, sigma_c=sig_c,
                       sigma_c_raw=sig_c_raw, seed_sd_used=ssd, band_total=tot,
                       between_draw_sd=dr.get("between_sd"), between_draw_n=len(dr.get("values", [])))
        bd = dr.get("between_sd")
        band_s = "" if band is None else f"{band:.5f}"
        tot_s = "" if tot is None else f"{tot:.5f}"
        z_s = "" if not tot else f"{o['mean'] / tot:+.1f}"
        bd_s = "" if bd is None else f"{bd:.5f}"
        rat_s = "" if (bd is None or not tot) else f"{bd / tot:.2f}"
        print(f"| {k_} | {o['mean']:+.5f} | {band_s} | {sig_c:.5f} / {sig_c_raw:.5f} | {tot_s} | {z_s} | {bd_s} | {rat_s} |")
    (W / f"band553_{base}{tag}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
