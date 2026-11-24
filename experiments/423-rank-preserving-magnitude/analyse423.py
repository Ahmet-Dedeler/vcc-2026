#!/usr/bin/env python3
"""423 table: every post-pass variant vs its same-seed base (same emitted blocks), mean +- sd over seeds.

usage: uv run python analyse423.py [--bulk base|ps] [--out work/analysis423]
Board slopes per raw unit (416's mandate): pds +2.2138, mse -1.051, nmae -1.6532, fid +3.3645, reach +1.1292,
jac +2.7263; overall = sum(slope * delta) / 6.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = HERE / "work" / "h1"
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
MEM = ("pds", "mse", "nmae", "fid", "reach", "jac")
SUMK = ("frac_massw_mean", "frac_massw_up_mean", "frac_massw_dn_mean", "frac_massw_zs_mean", "frac_full_mean",
        "n_moved_mean", "n_touch_mean", "counts_moved_mean", "rolled_mean", "n_dU_ball_mean", "dU_ball_abs_max_mean",
        "n_genes_U_changed_mean", "n_touched_U_changed_mean", "n_untouched_U_changed_mean", "U_abs_max_mean",
        "call_flips_mean", "sign_flips_mean", "p_logmax_mean", "lim_ballast_subgate_mean")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", choices=["base", "ps"], default="base")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    sfx = "_bps" if a.bulk == "ps" else ""
    pat = re.compile(r"metrics-(.+)_s(\d+)" + sfx + r"\.json$")
    runs = defaultdict(dict)
    for f in sorted(W.glob("metrics-*.json")):
        m = pat.match(f.name)
        if not m or (not sfx and f.name.endswith("_bps.json")):
            continue
        runs[m.group(1)][int(m.group(2))] = json.loads(f.read_text())
    base = runs.pop("base", {})
    rows = []
    for v in sorted(runs, key=lambda x: (x.split("_")[-1] if x.count("_") > 1 else "pc", x)):
        seeds = sorted(s for s in runs[v] if s in base)
        if not seeds:
            continue
        d = {m: np.array([runs[v][s]["raw"][m] - base[s]["raw"][m] for s in seeds]) for m in MEM}
        ovr = np.array([sum(SLOPE[m] * d[m][i] for m in MEM) / 6 for i in range(len(seeds))])
        de = np.array([sum(SLOPE[m] * d[m][i] for m in ("nmae", "fid", "reach", "jac")) / 6 for i in range(len(seeds))])
        summ = {k: float(np.mean([runs[v][s]["summary"].get(k, np.nan) for s in seeds])) for k in SUMK}
        rows.append(dict(variant=v, seeds=seeds, **{f"d_{m}": float(d[m].mean()) for m in MEM},
                         **{f"sd_{m}": float(d[m].std(ddof=1)) if len(seeds) > 1 else None for m in MEM},
                         ovr=float(ovr.mean()), ovr_sd=float(ovr.std(ddof=1)) if len(seeds) > 1 else None,
                         nmae_ovr=float(-1.6532 * d["nmae"].mean() / 6), de_ovr=float(de.mean()), **summ))
    bl = {m: [base[s]["raw"][m] for s in sorted(base)] for m in MEM}
    lines = [f"# 423 table ({a.bulk} bulk); base seeds {sorted(base)}", "",
             "base: " + ", ".join(f"{m} {np.mean(bl[m]):.5f}" for m in MEM), "",
             "| variant | seeds | Δnmae (±sd) | Δpds | Δmse | Δfid | Δreach | Δjac | nmae→ovr | total ovr (±sd) | feasible (mass-w) up/dn | zs share | genes moved | ΔU genes (ballast) | call flips | sign flips |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        sd = lambda m: (f" ± {r['sd_' + m]:.5f}" if r["sd_" + m] is not None else "")
        osd = "" if r["ovr_sd"] is None else f" ± {r['ovr_sd']:.5f}"
        lines.append(
            f"| {r['variant']} | {len(r['seeds'])} | {r['d_nmae']:+.5f}{sd('nmae')} | {r['d_pds']:+.5f} | {r['d_mse']:+.6f} | "
            f"{r['d_fid']:+.5f} | {r['d_reach']:+.5f} | {r['d_jac']:+.6f} | {r['nmae_ovr']:+.5f} | {r['ovr']:+.5f}"
            f"{osd} | "
            f"{r['frac_massw_mean']:.3f} ({r['frac_massw_up_mean']:.3f}/{r['frac_massw_dn_mean']:.3f}) | "
            f"{r['frac_massw_zs_mean']:.3f} | {r['n_moved_mean']:.0f} | {r['n_genes_U_changed_mean']:.0f} | "
            f"{r['call_flips_mean']:.2f} | {r['sign_flips_mean']:.2f} |")
    out = Path(a.out) if a.out else HERE / "work" / f"analysis423{sfx}"
    out.with_suffix(".md").write_text("\n".join(lines) + "\n")
    out.with_suffix(".json").write_text(json.dumps(rows, indent=1) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
