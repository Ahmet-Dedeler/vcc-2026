#!/usr/bin/env python3
"""528 — markdown tables for RESULTS §4 from work/price528.json (price528.py's output): per arm (zs form = 446bj2z's DE
construction), the overall value at the four regimes with the seed-averaged target CI, the member deltas, fid + jac,
board columns, pds/mse pins and deliverability.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/summ528.py > experiments/528-callprop-block/work/summ528.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ORDER = ["lcap-kept", "lcap", "lcap-dec", "lcap-k1980", "lcap-k1650", "lcap-k1100", "lcap-k1980-kept", "lcap-k1650-kept", "ltv",
         "p2f", "shuf", "lcap+mzg", "lcap-kept+mzg", "lcap-n1200", "lcap-kept-n1200"]


def f(x, d=4):
    return f"{x:+.{d}f}"


def main():
    d = json.load(open(HERE / "work/price528.json"))
    arms = d["arms"]
    keys = [k for k in ORDER if f"{k} +zs" in arms] + sorted(k[:-4] for k in arms if k.endswith(" +zs") and k[:-4] not in ORDER)
    print("| arm (vs its base) | seeds | overall f 1.0 [95 % CI] | f 0.5 | f 0.35 | cmp [95 % CI] | Δfid / Δjac / Δreach / Δnmae raw (f 1) | fid+jac overall f 1 / cmp | board f 1 [reach ×1 \\| ×0 \\| ×−0.76] |")
    print("|---|---|---|---|---|---|---|---|---|")
    for k in keys:
        r = arms[f"{k} +zs"]
        a, b, c, m = r["1.0"], r["0.5"], r["0.35"], r["cmp"]
        dm = a["d_mean"]
        seeds = " / ".join(f(p["value"]) for p in a["per_seed"])
        print(f"| `{k}` (vs {r['base']}) | {len(r['seeds'])}: {seeds} | **{f(a['value_mean'])}** [{f(a['ci'][0])}, {f(a['ci'][1])}] | "
              f"{f(b['value_mean'])} | {f(c['value_mean'])} | {f(m['value_mean'])} [{f(m['ci'][0])}, {f(m['ci'][1])}] | "
              f"{f(dm['fid'])} / {f(dm['jac'])} / {f(dm['reach'])} / {f(dm['nmae'])} | {f(a['fid_jac_overall'])} / "
              f"{f(m['fid_jac_overall'])} | {f(a['board'][0])} \\| {f(a['board'][1])} \\| {f(a['board'][2])} |")
    print()
    print("| arm | pds / mse Δ (all seeds) | per-gene totals == base | swapped-in / target | swapped-in call rate | missed / target (base) | unsuppressible / target (base) | swapped-out now called | realized vs prescribed jac (f 1) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for k in keys:
        r = arms.get(k)
        if not r or "deliverability" not in r:
            continue
        pm = r["pds_mse"]
        pins = max(max(abs(v["pds"]), abs(v["mse"])) for v in pm.values())
        eq = all(v["sums_equal"] for v in pm.values())
        dl = r["deliverability"]
        A, B = dl["arm"], dl["base"]
        sw = A["block swapped-in"]
        rv = dl["realized_vs_prescribed"]["jac"]
        print(f"| `{k}` | max {pins:.1e} | {eq} | {sw['n_per_target']:.0f} | {sw['call_rate']:.4f} | {A['missed_per_target']:.1f} "
              f"({B['missed_per_target']:.1f}) | {A['unsuppressible_per_target']:.1f} ({B['unsuppressible_per_target']:.1f}) | "
              f"{A['swapped-out (now complement)']['call_rate']:.4f} | {rv[0]:.5f} / {rv[1]:.5f} |")


if __name__ == "__main__":
    main()
