"""212 — the paired target-cluster bootstrap on Δcos and Δpds, 204's protocol unchanged.

204 measured this instrument's resolution as +-0.000044 absolute on Δcos (4,000 draws,
tightest arms), i.e. **+-0.20% relative = +-0.00016 board cosine ~ +-0.00006 overall**
through 8.2's exchange rate. This script re-measures it on the gamma arms so every Δcos in
RESULTS.md carries its own interval rather than borrowing 204's.

Targets are resampled WITHIN each fold with replacement, and the SAME resampled index
vector is used for the arm and for the incumbent, so the interval is on the paired
difference and the fold-to-fold level cancels.

    uv run python boot212.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
NDRAW = 4000
SEED = 204


def macro(rows, key, idx=None):
    """Mean over folds of the mean per-target value, optionally on resampled indices."""
    out = []
    for f, v in rows.items():
        a = np.asarray(v[key])
        out.append(a[idx[f]].mean() if idx is not None else a.mean())
    return float(np.mean(out))


def main() -> None:
    src = json.loads((WORK / "plolo212.json").read_text())
    folds = sorted({r["fold"] for r in src})
    gammas = sorted({r["gamma"] for r in src})
    by = {(r["fold"], r["gamma"]): r for r in src}
    rng = np.random.default_rng(SEED)
    n = {f: len(by[(f, 0.65)]["_cos_per_target_a23"]) for f in folds}
    draws = [{f: rng.integers(0, n[f], n[f]) for f in folds} for _ in range(NDRAW)]

    base = {f: by[(f, 0.65)] for f in folds}
    out = []
    for g in gammas:
        arm = {f: by[(f, g)] for f in folds}
        row = dict(gamma=g)
        for key, tag in (("_cos_per_target_a23", "cos_a23"),
                         ("_cos_per_target_a15", "cos_a15"),
                         ("_pds_per_target_a23", "pds_a23")):
            point = macro(arm, key) - macro(base, key)
            d = np.asarray([macro(arm, key, idx) - macro(base, key, idx) for idx in draws])
            lo, hi = np.percentile(d, [2.5, 97.5])
            row[tag] = dict(delta=point, lo=float(lo), hi=float(hi),
                            half_width=float((hi - lo) / 2),
                            p_gt0=float((d > 0).mean()),
                            rel_half_width=float((hi - lo) / 2 / macro(base, key)))
        out.append(row)
        c = row["cos_a15"]
        print(f"g{g:5.3f}  dcos(a15) {c['delta']:+.7f} "
              f"[{c['lo']:+.7f}, {c['hi']:+.7f}]  half {c['half_width']:.7f} "
              f"= {c['rel_half_width']:.3%} of the incumbent  P(>0) {c['p_gt0']:.3f}   "
              f"dpds {row['pds_a23']['delta']:+.6f} "
              f"[{row['pds_a23']['lo']:+.6f}, {row['pds_a23']['hi']:+.6f}]")

    with open(WORK / "boot212.json", "w") as fh:
        json.dump(dict(ndraw=NDRAW, seed=SEED, n_targets=n, rows=out), fh, indent=1)

    # 204's own published resolution, re-derived on this sweep's tightest row
    tight = min(out, key=lambda r: r["cos_a23"]["half_width"] if r["gamma"] != 0.65 else 9)
    print(f"\n[boot] tightest non-identity arm: gamma {tight['gamma']}, Δcos half-width "
          f"{tight['cos_a23']['half_width']:.2e} absolute = "
          f"{tight['cos_a23']['rel_half_width']:.3%} relative "
          f"(204 published 4.3e-5 / 0.20%)")


if __name__ == "__main__":
    main()
