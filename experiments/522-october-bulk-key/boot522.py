"""522 — target bootstrap of the paired key-form differences (arms522's saved per-target arrays; the arms themselves are
fixed, targets are the resampling unit, as 518's robust518 / 438's CIs).

cos_agg(arm) on a resample = mean over frames of sum_p x_p / sqrt(sum_p a_p * sum_p d_p) (row_score's identity). For each
comparison A vs B, per fold: r = cos_agg(A) / cos_agg(B) - 1 on the same resample; a panel's statistic is the mean of its
cancer folds' r (each fold resampled independently). Reports the point value (== arms522's exact ratio), the 95 %
interval and P(r <= 0), 2,000 draws, seed 522.

    uv run python experiments/522-october-bulk-key/boot522.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
PANELS = ("abc", "efp", "efr", "real")
NB = 2000
PAIRS = (("pctPROD", "cap3"), ("cap6", "cap3"), ("pctPROD", "cap6"), ("S1", "cap3"), ("S1", "pctPROD"),
         ("pctPRODxS0.25", "pctPROD"), ("pctPRODxS0.5", "pctPROD"), ("cap6xS0.5", "cap6"), ("cap1e+09", "pctPROD"),
         ("gate", "pctPROD"))


def cos_agg(z, arm, idx):
    x, a, d = z[f"{arm}__x"][:, idx], z[f"{arm}__a"][:, idx], z[f"{arm}__d"][idx]
    return float(np.mean([x[j].sum() / np.sqrt(a[j].sum() * d.sum()) for j in range(x.shape[0])]))


def main():
    rng = np.random.default_rng(522)
    out = {}
    for p in PANELS:
        Z = {f: np.load(WORK / f"arms522-{p}-{f}.npz") for f in CANCER if (WORK / f"arms522-{p}-{f}.npz").exists()}
        if len(Z) < 3:
            continue
        n = {f: Z[f]["cap3__d"].size for f in Z}
        draws = {f: [rng.integers(0, n[f], n[f]) for _ in range(NB)] for f in Z}
        out[p] = {}
        for A, B in PAIRS:
            pt = float(np.mean([cos_agg(Z[f], A, np.arange(n[f])) / cos_agg(Z[f], B, np.arange(n[f])) - 1 for f in Z]))
            J = {f: json.loads((WORK / f"arms522-{p}-{f}.json").read_text())["arms"] for f in Z}
            ex = float(np.mean([J[f][A]["cos_agg"] / J[f][B]["cos_agg"] - 1 for f in Z]))
            assert abs(pt - ex) < 1e-12, f"KA arrays vs exact arms: {p} {A} vs {B}: {pt} {ex}"   # row_score identity
            bs = np.asarray([np.mean([cos_agg(Z[f], A, draws[f][b]) / cos_agg(Z[f], B, draws[f][b]) - 1 for f in Z])
                             for b in range(NB)])
            lo, hi = np.quantile(bs, [0.025, 0.975])
            out[p][f"{A} vs {B}"] = dict(r=pt, ci95=[float(lo), float(hi)], p_le0=float(np.mean(bs <= 0)))
            print(f"{p:5s} {A + ' vs ' + B:30s} {100 * pt:+7.2f} %  [{100 * lo:+6.2f}, {100 * hi:+6.2f}]  P(<=0) {np.mean(bs <= 0):.3f}",
                  flush=True)
    (WORK / "boot522.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
