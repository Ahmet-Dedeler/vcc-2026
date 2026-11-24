"""279 — POST-HOC descriptives and two stability checks on the matched effect.

Everything here is labelled POST-HOC in RESULTS.md: it was written after the preregistered
numbers were read.  Nothing here can change the verdict, which is fixed by PREREGISTRATION §4.

  1. the ZNF share of the bottom/top decile of the held map (276 §5's panel statistic);
  2. the STRATIFIED AUC — P(a ZNF gene sits below a non-ZNF gene from the SAME
     (expr-decile, mag-decile) cell) — the rank version of the matched effect, which does not
     depend on a matching draw;
  3. the 1:1 matching repeated 200 times, to show the whole draw distribution.

    PYTHONDONTWRITEBYTECODE=1 uv run python extra279.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/279-znf-state-dependence"
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/276-conservation-map"))
sys.path.insert(0, str(HERE))
import analyse279 as A  # noqa: E402

SEED = 20260906


def stratified_auc(s, cell, is_z, is_n):
    """P(ZNF < non-ZNF | same cell), pooled over cells."""
    num = den = 0.0
    for c in np.unique(cell[cell >= 0]):
        x = s[is_z & (cell == c)]
        y = s[is_n & (cell == c)]
        if x.size == 0 or y.size == 0:
            continue
        # U = #(x<y) + 0.5 ties, computed by ranks
        allv = np.concatenate([x, y])
        r = np.argsort(np.argsort(allv, kind="stable"), kind="stable").astype(float)
        # average ranks for ties
        order = np.argsort(allv, kind="stable")
        sv = allv[order]
        i = 0
        while i < sv.size:
            j = i
            while j + 1 < sv.size and sv[j + 1] == sv[i]:
                j += 1
            if j > i:
                r[order[i:j + 1]] = (i + j) / 2.0
            i = j + 1
        rx = r[:x.size].sum()
        u_x_greater = rx - x.size * (x.size - 1) / 2.0     # #(x>y) + 0.5 ties
        num += x.size * y.size - u_x_greater
        den += x.size * y.size
    return float(num / den), float(den)


def main() -> None:
    t0 = time.time()
    zh = np.load(WORK / "held279.npz", allow_pickle=True)
    genes = [str(g) for g in zh["genes"]]
    c = np.asarray(zh["c_held"], np.float64)
    expr = np.asarray(zh["expr"], np.float64)
    mag = np.asarray(zh["mag"], np.float64)
    cls, groups = A.classify(genes)
    is_z, is_n = cls == "ZNF", cls == "NONZNF"
    de, dm = A.deciles(expr), A.deciles(mag)
    cell = np.where((de >= 0) & (dm >= 0), de * 10 + dm, -1)
    out: dict = {}

    # 1. decile shares (276 §5's panel statistic, on the held map) --------------------------------
    f = np.isfinite(c) & (is_z | is_n)
    q = np.quantile(c[f], [0.1, 0.9])
    bot = f & (c <= q[0]); top = f & (c >= q[1])
    out["deciles"] = dict(
        cut_bottom=float(q[0]), cut_top=float(q[1]),
        znf_share_bottom=float(is_z[bot].mean()), znf_share_top=float(is_z[top].mean()),
        znf_share_overall=float(is_z[f].mean()),
        share_le0_znf=float((c[is_z & f] <= 0).mean()),
        share_le0_nonznf=float((c[is_n & f] <= 0).mean()))
    d = out["deciles"]
    print(f"decile shares: ZNF is {d['znf_share_bottom']*100:.1f}% of the bottom decile, "
          f"{d['znf_share_top']*100:.1f}% of the top, {d['znf_share_overall']*100:.1f}% overall; "
          f"C<=0 for {d['share_le0_znf']*100:.1f}% of ZNF vs {d['share_le0_nonznf']*100:.1f}% "
          f"of non-ZNF", flush=True)

    # 1b. the same on the BALANCED CORE (identical 3 pairs for every gene): the held map's tails
    #     are dominated by genes with few pairs, which is why 1's decile shares read low.
    import csv as _csv
    _rows = list(_csv.reader(open(WORK / "znf279.csv")))[1:]
    s3d = np.asarray([float(r[14]) if r[14] else np.nan for r in _rows])
    f3 = np.isfinite(s3d) & (is_z | is_n)
    q3 = np.quantile(s3d[f3], [0.1, 0.9])
    b3 = f3 & (s3d <= q3[0]); t3 = f3 & (s3d >= q3[1])
    out["deciles_balanced_core"] = dict(
        znf_share_bottom=float(is_z[b3].mean()), znf_share_top=float(is_z[t3].mean()),
        znf_share_overall=float(is_z[f3].mean()),
        share_le0_znf=float((s3d[is_z & f3] <= 0).mean()),
        share_le0_nonznf=float((s3d[is_n & f3] <= 0).mean()))
    d3 = out["deciles_balanced_core"]
    print(f"decile shares (balanced core): ZNF is {d3['znf_share_bottom']*100:.1f}% of the bottom "
          f"decile, {d3['znf_share_top']*100:.1f}% of the top, {d3['znf_share_overall']*100:.1f}% "
          f"overall; C<=0 for {d3['share_le0_znf']*100:.1f}% vs {d3['share_le0_nonznf']*100:.1f}%",
          flush=True)

    # 2. stratified AUC, with a bootstrap over genes ----------------------------------------------
    out["stratified_auc"] = {}
    scores = {"C_held": c}
    s3 = np.asarray([float(r) if r else np.nan for r in
                     (row[14] for row in list(__import__("csv").reader(open(WORK / "znf279.csv")))[1:])])
    scores["c3_balanced_core"] = s3
    rng = np.random.default_rng(SEED)
    for nm, s in scores.items():
        ff = np.isfinite(s)
        a, den = stratified_auc(s, cell, is_z & ff, is_n & ff)
        iz = np.flatnonzero(is_z & ff); inn = np.flatnonzero(is_n & ff)
        bs = []
        for _ in range(500):
            zz = iz[rng.integers(0, iz.size, iz.size)]
            nn = inn[rng.integers(0, inn.size, inn.size)]
            sb, cb, lb = [], [], []
            sb.append(s[zz]); cb.append(cell[zz]); lb.append(np.ones(zz.size, bool))
            sb.append(s[nn]); cb.append(cell[nn]); lb.append(np.zeros(nn.size, bool))
            sv = np.concatenate(sb); cv = np.concatenate(cb); lv = np.concatenate(lb)
            aa, _ = stratified_auc(sv, cv, lv, ~lv)
            bs.append(aa)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        out["stratified_auc"][nm] = dict(auc=a, ci=[float(lo), float(hi)], n_pairs=den)
        print(f"stratified AUC ({nm}): {a:.4f} [{lo:.4f}, {hi:.4f}] over {den:.3g} within-cell "
              f"gene pairs", flush=True)

    # 3. the matching draw distribution ------------------------------------------------------------
    mt = A.matched_effect(c, cell, is_z & np.isfinite(c), is_n & np.isfinite(c), n_draw=200, b=2000)
    out["matched_200_draws"] = mt
    print(f"matched (200 draws): Δ {mt['d_median']:+.5f} [{mt['ci'][0]:+.5f}, {mt['ci'][1]:+.5f}], "
          f"per-draw [{mt['per_draw_min']:+.5f}, {mt['per_draw_max']:+.5f}]", flush=True)

    out["seconds"] = time.time() - t0
    (WORK / "extra279.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
