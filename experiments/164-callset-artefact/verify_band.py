"""164 — is the h1-test prediction's 5-10 CPM band ZEROED by a --min-cpm gate, or not?

`decomp2.py` compared median |log2(pred_pb / ctrl_pb)| in the band (0.173) against the 10-20
CPM slice (0.123) and concluded "not zeroed".  That comparison is NOT safe on its own: the
emitter floors each gene's panel total, so a genuinely zeroed low-count gene still shows a
nonzero ratio, and the size of that floor artefact grows as counts fall -- exactly where the
band lives.  A detector whose noise grows into the signal it is looking for cannot separate
"small effect" from "zero effect plus quantisation".

This settles it with a statistic quantisation cannot fake: the BETWEEN-TARGET spread of a
gene's predicted pseudobulk.

  * If the gene is zeroed, its profile is `expm1(b_ctrl)/Z` for EVERY target, so its total is
    `floor(N*L_t*p_g)` and `pred_pb` is the same number for every target up to one count --
    a relative spread bounded by ~1/T_g.
  * If the gene carries a real donor effect, `pred_pb` moves from target to target by the
    effect itself, which is shared across all 400 cells and is not a rounding artefact.

So: per gene, (max - min) of log2(pred_pb) across targets, against the flooring bound 2/T_g.
FAILURE MODE, so the check is not vacuous: if the band IS zeroed, its between-target spread
sits at the flooring bound and the ratio below is ~1 for the band while being >>1 for the
gated slice.
"""
from __future__ import annotations

import anndata as ad
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = (f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad")
CTRL, COL = "non-targeting", "target_gene"
NTGT = 8


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def main() -> None:
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    gfull = np.asarray([str(g) for g in R.var_names])
    ci = np.flatnonzero(rlab == CTRL)
    csum = np.zeros(gfull.size)
    cellcpm = np.zeros(gfull.size)
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        csum += m.sum(0)
        cellcpm += (m * (1e6 / s[:, None])).sum(0)
    cellcpm /= ci.size
    ctrl_pb = 1e6 * csum / csum.sum()
    del R

    P = ad.read_h5ad(PRED)
    plab = np.asarray(P.obs[COL].astype(str).values, dtype="U64")
    pg = np.asarray([str(g) for g in P.var_names])
    tg = sorted(set(plab) - {CTRL})[:NTGT]
    pos = {g: i for i, g in enumerate(pg)}
    cols = np.array([pos.get(g, -1) for g in gfull])
    ok = cols >= 0
    tot, pbs = [], []
    for t in tg:
        v = dense(P, np.flatnonzero(plab == t)).sum(0)
        T = np.zeros(gfull.size)
        T[ok] = v[cols[ok]]
        tot.append(T)
        pbs.append(1e6 * T / T.sum())
    del P
    tot = np.array(tot)
    pbs = np.array(pbs)
    print(f"targets: {len(tg)}   median panel total for a 5-CPM gene: "
          f"{np.median(tot[:, (ctrl_pb > 4.5) & (ctrl_pb < 5.5)]):.0f} counts")

    lg = np.log2(np.where(pbs > 0, pbs, np.nan))
    spread = np.nanmax(lg, 0) - np.nanmin(lg, 0)          # between-target spread
    Tmin = tot.min(0)
    floor_bound = np.log2(1.0 + 2.0 / np.maximum(Tmin, 1))  # 2 counts of flooring slack

    scored = cellcpm > 5.0
    band = scored & (ctrl_pb <= 10.0)
    print(f"\n{'block':<34}{'n':>7}{'med between-tgt spread':>24}"
          f"{'med flooring bound':>20}{'ratio':>9}")
    for nm, m in (("sub-gate band (cell>5, pb<=10)", band),
                  ("pb-CPM 10-20", scored & (ctrl_pb > 10) & (ctrl_pb <= 20)),
                  ("pb-CPM 20-50", scored & (ctrl_pb > 20) & (ctrl_pb <= 50)),
                  ("pb-CPM > 50", scored & (ctrl_pb > 50))):
        if not m.any():
            continue
        r = np.nanmedian(spread[m]) / np.nanmedian(floor_bound[m])
        print(f"{nm:<34}{int(m.sum()):>7}{np.nanmedian(spread[m]):>24.4f}"
              f"{np.nanmedian(floor_bound[m]):>20.4f}{r:>9.2f}")

    r_band = np.nanmedian(spread[band]) / np.nanmedian(floor_bound[band])
    print(f"\nVERDICT: band ratio {r_band:.2f} -> "
          + ("ZEROED (between-target spread is pure flooring)" if r_band < 1.5 else
             "NOT ZEROED (band genes carry real, target-specific predicted effects); "
             "the 095 h1-test artifact was not built at --min-cpm 10, so the challenge's "
             "794/958/997-gene sub-gate block has NO analogue on this surface"))


if __name__ == "__main__":
    main()
