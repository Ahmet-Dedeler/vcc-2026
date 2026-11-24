"""The questions the depth-confounded harness answered wrongly, asked again on a clean one.

Every number in experiments 092/093/094 that touched direction, `fid`, `reach`, `jac` or
`n_real` was measured on `data/harness/orion-hct116.h5ad`, whose control cells are 38%
deeper than its perturbation cells (docs/HARNESS-DEPTH-CONFOUND.md). This re-asks them
against a repaired harness, scored on the anchors fitted from 780 live board entries
rather than the harness's own.

    uv run python experiments/095-clean-harness/analyse.py <de_tables.npz>
"""
from __future__ import annotations

import sys

import numpy as np

ALPHA = 0.05
# fitted from 780 leaderboard entries, R^2 > 0.9999 (see docs/STRATEGY.md addendum 6)
ANCH = {"fid": (0.51239, 0.80964), "reach": (0.07902, 0.96414), "jac": (0.03038, 0.39727)}


def scaled(k, v):
    b, r = ANCH[k]
    return (v - b) / (r - b)


def main(path: str) -> None:
    d = np.load(path, allow_pickle=True)
    genes, perts = d["genes"].astype(str), d["perts"].astype(str)
    rq, rl = d["real_padj"], d["real_lfc"]
    pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
    G = genes.size
    gpos = {g: i for i, g in enumerate(genes)}
    fr = np.empty(G, np.int64)
    fr[np.argsort(genes)] = np.arange(G)

    def reach_of(order, in_denom, mt, nc):
        nd, nmt = np.cumsum(in_denom[order]), np.cumsum(mt[order])
        with np.errstate(invalid="ignore", divide="ignore"):
            pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
        good = np.flatnonzero((nd > 0) & (pur >= 0.9))
        return (nd[good[-1]] / nc) if good.size else 0.0

    n_real, n_pred, acc, negfrac = [], [], [], []
    R = {k: [] for k in ("prod", "lfc", "oracle", "random")}
    fid, jac = [], []
    rng = np.random.default_rng(0)
    for i, t in enumerate(perts):
        ok = np.ones(G, bool)
        j = gpos.get(t)
        if j is not None:
            ok[j] = False
        rs = ok & (rq[i] < ALPHA)
        ps = ok & (pq[i] < ALPHA)
        ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
        scored = ps & ind
        match = scored & (pl[i] != 0) & ~np.isnan(pl[i]) & (np.sign(pl[i]) == np.sign(rl[i]))
        nr, npd, k = int(rs.sum()), int(scored.sum()), int(match.sum())
        n_real.append(nr)
        n_pred.append(npd)
        if max(npd, nr) > 0:
            fid.append(k / max(npd, nr))
        u = int((rs | ps).sum())
        jac.append(1.0 if u == 0 else int((rs & ps).sum()) / u)
        gs = np.flatnonzero(rs)
        if gs.size == 0:
            continue
        rlv, plv = rl[i][gs], pl[i][gs]
        idn = (~np.isnan(rlv)) & (rlv != 0)
        mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
        if idn.sum():
            acc.append(mt[idn].sum() / idn.sum())
            negfrac.append(float((rlv[idn] < 0).mean()))
        absl = np.abs(np.nan_to_num(plv, nan=-np.inf))
        sig = ps[gs]
        R["prod"].append(reach_of(np.lexsort((fr[gs], -absl, pp[i][gs], pq[i][gs], ~sig)), idn, mt, gs.size))
        R["lfc"].append(reach_of(np.lexsort((fr[gs], -absl)), idn, mt, gs.size))
        R["oracle"].append(reach_of(np.lexsort((fr[gs], ~mt)), idn, mt, gs.size))
        R["random"].append(reach_of(rng.permutation(gs.size), idn, mt, gs.size))

    print(f"perturbations {len(n_real)}  genes {G}\n")
    print("=== the reference ===")
    print(f"  n_real      median {np.median(n_real):7.0f}  mean {np.mean(n_real):7.0f}  "
          f"range [{min(n_real)}, {max(n_real)}]")
    print(f"  fraction of reference-DE genes NEGATIVE: {np.mean(negfrac):.4f}")
    print("     (0.908 on the broken harness, ~0.50 on h1-test / Replogle)\n")
    print("=== us ===")
    print(f"  n_pred      median {np.median(n_pred):7.0f}   n_pred/n_real median "
          f"{np.median(np.array(n_pred)/np.maximum(np.array(n_real),1)):.1f}x")
    print(f"  sign accuracy on the reference's DE set: {np.mean(acc):.4f}")
    print(f"  raw fid {np.nanmean(fid):.4f} (scaled {scaled('fid', np.nanmean(fid)):+.4f})   "
          f"raw jac {np.mean(jac):.4f} (scaled {scaled('jac', np.mean(jac)):+.4f})\n")
    print("=== reach by ordering (the only free variable) ===")
    base = np.mean(R["prod"])
    print(f"  {'ordering':<24}{'raw':>10}{'scaled':>10}{'Δoverall vs prod':>20}")
    for k in ("prod", "random", "lfc", "oracle"):
        v = np.mean(R[k])
        print(f"  {k:<24}{v:>10.5f}{scaled('reach', v):>10.4f}{(v-base)/(ANCH['reach'][1]-ANCH['reach'][0])/6:>+20.5f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
