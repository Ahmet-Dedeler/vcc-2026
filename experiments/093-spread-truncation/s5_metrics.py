"""Step 5: score every emitted policy on all six metrics.

`fid` / `jac` follow experiments/092/step3_metrics.py, which reproduces the scorer to
2.5e-6; `reach` / `nmae` follow experiments/092/step12_reach_exact.py, which is transcribed
from cell_eval2 0.16.0 `metrics/direction.py` and `metrics/de.py`.

`pds` and `mse` are read off the pseudobulk. cell_eval2 builds that with
`prep._grouped_sums` -> `prep.bulk_lognorm_means`, i.e. from RAW per-group column sums, so
the emitter's exact-total repair makes both bit-identical across every policy; the sweep
asserts that equality gene-by-gene as an integer at emission time. (`expr_mse_unbiased`
additionally SUBTRACTS a delete-1 jackknife `C_p` that does read cell-level spread, so more
spread can only lower the raw mse -- and the official mse is already clamped at 0 from far
above its baseline, 2.175 against 0.989, so it cannot move either way.)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

A = L.ALPHA
d = np.load(os.path.join(L.E092, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str)
perts = d["perts"].astype(str)
rp, rq, rl = d["real_p"], d["real_padj"], d["real_lfc"]
G, P = genes.size, perts.size
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64)
fr[np.argsort(genes)] = np.arange(G)
OKS = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None:
        ok[j] = False
    OKS.append(ok)


def score(pp, pq, pl):
    """(fid, jac, reach, nmae, n_pred array) for one emitted prediction."""
    fid = np.empty(P)
    jac = np.empty(P)
    rch = np.empty(P)
    npred = np.empty(P, np.int64)
    nm = []
    for i in range(P):
        ok = OKS[i]
        rs = ok & (rq[i] < A)                       # reference significant set
        ps = ok & (pq[i] < A)                       # our significant set
        ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])  # adjudicable
        scored = ps & ind
        match = scored & (pl[i] != 0) & ~np.isnan(pl[i]) & (np.sign(pl[i]) == np.sign(rl[i]))
        n_real, n_pred, k = int(rs.sum()), int(scored.sum()), int(match.sum())
        npred[i] = n_pred
        fid[i] = k / max(n_pred, n_real) if max(n_pred, n_real) > 0 else np.nan
        union = int((rs | ps).sum())
        jac[i] = 1.0 if union == 0 else int((rs & ps).sum()) / union

        # ---- reach, exactly per cell_eval2 metrics/direction.py ----
        gsel = np.flatnonzero(rs)
        nc = gsel.size
        if nc == 0:
            rch[i] = np.nan
        else:
            rlv, plv = rl[i][gsel], pl[i][gsel]
            in_denom = (~np.isnan(rlv)) & (rlv != 0)
            committed = (~np.isnan(plv)) & (plv != 0)
            mt = in_denom & committed & (np.sign(plv) == np.sign(rlv))
            sig_pred = ps[gsel]
            order = np.lexsort((fr[gsel], -np.abs(np.nan_to_num(plv, nan=-np.inf)),
                                pp[i][gsel], pq[i][gsel], ~sig_pred))
            nd = np.cumsum(in_denom[order])
            nmt = np.cumsum(mt[order])
            with np.errstate(invalid="ignore", divide="ignore"):
                pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
            good = np.flatnonzero((nd > 0) & (pur >= 0.9))
            rch[i] = (nd[good[-1]] / nc) if good.size else 0.0

            # ---- nmae, exactly per cell_eval2 metrics/de.py ----
            gate_ok = np.isfinite(rlv)
            if gate_ok.sum() >= 10:
                pv = np.where(np.isfinite(plv), plv, 0.0)[gate_ok]
                rv = rlv[gate_ok]
                den = np.abs(rv).mean()
                if den > 0 and np.isfinite(den):
                    nm.append(np.abs(pv - rv).mean() / den)
    return (float(np.nanmean(fid)), float(np.nanmean(jac)), float(np.nanmean(rch)),
            float(np.mean(nm)), npred)


SWEEP = os.path.join(L.SCRATCH, "sweep093")
base = score(d["pred_p"], d["pred_padj"], d["pred_lfc"])
bs = {k: L.scaled(k, v) for k, v in zip(("fid", "jac", "reach", "nmae"), base[:4])}
print(f"PRODUCTION (the shipped emission, spread 0.25)")
print(f"  fid {base[0]:.6f}  jac {base[1]:.6f}  reach {base[2]:.6f}  nmae {base[3]:.6f}")
print(f"  scorer values for this run: fid 0.5312814007682398  jac 0.0340283619023056")
print(f"  n_pred median {np.median(base[4]):.0f}  [{base[4].min()}, {base[4].max()}]   "
      f"over-calling perts {int((base[4] > np.array([int((OKS[i] & (rq[i] < A)).sum()) for i in range(P)])).sum())}/{P}")
print()
hdr = (f"{'policy':<16}{'n_pred med':>11}{'n_pred p10':>11}{'n_pred p90':>11}"
       f"{'fid':>10}{'jac':>9}{'reach':>9}{'nmae':>9}{'d_fid':>9}{'d_jac':>9}"
       f"{'d_reach':>9}{'d_nmae':>9}{'d_avg6':>9}")
print(hdr)
print("=" * len(hdr))
rows = []
names = sorted(os.listdir(SWEEP))
order = ([f"global_s{s}.npz" for s in (0.25, 0.4, 0.55, 0.7, 0.85, 1.0)]
         + [n for n in names if n.startswith("eff")]
         + [n for n in names if n.startswith("padj") or n.startswith("tight")])
for fn in order:
    if fn not in names:
        continue
    z = np.load(os.path.join(SWEEP, fn))
    f, j, r, n, npred = score(z["pred_p"], z["pred_padj"], z["pred_lfc"])
    sc = {k: L.scaled(k, v) for k, v in zip(("fid", "jac", "reach", "nmae"), (f, j, r, n))}
    dd = {k: sc[k] - bs[k] for k in sc}
    d6 = sum(dd.values()) / 6
    rows.append((fn[:-4], np.median(npred), f, j, r, n, dd, d6))
    print(f"{fn[:-4]:<16}{np.median(npred):>11.0f}{np.percentile(npred, 10):>11.0f}"
          f"{np.percentile(npred, 90):>11.0f}{f:>10.6f}{j:>9.6f}{r:>9.6f}{n:>9.6f}"
          f"{dd['fid']:>+9.4f}{dd['jac']:>+9.4f}{dd['reach']:>+9.4f}{dd['nmae']:>+9.4f}"
          f"{d6:>+9.4f}")
best = max(rows, key=lambda x: x[-1])
print(f"\nbest: {best[0]}  d_avg6 {best[-1]:+.4f}  (n_pred median {best[1]:.0f})")
np.save(os.path.join(L.SCRATCH, "s5_rows.npy"),
        np.array([(r[0], r[1], r[2], r[3], r[4], r[5], r[7]) for r in rows], dtype=object),
        allow_pickle=True)
