"""Step 2 (no emission): a ranking key that actually tracks the emitted p-value order.

s1 measured that |predicted effect| recovers only +0.0041 of the +0.0299 counterfactual --
it is not the ordering the truncation needs. The reason is structural, not incidental.

Under a low-variance emitter every one of the 400 predicted cells carries essentially the
same value `v_g = log1p(CPM_pred_g)`, so the Wilcoxon rank-sum against the m real control
cells reduces to *where that single value sits inside the control's own empirical
distribution*:

    R1  ~  n * (m * F_ctrl(v_g) + (n + 1) / 2),      z  ~  sqrt(12 n m / (N + 1)) * (F - 1/2)

i.e. the p-value order is the order of |F_ctrl(v_g) - 1/2| -- a *quantile displacement*,
not an effect size. A large predicted log-fold change on a gene whose control distribution
is broad moves F very little; a small one on a tight, high-expressed gene saturates it.
That is exactly the gene class the reference also calls, which is why the p-adj order beats
|effect| by 7x in the counterfactual.

Everything F needs -- the control cells and our own profile -- is in hand inside the
submission builder before a single cell is emitted, so this key is as production-available
as |effect| is. Measured here against the same 092 counterfactual.
"""
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402

sys.path.insert(0, L.E092)
import trunc_lib as T  # noqa: E402

gate = L.Gate()
tot = np.load(os.path.join(L.SCRATCH, "totals093.npy"))
d = np.load(os.path.join(L.E092, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str)
perts = d["perts"].astype(str)
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
rq, rl = d["real_padj"], d["real_lfc"]
G = genes.size
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64)
fr[np.argsort(genes)] = np.arange(G)

ctrl_sorted = np.load(os.path.join(L.SCRATCH, "ctrl_log_sorted.npy"), mmap_mode="r")
m = ctrl_sorted.shape[1]
print(f"control block {ctrl_sorted.shape}")


def quantile_key(totals):
    """|F_ctrl(v_g) - 1/2| for our predicted per-cell value v_g."""
    t = np.asarray(totals, dtype=np.float64)
    cpm = 1e6 * t / t.sum()
    v = np.log1p(cpm[gate.keep_idx])
    lo = np.empty(G)
    hi = np.empty(G)
    for g in range(G):
        row = ctrl_sorted[g]
        lo[g] = np.searchsorted(row, v[g], "left")
        hi[g] = np.searchsorted(row, v[g], "right")
    F = 0.5 * (lo + hi) / m          # mid-rank quantile, the tie convention the test uses
    return np.abs(F - 0.5)


KEYS = {}
for i in range(len(perts)):
    KEYS[i] = quantile_key(tot[i])
    if i % 40 == 0:
        print(f"  key {i + 1}/{len(perts)}", flush=True)

print("\n" + "=" * 92)
print("ordering diagnostics against the emitted p_adj order (over our significant set)")
print("=" * 92)
sp_e, sp_q, ov_e, ov_q = [], [], [], []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None:
        ok[j] = False
    idx = np.flatnonzero(ok & (pq[i] < L.ALPHA))
    order = idx[np.lexsort((fr[idx], -np.abs(pl[i][idx]), pp[i][idx], pq[i][idx]))]
    ke = L.rank_key(tot[i], gate)
    kq = KEYS[i]
    sp_e.append(stats.spearmanr(ke[idx], -pp[i][idx]).statistic)
    sp_q.append(stats.spearmanr(kq[idx], -pp[i][idx]).statistic)
    top = set(order[:850].tolist())
    ov_e.append(len(top & set(idx[np.argsort(-ke[idx], kind="stable")][:850])) / 850)
    ov_q.append(len(top & set(idx[np.argsort(-kq[idx], kind="stable")][:850])) / 850)
print(f"spearman vs -p_value :  |effect| {np.mean(sp_e):+.4f}   quantile {np.mean(sp_q):+.4f}")
print(f"top-850 overlap      :  |effect| {np.mean(ov_e):.4f}    quantile {np.mean(ov_q):.4f}")


def build_pre(keyfn):
    out = []
    for i, t in enumerate(perts):
        ok = np.ones(G, bool)
        j = gpos.get(t)
        if j is not None:
            ok[j] = False
        rs = ok & (rq[i] < L.ALPHA)
        ps = ok & (pq[i] < L.ALPHA)
        ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
        idx = np.flatnonzero(ps)
        k = keyfn(i)
        r = idx[np.lexsort((fr[idx], -k[idx]))]
        adj = ind[r]
        hit = adj & (pl[i][r] != 0) & (np.sign(pl[i][r]) == np.sign(rl[i][r]))
        out.append(dict(t=t, n_real=int(rs.sum()), n_all=int(r.size),
                        cum_adj=np.cumsum(adj), cum_hit=np.cumsum(hit),
                        cum_in=np.cumsum(rs[r])))
    return out


pre_p = T.precompute()
pre_e = build_pre(lambda i: L.rank_key(tot[i], gate))
pre_q = build_pre(lambda i: KEYS[i])
P = len(pre_p)
bf, bj, _ = T.fid_jac(pre_p, [None] * P)
b_of, b_oj = T.scaled(float(np.nanmean(bf)), float(np.nanmean(bj)))


def row(lab, mv_of):
    cells = []
    for pre in (pre_p, pre_e, pre_q):
        f, j, _ = T.fid_jac(pre, mv_of(pre))
        o = T.scaled(float(np.nanmean(f)), float(np.nanmean(j)))
        cells.append((float(np.nanmean(f)), ((o[0] - b_of) + (o[1] - b_oj)) / 6))
    print(f"{lab:>9}" + "".join(f"{c[0]:>11.6f}{c[1]:>+10.4f}" for c in cells))


print("\n" + "=" * 92)
print("truncation counterfactual (fid, d_avg6) by ORDERING KEY")
print(f"{'m':>9}{'fid p_adj':>11}{'d6':>10}{'fid |eff|':>11}{'d6':>10}{'fid quant':>11}{'d6':>10}")
print("=" * 92)
for mm in (400, 600, 850, 1200, 2000, 3000, 4000):
    row(str(mm), lambda pre, mm=mm: [float(mm)] * P)
for frac in (0.10, 0.15, 0.20, 0.30):
    row(f"top {frac:.0%}", lambda pre, fr_=frac: [fr_ * s["n_all"] for s in pre])
