"""Step 1 (no emission): is |predicted effect| a usable stand-in for our own p-value rank?

Two questions, both answered off the tables experiment 092 already produced:

  (a) how monotone is the emitted p-value ordering in |predicted effect|, restricted to
      the genes we actually call significant -- exp013's claim, re-measured here;
  (b) what does the 092 truncation counterfactual pay if the call set is truncated by
      |predicted effect| rank instead of by our own adjusted-p rank? The spread policy can
      only rank by (b)'s key, so (b) is the ceiling the real mechanism is chasing.
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
rq = d["real_padj"]
G = genes.size
assert list(genes) == list(gate.genes)

gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64)
fr[np.argsort(genes)] = np.arange(G)

print("=" * 90)
print("(a) rank key |log2 pred_cpm/ctrl_cpm| from TOTALS  vs  the emitted p-value order")
print("=" * 90)
sp_all, sp_sig, ov = [], [], []
for i, t in enumerate(perts):
    key = L.rank_key(tot[i], gate)
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None:
        ok[j] = False
    sp_all.append(stats.spearmanr(key[ok], -pq[i][ok]).statistic)
    sig = ok & (pq[i] < L.ALPHA)
    idx = np.flatnonzero(sig)
    # our call list, cell_eval2 ranking key: p_adj asc, p asc, |lfc| desc, feature asc
    order = idx[np.lexsort((fr[idx], -np.abs(pl[i][idx]), pp[i][idx], pq[i][idx]))]
    sp_sig.append(stats.spearmanr(key[idx], -pq[i][idx]).statistic)
    for m in (850,):
        a = set(order[:m].tolist())
        b = set(idx[np.argsort(-key[idx], kind="stable")][:m].tolist())
        ov.append(len(a & b) / m)
print(f"spearman(|effect|, -p_adj) over all gated genes : {np.mean(sp_all):+.4f} "
      f"(min {np.min(sp_all):+.4f})")
print(f"    ... restricted to our significant call set   : {np.mean(sp_sig):+.4f} "
      f"(min {np.min(sp_sig):+.4f})")
print(f"top-850 set overlap, |effect| rank vs p_adj rank : {np.mean(ov):.4f} "
      f"(min {np.min(ov):.4f})")

print()
print("=" * 90)
print("(b) truncation counterfactual: p_adj rank order vs |effect| rank order")
print("=" * 90)
pre = T.precompute()
P = len(pre)

# Rebuild the same per-pert structures but with the call list ordered by |effect|.
pre_eff = []
for i, t in enumerate(perts):
    key = L.rank_key(tot[i], gate)
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None:
        ok[j] = False
    rs = ok & (rq[i] < L.ALPHA)
    ps = ok & (pq[i] < L.ALPHA)
    ind = ok & (d["real_lfc"][i] != 0) & ~np.isnan(d["real_lfc"][i])
    idx = np.flatnonzero(ps)
    r = idx[np.lexsort((fr[idx], -key[idx]))]
    adj = ind[r]
    hit = adj & (pl[i][r] != 0) & (np.sign(pl[i][r]) == np.sign(d["real_lfc"][i][r]))
    pre_eff.append(dict(t=t, n_real=int(rs.sum()), n_all=int(r.size),
                        cum_adj=np.cumsum(adj), cum_hit=np.cumsum(hit),
                        cum_in=np.cumsum(rs[r])))

bf, bj, _ = T.fid_jac(pre, [None] * P)
b_of, b_oj = T.scaled(float(np.nanmean(bf)), float(np.nanmean(bj)))
print(f"baseline           fid {np.nanmean(bf):.6f}  jac {np.nanmean(bj):.6f}")
print(f"{'m':>7}{'fid(padj)':>12}{'jac(padj)':>12}{'d6(padj)':>11}"
      f"{'fid(|eff|)':>13}{'jac(|eff|)':>12}{'d6(|eff|)':>11}")
for m in (400, 600, 850, 1200, 2000, 3000):
    f1, j1, _ = T.fid_jac(pre, [float(m)] * P)
    f2, j2, _ = T.fid_jac(pre_eff, [float(m)] * P)
    o1 = T.scaled(float(np.nanmean(f1)), float(np.nanmean(j1)))
    o2 = T.scaled(float(np.nanmean(f2)), float(np.nanmean(j2)))
    print(f"{m:>7}{np.nanmean(f1):>12.6f}{np.nanmean(j1):>12.6f}"
          f"{((o1[0] - b_of) + (o1[1] - b_oj)) / 6:>+11.4f}"
          f"{np.nanmean(f2):>13.6f}{np.nanmean(j2):>12.6f}"
          f"{((o2[0] - b_of) + (o2[1] - b_oj)) / 6:>+11.4f}")

# scale-adaptive form: top 15% of our own call list
for frac in (0.10, 0.15, 0.20):
    mv = [frac * s["n_all"] for s in pre]
    f1, j1, mm = T.fid_jac(pre, mv)
    f2, j2, _ = T.fid_jac(pre_eff, [frac * s["n_all"] for s in pre_eff])
    o1 = T.scaled(float(np.nanmean(f1)), float(np.nanmean(j1)))
    o2 = T.scaled(float(np.nanmean(f2)), float(np.nanmean(j2)))
    print(f"{f'top {frac:.0%}':>7}{np.nanmean(f1):>12.6f}{np.nanmean(j1):>12.6f}"
          f"{((o1[0] - b_of) + (o1[1] - b_oj)) / 6:>+11.4f}"
          f"{np.nanmean(f2):>13.6f}{np.nanmean(j2):>12.6f}"
          f"{((o2[0] - b_of) + (o2[1] - b_oj)) / 6:>+11.4f}   mean m={mm.mean():.0f}")

n_all = np.array([s["n_all"] for s in pre])
n_real = np.array([s["n_real"] for s in pre])
print(f"\nn_pred(sig, pre-truncation) median {np.median(n_all):.0f} "
      f"[{n_all.min()}, {n_all.max()}]   n_real median {np.median(n_real):.0f}")
