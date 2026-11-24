"""Which K does `--keep-top K --keep-spread 0` actually want?

`--keep-spread` below `--tail-spread` lowers the p-values of the top-K genes by
|predicted effect|, which moves them toward the FRONT of the confidence order the scorer
walks. It is a two-level approximation of a full |lfc| ordering, which experiment 095
measured at +0.0045 overall. This models the two-level version directly: promote the top-K
ahead of the rest of our significant calls, leave everything else in production order.
"""
import sys
import numpy as np

d = np.load(sys.argv[1] if len(sys.argv) > 1 else
            "experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
genes, perts = d["genes"].astype(str), d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
G, ALPHA = genes.size, 0.05
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)
B, R = 0.079022, 0.964145

def reach_of(order, idn, mt, nc):
    nd, nmt = np.cumsum(idn[order]), np.cumsum(mt[order])
    with np.errstate(invalid="ignore", divide="ignore"):
        pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
    good = np.flatnonzero((nd > 0) & (pur >= 0.9))
    return (nd[good[-1]] / nc) if good.size else 0.0

KS = [0, 250, 500, 1000, 2000, 3000, 4000, 6000]
out = {k: [] for k in KS}
out["full_lfc"] = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA); ps = ok & (pq[i] < ALPHA)
    gs = np.flatnonzero(rs); nc = gs.size
    if nc == 0: continue
    rlv, plv = rl[i][gs], pl[i][gs]
    idn = (~np.isnan(rlv)) & (rlv != 0)
    mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
    absl_all = np.abs(np.nan_to_num(pl[i], nan=-np.inf))
    absl = absl_all[gs]; sig = ps[gs]
    # rank of every gene by |predicted effect| over the whole gated axis, as the builder does
    rank_all = np.argsort(np.argsort(-absl_all))
    for K in KS:
        promoted = (rank_all[gs] < K) if K else np.zeros(nc, bool)
        out[K].append(reach_of(
            np.lexsort((fr[gs], -absl, pp[i][gs], pq[i][gs], ~promoted, ~sig)), idn, mt, nc))
    out["full_lfc"].append(reach_of(np.lexsort((fr[gs], -absl)), idn, mt, nc))

base = np.mean(out[0])
print(f"{'policy':<28}{'raw reach':>11}{'scaled':>10}{'Δoverall':>12}")
for K in KS:
    v = np.mean(out[K])
    tag = "production (K=0)" if K == 0 else f"promote top-{K} by |effect|"
    print(f"{tag:<28}{v:>11.5f}{(v-B)/(R-B):>10.4f}{(v-base)/(R-B)/6:>+12.5f}")
v = np.mean(out["full_lfc"])
print(f"{'full |lfc| order (ceiling)':<28}{v:>11.5f}{(v-B)/(R-B):>10.4f}{(v-base)/(R-B)/6:>+12.5f}")
