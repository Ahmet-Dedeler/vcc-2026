"""The buildable form: promote the top-K by a candidate key, keep the emitted order inside.

`--keep-top K --keep-spread 0` ranks by |effect|. This asks what the same policy is worth
if it ranked by cross-donor agreement instead -- the key experiment 013 nominated and
nobody ever tested against the real metric.
"""
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
d = np.load(f"{ROOT}/experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
genes, perts = d["genes"].astype(str), d["perts"].astype(str)
rq, rl, pp, pq, pl = d["real_padj"], d["real_lfc"], d["pred_p"], d["pred_padj"], d["pred_lfc"]
G, ALPHA = genes.size, 0.05
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)
B, R = 0.079022, 0.964145
DON = {}
for name in ("k562_gwps", "orion_hek293t"):
    root = f"{ROOT}/data/donors/{name}"
    dg = np.array(open(f"{root}/genes.txt").read().split())
    dp = np.array(open(f"{root}/perts.txt").read().split())
    gi = {g: i for i, g in enumerate(dg)}
    DON[name] = (np.load(f"{root}/delta.npy", mmap_mode="r"),
                 np.array([gi.get(g, -1) for g in genes]), {p: i for i, p in enumerate(dp)})
def dv(name, t):
    delta, cols, pi = DON[name]
    if t not in pi: return None
    row = np.asarray(delta[pi[t]], float); out = np.zeros(G); ok = cols >= 0
    out[ok] = row[cols[ok]]; return out
def reach_of(o, idn, mt, nc):
    nd, nmt = np.cumsum(idn[o]), np.cumsum(mt[o])
    with np.errstate(invalid="ignore", divide="ignore"):
        pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
    g = np.flatnonzero((nd > 0) & (pur >= 0.9))
    return (nd[g[-1]] / nc) if g.size else 0.0

KS = [500, 1000, 2000, 3000, 4000]
res = {("effect", k): [] for k in KS}
res.update({("agree", k): [] for k in KS})
base = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool); j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA); ps = ok & (pq[i] < ALPHA)
    gs = np.flatnonzero(rs); nc = gs.size
    if nc == 0: continue
    rlv, plv = rl[i][gs], pl[i][gs]
    idn = (~np.isnan(rlv)) & (rlv != 0)
    mt = idn & (~np.isnan(plv)) & (plv != 0) & (np.sign(plv) == np.sign(rlv))
    absl = np.abs(np.nan_to_num(plv, nan=-np.inf)); sig = ps[gs]
    base.append(reach_of(np.lexsort((fr[gs], -absl, pp[i][gs], pq[i][gs], ~sig)), idn, mt, nc))
    eff_all = np.abs(np.nan_to_num(pl[i], nan=-np.inf))
    A, Bv = dv("k562_gwps", t), dv("orion_hek293t", t)
    ag_all = np.sqrt(np.abs(A) * np.abs(Bv)) * ((np.sign(A) == np.sign(Bv)) & (A != 0) & (Bv != 0))
    for tag, score in (("effect", eff_all), ("agree", ag_all)):
        rank = np.argsort(np.argsort(-score))
        for K in KS:
            pr = rank[gs] < K
            res[(tag, K)].append(reach_of(
                np.lexsort((fr[gs], -absl, pp[i][gs], pq[i][gs], ~pr, ~sig)), idn, mt, nc))
b = np.mean(base)
print(f"production baseline raw reach {b:.5f}\n")
print(f"{'K':>6}{'promote by |effect|':>22}{'promote by agreement':>24}")
for K in KS:
    e = np.mean(res[("effect", K)]); a = np.mean(res[("agree", K)])
    print(f"{K:>6}{(e-b)/(R-B)/6:>+22.5f}{(a-b)/(R-B)/6:>+24.5f}")
