"""094 step 1 — the anatomy of `reach`.

`reach` = mean over perturbations of k*/n_real, where the pool is the REFERENCE's
significant set, ordered by OUR confidence, and k* is the deepest prefix with
sign purity >= 0.9.

The pool and the sign-correctness of each gene are fixed by the model. The ONLY
free variable is the ORDER. This asks how much of our reach deficit is
sign accuracy (a modelling problem) and how much is ordering (a free channel,
because per-gene emission `spread` moves p-values without touching pseudobulk).
"""
import os, sys
import numpy as np

E092 = "~/Code/virtual-cell-challenge-2026/experiments/092-call-set-truncation"
d = np.load(os.path.join(E092, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rp, rq, rl = d["real_p"], d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
G, P = genes.size, perts.size
ALPHA = 0.05
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, np.int64); fr[np.argsort(genes)] = np.arange(G)

def reach_from_order(order, in_denom, mt, nc):
    nd = np.cumsum(in_denom[order]); nmt = np.cumsum(mt[order])
    with np.errstate(invalid="ignore", divide="ignore"):
        pur = np.where(nd > 0, nmt / np.maximum(nd, 1), np.nan)
    good = np.flatnonzero((nd > 0) & (pur >= 0.9))
    return (nd[good[-1]] / nc) if good.size else 0.0

rows = {k: [] for k in ("prod","lfc_only","p_only","oracle","antioracle","random","sig_last")}
acc = []; nreals = []; nsigs = []
rng = np.random.default_rng(0)

for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA)
    ps = ok & (pq[i] < ALPHA)
    gsel = np.flatnonzero(rs); nc = gsel.size
    if nc == 0: continue
    rlv, plv = rl[i][gsel], pl[i][gsel]
    in_denom = (~np.isnan(rlv)) & (rlv != 0)
    committed = (~np.isnan(plv)) & (plv != 0)
    mt = in_denom & committed & (np.sign(plv) == np.sign(rlv))
    sig_pred = ps[gsel]
    absl = np.abs(np.nan_to_num(plv, nan=-np.inf))
    nreals.append(nc); nsigs.append(int(sig_pred.sum()))
    if in_denom.sum(): acc.append(mt[in_denom].sum() / in_denom.sum())

    # production ordering (exactly the scorer's)
    o_prod = np.lexsort((fr[gsel], -absl, pp[i][gsel], pq[i][gsel], ~sig_pred))
    # |lfc| only, ignoring significance and p-value
    o_lfc  = np.lexsort((fr[gsel], -absl))
    # p-value only (drop the significant-first primary key)
    o_p    = np.lexsort((fr[gsel], -absl, pp[i][gsel], pq[i][gsel]))
    # oracle: correct signs first, then the rest
    o_orc  = np.lexsort((fr[gsel], ~mt))
    o_anti = np.lexsort((fr[gsel], mt))
    o_rnd  = rng.permutation(nc)
    # significant LAST (isolates what the significant-first key is worth)
    o_sl   = np.lexsort((fr[gsel], -absl, pp[i][gsel], pq[i][gsel], sig_pred))

    for k, o in (("prod",o_prod),("lfc_only",o_lfc),("p_only",o_p),("oracle",o_orc),
                 ("antioracle",o_anti),("random",o_rnd),("sig_last",o_sl)):
        rows[k].append(reach_from_order(o, in_denom, mt, nc))

print(f"perturbations scored: {len(rows['prod'])}/{P}")
print(f"n_real  : median {np.median(nreals):.0f}  mean {np.mean(nreals):.0f}  range [{min(nreals)}, {max(nreals)}]")
print(f"of the reference's DE genes, how many WE also call significant: median {np.median(nsigs):.0f} "
      f"({np.median(np.array(nsigs)/np.array(nreals))*100:.1f}% of the pool)")
print(f"\nSIGN ACCURACY on the reference's DE set: mean {np.mean(acc):.4f}  median {np.median(acc):.4f}")
print(f"  (0.9 is the purity pass mark. Genes correct on sign are the entire budget for k*.)\n")

BASE, SPAN = 0.07902, 0.8851   # fitted from the live leaderboard
print(f"{'ordering':<12}{'raw reach':>11}{'scaled':>10}{'vs prod (overall pts)':>24}")
base = np.mean(rows["prod"])
for k in ("prod","p_only","lfc_only","sig_last","random","oracle","antioracle"):
    v = np.mean(rows[k]); sc = (v - BASE)/SPAN
    print(f"{k:<12}{v:>11.5f}{sc:>10.4f}{((v-base)/SPAN/6):>+24.5f}")
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r1.npz"),
         **{k: np.array(v) for k, v in rows.items()},
         acc=np.array(acc), nreal=np.array(nreals), nsig=np.array(nsigs))
