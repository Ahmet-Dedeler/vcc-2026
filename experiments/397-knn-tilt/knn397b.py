#!/usr/bin/env python3
"""397b -- what else the neighbourhood confidence buys on H1 (analytic, 397's frame): (a) reach order that promotes kNN-confident
block genes; (b) head genes: does design/kNN sign agreement predict correctness; (c) nmae from confidence-targeted magnitudes.
kNN25 shrink50 over GWPS + K562-essential + Jurkat (397's construction)."""
import numpy as np, json
exec(open("~/Code/virtual-cell-challenge-2026/experiments/397-knn-tilt/knn397.py").read().split("# ---- kNN tilts")[0].replace('print("KNOWN ANSWER', 'pass  # print("KNOWN ANSWER'))
has = np.array([t in gwp for t in targets]); K = 25; m = 50.0
PK = np.full((nT, nG), np.nan)
for i, t in enumerate(targets):
    if not has[i]: continue
    v = GWn[gwp[t]]; s_ = GWn @ v; s_[gwp[t]] = -np.inf; nb = np.argsort(-s_)[:K]
    P, nn = votes(nb, ["k562_essential", "nadig_jurkat"])
    PK[i] = np.where(np.isfinite(P) & np.isfinite(P_marg_3), (nn * np.nan_to_num(P) + m * np.nan_to_num(P_marg_3)) / (nn + m), np.where(np.isfinite(P), P, P_marg_3))
conf = np.abs(PK - 0.5); sgk = np.where(np.isfinite(PK), np.where(PK > 0.5, -1.0, 1.0), sg4[None, :])
SGK = np.where(cb, sgk, SG0)   # 399's signs: kNN on the called block, delivered elsewhere
truth = np.sign(lrm)
# (b) head genes: correctness of the delivered sign by kNN agreement / confidence
hc = head & sig & adjud & ~own
agree = (SG0 == sgk); right = (SG0 == truth)
for lo, hi in ((0, .05), (.05, .1), (.1, .2), (.2, .3), (.3, .5)):
    mk = hc & np.isfinite(PK) & (conf >= lo) & (conf < hi)
    if mk.sum() == 0: continue
    print(f"HEAD conf[{lo},{hi}): n={mk.sum():6d}  P(right|agree)={right[mk & agree].mean():.3f} (n={int((mk & agree).sum())})  P(right|disagree)={right[mk & ~agree].mean():.3f} (n={int((mk & ~agree).sum())})")
bk = cb & adjud & np.isfinite(PK)
for lo, hi in ((0, .05), (.05, .1), (.1, .2), (.2, .3), (.3, .5)):
    mk = bk & (conf >= lo) & (conf < hi)
    print(f"BLOCK conf[{lo},{hi}): n={mk.sum():6d}  kNN sign precision={(sgk == truth)[mk].mean():.3f}  delivered={(SG0 == truth)[mk].mean():.3f}  in ref gate={pool[mk].mean():.3f}")
# (a) reach: promote confident block genes to the front of the order
def order_with_front(front_mask, key):
    cur2 = np.empty_like(cur)
    for i in range(nT):
        k = np.where(front_mask[i], -1e6 - key[i], 0.0)   # most confident first, then the delivered order
        cur2[i] = np.lexsort((cur[i].argsort(), k))  # stable: ties keep delivered rank
    return cur2
base_reach = walk_reach(cur, adjud & (SGK == truth) & pool, adjud & pool, n_conf)
print("reach with 399 signs, delivered order:", round(float(np.nanmean(base_reach)), 5), " (E0 signs:", round(float(np.nanmean(M0['reach'])), 5), ")")
rank_cur = np.argsort(cur, axis=1)
for tau in (.25, .3, .35, .4, .45):
    fm = cb & np.isfinite(PK) & (conf >= tau)
    cur2 = np.empty_like(cur)
    for i in range(nT):
        k = np.where(fm[i], -1e6 + (0.5 - conf[i]), rank_cur[i].astype(float))
        cur2[i] = np.argsort(k, kind="stable")
    R = walk_reach(cur2, adjud & (SGK == truth) & pool, adjud & pool, n_conf); d = R - base_reach; okk = np.isfinite(d)
    D = rng.integers(0, okk.sum(), (4000, okk.sum())); bs = d[okk][D].mean(1)
    print(f"front tau={tau}: genes promoted/target median {int(np.median(fm.sum(1)))}  reach {np.nanmean(R):.5f}  Δ {np.nanmean(d):+.5f} [{np.quantile(bs,.025):+.5f}, {np.quantile(bs,.975):+.5f}]")
# (c) nmae from confidence-targeted block magnitudes (signs = 399's)
def nmae_of(SG, MAG):
    X = SG * MAG; return np.array([np.abs(X[k, GIDX[j]] - LREF[j]).sum() / np.abs(LREF[j]).sum() for j, k in enumerate(elig)])
base_n = nmae_of(SGK, np.abs(E0)); print("nmae 399 signs |E0|:", round(float(base_n.mean()), 5), " E0:", round(float(M0['nmae'].mean()), 5))
for tau in (.2, .3, .4):
    for a in (2.0, 4.0, 8.0):
        MAG = np.where(cb & np.isfinite(PK) & (conf >= tau), a * np.abs(E0), np.abs(E0))
        n = nmae_of(SGK, MAG); d = n - base_n; D = rng.integers(0, len(d), (4000, len(d))); bs = d[D].mean(1)
        print(f"block mag x{a} where conf>={tau}: nmae {n.mean():.5f}  Δ {d.mean():+.5f} [{np.quantile(bs,.025):+.5f}, {np.quantile(bs,.975):+.5f}]  (median |E0| on those genes {np.median(np.abs(E0)[cb & np.isfinite(PK) & (conf >= tau)]):.3f} log2)")
