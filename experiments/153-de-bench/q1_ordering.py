"""Q1 — is there any subset of genes where our sign is reliable enough for `reach`?

`reach` needs a prefix that is >=90% sign-pure.  Our aggregate sign accuracy on the
reference's DE set is 0.573.  If no stratum reaches 0.90, reach is structurally capped no
matter how good the ordering key is, and the whole ordering channel is a mirage.
"""
import numpy as np
from debench import Surface, overall_delta

s = Surface("../095-clean-harness/de_thinned.npz")
P = len(s.perts)

# ---- sign accuracy stratified by |predicted lfc| percentile, within each pert's DE pool --
print("sign accuracy on the reference's DE pool, by |predicted lfc| percentile within pert")
print(f"{'bucket':>12} {'n pairs':>9} {'sign acc':>9}")
buckets = [(90, 100), (75, 90), (50, 75), (25, 50), (0, 25)]
for lo, hi in buckets:
    num = den = 0
    for p in range(P):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10:
            continue
        a = np.abs(s.pred_lfc[p, i])
        l, h = np.percentile(a, lo), np.percentile(a, hi)
        m = i[(a >= l) & (a <= h)] if hi == 100 else i[(a >= l) & (a < h)]
        num += int((np.sign(s.pred_lfc[p, m]) == np.sign(s.real_lfc[p, m])).sum())
        den += m.size
    print(f"{f'p{lo}-{hi}':>12} {den:9d} {num/max(den,1):9.4f}")

# ---- the same, but the strongest key 095 §6 found: cross-donor agreement is unavailable
# here, so use the two computable proxies the table carries ------------------------------
print("\nreach under different ordering keys (deterministic, no emitter, no seed)")
keys = {
    "production (emitted p-order)": None,
    "|predicted lfc|":             [np.abs(s.pred_lfc[p]) for p in range(P)],
    "-predicted p_adj":            [-s.pred_padj[p] for p in range(P)],
    "random":                      [np.random.default_rng(0).random(s.genes.size) for _ in range(P)],
    "ORACLE (correct signs first)":[(np.sign(s.pred_lfc[p]) == np.sign(s.real_lfc[p])).astype(float)
                                    for p in range(P)],
}
base = None
print(f"{'key':32} {'raw reach':>10} {'Δoverall':>10}")
for name, k in keys.items():
    r = s.reach(key=k)
    if base is None:
        base = r
    print(f"{name:32} {r:10.5f} {overall_delta(d_reach=r-base):+10.5f}")

# ---- how deep does purity actually survive? --------------------------------------------
print("\nwhere the 90% purity bar bites, under the production order")
depths, pools = [], []
for p in range(P):
    pool = np.flatnonzero(s.real_sig[p])
    if pool.size == 0:
        continue
    o = s._order(p, pool, None, None)
    pool = pool[o]
    pool = pool[s.adj[p, pool]]
    if pool.size == 0:
        depths.append(0); pools.append(0); continue
    ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
    pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
    d = np.flatnonzero(pur >= 0.9)
    depths.append(int(d[-1] + 1) if d.size else 0)
    pools.append(pool.size)
depths, pools = np.array(depths), np.array(pools)
print(f"  perturbations with k* == 0 (first call wrong, never recovers): "
      f"{int((depths == 0).sum())}/{depths.size} = {(depths==0).mean():.3f}")
print(f"  median k* {np.median(depths):.1f} of median pool {np.median(pools):.1f}")
print(f"  mean k*/n_real (= raw reach) {(depths/np.maximum(pools,1)).mean():.5f}")
