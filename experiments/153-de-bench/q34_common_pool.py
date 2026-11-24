"""Q34 — Q33's 0.57-vs-0.71 gap on a COMMON gene pool.  Q33 compared different pools.

Q32 restricted to genes passing HEK293T's own 5-CPM gate (9,645 pairs).  Q33's donor-table run
admitted every gene the 38,606-column table carries (11,279 pairs), which pulls in low-expression
coordinates where sign agreement is much worse.  Comparing those two numbers was my error.

Here every source is scored on the SAME (target, gene) pairs: HCT116's DE pool, intersected with
HEK293T's gated axis and with both donor tables' axes.
"""
import numpy as np, csv
from debench import Surface

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
DD = "../../data/donors"
gh = {g: i for i, g in enumerate(hek.genes)}; th = {t: i for i, t in enumerate(hek.perts)}
src = {}
for name in ("orion_hek293t", "k562_gwps"):
    gs = [l.strip() for l in open(f"{DD}/{name}/genes.txt") if l.strip()]
    ps = [l.strip() for l in open(f"{DD}/{name}/perts.txt") if l.strip()]
    src[name] = ({g: i for i, g in enumerate(gs)}, {t: i for i, t in enumerate(ps)},
                 np.load(f"{DD}/{name}/delta.npy", mmap_mode="r"))

common_g = [i for i, g in enumerate(hct.genes)
            if g in gh and all(g in src[n][0] for n in src)]
tgt = [p for p, t in enumerate(hct.perts)
       if t in th and all(t in src[n][1] for n in src)]
print(f"common pool: {len(tgt)} targets, {len(common_g)} genes\n")

acc = {k: [0, 0] for k in ("HEK measured (harness truth)", "donor orion_hek293t",
                           "donor k562_gwps", "donor mean (equal weights)",
                           "our shipped prediction")}
for p in tgt:
    pool = np.array([g for g in common_g if hct.real_sig[p, g] and hct.adj[p, g]])
    if pool.size < 10: continue
    a = np.sign(hct.real_lfc[p, pool])
    cand = {}
    cand["HEK measured (harness truth)"] = np.sign(
        hek.real_lfc[th[hct.perts[p]], [gh[hct.genes[g]] for g in pool]])
    ds = {}
    for n in src:
        gi, tiy, dl = src[n]
        ds[n] = np.asarray(dl[tiy[hct.perts[p]], [gi[hct.genes[g]] for g in pool]], dtype=np.float64)
        cand["donor " + n] = np.sign(ds[n])
    z = {n: ds[n] / max(np.linalg.norm(ds[n]), 1e-12) for n in ds}   # equal-weight after L2
    cand["donor mean (equal weights)"] = np.sign(sum(z.values()))
    cand["our shipped prediction"] = np.sign(hct.pred_lfc[p, pool])
    for k, b in cand.items():
        acc[k][0] += int(((a == b) & (b != 0)).sum()); acc[k][1] += int((b != 0).sum())

print(f"{'source':34} {'sign acc':>9} {'n pairs':>10}")
for k, (n_, d_) in acc.items():
    print(f"{k:34} {n_/max(d_,1):9.4f} {d_:10,}")
