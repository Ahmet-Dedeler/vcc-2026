"""Q36 — we score WORSE on per-gene sign than a target-INDEPENDENT shared direction.

Q35's control: HEK293T's effect for a RANDOM DIFFERENT target predicts HCT116's signs at
**0.6836**; matched targets give 0.7051.  So target-specific transfer is worth only +0.0215 and
the bulk of cross-context sign agreement is a shared, perturbation-independent program.

Our shipped prediction scores **0.5681** on the same pairs -- *worse than the shared component
alone*.  We remove that component deliberately: whitening and row-centering exist to strip it,
because it is non-discriminative and `pds` punishes it (129: "adding broad direction can make
each query more similar to decoys faster than to its own truth").

But `pds` is the member we lead and the DE members carry 96% of the gap.  So: add the shared
direction back at strength beta and price BOTH sides on the same surface.

The shared vector is the donor mean effect -- deployable, no target-context truth involved.
"""
import numpy as np, csv
from debench import Surface, overall_delta
from exact_pds import pds

hct = Surface("../095-clean-harness/de_thinned.npz")
canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
DS = "../104-context-heldout-credit-adapter/work/donor-stats"
mean_vec = np.load(f"{DS}/orion_hek293t-mean.npy")
print(f"donor shared-mean vector: shape {mean_vec.shape}")
# map the donor mean onto the harness gene axis via the DONOR's own gene list -- it is on the
# 38,606-column Orion axis, not the canonical 18,533.  A first version guarded on
# `mean_vec.size == len(canon)` and silently produced an all-zero vector, which read as a
# perfectly flat sweep.  Known-answer check below: beta must change sign accuracy.
dg = [l.strip() for l in open("../../data/donors/orion_hek293t/genes.txt") if l.strip()]
dpos = {g: i for i, g in enumerate(dg)}
shared = np.array([mean_vec[dpos[g]] if g in dpos else 0.0 for g in hct.genes], dtype=np.float64)
assert np.linalg.norm(shared) > 0, "shared vector is empty -- gene mapping failed"
print(f"shared vector: {int((shared != 0).sum())}/{shared.size} genes non-zero, "
      f"|.| = {np.linalg.norm(shared):.4f}")
shared = shared / np.linalg.norm(shared)

R = np.load("real_delta_hct116.npz", allow_pickle=True)
rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
pg = np.array([i for i, g in enumerate(hct.genes) if g in rg and g not in panel])
pcol = np.array([rg[hct.genes[i]] for i in pg])
prow = np.array([rt.get(t, -1) for t in hct.perts]); keep = prow >= 0
Dreal = R["delta"][np.ix_(prow[keep], pcol)].astype(np.float64)

P0 = hct.pred_lfc.copy()
scale = np.median(np.abs(P0[P0 != 0]))
print(f"{'beta':>7} {'sign acc':>9} {'head pur':>9} {'reach':>9} {'nmae':>9} {'pds':>9} "
      f"{'Δoverall(reach)':>16}")
base = None
for beta in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
    P = P0 + beta * scale * shared[None, :]
    hct.pred_lfc = P
    num = den = hnum = hden = 0
    for p in range(hct.perts.size):
        i = np.flatnonzero(hct.real_sig[p] & hct.adj[p])
        if i.size < 10: continue
        num += int((np.sign(P[p, i]) == np.sign(hct.real_lfc[p, i])).sum()); den += i.size
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        hnum += int((np.sign(P[p, m]) == np.sign(hct.real_lfc[p, m])).sum()); hden += m.size
    r = hct.reach(key=[np.abs(P[p]) for p in range(hct.perts.size)])
    n = hct.nmae(); pd_ = pds(P[np.ix_(keep, pg)].astype(np.float64), Dreal)
    if base is None: base = (r, n, pd_)
    print(f"{beta:7.2f} {num/den:9.4f} {hnum/max(hden,1):9.4f} {r:9.5f} {n:9.5f} {pd_:9.5f} "
          f"{overall_delta(d_reach=r-base[0]):16.5f}")
