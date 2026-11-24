"""Q7 — re-screen experiment 108 (signed regulon) on DE endpoints instead of PDS.

108 killed the DoRothEA signed-regulon family on a macro raw-PDS delta of -0.00035, with
"the network subspace itself is strongly lossy for this retrieval problem."  Retrieval is
`pds`, where we are 98th percentile on the board.  Live right now, `tydymy`'s
`regprop-collectri-v1` runs reach 0.2549 at pds 0.7732 -- the same family, PDS-negative and
reach-positive.  So the question 108 never asked is whether the transform helps the members
that carry 96% of our gap.

Compares 108's own three frozen arms against each other on the shared gene axis, so the
comparison is like-for-like.  No new fitting, no tuning -- the arms are frozen.
"""
import numpy as np
from debench import Surface, overall_delta

s = Surface("../095-clean-harness/de_thinned.npz")
z = np.load("../108-signed-regulon/work/frozen.npz", allow_pickle=True)
g108 = np.asarray([str(x) for x in z["genes"]])
t108 = np.asarray([str(x) for x in z["targets"]])

gi = {g: i for i, g in enumerate(g108)}
ti = {t: i for i, t in enumerate(t108)}
gcol = np.array([gi.get(g, -1) for g in s.genes])
prow = np.array([ti.get(p, -1) for p in s.perts])
gok, pok = gcol >= 0, prow >= 0
print(f"shared axis: {gok.sum()}/{s.genes.size} genes, {pok.sum()}/{s.perts.size} perturbations")

# Restrict the reference pool to the shared axis so every arm is scored on the same genes.
keep = np.zeros_like(s.real_sig)
keep[np.ix_(pok, gok)] = True
s.real_sig &= keep
s.adj &= keep

def arm(name, M):
    """M: (300, 7091) effect matrix -> per-(pert,gene) predicted direction+magnitude."""
    P = np.zeros_like(s.pred_lfc)
    P[np.ix_(pok, gok)] = M[np.ix_(prow[pok], gcol[gok])]
    s.pred_lfc = P
    key = [np.abs(P[p]) for p in range(len(s.perts))]
    acc = s.sign_accuracy()
    r = s.reach(key=key)
    return name, acc, r

print(f"\n{'arm':26} {'sign acc':>9} {'reach (|effect| key)':>21} {'Δreach':>9} {'Δoverall':>10}")
base_acc = base_r = None
for name, col in (("base (incumbent)", "base"), ("regulon reweight", "reweight"),
                  ("10% direction blend", "blend10"), ("regulon only", "reconstruction")):
    n, a, r = arm(name, z[col])
    if base_r is None:
        base_acc, base_r = a, r
    print(f"{n:26} {a:9.4f} {r:21.5f} {r-base_r:+9.5f} {overall_delta(d_reach=r-base_r):+10.5f}")

# 108's own endpoint, for reference: macro reweight-minus-base PDS was -0.00035.
print("\n108's frozen PDS endpoint (from work/results.json): macro reweight - base = -0.00035")
print("Exchange rate: +0.01 raw reach = +0.00188 overall; +0.01 raw pds = +0.00369 overall.")
# ERRATUM 2026-08-31: this line read +0.00019, exactly 10x low (reach span is 0.88512, so
# 0.01/0.88512/6 = 0.00188).  Caught by the peer session.  It was a bare print and every
# COMPUTED figure in this directory goes through overall_delta(), which always used the
# correct span -- so no reported result was affected.  This script also cannot run: 108's
# gallery has 0/132 panel overlap (see RESULTS.md section 7), so pok.sum() == 0 by
# construction.  Kept as the record of a blocked re-screen, not as a live script.
