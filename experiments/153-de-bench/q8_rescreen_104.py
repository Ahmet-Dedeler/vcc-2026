"""Q8 — re-screen experiment 104 (context-heldout credit adapter) on DE endpoints.

104's verdict: "The full STRING/context ridge improved continuous signed credit and PDS on
all three held surfaces, but it reduced H1 top-20% sign precision by 0.0414 and failed the
frozen permutation-null gate on H1 and HCT116.  The signal is real but not safe enough to
ship."  Its endpoints were signed credit, PDS, and a top-20% sign-precision gate -- none of
which is `reach`, `nmae`, `fid` or `jac`.

104's frozen HCT116 artifact carries `baseline`, `agreement` and `full` effect matrices on
the SAME 132 panel targets and the full 18,533-gene axis, so the arms can be scored against
the depth-repaired harness's own reference DE table.

CAVEAT: the reference here is orion-hct116-thinned (150 cells/perturbation, depth ratio
1.0004), not the full held HCT116 that 104 opened.  Same cell line, same targets, clean
depth -- adequate for a directional screen, not a promotion gate.  Nothing here authorizes a
build; this is a re-analysis of a closed experiment on a different endpoint.
"""
import numpy as np
from debench import Surface, overall_delta

s = Surface("../095-clean-harness/de_thinned.npz")
z = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = np.asarray([str(x) for x in z["genes"]]); t104 = np.asarray([str(x) for x in z["targets"]])
gi = {g: i for i, g in enumerate(g104)}; ti = {t: i for i, t in enumerate(t104)}
gcol = np.array([gi.get(g, -1) for g in s.genes]); prow = np.array([ti.get(p, -1) for p in s.perts])
gok, pok = gcol >= 0, prow >= 0
print(f"shared axis: {gok.sum()}/{s.genes.size} genes, {pok.sum()}/{s.perts.size} perturbations")

keep = np.zeros_like(s.real_sig); keep[np.ix_(pok, gok)] = True
s.real_sig &= keep; s.adj &= keep

def arm(M):
    P = np.zeros_like(s.pred_lfc)
    P[np.ix_(pok, gok)] = M[np.ix_(prow[pok], gcol[gok])]
    s.pred_lfc = P
    return s.sign_accuracy(), s.reach(key=[np.abs(P[p]) for p in range(len(s.perts))]), P

print(f"\n{'arm':30} {'sign acc':>9} {'reach':>9} {'Δreach':>9} {'Δoverall':>10}")
b_acc = b_r = None
for name, col in (("baseline (incumbent)", "baseline"),
                  ("agreement-only ridge", "agreement"),
                  ("full STRING/context ridge", "full")):
    a, r, _ = arm(z[col])
    if b_r is None: b_acc, b_r = a, r
    print(f"{name:30} {a:9.4f} {r:9.5f} {r-b_r:+9.5f} {overall_delta(d_reach=r-b_r):+10.5f}")

# Top-decile purity is what `reach` actually consumes -- 153 §1 showed the aggregate hides it.
print(f"\ntop-decile (by |effect|) sign accuracy on the reference DE pool")
for name, col in (("baseline", "baseline"), ("agreement", "agreement"), ("full", "full")):
    P = np.zeros_like(s.pred_lfc); P[np.ix_(pok, gok)] = z[col][np.ix_(prow[pok], gcol[gok])]
    num = den = 0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & s.adj[p])
        if i.size < 10: continue
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        num += int((np.sign(P[p, m]) == np.sign(s.real_lfc[p, m])).sum()); den += m.size
    print(f"  {name:12} {num/max(den,1):.4f}  (n={den})")
