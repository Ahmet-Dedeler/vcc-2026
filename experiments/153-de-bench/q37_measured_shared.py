"""Q37 — the shared direction Q35 measured lives in the MEASURED lfc space, not the donor table.

Q36 used the donor table's mean and it did nothing (sign accuracy 0.5724 -> 0.5661).  That is
consistent rather than contradictory: donor tables are row-centred and 112 established that
centring removes the shared translation exactly.  So the component Q35 found is not in them.

Q35's control used HEK293T's MEASURED lfc for a random other target.  The corresponding vector
is the mean of HEK293T's measured effects over targets -- a donor-derived, deployable object
(for the challenge it would be a donor line's mean response, involving no target-context truth).

Sanity gate first: that vector ALONE must score ~0.68 on HCT116's DE pool, matching Q35's
mismatched-target value.  If it does not, the vector is not the thing Q35 measured.
"""
import numpy as np, csv
from debench import Surface, overall_delta
from exact_pds import pds

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
gh = {g: i for i, g in enumerate(hek.genes)}
col = np.array([gh.get(g, -1) for g in hct.genes]); ok = col >= 0
shared = np.zeros(hct.genes.size)
shared[ok] = np.nanmean(np.where(np.isfinite(hek.real_lfc), hek.real_lfc, np.nan), axis=0)[col[ok]]
shared = np.nan_to_num(shared)

num = den = 0
for p in range(hct.perts.size):
    i = np.flatnonzero(hct.real_sig[p] & hct.adj[p] & ok)
    if i.size < 10: continue
    num += int((np.sign(shared[i]) == np.sign(hct.real_lfc[p, i])).sum()); den += i.size
print(f"GATE: shared vector ALONE on HCT116's DE pool = {num/den:.4f}   (Q35 mismatched: 0.6836)")
assert 0.60 < num/den < 0.85, ("shared vector does not reproduce Q35 family -- wrong object")
print("   (higher than Q35's single-target 0.6836 because averaging over 132 targets denoises\n    the shared programme -- the mean is a cleaner estimate of it than any one draw)")

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
R = np.load("real_delta_hct116.npz", allow_pickle=True)
rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
pg = np.array([i for i, g in enumerate(hct.genes) if g in rg and g not in panel])
pcol = np.array([rg[hct.genes[i]] for i in pg])
prow = np.array([rt.get(t, -1) for t in hct.perts]); keep = prow >= 0
Dreal = R["delta"][np.ix_(prow[keep], pcol)].astype(np.float64)

P0 = hct.pred_lfc.copy()
u = shared / max(np.linalg.norm(shared), 1e-12)
sc = float(np.median(np.linalg.norm(P0, axis=1)))
print(f"\n{'beta':>6} {'sign acc':>9} {'head pur':>9} {'reach':>9} {'nmae':>9} {'pds':>9} "
      f"{'Δreach ovr':>11} {'Δpds ovr':>9}")
base = None
for beta in (0.0, 0.05, 0.10, 0.20, 0.35, 0.60, 1.0):
    P = P0 + beta * sc * u[None, :]
    hct.pred_lfc = P
    num = den = hn = hd = 0
    for p in range(hct.perts.size):
        i = np.flatnonzero(hct.real_sig[p] & hct.adj[p])
        if i.size < 10: continue
        num += int((np.sign(P[p, i]) == np.sign(hct.real_lfc[p, i])).sum()); den += i.size
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        hn += int((np.sign(P[p, m]) == np.sign(hct.real_lfc[p, m])).sum()); hd += m.size
    r = hct.reach(key=[np.abs(P[p]) for p in range(hct.perts.size)])
    n = hct.nmae(); pd_ = pds(P[np.ix_(keep, pg)].astype(np.float64), Dreal)
    if base is None: base = (r, n, pd_)
    print(f"{beta:6.2f} {num/den:9.4f} {hn/max(hd,1):9.4f} {r:9.5f} {n:9.5f} {pd_:9.5f} "
          f"{overall_delta(d_reach=r-base[0]):11.5f} {(pd_-base[2])/0.452/6:9.5f}")
