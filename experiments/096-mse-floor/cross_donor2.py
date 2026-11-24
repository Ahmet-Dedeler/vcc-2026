"""Cross-donor transferability, measured where there is signal to transfer.

v1 correlated each perturbation's effect across ALL ~8,000 shared genes and got rho=0.035.
That is almost certainly dilution, not absence: if 95% of genes carry no true response, a
Pearson over all of them is near zero even when the responding genes agree perfectly. Our
`pds_cosine` is the best on the leaderboard, which is impossible if nothing transfers.

So measure it three ways instead, each answering a different question:

  A. ON-TARGET. CRISPRi knocks down the target gene itself. Does the knockdown reproduce?
     This is the floor of transferability and should be strongly positive.
  B. TOP-K RESPONDERS. Take the union of each donor's top-K genes by |effect| for that
     perturbation and correlate there. This asks whether the genes that DO respond agree.
  C. MAGNITUDE-WEIGHTED. Correlate over all genes but weight by the geometric mean of the
     two donors' |effect|, so noise genes contribute in proportion to how little they move.
"""
import numpy as np, sys
sys.path.insert(0, "src")
from vcc2026.donor_transfer import Donor

A, B = sys.argv[1], sys.argv[2]
da, db = Donor.load(f"data/donors/{A}"), Donor.load(f"data/donors/{B}")
gb = {g: i for i, g in enumerate(db.genes)}
shared_g = [g for g in da.genes if g in gb]
ia = np.array([{g: i for i, g in enumerate(da.genes)}[g] for g in shared_g])
ib = np.array([gb[g] for g in shared_g])
pb = {p: i for i, p in enumerate(db.perts)}
pa_i = {p: i for i, p in enumerate(da.perts)}
shared_p = [p for p in da.perts if p in pb and p]
Xa = da.delta[np.array([pa_i[p] for p in shared_p])][:, ia]
Xb = db.delta[np.array([pb[p] for p in shared_p])][:, ib]
del da, db
Xa -= Xa.mean(0); Xb -= Xb.mean(0)          # specific part
gpos = {g: i for i, g in enumerate(shared_g)}
print(f"{len(shared_p)} shared perturbations, {len(shared_g)} shared genes\n")

# --- A. on-target knockdown -------------------------------------------------------------
ok = [(i, gpos[p]) for i, p in enumerate(shared_p) if p in gpos]
ra = np.array([Xa[i, j] for i, j in ok]); rb = np.array([Xb[i, j] for i, j in ok])
print(f"A. ON-TARGET ({len(ok)} perturbations whose target is on the shared axis)")
print(f"   {A}: mean self-effect {ra.mean():+.4f}   {B}: {rb.mean():+.4f}")
print(f"   correlation of the two knockdown depths: {np.corrcoef(ra, rb)[0,1]:+.4f}")
print(f"   both negative (real knockdown) in {100*np.mean((ra<0)&(rb<0)):.1f}% of cases\n")

# --- B. top-K responders ----------------------------------------------------------------
for K in (25, 100, 400):
    rs = []
    for i in range(len(shared_p)):
        a, b = Xa[i], Xb[i]
        idx = np.union1d(np.argpartition(-np.abs(a), K)[:K], np.argpartition(-np.abs(b), K)[:K])
        u, v = a[idx], b[idx]
        u = u - u.mean(); v = v - v.mean()
        d = np.sqrt((u@u)*(v@v))
        if d > 0: rs.append((u@v)/d)
    rs = np.array(rs)
    print(f"B. TOP-{K} UNION: mean rho {rs.mean():+.4f}  median {np.median(rs):+.4f}  "
          f"mean rho^2 {(rs**2).mean():.4f}  frac>0 {100*np.mean(rs>0):.1f}%")

# --- C. magnitude-weighted --------------------------------------------------------------
W = np.sqrt(np.abs(Xa) * np.abs(Xb))
num = (W * Xa * Xb).sum(1)
den = np.sqrt((W * Xa**2).sum(1) * (W * Xb**2).sum(1))
rw = num / np.maximum(den, 1e-30)
print(f"\nC. MAGNITUDE-WEIGHTED: mean rho {rw.mean():+.4f}  median {np.median(rw):+.4f}  "
      f"mean rho^2 {(rw**2).mean():.4f}")

# --- pseudobulk-level: does the AVERAGE over perturbations transfer? ---------------------
print(f"\nD. How much signal is there at all? per-perturbation effect norms:")
na, nb = np.linalg.norm(Xa, axis=1), np.linalg.norm(Xb, axis=1)
print(f"   ||s|| {A}: median {np.median(na):.3f}   {B}: {np.median(nb):.3f}")
print(f"   correlation of the two NORMS (does 'this perturbation is strong' transfer?): "
      f"{np.corrcoef(na, nb)[0,1]:+.4f}")
