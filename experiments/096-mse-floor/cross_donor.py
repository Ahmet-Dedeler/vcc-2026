"""How well does a perturbation's effect in one cell line predict it in another?

This sets the CEILING on donor transfer as an approach, and it is the quantity the mse
member actually reads.

The mse baseline `b` is *the context's mean perturbation response assigned identically to
every perturbation*. So beating the baseline does NOT require getting the mean response
right -- it requires getting the perturbation-SPECIFIC deviation right. We therefore split
each donor's effect matrix into

    delta[p] = m + s[p],        m = mean over perturbations,  s[p] = the specific part

and correlate the *specific* parts across donors. Correlating raw `delta` instead would
report the shared mean response as if it were transferable signal, which is exactly the
error that makes a model look better than the baseline it cannot beat.
"""
import numpy as np, sys, json
from pathlib import Path
sys.path.insert(0, "src")
from vcc2026.donor_transfer import Donor

A, B = sys.argv[1], sys.argv[2]
da, db = Donor.load(f"data/donors/{A}"), Donor.load(f"data/donors/{B}")
print(f"{A}: delta{da.delta.shape}   {B}: delta{db.delta.shape}")

# shared gene axis
ga = {g: i for i, g in enumerate(da.genes)}
gb = {g: i for i, g in enumerate(db.genes)}
shared_g = [g for g in da.genes if g in gb]
ia = np.array([ga[g] for g in shared_g]); ib = np.array([gb[g] for g in shared_g])
# shared perturbations
pa = {p: i for i, p in enumerate(da.perts)}
pb = {p: i for i, p in enumerate(db.perts)}
shared_p = [p for p in da.perts if p in pb and p]
print(f"shared genes {len(shared_g)}   shared perturbations {len(shared_p)}")

Xa = da.delta[np.array([pa[p] for p in shared_p])][:, ia]
Xb = db.delta[np.array([pb[p] for p in shared_p])][:, ib]
del da, db

def stats(Xa, Xb, label):
    # per-perturbation Pearson across genes
    A_ = Xa - Xa.mean(1, keepdims=True); B_ = Xb - Xb.mean(1, keepdims=True)
    num = (A_ * B_).sum(1); den = np.sqrt((A_**2).sum(1) * (B_**2).sum(1))
    r = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    r = r[np.isfinite(r)]
    print(f"\n  {label}: n={r.size}")
    print(f"    mean rho = {r.mean():+.4f}   median = {np.median(r):+.4f}   "
          f"mean rho^2 = {(r**2).mean():.4f}")
    for q in (10, 25, 50, 75, 90):
        print(f"      p{q:<3} {np.percentile(r, q):+.4f}", end="")
    print()
    return r

print("\n=== RAW delta (includes the shared mean response -- NOT what the metric rewards) ===")
stats(Xa, Xb, "raw")

ma, mb = Xa.mean(0), Xb.mean(0)
Sa, Sb = Xa - ma, Xb - mb
print("\n=== SPECIFIC part s = delta - mean (this is what beats the mse baseline) ===")
r = stats(Sa, Sb, "specific")

# what mse would a perfectly-calibrated transfer achieve?
# best global alpha on the specific part: mse_ratio ~ 1 - rho^2 (see 096/PREDICTION.md)
print(f"\n  implied best-case mse_raw with a global amplitude on the specific part:")
print(f"    using mean rho^2      : {1-(r**2).mean():.4f}")
print(f"    using (mean rho)^2    : {1-r.mean()**2:.4f}")
print(f"  the board baseline is ~0.993 and the leader scores 0.827")

# how much of the norm is the shared mean?
print(f"\n  ||mean response|| / mean ||delta||  =  "
      f"{np.linalg.norm(ma)/np.mean(np.linalg.norm(Xa,axis=1)):.4f} ({A})")
print(f"  fraction of {A}'s effect variance that is the shared mean: "
      f"{(np.linalg.norm(ma)**2*Xa.shape[0])/ (Xa**2).sum():.4f}")
