"""Q18 — the emission artefact, quantified in closed form on the REAL challenge controls.

The scorer's DE is a Mann-Whitney U of our 400 emitted cells against the reference's 18,400
real control cells, "p-value from the normal approximation to the rank sum" (BRIEF section 3).
So for ANY emission law the whole call set is predictable analytically -- no emitting, no
seeds:

    Fbar = P(emitted > control) + 0.5 P(emitted == control)
    z    = (Fbar - 0.5) * sqrt(12 n m / (N + 1)),   n = 400, m = 18400

At n=400 / m=18400 that factor is 68.5, and clearing BH-0.05 over ~10k genes needs |z| ~ 2.8.
So **|Fbar - 0.5| > 0.041 makes a gene significant AT ZERO PREDICTED EFFECT.**

A minimum-variance emitter puts all 400 cells at the predicted MEAN.  For a right-skewed count
distribution the mean sits well above the median, so Fbar is far from 0.5 and the gene is
called significant with no effect at all.  A resampling emitter draws from the control's own
law, giving Fbar = 0.5 exactly by exchangeability -- the artefact is removed at its root rather
than damped.  This measures how big that is on the actual challenge controls.
"""
import anndata as ad, numpy as np, scipy.sparse as sp

n_emit, TSCALE = 400, 1e6
for ctx in "ABC":
    A = ad.read_h5ad(f"../../data/challenge/context_{ctx}.h5ad")
    X = A.X.tocsc().astype(np.float64)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
    D = (sp.diags(TSCALE / lib) @ X).tocsc()          # per-cell CPM, the values the scorer ranks
    m = D.shape[0]
    cpm_mean = np.asarray(D.mean(0)).ravel()
    gated = np.flatnonzero(cpm_mean > 5.0)
    fac = np.sqrt(12.0 * n_emit * m / (m + n_emit + 1.0))
    fbar = np.empty(gated.size)
    for k, g in enumerate(gated):
        col = D.data[D.indptr[g]:D.indptr[g + 1]]      # nonzero values only
        mu = cpm_mean[g]
        nz_below = int((col < mu).sum()); nz_eq = int((col == mu).sum())
        below = (m - col.size) + nz_below              # every implicit zero is below the mean
        fbar[k] = (below + 0.5 * nz_eq) / m
    z = (fbar - 0.5) * fac
    sig = np.abs(z) > 2.8
    print(f"context {ctx}: {gated.size} gated genes, mean-emitter null behaviour")
    print(f"   Fbar   median {np.median(fbar):.4f}   p10 {np.percentile(fbar,10):.4f}   "
          f"p90 {np.percentile(fbar,90):.4f}   (0.5 = correctly null)")
    print(f"   |z|    median {np.median(np.abs(z)):.1f}")
    print(f"   genes called significant AT ZERO EFFECT: {sig.sum()} / {gated.size} "
          f"= {sig.mean():.1%}   (of which {(z[sig]>0).mean():.1%} called UP)")
    del A, X, D
