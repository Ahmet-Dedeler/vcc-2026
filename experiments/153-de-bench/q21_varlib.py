"""Q21 — does NOT collapsing lib_sizes to their median dissolve the residual artefact?

`emit_counts` already takes per-cell `lib_sizes`; quasi collapses them to the median
(`donor_transfer.py:299`).  That collapse puts every emitted cell on one CPM lattice while real
cells span 16-fold, concentrating ties at zero.  If the collapse is the cause, allocating totals
proportional to the REAL library sizes should move Fbar toward 0.5 on its own -- a change to an
existing function rather than a resampling rebuild.

THE ALLOCATION MUST PRESERVE THE GENE TOTAL.  A first version of this file floored the expected
count per cell independently, which silently destroys up to one count per cell, drove Fbar to
0.39 and inverted the call direction.  It disagreed with q20 on the IDENTICAL fixed-L
configuration, which is how it was caught -- so the fixed-L row is retained below purely as a
self-check and must reproduce q20's 19.0%.
"""
import anndata as ad, numpy as np, scipy.sparse as sp

N, PHI, SPREAD = 400, 4.0, 0.25
rng = np.random.default_rng(0)

def largest_remainder(w, T):
    """Split integer total T across cells proportional to w, exactly."""
    if T <= 0: return np.zeros(w.size, dtype=np.int64)
    q = w / w.sum() * T
    fl = np.floor(q).astype(np.int64)
    r = int(T - fl.sum())
    if r > 0: fl[np.argsort(-(q - fl))[:r]] += 1
    return fl

def fbar(vals, cts, nz, n_zero, m):
    tot = 0.0
    for e, c in zip(vals, cts):
        if c == 0: continue
        if e == 0.0: below, eq = 0, n_zero
        else:
            lo = np.searchsorted(nz, e, "left"); hi = np.searchsorted(nz, e, "right")
            below, eq = n_zero + lo, hi - lo
        tot += c * (below + 0.5 * eq) / m
    return tot / N

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsc().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
L = int(np.median(lib))
D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
colsum = np.asarray(X.sum(0)).ravel(); p = colsum / colsum.sum()
fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
libs_v = rng.choice(lib, N, replace=False)
print(f"control libraries: min {lib.min():.0f} median {L} max {lib.max():.0f} "
      f"({lib.max()/lib.min():.1f}-fold)\n")

for label, varlib, disp in (("min-variance, fixed L", False, False),
                            ("production quasi .25/phi4, fixed L  [must match q20: 19.0%]", False, True),
                            ("quasi .25/phi4, VARYING libraries", True, True),
                            ("min-variance, VARYING libraries", True, False)):
    fb = np.empty(gated.size)
    for k, g in enumerate(gated):
        nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
        Ls = libs_v if varlib else np.full(N, float(L))
        T = int(np.floor(Ls.sum() * p[g]))
        if disp:
            t_r = int(np.floor(T * SPREAD)); t_e = T - t_r
            c = largest_remainder(Ls, t_e).astype(np.float64)
            if t_r > 0:
                mu_r = Ls / Ls.sum() * t_r
                lam = rng.gamma(1.0 / PHI, PHI * np.maximum(mu_r, 1e-12))
                c = c + rng.poisson(lam)
        else:
            c = largest_remainder(Ls, T).astype(np.float64)
        cpm = c * 1e6 / Ls
        vals, cts = np.unique(cpm, return_counts=True)
        fb[k] = fbar(vals, cts, nz, n_zero, m)
    z = (fb - 0.5) * fac; sig = np.abs(z) > 2.8
    print(f"{label:56} Fbar med {np.median(fb):.4f}  sig@0 {sig.sum():5d}/{gated.size} "
          f"= {sig.mean():5.1%}  up {(z[sig]>0).mean() if sig.any() else 0:5.1%}")
