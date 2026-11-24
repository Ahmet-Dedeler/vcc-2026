"""Q20 — the artefact with the emitter's ACTUAL integer output law and exact tie handling.

Both earlier estimates modelled a float point mass at the mean.  The real emitter does not do
that.  `emit_counts` quasi (donor_transfer.py:280+):

    L      = int(min(median(lib_sizes), max_counts))     # ONE fixed library for all cells
    totals = floor(n * L * profile_p)                    # integer count per gene
    even part: base, rem = divmod(t_even, n) -> (n-rem) cells at `base`, rem at `base+1`

So emitted CPM sits on a lattice of 1e6/L, while control CPM = count * 1e6 / lib_i with lib_i
varying 3,275-52,420 does not.  Exact ties are therefore concentrated at **zero**, where every
zero-count emitted cell ties every zero-count control cell -- and for a low-expression gene that
is most of both samples.  Mann-Whitney scores ties at 0.5, which pulls Fbar toward 0.5 and can
push a gene back inside the null band.  That is the fat low-|z| tail the two earlier runs
disagreed about.

Reports the integerised min-variance emitter (artefact's full size, computed consistently) and
production's quasi/phi=4 (what we carry).  The gap is the recoverable amount.
"""
import anndata as ad, numpy as np, scipy.sparse as sp

N, PHI, SPREAD, NDRAW = 400, 4.0, 0.25, 400
rng = np.random.default_rng(0)

def fbar_for(values, counts, nz_sorted, n_zero, m):
    """values: emitted CPM levels; counts: how many of the 400 cells sit at each."""
    tot = 0.0
    for e, c in zip(values, counts):
        if c == 0:
            continue
        if e == 0.0:
            below, eq = 0, n_zero                      # ties with every zero-count control
        else:
            lo = np.searchsorted(nz_sorted, e, "left")
            hi = np.searchsorted(nz_sorted, e, "right")
            below, eq = n_zero + lo, hi - lo
        tot += c * (below + 0.5 * eq) / m
    return tot / N

for ctx in "ABC":
    A = ad.read_h5ad(f"../../data/challenge/context_{ctx}.h5ad")
    X = A.X.tocsc().astype(np.float64)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
    L = int(np.median(lib))                            # what the emitter picks
    D = (sp.diags(1e6 / lib) @ X).tocsc()
    m = D.shape[0]
    cpm_mean = np.asarray(D.mean(0)).ravel()
    gated = np.flatnonzero(cpm_mean > 5.0)
    colsum = np.asarray(X.sum(0)).ravel()
    p = colsum / colsum.sum()                          # zero-effect profile = control pseudobulk
    fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
    step = 1e6 / L
    res = {}
    for label in ("integerised min-variance", "production quasi .25 / phi 4"):
        fb = np.empty(gated.size)
        for k, g in enumerate(gated):
            nz = np.sort(D.data[D.indptr[g]:D.indptr[g + 1]]); n_zero = m - nz.size
            T = int(np.floor(N * L * p[g]))
            if label.startswith("integer"):
                base, rem = divmod(T, N)
                vals = np.array([base, base + 1]) * step
                cts = np.array([N - rem, rem])
            else:
                t_rand = int(np.floor(T * SPREAD)); t_even = T - t_rand
                b, r = divmod(t_even, N)
                even = np.full(N, b); even[:r] += 1
                if t_rand > 0:
                    lam = rng.gamma(1.0 / PHI, PHI * (t_rand / N), N)
                    draw = rng.poisson(lam)
                else:
                    draw = np.zeros(N, dtype=int)
                c = even + draw
                vals, cts = np.unique(c * step, return_counts=True)
            fb[k] = fbar_for(vals, cts, nz, n_zero, m)
        z = (fb - 0.5) * fac; sig = np.abs(z) > 2.8
        res[label] = (fb, z, sig)
        print(f"ctx {ctx} L={L:6d}  {label:30} Fbar med {np.median(fb):.4f}  |z| med {np.median(np.abs(z)):5.1f}  "
              f"sig@0 {sig.sum():5d}/{gated.size} = {sig.mean():5.1%}  up {(z[sig]>0).mean() if sig.any() else 0:5.1%}")
    del A, X, D
