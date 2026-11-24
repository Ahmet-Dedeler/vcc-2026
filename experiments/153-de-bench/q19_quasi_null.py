"""Q19 — the SAME null calculation under production's actual emitter (quasi, spread .25, phi 4).

emit_counts(mode="quasi"): per-cell count = (1-spread)*T/n  +  NB(mean=spread*T/n, phi), so the
count's mean is the full T/n and only `spread` of it is scattered.  Production emits at a fixed
library size EMIT_L = 100_000, so emitted CPM = count * 10.

q18 priced the minimum-variance emitter (95-97% of genes significant at zero effect).  This
prices what we ACTUALLY ship, which is the residual a resampling emitter would recover.
"""
import anndata as ad, numpy as np, scipy.sparse as sp

N_EMIT, EMIT_L, SPREAD, PHI, NDRAW = 400, 100_000, 0.25, 4.0, 4000
rng = np.random.default_rng(0)
for ctx in "A":
    A = ad.read_h5ad(f"../../data/challenge/context_{ctx}.h5ad")
    X = A.X.tocsc().astype(np.float64)
    lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
    D = (sp.diags(1e6 / lib) @ X).tocsc()
    m = D.shape[0]
    cpm = np.asarray(D.mean(0)).ravel()
    gated = np.flatnonzero(cpm > 5.0)
    fac = np.sqrt(12.0 * N_EMIT * m / (m + N_EMIT + 1.0))
    out = {}
    for label, mode in (("min-variance (all cells at the mean)", "even"),
                        ("production quasi spread .25 phi 4", "quasi")):
        fbar = np.empty(gated.size)
        for k, g in enumerate(gated):
            nz = np.sort(D.data[D.indptr[g]:D.indptr[g + 1]])
            nzero = m - nz.size
            mu_ct = cpm[g] / 10.0                      # expected emitted COUNT at EMIT_L
            if mode == "even":
                e = np.full(1, cpm[g])
            else:
                mu_s = SPREAD * mu_ct
                # gamma-Poisson: Poisson(Gamma(shape=1/phi, scale=phi*mu_s))
                lam = rng.gamma(1.0 / PHI, PHI * mu_s, NDRAW) if mu_s > 0 else np.zeros(NDRAW)
                e = ((1 - SPREAD) * mu_ct + rng.poisson(lam)) * 10.0   # -> CPM
            below = np.where(e > 0, nzero, 0) + np.searchsorted(nz, e, "left")
            eq = np.searchsorted(nz, e, "right") - np.searchsorted(nz, e, "left")
            eq = eq + np.where(e == 0, nzero, 0)
            fbar[k] = np.mean((below + 0.5 * eq) / m)
        z = (fbar - 0.5) * fac
        sig = np.abs(z) > 2.8
        out[label] = (fbar, z, sig)
        print(f"{label:38} Fbar med {np.median(fbar):.4f}  |z| med {np.median(np.abs(z)):5.1f}  "
              f"sig-at-zero {sig.sum():5d}/{gated.size} = {sig.mean():5.1%}  up {(z[sig]>0).mean():5.1%}")
    print(f"\nboard n_pred for a15 is ~4,400 of ~9,900 gated genes = ~44%.")
    print(f"quasi model predicts {out['production quasi spread .25 phi 4'][2].mean():.1%} "
          f"significant with NO effect at all.")
