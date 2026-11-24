"""157 Q1 — the zero-effect call artefact under the SHIPPED tier spread, on the real panel.

153 §15 (`q20_integer_null.py`) measured how many genes the emitter calls significant at ZERO
predicted effect, with the exact integer output law and exact Mann-Whitney tie handling, on
`data/challenge/context_{A,B,C}.h5ad`.  It reported the minimum-variance emitter and a FLAT
`spread 0.25 / phi 4` law -- i.e. a8's.  a12/a14/a15 do not ship that: they ship
`--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25`, so roughly 8,000 gated genes sit at
0.25 and the remaining ~1,900 sit at **1.0**.

157's HCT116 sweep found that the spread-1.0 tail block supplies ~80% of the first ten genes
in `reach`'s order.  This checks that mechanism in the challenge's own regime, which no
harness can represent, by computing the artefact rate at each of the two spreads and
reporting the shipped mixture.

Closed form, exact ties, one gamma-Poisson draw per gene.  Reads only the CONTROL cells the
challenge gives us -- no reference truth exists for these contexts.  No cell-eval2, no build,
no submission.
"""
from __future__ import annotations
import anndata as ad
import numpy as np
import scipy.sparse as sp

N, PHI = 400, 4.0
KEEP_TOP = 8000          # production's --keep-top
rng = np.random.default_rng(0)


def fbar_for(values, counts, nz_sorted, n_zero, m):
    tot = 0.0
    for e, c in zip(values, counts):
        if c == 0:
            continue
        if e == 0.0:
            below, eq = 0, n_zero
        else:
            lo = np.searchsorted(nz_sorted, e, "left")
            hi = np.searchsorted(nz_sorted, e, "right")
            below, eq = n_zero + lo, hi - lo
        tot += c * (below + 0.5 * eq) / m
    return tot / N


def main() -> None:
    print(f"{'ctx':>4}{'gated':>8}{'spread/phi':>11}{'Fbar med':>10}{'|z| med':>9}"
          f"{'sig@0':>9}{'rate':>8}{'up':>8}")
    for ctx in "ABC":
        A = ad.read_h5ad(f"../../data/challenge/context_{ctx}.h5ad")
        X = A.X.tocsc().astype(np.float64)
        lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
        L = int(np.median(lib))
        D = (sp.diags(1e6 / lib) @ X).tocsc()
        m = D.shape[0]
        gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
        colsum = np.asarray(X.sum(0)).ravel()
        p = colsum / colsum.sum()          # zero-effect profile = the control pseudobulk
        fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
        step = 1e6 / L
        rate = {}
        # (spread, phi). The tail block exists to raise emitted dispersion for `mse`; the
        # (1.0, 0) row asks whether the artefact comes from the spread or from the phi, i.e.
        # whether the dispersion can be kept while the hyper-significance is dropped.
        for spread, phi in ((0.25, PHI), (1.0, PHI), (1.0, 0.0), (0.25, 0.0)):
            fb = np.empty(gated.size)
            for k, g in enumerate(gated):
                nz = np.sort(D.data[D.indptr[g]:D.indptr[g + 1]])
                n_zero = m - nz.size
                T = int(np.floor(N * L * p[g]))
                t_rand = int(np.floor(T * spread)); t_even = T - t_rand
                b, r = divmod(t_even, N)
                even = np.full(N, b); even[:r] += 1
                if t_rand > 0 and phi > 0:
                    draw = rng.poisson(rng.gamma(1.0 / phi, phi * (t_rand / N), N))
                elif t_rand > 0:
                    draw = rng.poisson(np.full(N, t_rand / N))
                else:
                    draw = np.zeros(N, dtype=int)
                vals, cts = np.unique((even + draw) * step, return_counts=True)
                fb[k] = fbar_for(vals, cts, nz, n_zero, m)
            z = (fb - 0.5) * fac
            sig = np.abs(z) > 2.8          # BH-0.05 over ~10k genes needs |z| ~ 2.8
            rate[(spread, phi)] = sig.mean()
            print(f"{ctx:>4}{gated.size:>8}{f'{spread:.2f}/phi{phi:g}':>11}{np.median(fb):>10.4f}"
                  f"{np.median(np.abs(z)):>9.1f}{int(sig.sum()):>9d}{sig.mean():>8.1%}"
                  f"{((z[sig] > 0).mean() if sig.any() else 0):>8.1%}")
        ntail = max(gated.size - KEEP_TOP, 0)
        mix = (KEEP_TOP * rate[(0.25, PHI)] + ntail * rate[(1.0, PHI)]) / gated.size
        alt = (KEEP_TOP * rate[(0.25, PHI)] + ntail * rate[(1.0, 0.0)]) / gated.size
        print(f"{ctx:>4}{gated.size:>8}{'SHIPPED':>11}{'':>10}{'':>9}"
              f"{'':>9}{mix:>8.1%}        (8,000 at 0.25/phi4 + {ntail} at 1.0/phi4)")
        print(f"{ctx:>4}{gated.size:>8}{'tail phi0':>11}{'':>10}{'':>9}"
              f"{'':>9}{alt:>8.1%}        (8,000 at 0.25/phi4 + {ntail} at 1.0/phi0)")
        del A, X, D


if __name__ == "__main__":
    main()
