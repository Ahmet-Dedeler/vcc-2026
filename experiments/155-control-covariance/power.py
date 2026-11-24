"""155 — is there any REPRODUCIBLE covariance signal in the controls at all?

The scored results can fail for two very different reasons and they need separating:

  (i)  the estimator recovers nothing -- the control covariance at this depth is noise;
  (ii) the estimator recovers something real that simply is not the knockdown direction.

Both diagnostics here use CONTROLS ONLY and never touch the reference:

  A. INFLATION.  Under the null of zero partial correlation, r ~ N(0, 1/(N-k-2)).  The ratio
     var(observed r) / (1/(N-k-2)) is the fraction of coefficient variance that is real.
     1.00 means pure noise.
  B. SPLIT-HALF.  Fit beta on two disjoint halves of the 18,400 control cells and correlate
     the two coefficient matrices.  This is the estimator's own reproducibility and it is an
     upper bound on anything it can transfer.
"""
from __future__ import annotations

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = ("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
         "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/ctrl155.npz")


def coefs(X, lib, tcol, k, pcs):
    cols = [np.ones(X.shape[0]), (np.log(lib) - np.log(lib).mean()) / np.log(lib).std()]
    if k:
        cols += [pcs[:, j] / pcs[:, j].std() for j in range(k)]
    Q = np.linalg.qr(np.column_stack(cols))[0].astype(np.float32)
    R = X.copy()
    R -= Q @ (Q.T @ R)
    T = R[:, tcol].astype(np.float64)
    tss = (T ** 2).sum(0)
    corr = np.empty((tcol.size, R.shape[1]))
    for c0 in range(0, R.shape[1], 1000):
        B = R[:, c0:c0 + 1000].astype(np.float64)
        corr[:, c0:c0 + B.shape[1]] = (T.T @ B) / np.sqrt(np.outer(tss, (B ** 2).sum(0)))
    return corr


def main() -> None:
    d = np.load(CACHE, allow_pickle=True)
    X, lib = d["X"], d["lib"]
    axis = np.asarray([str(g) for g in d["axis"]])
    targets = np.asarray([str(t) for t in d["targets"]])
    n_gated = int((np.asarray([str(g) for g in d["genes"]])).size)
    apos = {g: i for i, g in enumerate(axis)}
    tcol = np.array([apos[t] for t in targets])
    N = X.shape[0]

    from sklearn.decomposition import PCA
    pcs = PCA(n_components=50, svd_solver="randomized",
              random_state=0).fit_transform(X).astype(np.float64)

    print("A. coefficient-variance inflation over the pure-noise expectation")
    print(f"{'k PCs':>6} {'sd(r) observed':>16} {'sd(r) null':>12} {'inflation':>11} "
          f"{'real var frac':>14}")
    for k in (0, 5, 10, 20, 50):
        C = coefs(X, lib, tcol, k, pcs)[:, :n_gated]
        # drop each target's own column (r = 1 by construction)
        for i, c in enumerate(tcol):
            if c < n_gated:
                C[i, c] = np.nan
        obs = np.nanstd(C)
        null = 1.0 / np.sqrt(N - k - 2)
        print(f"{k:>6} {obs:16.5f} {null:12.5f} {(obs/null)**2:11.3f} "
              f"{1 - (null/obs)**2:14.3f}")
        del C

    print()
    print("B. split-half reproducibility of the coefficient matrix (controls only)")
    rng = np.random.default_rng(0)
    perm = rng.permutation(N)
    h1, h2 = perm[:N // 2], perm[N // 2:]
    print(f"{'k PCs':>6} {'r(half1, half2)':>17} {'r after double-centering':>26}")
    for k in (0, 5, 20, 50):
        A = coefs(X[h1], lib[h1], tcol, k, pcs[h1])[:, :n_gated]
        B = coefs(X[h2], lib[h2], tcol, k, pcs[h2])[:, :n_gated]
        m = np.isfinite(A) & np.isfinite(B)
        for i, c in enumerate(tcol):
            if c < n_gated:
                m[i, c] = False
        raw = np.corrcoef(A[m], B[m])[0, 1]

        def dbl(M):
            M = M - M.mean(0, keepdims=True)
            return M - M.mean(1, keepdims=True)

        ra, rb = dbl(A), dbl(B)
        cen = np.corrcoef(ra[m], rb[m])[0, 1]
        print(f"{k:>6} {raw:17.4f} {cen:26.4f}")
        del A, B, ra, rb


if __name__ == "__main__":
    main()
