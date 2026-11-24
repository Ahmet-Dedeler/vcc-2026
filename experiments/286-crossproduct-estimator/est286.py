"""286 -- the CROSS-PRODUCT response estimator: shared response DIRECTIONS learned from cross-products of two
INDEPENDENT views, with a full-rank residual path.  FROZEN before any held number (this file's sha256 is recorded in
every output; every consumer asserts it did not change between folds / stages).

WHY THIS OBJECT IS DIFFERENT FROM 271 AND 272
---------------------------------------------
271 tested SCALAR reliability weights and a DIAGONAL per-gene Wiener on the donor rows (0/10 cells).  272 re-weighted
the EIGENVALUES of the FITTED whitener's eigenbasis (0/11) and named the object it did not measure: "a noise correction
that re-estimates the eigenvectors ... is a different object" (R272 section 5).  272's mechanism is the reason to
expect that object to differ: on every operator on disk the fitted eigenvectors are 77-79 % noise directions carrying
87-92 % of the whitened prediction energy, and `pds` NEEDS those directions because the target-specific signal is
spread over them (R272 section 0.10).

A cross-product between two INDEPENDENT views a, b of the same perturbation has expectation Sigma_signal -- the
sampling noise is independent across views and averages out -- so the DIRECTIONS of the shared response can be learned
without the within-view noise inflating them, and each target's response is then estimated by a FULL-MATRIX
combination (not a diagonal one) with a FLOOR, so that every direction keeps a path.

THE MOMENTS (`accumulate`, `moments`)
-------------------------------------
For a bank with two independent views A, B of each row (each view centred on its OWN mean over the bank's rows, in the
production row's units), over n rows and m columns, write U = (A + B)/2 and D = A - B.  Then

    S-hat = sym(A^T B) / n              = (U^T U - D^T D / 4) / n        the SIGNAL covariance
    N-hat = (A - B)^T (A - B) / (4 n)   = D^T D / (4 n)                  the noise covariance OF THE MEAN ROW
    C-hat = (A^T A + B^T B) / (2 n)     = (U^T U + D^T D / 4) / n        the half row's covariance

(the second form of each is the identity used to accumulate: two PSD Grams, G_U = U^T U and G_D = D^T D, so S-hat is
exactly symmetric by construction and no `sym()` step is needed).  With A = s + e_A, B = s + e_B and e independent,
E[S-hat] = Sigma_s and E[N-hat] = Var(noise of (A+B)/2) = 271's identity Var((A+B)/2) = Var(A-B)/4.  The EXACT
algebraic identity is therefore

    C-hat = S-hat + 2 N-hat          (KA-M, machine precision; it holds row by row, no statistics involved)

and the mandate's literal moment residual ||C-hat - S-hat - N-hat||_F / ||C-hat||_F is reported alongside it: it is
NOT zero, it equals ||N-hat||_F / ||C-hat||_F, i.e. the noise share of the half row's covariance.  Both are stated.

S-hat is symmetric but NOT PSD (the cross-product of two noisy views has negative eigenvalues wherever the shared
signal is below the sampling error).  `psd_eigh` projects it: S+ = V diag(max(lam, 0)) V^T, and reports the clipped
energy sum|lam<0| / sum|lam|.

THE ESTIMATOR FAMILY (frozen; every row arm is scored BOTH norm-matched `:nm` and literal `:lit`)
------------------------------------------------------------------------------------------------
E0  `mean`       the equal-weight mean of the donor rows = the incumbent (known answer at 0.0e+00).

E1  `e1:phi`     FULL-MATRIX Wiener on the row with a FLOOR.  With N' = N-hat / J + eps I (J = the number of
                 independent views already averaged into the row being filtered; J = 1 for a single donor's
                 production row, which IS the mean of its two halves, so N-hat is already that row's noise
                 covariance), G = S+ (S+ + N')^-1.  Factorise N' = L L^T (Cholesky) and K = L^-1 S+ L^-T = Q diag(k)
                 Q^T (symmetric eigenproblem); then EXACTLY

                     G = L Q diag(k / (1 + k)) (L Q)^-1

                 so G's eigenvalues are the Wiener gains g_i = k_i / (1 + k_i) in [0, 1) and its eigenbasis is L Q
                 (the basis that simultaneously diagonalises S+ and N').  The FLOOR raises every gain:
                 g_i -> max(g_i, phi), phi in PHIS -- the full-rank residual path, NO direction is zeroed.
                 phi = 0 is the un-floored Wiener, reported as the comparison.
                 x-hat = G_phi x is applied by three cheap steps (triangular solve, diagonal, triangular multiply);
                 the m x m matrix G is never formed.

E2  `e2:r`       the LOW-RANK version the review warned against: S-hat truncated to rank r (top-r eigenpairs of S+),
                 NO floor.  G_r = S_r (S_r + N')^-1 by Woodbury (exact, O(m r) after one Cholesky solve).

E3  `e3:c`       the SIGNAL-COVARIANCE WHITENER: 062's operator recipe (op272.fit_operator's arithmetic) with the
                 full-row covariance replaced by S-hat + a floor, used IN THE OPERATOR'S PLACE on the incumbent rows:
                     W = V diag( (max(lam_i, 0) + c * nu_i + eps)^(-gamma/2) ) V^T,   nu_i = v_i^T N-hat v_i
                 with production's gamma = 0.65 and production's on-axis preserve (ONE global scalar), so the arm is
                 shape-only up to the render.  c = 1 is the measured noise level (272's A' floor, but in S-hat's OWN
                 eigenbasis rather than the fitted operator's -- the corrected-DIRECTIONS object); c = 0 floors at
                 E3_LMAX_FLOOR * lam_max instead.

E4  `e4`         271's DIAGONAL per-gene Wiener (its (c) `wiener_own`, est271.wiener_own): the ALREADY-FAILED
                 baseline (271: -0.013 pessimistic).  Re-run through this experiment's plumbing.
E5  `e5`         272's eigenvalue floor in the FITTED operator's eigenbasis at c = 1 (op272.DampedOperator with
                 floor_c = 1): the second ALREADY-FAILED baseline (272: -0.0022 [-0.0040, -0.0004]).
E6  `e6`         272's properly cross-fitted operator (its construction B, `xfit_avg`) where 272's refit arrays are
                 on disk; skipped and reported as skipped where they are not.

CONTROLS (frozen with the estimators; every one is built by the SAME code path as the arm it controls)
-----------------------------------------------------------------------------------------------------
Cs `:shufpair`   SHUFFLED TARGET PAIRING: S-hat built from view a of row t and view b of row pi(t) for a fixed
                 derangement pi (expectation zero signal).  N-hat is rebuilt from the same shuffled pairing, so the
                 control differs from the arm ONLY in the pairing.  An arm on S_shuf must not gain.
Cw `:wrong`      WRONG-SOURCE S-hat: another donor line's S-hat on the shared columns (its N-hat kept), i.e. the
                 directions are response directions but of the wrong context.  A gain there means the covariance is
                 generic, not response-specific.
Cr `:rand`       an energy-matched random PSD covariance: W W^T scaled so that trace and Frobenius norm match S+
                 (W is n_rand x m standard normal, seed RAND_SEED, n_rand = the arm's n).
Cm `:mirror`     a mirrored gain G' = 2 I - G_phi (271/272's mirror: what the arm does to a direction is undone).

Everything downstream of the row is cd4svd5 verbatim (lib263: whitening through plolo.Operator.apply, CD4 deflated
k = 5 re-solved at f = 0.34, the H1-pair union re-solved by the carrying-row rule at g = 0.28, one add, one preserve).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import scipy.linalg as sla

# ---------------------------------------------------------------------------------------------- frozen constants
EPS_REL = 1e-4                 # eps I with eps = EPS_REL * trace(N-hat) / m  (a relative jitter on the noise metric)
PHIS = (0.0, 0.10, 0.25)       # the gain floor; 0.0 = the un-floored Wiener (the comparison), 0.10 / 0.25 = the arms
RANKS = (20, 100)              # E2's truncation ranks
E3_CS = (0.0, 1.0)             # E3's noise-floor multiple
E3_LMAX_FLOOR = 1e-6           # E3's floor at c = 0: E3_LMAX_FLOOR * lam_max
GAMMA = 0.65                   # production's whitening exponent
SHUF_SEED = 20260906
RAND_SEED = 20260906
DERANGE_SEED = 20260906
EIGH_DTYPE = np.float32        # scipy driver 'evd' in float32: ~3 min and ~5 GiB at m = 17.8k (f64 needs ~10 GiB)
EIGH_DRIVER = "evd"


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def derangement(n: int, seed: int = DERANGE_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if n < 2 or not np.any(p == np.arange(n)):
            return p


# ------------------------------------------------------------------------------------------------- the moments
def accumulate(rows_A, rows_B, m: int, *, perm: np.ndarray | None = None, block: int = 512, log=None):
    """G_U = U^T U and G_D = D^T D accumulated in float64 row blocks, U = (A + B)/2, D = A - B.

    `rows_A(lo, hi)` / `rows_B(lo, hi)` return (hi - lo, m) float64 blocks of the two views in the production row's
    units.  `perm` (the Cs control) pairs view a of row t with view b of row perm[t]: B is then gathered row by row.
    Returns (G_U, G_D, n).
    """
    G_U = np.zeros((m, m))
    G_D = np.zeros((m, m))
    n = 0
    lo = 0
    while True:
        A = rows_A(lo, lo + block)
        if A is None or A.shape[0] == 0:
            break
        hi = lo + A.shape[0]
        B = rows_B(lo, hi) if perm is None else rows_B(None, perm[lo:hi])
        U = 0.5 * (A + B)
        D = A - B
        G_U += U.T @ U
        G_D += D.T @ D
        n += A.shape[0]
        if log is not None:
            log(f"    block {lo}:{hi} ({n} rows)")
        lo = hi
    return G_U, G_D, n


def moments(G_U: np.ndarray, G_D: np.ndarray, n: int):
    """(S-hat, N-hat) in float64 from the two Grams; C-hat = S-hat + 2 N-hat exactly."""
    S = (G_U - 0.25 * G_D) / n
    N = G_D / (4.0 * n)
    return S, N


def moment_ka(S: np.ndarray, N: np.ndarray, G_U: np.ndarray, G_D: np.ndarray, n: int) -> dict:
    """KA-M: the exact identity C = S + 2N at machine precision, and the mandate's literal residual
    ||C - S - N||_F / ||C||_F (= ||N||_F / ||C||_F, the noise share of the half row's covariance)."""
    C = (G_U + 0.25 * G_D) / n
    fC = float(np.linalg.norm(C))
    exact = float(np.linalg.norm(C - S - 2.0 * N)) / fC
    literal = float(np.linalg.norm(C - S - N)) / fC
    return dict(frob_C=fC, frob_S=float(np.linalg.norm(S)), frob_N=float(np.linalg.norm(N)),
                ka_exact_C_eq_S_plus_2N=exact, literal_C_minus_S_minus_N=literal,
                trace_C=float(np.trace(C)), trace_S=float(np.trace(S)), trace_N=float(np.trace(N)),
                signal_share_of_trace=float(np.trace(S) / np.trace(C)) if np.trace(C) else float("nan"))


# ------------------------------------------------------------------------------------------------- eigen tools
def psd_eigh(S: np.ndarray, *, log=None):
    """Symmetric eigendecomposition of S-hat (float32, driver 'evd'), the PSD projection and the clipped energy.

    Returns (lam, V, audit) with lam ASCENDING (LAPACK order) and V the float32 eigenvectors (columns).
    """
    A = np.asarray(S, dtype=EIGH_DTYPE, order="F")
    lam, V = sla.eigh(A, driver=EIGH_DRIVER, overwrite_a=True, check_finite=False)
    del A
    lam = np.asarray(lam, dtype=np.float64)
    neg = lam < 0
    aud = dict(m=int(lam.size), lam_max=float(lam.max()), lam_min=float(lam.min()),
               n_negative=int(neg.sum()), share_negative=float(neg.mean()),
               clipped_energy=float(np.abs(lam[neg]).sum() / np.abs(lam).sum()) if lam.size else float("nan"),
               clipped_energy_sq=float((lam[neg] ** 2).sum() / (lam ** 2).sum()) if lam.size else float("nan"),
               trace_pos=float(np.clip(lam, 0, None).sum()), trace_all=float(lam.sum()))
    if log is not None:
        log(f"    psd_eigh m={aud['m']}: lam in [{aud['lam_min']:.4e}, {aud['lam_max']:.4e}], "
            f"{aud['n_negative']} ({100*aud['share_negative']:.1f} %) negative, clipped energy {aud['clipped_energy']:.4f}")
    return lam, V, aud


def s_plus(lam: np.ndarray, V: np.ndarray) -> np.ndarray:
    """S+ = V diag(max(lam, 0)) V^T (float32)."""
    w = np.clip(lam, 0.0, None).astype(EIGH_DTYPE)
    return (V * w) @ V.T


def wiener_factors(Sp: np.ndarray, N: np.ndarray, J: float, *, eps_rel: float = EPS_REL, log=None):
    """The exact factorisation of G = S+ (S+ + N/J + eps I)^-1.

    Returns (L, Q, k, audit) with N' = N/J + eps I = L L^T, K = L^-1 S+ L^-T = Q diag(k) Q^T, so
    G = L Q diag(k/(1+k)) (L Q)^-1 and the Wiener gains are g = k / (1 + k).
    """
    m = N.shape[0]
    eps = eps_rel * float(np.trace(N)) / m
    Np = (np.asarray(N, dtype=EIGH_DTYPE) / EIGH_DTYPE(J)).copy(order="F")
    Np[np.diag_indices(m)] += EIGH_DTYPE(eps)
    L = sla.cholesky(Np, lower=True, overwrite_a=True, check_finite=False)
    del Np
    # K = L^-1 S+ L^-T   (two triangular solves; solve_triangular works on the left)
    K = sla.solve_triangular(L, np.asarray(Sp, dtype=EIGH_DTYPE), lower=True, check_finite=False)
    K = sla.solve_triangular(L, K.T, lower=True, check_finite=False)
    K = 0.5 * (K + K.T)                                    # exactly symmetric (K is symmetric in exact arithmetic)
    k, Q = sla.eigh(np.asfortranarray(K), driver=EIGH_DRIVER, overwrite_a=True, check_finite=False)
    del K
    k = np.clip(np.asarray(k, dtype=np.float64), 0.0, None)
    g = k / (1.0 + k)
    aud = dict(eps=float(eps), eps_rel=float(eps_rel), J=float(J), m=int(m),
               gain_median=float(np.median(g)), gain_mean=float(g.mean()), gain_p10=float(np.percentile(g, 10)),
               gain_p90=float(np.percentile(g, 90)), gain_max=float(g.max()),
               share_gain_lt_010=float((g < 0.10).mean()), share_gain_lt_025=float((g < 0.25).mean()),
               k_max=float(k.max()))
    if log is not None:
        log(f"    wiener_factors: eps {eps:.4e}; Wiener gain median {aud['gain_median']:.4f} "
            f"(p10 {aud['gain_p10']:.4f}, p90 {aud['gain_p90']:.4f}, max {aud['gain_max']:.4f}); "
            f"share < 0.10 {aud['share_gain_lt_010']:.3f}, < 0.25 {aud['share_gain_lt_025']:.3f}")
    return L, Q, k, aud


def apply_e1(X: np.ndarray, L: np.ndarray, Q: np.ndarray, k: np.ndarray, phi: float, *, mirror: bool = False) -> np.ndarray:
    """x-hat = G_phi x with G_phi = L Q diag(max(k/(1+k), phi)) (L Q)^-1, rows of X in the columns' space.

    (L Q)^-1 x = Q^T (L^-1 x);  then the diagonal;  then L (Q y).  `mirror` (control Cm) uses 2 I - G_phi, i.e. the
    gains 2 - g_phi.
    """
    g = np.clip(k / (1.0 + k), phi, None)
    if mirror:
        g = 2.0 - g
    Xf = np.asarray(X, dtype=EIGH_DTYPE)
    Y = sla.solve_triangular(L, Xf.T, lower=True, check_finite=False)      # L^-1 x   (m, nrows)
    Y = Q.T @ Y
    Y *= g.astype(EIGH_DTYPE)[:, None]
    Y = Q @ Y
    Y = L @ Y
    return np.asarray(Y.T, dtype=np.float64)


def apply_e2(X: np.ndarray, L: np.ndarray, lam: np.ndarray, V: np.ndarray, r: int) -> np.ndarray:
    """x-hat = S_r (S_r + N')^-1 x with S_r = V_r diag(lam_r) V_r^T (top r POSITIVE eigenpairs), by Woodbury:
        (S_r + N')^-1 = Ninv - Bt (diag(1/lam_r) + V_r^T Bt)^-1 Bt^T,   Bt = N'^-1 V_r,
    then multiply by S_r.  N' enters only through L (its Cholesky factor)."""
    order = np.argsort(lam)[::-1][:r]
    lr = lam[order].astype(np.float64)
    assert np.all(lr > 0), "E2: fewer than r positive eigenvalues"
    Vr = np.asarray(V[:, order], dtype=EIGH_DTYPE)                           # (m, r)
    Xf = np.asarray(X, dtype=EIGH_DTYPE).T                                   # (m, nrows)
    Bt = sla.cho_solve((L, True), Vr, check_finite=False)                    # N'^-1 V_r   (m, r)
    Ninv_x = sla.cho_solve((L, True), Xf, check_finite=False)                # N'^-1 x     (m, nrows)
    Mid = np.diag(1.0 / lr).astype(EIGH_DTYPE) + Vr.T @ Bt                   # (r, r)
    Z = Ninv_x - Bt @ np.linalg.solve(Mid.astype(np.float64), (Bt.T @ Xf).astype(np.float64)).astype(EIGH_DTYPE)
    Y = Vr @ (lr.astype(EIGH_DTYPE)[:, None] * (Vr.T @ Z))                   # S_r Z
    return np.asarray(Y.T, dtype=np.float64)


def e3_weights(lam: np.ndarray, nu: np.ndarray, c: float, *, gamma: float = GAMMA,
               lmax_floor: float = E3_LMAX_FLOOR) -> np.ndarray:
    """W = V diag(w) V^T with w_i = (max(lam_i, 0) + floor_i)^(-gamma/2); floor_i = c * nu_i (the noise level in
    direction i) when c > 0, else lmax_floor * lam_max."""
    lp = np.clip(np.asarray(lam, dtype=np.float64), 0.0, None)
    floor = c * np.asarray(nu, dtype=np.float64) if c > 0 else np.full(lp.shape, lmax_floor * float(lp.max()))
    return (lp + floor) ** (-gamma / 2.0)


def apply_e3(X: np.ndarray, V: np.ndarray, w: np.ndarray) -> np.ndarray:
    """The un-preserved whitened rows: ((X V) * w) V^T."""
    Xf = np.asarray(X, dtype=EIGH_DTYPE)
    return np.asarray(((Xf @ V) * w.astype(EIGH_DTYPE)) @ V.T, dtype=np.float64)


def direction_noise(N: np.ndarray, V: np.ndarray) -> np.ndarray:
    """nu_i = v_i^T N v_i for every column of V (the noise variance in each of S-hat's directions)."""
    Z = np.asarray(N, dtype=EIGH_DTYPE) @ V
    return np.asarray((V * Z).sum(0), dtype=np.float64)


# ------------------------------------------------------------------------------------------------- controls
def random_psd(m: int, n_rand: int, target_trace: float, seed: int = RAND_SEED) -> np.ndarray:
    """Cr: an energy-matched random PSD covariance -- a Wishart with n_rand degrees of freedom, rescaled so its
    trace equals `target_trace` (the same total energy as S+, spread isotropically)."""
    rng = np.random.default_rng(seed)
    W = rng.standard_normal((n_rand, m)).astype(EIGH_DTYPE)
    R = (W.T @ W) / EIGH_DTYPE(n_rand)
    del W
    R *= EIGH_DTYPE(target_trace / float(np.trace(R)))
    return R


# ------------------------------------------------------------------------------------------------- row helpers
def norm_match(new: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """263's `:nm`: rescale every row of `new` to the corresponding row norm of `ref` (shape-only)."""
    out = np.array(new, dtype=np.float64, copy=True)
    na = np.linalg.norm(out, axis=1)
    nb = np.linalg.norm(ref, axis=1)
    ok = (na > 0) & (nb > 0)
    out[ok] *= (nb[ok] / na[ok])[:, None]
    return out


def row_cos(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


def principal_angles(V1: np.ndarray, V2: np.ndarray, k: int) -> np.ndarray:
    """The cosines of the principal angles between the top-k columns of two orthonormal bases (descending)."""
    A = np.asarray(V1[:, :k], dtype=np.float64)
    B = np.asarray(V2[:, :k], dtype=np.float64)
    s = np.linalg.svd(A.T @ B, compute_uv=False)
    return np.clip(s, -1.0, 1.0)
