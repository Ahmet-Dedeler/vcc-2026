"""Q22 — does varying the emitted library sizes move the `mse` jackknife credit?

This is the blocking question for the lib_sizes fix.  BRIEF section 2:

    C_p = (n-1)/n * sum_g sum_i (v_ig - vbar_g)^2,
    v_ig = log(1 + TS * (P_pg - y_ig) / (S_p - l_i)),   TS = 5e4

and the credited correction is `rho * min(C_hat_p, C_p)` — capped by the REFERENCE's own
correction.  a15's entire `mse` member (+0.120 scaled, 6th best on the board) comes from having
raised emitted dispersion to that cap.  So the question is not "does C_hat move" but "does it
move toward the cap or past it".

The right comparator is **400 real control cells**: that is exactly what C_hat would be if we
emitted real cells, so it stands in for the cap's scale.  Under-shooting it forfeits credit;
over-shooting wastes the excess (min()) while `||b_hat - b||^2` still grows with the emitted
variance — which section 4.1b of AGENTS.md measured as strictly harmful at rho < 1.
"""
import anndata as ad, numpy as np, scipy.sparse as sp

TS, N, PHI, SPREAD = 5e4, 400, 4.0, 0.25
rng = np.random.default_rng(0)

def jackknife_C(Y):
    """BRIEF section 2's delete-one jackknife, exactly."""
    Y = np.asarray(Y, dtype=np.float64)
    n = Y.shape[0]
    P = Y.sum(0); S = Y.sum(); l = Y.sum(1)
    denom = (S - l)[:, None]
    V = np.log1p(TS * (P[None, :] - Y) / np.maximum(denom, 1e-9))
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())

def largest_remainder(w, T):
    if T <= 0: return np.zeros(w.size, dtype=np.int64)
    q = w / w.sum() * T; fl = np.floor(q).astype(np.int64)
    r = int(T - fl.sum())
    if r > 0: fl[np.argsort(-(q - fl))[:r]] += 1
    return fl

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsr().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
L = int(np.median(lib))
colsum = np.asarray(X.sum(0)).ravel(); p = colsum / colsum.sum()
G = X.shape[1]
picks = rng.choice(X.shape[0], N, replace=False)          # what the builder does
libs_v = lib[picks]

def emit(varlib, disp):
    Ls = libs_v if varlib else np.full(N, float(L))
    Y = np.zeros((N, G))
    T_all = np.floor(Ls.sum() * p).astype(np.int64)
    for g in np.flatnonzero(T_all > 0):
        T = int(T_all[g])
        if disp:
            t_r = int(np.floor(T * SPREAD)); t_e = T - t_r
            c = largest_remainder(Ls, t_e).astype(np.float64)
            if t_r > 0:
                mu_r = Ls / Ls.sum() * t_r
                c = c + rng.poisson(rng.gamma(1.0 / PHI, PHI * np.maximum(mu_r, 1e-12)))
        else:
            c = largest_remainder(Ls, T).astype(np.float64)
        Y[:, g] = c
    return Y

real400 = np.asarray(X[picks].todense())
ref = jackknife_C(real400)
print(f"context A | emitted cells n={N}, genes={G}, TS={TS:.0f}")
print(f"{'configuration':44} {'C_hat':>12} {'vs 400 real cells':>18}")
print(f"{'400 REAL control cells (the cap analogue)':44} {ref:12.1f} {1.0:18.3f}")
for label, (vl, dp) in (("min-variance, fixed L", (False, False)),
                        ("production quasi .25 / phi 4, fixed L", (False, True)),
                        ("quasi .25 / phi 4, VARYING libraries", (True, True)),
                        ("min-variance, VARYING libraries", (True, False))):
    c = jackknife_C(emit(vl, dp))
    print(f"{label:44} {c:12.1f} {c/ref:18.3f}")
