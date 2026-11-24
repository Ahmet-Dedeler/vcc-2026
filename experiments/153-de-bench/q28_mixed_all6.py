"""Q28 — the x=0.256 mixed law, priced on the members a19 was blind to.

a19 lost -0.0959 on `mse` alone, a member my (gamma, lam) screen never measured on any operator
at any cell.  So before the mixed arm: what does it do to `mse`?

`mse` numerator N_p = (1/G_p)( ||b_hat - b||^2 - rho*min(C_hat, C_p) - C_p ).  An EMISSION change
reaches it through exactly two terms:
  * ||b_hat - b||^2 -- only if the emitted pseudobulk drifts off the predicted profile.  Variant
    B's total repair pins it (rel L2 0.0013), so this term should be inert.
  * C_hat -- the credit.  The cap fires on 0/900 rows (COMPETITIVE-INTEL:162), so raising C_hat
    toward C_p is credited rather than wasted.

That is the structural reason this arm is not a19: a19 changed the OPERATOR, which changes the
effect vector's shape and therefore ||b_hat - b||^2 directly, with alpha still solved at the old
gamma.  The mixed law changes only which law emits each gene.
"""
import anndata as ad, numpy as np, scipy.sparse as sp, csv

TS, N, PHI, SPREAD = 5e4, 400, 4.0, 0.25
rng = np.random.default_rng(0)

def jackknife_C(Y):
    Y = np.asarray(Y, float); n = Y.shape[0]
    P = Y.sum(0); S = Y.sum(); l = Y.sum(1)
    V = np.log1p(TS * (P[None, :] - Y) / np.maximum((S - l)[:, None], 1e-9))
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())

def lr(w, T):
    if T <= 0: return np.zeros(w.size)
    s = w.sum(); q = (np.full(w.size, T / w.size) if s <= 0 else w / s * T)
    fl = np.floor(q); r = int(round(T - fl.sum()))
    if r > 0: fl[np.argsort(-(q - fl))[:r]] += 1
    return fl

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsr().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
L = int(np.median(lib)); G = X.shape[1]
D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
colsum = np.asarray(X.sum(0)).ravel(); p0 = colsum / colsum.sum()
picks = rng.choice(X.shape[0], N, replace=False)
Y0 = np.asarray(X[picks].todense()); libs_v = lib[picks]
refC = jackknife_C(Y0); fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
eff = z104["baseline"][0][[g104[g] for g in canon]].astype(np.float64)
f = np.exp(eff); p_eff = p0 * f; p_eff /= p_eff.sum()
Ls = np.full(N, float(L))

# resampled block (variant B: fold change, per-cell library preserved, then total repair)
Yf = Y0 * f[None, :]; Yf *= (libs_v / np.maximum(Yf.sum(1), 1e-9))[:, None]
YR = np.floor(Yf)
for i in range(N):
    r = int(libs_v[i] - YR[i].sum())
    if r > 0: YR[i, np.argsort(-(Yf[i] - YR[i]))[:r]] += 1
T_tgt = np.floor(N * libs_v.mean() * p_eff)
for g in np.flatnonzero(T_tgt > 0):
    YR[:, g] = lr(YR[:, g], int(T_tgt[g]))
# parametric block (a15 keep level)
YP = np.zeros((N, G))
T_all = np.floor(N * L * p_eff).astype(np.int64)
for g in np.flatnonzero(T_all > 0):
    T = int(T_all[g]); t_r = int(np.floor(T * SPREAD)); t_e = T - t_r
    c = lr(Ls, t_e)
    if t_r > 0: c = c + rng.poisson(rng.gamma(1.0/PHI, PHI*np.full(N, t_r/N)))
    YP[:, g] = c

order = gated[np.argsort(-np.abs(eff[gated]))]
tgt = np.log1p(TS * p_eff)
print(f"{'x (frac parametric)':>20} {'genes param':>12} {'C_hat':>7} {'vs real':>8} {'pseudobulk relL2':>17} {'sig%':>7}")
for x in (0.0, 0.164, 0.256, 0.411, 1.0):
    kp = set(order[:int(round(x * gated.size))].tolist())
    Y = YR.copy()
    for g in kp: Y[:, g] = YP[:, g]
    # POST-MIX TOTAL REPAIR.  The two blocks carry different per-cell library structures
    # (parametric at fixed L, resampled at real libraries), so splicing genes between them
    # breaks the pseudobulk pinning -- measured at rel L2 0.0155 without this, against 0.0015
    # for pure resampling.  mse reads ||b_hat - b||^2 directly, so that drift is exactly the
    # channel that cost a19 its mse member.  Repair AFTER the mix, not before.
    for g in np.flatnonzero(T_tgt > 0):
        Y[:, g] = lr(Y[:, g], int(T_tgt[g]))
    C = jackknife_C(Y)
    pb = np.log1p(TS * Y.sum(0) / max(Y.sum(), 1))
    rl2 = float(np.linalg.norm(pb - tgt) / max(np.linalg.norm(tgt), 1e-12))
    Dp = (sp.diags(1e6/np.maximum(Y.sum(1),1)) @ sp.csr_matrix(Y)).tocsc(); sig = 0
    for g in gated:
        nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
        e = np.asarray(Dp[:, g].todense()).ravel()
        vals, cts = np.unique(e, return_counts=True); tot = 0.0
        for v, ct in zip(vals, cts):
            if v == 0.0: below, eq = 0, n_zero
            else:
                lo = np.searchsorted(nz, v, "left"); hi = np.searchsorted(nz, v, "right")
                below, eq = n_zero + lo, hi - lo
            tot += ct * (below + 0.5*eq) / m
        sig += abs((tot/N - 0.5)*fac) > 2.8
    print(f"{x:20.3f} {len(kp):12d} {C:7.1f} {C/refC:8.3f} {rl2:17.4f} {sig/gated.size:7.1%}")
print(f"\na15 reference: C_hat 9.4 (0.656x).  Higher C_hat = MORE mse credit (cap fires 0/900).")
