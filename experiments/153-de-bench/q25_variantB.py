"""Q25 — variant B (resample + repair gene totals) vs variant A (resample, no repair).

The repair pins each gene's total to n*Lbar*profile_p by redistributing counts BETWEEN cells,
using the resampled values as allocation weights -- so the cell-to-cell pattern is inherited
from real cells while the pseudobulk is exact.  That is the property every prior emission arm
rested on (donor_transfer.py:310-320): repaired totals make the mean profile bit-identical, so
only the DE members can move.

Also measures what variant A actually risks: cos(emitted pseudobulk, predicted profile) and
relative L2, in the log1p-CPM space the scorer reads.
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
    s = w.sum()
    q = (np.full(w.size, T / w.size) if s <= 0 else w / s * T)
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

def spurious(Y):
    Dp = (sp.diags(1e6 / np.maximum(Y.sum(1), 1)) @ sp.csr_matrix(Y)).tocsc()
    fb = np.empty(gated.size)
    for k, g in enumerate(gated):
        nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
        e = np.asarray(Dp[:, g].todense()).ravel()
        vals, cts = np.unique(e, return_counts=True); tot = 0.0
        for v, c in zip(vals, cts):
            if v == 0.0: below, eq = 0, n_zero
            else:
                lo = np.searchsorted(nz, v, "left"); hi = np.searchsorted(nz, v, "right")
                below, eq = n_zero + lo, hi - lo
            tot += c * (below + 0.5 * eq) / m
        fb[k] = tot / N
    return fb, (np.abs((fb - 0.5) * fac) > 2.8)

print(f"{'scale':>6} {'variant':>26} {'C_hat':>7} {'vs real':>8} {'sig%':>7} {'pseudobulk cos':>15} {'rel L2':>8}")
for s in (0.0, 1.0, 2.0):
    f = np.exp(s * eff); p_eff = p0 * f; p_eff /= p_eff.sum()
    # ---- A: resample, fold change, per-cell library preserved, no repair ----
    Yf = Y0 * f[None, :]; Yf *= (libs_v / np.maximum(Yf.sum(1), 1e-9))[:, None]
    YA = np.floor(Yf)
    for i in range(N):
        r = int(libs_v[i] - YA[i].sum())
        if r > 0: YA[i, np.argsort(-(Yf[i] - YA[i]))[:r]] += 1
    # ---- B: same, then repair each gene's total to n*Lbar*p_eff ----
    T_tgt = np.floor(N * libs_v.mean() * p_eff)
    YB = np.zeros_like(YA)
    for g in np.flatnonzero(T_tgt > 0):
        YB[:, g] = lr(YA[:, g], int(T_tgt[g]))
    for lbl, Y in (("A resample (no repair)", YA), ("B resample + total repair", YB)):
        pb = np.log1p(TS * Y.sum(0) / max(Y.sum(), 1))
        tgt = np.log1p(TS * p_eff)
        cos = float(pb @ tgt / max(np.linalg.norm(pb) * np.linalg.norm(tgt), 1e-12))
        rl2 = float(np.linalg.norm(pb - tgt) / max(np.linalg.norm(tgt), 1e-12))
        C = jackknife_C(Y); _, sg = spurious(Y)
        print(f"{s:6.1f} {lbl:>26} {C:7.1f} {C/refC:8.3f} {sg.mean():7.1%} {cos:15.6f} {rl2:8.4f}")
