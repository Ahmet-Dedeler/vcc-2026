"""Q24 — Fbar stratified by effect decile, both emitters carrying the SAME predicted effect.

The peer's caveat: measuring the null on genes with own-effect < 1e-4 selects coordinates the
donor doesn't cover -- the easy half.  Most genes live at |eff| ~ 0.0152, two orders of magnitude
above that threshold, and the question is whether the null degrades there.

For a gene WITH a real effect, Fbar SHOULD move off 0.5 -- that is signal, not artefact.  So the
diagnostic that separates them is **directional concordance**: of the genes an emitter calls
significant, what fraction move in the direction of our own predicted effect?  Production's calls
are ~88% UP regardless of the effect's sign, which is the artefact by definition.  A clean emitter's
calls should follow the model.

Both emitters carry the identical effect-bearing profile, so this is a like-for-like comparison.
"""
import anndata as ad, numpy as np, scipy.sparse as sp, csv

TS, N, PHI, SPREAD = 5e4, 400, 4.0, 0.25
rng = np.random.default_rng(0)

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
D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
colsum = np.asarray(X.sum(0)).ravel(); p0 = colsum / colsum.sum()
picks = rng.choice(X.shape[0], N, replace=False)
Y0 = np.asarray(X[picks].todense()); libs_v = lib[picks]
fac = np.sqrt(12.0 * N * m / (m + N + 1.0))

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
eff = z104["baseline"][0][[g104[g] for g in canon]].astype(np.float64)
f = np.exp(eff)
p_eff = p0 * f; p_eff /= p_eff.sum()                      # the effect-bearing profile

# --- emitter 1: production quasi, fixed L, carrying the effect ---------------------------
Yq = np.zeros((N, X.shape[1]))
T_all = np.floor(N * L * p_eff).astype(np.int64)
Ls_fixed = np.full(N, float(L))
for g in np.flatnonzero(T_all > 0):
    T = int(T_all[g]); t_r = int(np.floor(T * SPREAD)); t_e = T - t_r
    c = largest_remainder(Ls_fixed, t_e).astype(np.float64)
    if t_r > 0:
        mu_r = np.full(N, t_r / N)
        c = c + rng.poisson(rng.gamma(1.0 / PHI, PHI * np.maximum(mu_r, 1e-12)))
    Yq[:, g] = c
# --- emitter 2: resampling, carrying the same effect -------------------------------------
Yf = Y0 * f[None, :]
Yf *= (libs_v / np.maximum(Yf.sum(1), 1e-9))[:, None]
Yr = np.floor(Yf)
for i in range(N):
    r = int(libs_v[i] - Yr[i].sum())
    if r > 0: Yr[i, np.argsort(-(Yf[i] - Yr[i]))[:r]] += 1

def fbar_all(Y, libs):
    Dp = (sp.diags(1e6 / np.maximum(Y.sum(1), 1)) @ sp.csr_matrix(Y)).tocsc()
    out = np.empty(gated.size)
    for k, g in enumerate(gated):
        nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
        e = np.asarray(Dp[:, g].todense()).ravel()
        vals, cts = np.unique(e, return_counts=True)
        tot = 0.0
        for v, c in zip(vals, cts):
            if v == 0.0: below, eq = 0, n_zero
            else:
                lo = np.searchsorted(nz, v, "left"); hi = np.searchsorted(nz, v, "right")
                below, eq = n_zero + lo, hi - lo
            tot += c * (below + 0.5 * eq) / m
        out[k] = tot / N
    return out

fb_q, fb_r = fbar_all(Yq, Ls_fixed), fbar_all(Yr, libs_v)
zq, zr = (fb_q - 0.5) * fac, (fb_r - 0.5) * fac
e_g = eff[gated]; ae = np.abs(e_g)
edges = np.percentile(ae, np.arange(0, 101, 10))
print(f"{'|eff| decile':>13} {'median |eff|':>12} | {'PRODUCTION quasi fixed-L':^32} | {'RESAMPLING':^32}")
print(f"{'':13} {'':12} | {'Fbar med':>9} {'sig%':>7} {'concord':>7} {'up%':>6} | {'Fbar med':>9} {'sig%':>7} {'concord':>7} {'up%':>6}")
for d in range(10):
    lo, hi = edges[d], edges[d+1]
    msk = (ae >= lo) & (ae <= hi) if d == 9 else (ae >= lo) & (ae < hi)
    if msk.sum() == 0: continue
    row = f"{d+1:>13} {np.median(ae[msk]):12.4f} |"
    for fb, z in ((fb_q, zq), (fb_r, zr)):
        s = np.abs(z[msk]) > 2.8
        conc = (np.sign(z[msk][s]) == np.sign(e_g[msk][s])).mean() if s.any() else float("nan")
        up = (z[msk][s] > 0).mean() if s.any() else float("nan")
        row += f" {np.median(fb[msk]):9.4f} {s.mean():7.1%} {conc:7.1%} {up:6.1%} |"
    print(row)
for lbl, z in (("PRODUCTION", zq), ("RESAMPLING", zr)):
    s = np.abs(z) > 2.8
    print(f"{lbl}: overall sig {s.mean():.1%}, concordance with our own predicted sign "
          f"{(np.sign(z[s])==np.sign(e_g[s])).mean():.1%}")
