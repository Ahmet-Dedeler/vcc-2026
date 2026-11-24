"""Q23 — C_hat and Fbar for a TRUE resampling law, with the fold change actually applied.

The law under test, as the peer specified it:
  1. draw 400 real control cells
  2. apply the predicted per-gene fold change multiplicatively
  3. renormalise each cell back to ITS OWN original library size
  4. integerise, preserving each cell's library exactly (largest remainder)

At zero effect this is the identity on the drawn cells, so C_hat = the real-cell value and
Fbar = 0.5 by exchangeability.  The question is what the fold-change step costs, since step 3
redistributes mass across genes within each cell and step 4 quantises it.

Fold changes come from 104's frozen `baseline` -- our own shipped effect magnitudes, not a
synthetic -- swept over a scale so the answer is read as a function of effect size rather than
resting on one assumed log convention.  Null-gene Fbar is reported on the genes whose OWN effect
is ~0, which is where the artefact lived.
"""
import anndata as ad, numpy as np, scipy.sparse as sp, csv

TS, N, PHI, SPREAD = 5e4, 400, 4.0, 0.25
rng = np.random.default_rng(0)

def jackknife_C(Y):
    Y = np.asarray(Y, dtype=np.float64); n = Y.shape[0]
    P = Y.sum(0); S = Y.sum(); l = Y.sum(1)
    V = np.log1p(TS * (P[None, :] - Y) / np.maximum((S - l)[:, None], 1e-9))
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())

def integerise_rows(Yf, targets):
    """Round each row to integers preserving that row's total exactly."""
    out = np.floor(Yf)
    for i in range(Yf.shape[0]):
        r = int(targets[i] - out[i].sum())
        if r > 0:
            out[i, np.argsort(-(Yf[i] - out[i]))[:r]] += 1
        elif r < 0:
            cand = np.flatnonzero(out[i] > 0)
            out[i, cand[np.argsort(Yf[i, cand] - out[i, cand])[:(-r)]]] -= 1
    return out

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsr().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
picks = rng.choice(X.shape[0], N, replace=False)
Y0 = np.asarray(X[picks].todense()); libs_v = lib[picks]
ref_C = jackknife_C(Y0)
fac = np.sqrt(12.0 * N * m / (m + N + 1.0))

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
eff = z104["baseline"][0][[g104[g] for g in canon]].astype(np.float64)   # one target's effect
print(f"400 real control cells: C_hat = {ref_C:.1f}   (production quasi fixed-L = 9.4, 0.656x)")
print(f"effect vector: |eff| median {np.median(np.abs(eff)):.4f}  p99 {np.percentile(np.abs(eff),99):.4f}\n")
print(f"{'fold-change scale':>18} {'C_hat':>8} {'vs real':>8} {'null-gene Fbar med':>19} {'sig@null':>10}")
for s in (0.0, 0.5, 1.0, 2.0):
    f = np.exp(s * eff)
    Yf = Y0 * f[None, :]
    Yf *= (libs_v / np.maximum(Yf.sum(1), 1e-9))[:, None]      # preserve each cell's library
    Y = integerise_rows(Yf, libs_v)
    C = jackknife_C(Y)
    Dp = (sp.diags(1e6 / np.maximum(Y.sum(1), 1)) @ sp.csr_matrix(Y)).tocsc()
    null_g = gated[np.abs(eff[gated]) < 1e-4]                  # genes with ~no effect of their own
    fb = np.empty(null_g.size)
    for k, g in enumerate(null_g):
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
        fb[k] = tot / N
    zz = (fb - 0.5) * fac
    print(f"{s:18.1f} {C:8.1f} {C/ref_C:8.3f} {np.median(fb):19.4f} {np.mean(np.abs(zz)>2.8):10.1%}"
          + ("   <- zero effect" if s == 0 else "   <- our shipped magnitude" if s == 1 else ""))
