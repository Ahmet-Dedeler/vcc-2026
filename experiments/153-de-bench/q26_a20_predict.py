"""Q26 — pre-register what a20 should read under H_A vs H_B, before it scores.

a15 ships the a12 tiered spread policy: `--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0`,
so the top 8,000 genes by |effect| are emitted at spread 0.25 and the remaining gated genes at
spread 1.0.  a20 is `--keep-top 14000`, which exceeds the gated count, so EVERY gene goes to
spread 0.25.  Lower spread -> tighter emitted law -> more genes clear significance, so a20's
only effect is to raise n_pred.

That makes a20 a clean one-variable coverage probe, and the two hypotheses predict opposite
signs on `fid`:

    fid = precision * min(1, n_pred / n_real)

  H_B (n_real ~240):   coverage already 1 at a15's n_pred, so raising it only DILUTES
                       precision -> fid FALLS.
  H_A (n_real ~2,500): coverage is biting at a15's n_pred, so raising it RECOVERS coverage
                       -> fid RISES.

`jac` discriminates too and in the same direction under H_A.
"""
import anndata as ad, numpy as np, scipy.sparse as sp, csv

TS, N, PHI = 5e4, 400, 4.0
rng = np.random.default_rng(0)

def lr(w, T):
    if T <= 0: return np.zeros(w.size)
    q = w / w.sum() * T; fl = np.floor(q); r = int(round(T - fl.sum()))
    if r > 0: fl[np.argsort(-(q - fl))[:r]] += 1
    return fl

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsr().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
L = int(np.median(lib)); D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
colsum = np.asarray(X.sum(0)).ravel(); p0 = colsum / colsum.sum()
fac = np.sqrt(12.0 * N * m / (m + N + 1.0))
canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
eff = z104["baseline"][0][[g104[g] for g in canon]].astype(np.float64)
f = np.exp(eff); p_eff = p0 * f; p_eff /= p_eff.sum()

# a15's tiering: rank gated genes by |effect|, top 8000 -> spread .25, rest -> spread 1.0
rank = np.argsort(-np.abs(eff[gated]))
spread_a15 = np.full(gated.size, 1.0); spread_a15[rank[:8000]] = 0.25
spread_a20 = np.full(gated.size, 0.25)
Ls = np.full(N, float(L))

def n_pred_frac(spread_vec):
    sig = 0
    for k, g in enumerate(gated):
        T = int(np.floor(N * L * p_eff[g]))
        if T <= 0: continue
        s = spread_vec[k]; t_r = int(np.floor(T * s)); t_e = T - t_r
        c = lr(Ls, t_e)
        if t_r > 0:
            c = c + rng.poisson(rng.gamma(1.0 / PHI, PHI * np.full(N, t_r / N)))
        nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
        cpm = c * 1e6 / L
        vals, cts = np.unique(cpm, return_counts=True); tot = 0.0
        for v, ct in zip(vals, cts):
            if v == 0.0: below, eq = 0, n_zero
            else:
                lo = np.searchsorted(nz, v, "left"); hi = np.searchsorted(nz, v, "right")
                below, eq = n_zero + lo, hi - lo
            tot += ct * (below + 0.5 * eq) / m
        if abs((tot / N - 0.5) * fac) > 2.8: sig += 1
    return sig / gated.size

fa15, fa20 = n_pred_frac(spread_a15), n_pred_frac(spread_a20)
Ngat = gated.size
n15, n20 = fa15 * Ngat, fa20 * Ngat
scale = 4400.0 / max(n15, 1)          # calibrate the model to a15's observed board n_pred
print(f"gated {Ngat} | modelled n_pred: a15 {n15:.0f} ({fa15:.1%})  a20 {n20:.0f} ({fa20:.1%})")
print(f"ratio a20/a15 = {n20/max(n15,1):.3f}  ->  board a15 n_pred ~4,400 implies a20 ~{n20*scale:.0f}\n")
print("PRE-REGISTERED DISCRIMINATOR (fid = precision * min(1, n_pred/n_real), precision held flat):")
print(f"{'world':>22} {'n_real':>8} {'cov a15':>9} {'cov a20':>9} {'predicted fid ratio':>21} {'sign':>6}")
for lbl, nr in (("H_B  small n_real", 240.0), ("H_A  large n_real", 2500.0)):
    c15, c20 = min(1.0, 4400 / nr), min(1.0, n20 * scale / nr)
    print(f"{lbl:>22} {nr:8.0f} {c15:9.3f} {c20:9.3f} {c20/c15:21.3f} {'FALLS' if c20 <= c15 else 'RISES':>6}")
print("\nUnder H_B both coverages are saturated, so any fid change is precision dilution -> fid FALLS.")
print("Under H_A a15 is coverage-limited, so raising n_pred recovers coverage -> fid RISES.")
print("The SIGN of a20's fid delta is the readout. jac should move the same way under H_A.")
