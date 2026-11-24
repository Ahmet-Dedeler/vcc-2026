"""On-target share of the DE-NOISED signal.

`ontarget_share.py` measured shares of ||d_p - mean||^2 computed from finite samples, which
is inflated by sampling noise on every gene. The mse metric explicitly subtracts a jackknife
estimate of exactly that inflation, so the share that matters is of the de-noised quantity.

Unbiased estimator: split each perturbation's cells into two disjoint halves, form d1 and d2
independently, and use the CROSS product d1*d2 -- its expectation is the true squared signal
because the two halves' noise is independent. Same for the control.
"""
import numpy as np, anndata as ad, sys
TS = 10_000.0
REAL = sys.argv[1] if len(sys.argv) > 1 else "data/harness/orion-hct116-thinned.h5ad"
CTRL, COL = "non-targeting", "target_gene"
rng = np.random.default_rng(0)
R = ad.read_h5ad(REAL)
X = R.X; X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)
lab = np.asarray(R.obs[COL].astype(str).values)
genes = np.asarray([str(g) for g in R.var_names]); gpos = {g: i for i, g in enumerate(genes)}
def prof(M):
    s = M.sum(1, keepdims=True); s[s == 0] = 1.0
    return np.log1p(TS * (M / s)).mean(0)
ci = np.flatnonzero(lab == CTRL); rng.shuffle(ci)
cA, cB = prof(X[ci[:len(ci)//2]]), prof(X[ci[len(ci)//2:]])
targets = [t for t in sorted(set(lab)) if t != CTRL]
D1, D2, rows = [], [], []
for t in targets:
    idx = np.flatnonzero(lab == t)
    if idx.size < 4: continue
    rng.shuffle(idx); h = idx.size//2
    D1.append(prof(X[idx[:h]]) - cA); D2.append(prof(X[idx[h:]]) - cB); rows.append(t)
D1, D2 = np.array(D1), np.array(D2)
print(f"{REAL.split('/')[-1]}: {len(rows)} perturbations, split-half de-noised\n")
m1, m2 = D1.mean(0), D2.mean(0)
# de-noised squared norms via the cross product
tot_naive = ((D1+D2)/2)**2
sig_tot   = D1*D2                                   # unbiased for the true squared signal
res_naive = (((D1-m1)+(D2-m2))/2)**2
sig_res   = (D1-m1)*(D2-m2)
print(f"naive  sum ||d||^2            {tot_naive.sum():12.1f}")
print(f"DE-NOISED sum of true signal  {sig_tot.sum():12.1f}   "
      f"-> {100*sig_tot.sum()/tot_naive.sum():.1f}% of the naive number is real")
print(f"naive  sum ||d - mean||^2     {res_naive.sum():12.1f}")
print(f"DE-NOISED baseline residual   {sig_res.sum():12.1f}   "
      f"-> {100*sig_res.sum()/res_naive.sum():.1f}% real")
print(f"\nbaseline explains {100*(1-sig_res.sum()/sig_tot.sum()):.1f}% of the TRUE signal "
      f"(naive said {100*(1-res_naive.sum()/tot_naive.sum()):.1f}%)")
oi = [(i, gpos[t]) for i, t in enumerate(rows) if t in gpos]
on_naive = sum(res_naive[i, j] for i, j in oi)
on_sig   = sum(sig_res[i, j]   for i, j in oi)
print(f"\nON-TARGET share of the baseline residual:")
print(f"  naive      {100*on_naive/res_naive.sum():.2f}%")
print(f"  DE-NOISED  {100*on_sig/sig_res.sum():.2f}%   <-- the share that matters")
print(f"\n  oracle (mean response + exact on-target) would score mse_raw ~ "
      f"{0.993*(1-on_sig/sig_res.sum()):.4f}")
# concentration of the de-noised residual
per_gene = sig_res.sum(0); order = np.argsort(per_gene)[::-1]
cum = np.cumsum(per_gene[order])/per_gene.sum()
print(f"\nconcentration of the DE-NOISED baseline residual, by gene:")
for k in (1,10,50,200,1000,5000):
    if k <= cum.size: print(f"  top-{k:<5} genes carry {100*cum[k-1]:.1f}%")
