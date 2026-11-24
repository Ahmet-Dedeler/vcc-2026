"""What fraction of the mse distance is the perturbed gene itself?

The mse baseline is the context's MEAN perturbation response, applied identically to every
perturbation. Averaging over hundreds of different targets washes the on-target knockdown out
almost completely -- so the baseline captures essentially none of it. If the on-target gene is
a material share of ||b_p - b_ctrl||^2, then predicting JUST that one gene correctly, and the
mean response everywhere else, beats the baseline by that share.

Measured on real perturbation data, in the metric's own space: profiles are log1p(TS * count
fraction), the same transform cell_eval2 uses.
"""
import numpy as np, anndata as ad, sys
TS = 10_000.0
REAL = sys.argv[1] if len(sys.argv) > 1 else "data/harness/orion-hct116-thinned.h5ad"
CTRL, COL = "non-targeting", "target_gene"
R = ad.read_h5ad(REAL)
X = R.X; X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)
lab = np.asarray(R.obs[COL].astype(str).values)
genes = np.asarray([str(g) for g in R.var_names]); gpos = {g: i for i, g in enumerate(genes)}
def prof(M):
    s = M.sum(1, keepdims=True); s[s == 0] = 1.0
    return np.log1p(TS * (M / s)).mean(0)
b_ctrl = prof(X[lab == CTRL])
targets = [t for t in sorted(set(lab)) if t != CTRL]
rows, on_share, hit = [], [], 0
D = []
for t in targets:
    m = lab == t
    if m.sum() < 2: continue
    d = prof(X[m]) - b_ctrl
    D.append(d); rows.append(t)
D = np.array(D)
print(f"{REAL.split('/')[-1]}: {len(rows)} perturbations, {genes.size} genes")
# the metric's baseline: the mean response, applied to everyone
mean_resp = D.mean(0)
sq_total = (D**2).sum(1)
sq_base  = ((D - mean_resp)**2).sum(1)      # residual left by the baseline
print(f"\n||d_p||^2 median                       {np.median(sq_total):.4f}")
print(f"||d_p - mean||^2 median (the baseline) {np.median(sq_base):.4f}")
print(f"  -> baseline explains {100*(1-np.median(sq_base)/np.median(sq_total)):.1f}% of the distance")
# on-target share
idx = [(i, gpos[t]) for i, t in enumerate(rows) if t in gpos]
print(f"\n{len(idx)}/{len(rows)} targets are on the gene axis")
ot = np.array([D[i, j] for i, j in idx])
ot_base = np.array([(D[i, j] - mean_resp[j]) for i, j in idx])
sq_t = np.array([sq_total[i] for i, _ in idx]); sq_b = np.array([sq_base[i] for i, _ in idx])
print(f"  on-target effect: mean {ot.mean():+.4f}  median {np.median(ot):+.4f}  "
      f"{100*np.mean(ot<0):.1f}% negative")
print(f"  the baseline's on-target prediction: mean {np.median(mean_resp[[j for _,j in idx]]):+.5f} "
      f"(washed out, as expected)")
print(f"\n  SHARE OF THE BASELINE RESIDUAL THAT IS THE ON-TARGET GENE ALONE:")
print(f"    median (d_on - mean_on)^2 / ||d - mean||^2 = "
      f"{100*np.median(ot_base**2/sq_b):.2f}%")
print(f"    mean                                       = {100*np.mean(ot_base**2/sq_b):.2f}%")
print(f"    sum over perturbations                     = {100*(ot_base**2).sum()/sq_b.sum():.2f}%")
print(f"\n  => an oracle that predicts the mean response PLUS the exact on-target knockdown")
print(f"     would score mse_raw ~ {1 - (ot_base**2).sum()/sq_b.sum():.4f} x the baseline's,")
print(f"     i.e. about {0.993*(1-(ot_base**2).sum()/sq_b.sum()):.4f} on the board's scale.")
# how concentrated is the rest?
print(f"\n  concentration of the baseline residual (excluding on-target):")
Dr = D - mean_resp
for k in (1, 10, 50, 200, 1000):
    s = np.sort(Dr**2, axis=1)[:, ::-1][:, :k].sum(1) / (Dr**2).sum(1)
    print(f"    top-{k:<5} genes carry {100*np.median(s):.1f}% of it")
