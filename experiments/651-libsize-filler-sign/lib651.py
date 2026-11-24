"""651 step 1 (descriptive): per-cell library sizes in h1-test, controls vs each target's 400 cells.
Writes work/lib651.npz (per-cell L, target labels, batch suffix) and prints the per-target shift distribution.
Analysis only, < 3 GB, ~1 min. Usage: uv run python experiments/651-libsize-filler-sign/lib651.py"""
import numpy as np, anndata as ad, pandas as pd, scipy.sparse as sp
from pathlib import Path
W = Path(__file__).resolve().parent / "work"
a = ad.read_h5ad("data/harness/h1-test.h5ad", backed="r")
t = a.obs["target_gene"].astype(str).values
batch = pd.Series(a.obs_names).str.split("-", n=1).str[1].values
L = np.zeros(a.n_obs)
for s in range(0, a.n_obs, 4000):
    X = a.X[s:s + 4000]
    L[s:s + 4000] = np.asarray(X.sum(1)).ravel()
np.savez(W / "lib651.npz", L=L, t=t, batch=batch)
ctrl = t == "non-targeting"
lc = np.log(L[ctrl])
print("controls n %d  median L %.0f  mean logL %.4f  sd logL %.4f" % (ctrl.sum(), np.median(L[ctrl]), lc.mean(), lc.std()))
tg = sorted(set(t[~ctrl]))
rows = []
for g in tg:
    m = t == g
    rows.append((g, m.sum(), np.median(L[m]), np.log(L[m]).mean() - lc.mean()))
df = pd.DataFrame(rows, columns=["target", "n", "medL", "dlogL"]).sort_values("dlogL")
print(df.dlogL.describe())
print("se of a 400-cell mean logL: %.4f" % (lc.std() / np.sqrt(400)))
print("share of targets with |dlogL| > 3 se: %.2f" % (np.abs(df.dlogL) > 3 * lc.std() / np.sqrt(400)).mean())
print(df.head(8).to_string(index=False)); print(df.tail(8).to_string(index=False))
# batch balance: controls' batch mix vs perturbed cells' batch mix
cb = pd.Series(batch[ctrl]).value_counts(normalize=True)
pb = pd.Series(batch[~ctrl]).value_counts(normalize=True)
j = pd.concat([cb, pb], axis=1, keys=["ctrl", "pert"]).fillna(0)
print("batch mix TV distance ctrl vs pert: %.4f" % (0.5 * np.abs(j.ctrl - j.pert).sum()))
df.to_csv(W / "dlogL_h1.csv", index=False)
