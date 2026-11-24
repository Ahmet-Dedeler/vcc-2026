"""A/B/C jackknife margin: our C_hat (395's emitted blocks) vs the jackknife of 400 random real control cells (proxy for C)."""
import numpy as np, anndata as ad, h5py
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
TS = 5e4
def jk(X):
    X = np.asarray(X, dtype=np.float64); P = X.sum(0); S = P.sum(); L = X.sum(1); live = P > 0
    V = np.log1p(TS * (P[live][None, :] - X[:, live]) / (S - L)[:, None]); n = X.shape[0]
    return float((n - 1) / n * ((V - V.mean(0)) ** 2).sum())
H = ROOT / "experiments/395-protected-signs-arm/work/build/the-thermostat-has-a-union-rep.h5ad"
rng = np.random.default_rng(0)
with h5py.File(H, "r") as f:
    indptr = f["X/indptr"]; data = f["X/data"]; idx = f["X/indices"]
    ng = f["var/_index"].shape[0] if "_index" in f["var"] else 18533
    for ci, ctx in enumerate("ABC"):
        A = ad.read_h5ad(ROOT / f"data/challenge/context_{ctx}.h5ad"); Xc = A.X.tocsr()
        cref = [jk(Xc[rng.choice(Xc.shape[0], 400, replace=False)].toarray()) for _ in range(8)]
        ours = []
        for t in rng.choice(300, 12, replace=False):
            r0 = (ci * 300 + int(t)) * 400; r1 = r0 + 400
            a, b = int(indptr[r0]), int(indptr[r1])
            ip = np.asarray(indptr[r0:r1 + 1]) - a
            M = np.zeros((400, ng)); d = np.asarray(data[a:b]); j = np.asarray(idx[a:b])
            for k in range(400):
                M[k, j[ip[k]:ip[k + 1]]] = d[ip[k]:ip[k + 1]]
            ours.append(jk(M))
        print(f"{ctx}: C_ctrl400 mean {np.mean(cref):.4f} (sd {np.std(cref):.4f})  our C_hat q0/50/100 {np.quantile(ours,[0,.5,1]).round(4)}  ratio median {np.median(ours)/np.mean(cref):.3f} min {np.min(ours)/np.mean(cref):.3f}")
        del A, Xc
