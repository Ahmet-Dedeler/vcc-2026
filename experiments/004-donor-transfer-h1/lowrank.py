"""Does low-rank denoising of the donor effect matrix improve transfer?

K562-gwps has a median of 178 cells per perturbation, so each effect vector carries
substantial sampling noise. Perturbation responses are known to be low-dimensional --
a few hundred transcriptional programs -- so truncating the SVD should strip noise while
keeping signal. If the true shared structure is low-rank, this raises transfer cosine
for free, with no new data.

    uv run python experiments/004-donor-transfer-h1/lowrank.py
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad, scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TransferModel, TS

CONTROL, PERT_COL = "non-targeting", "target_gene"

a = ad.read_h5ad("data/harness/h1-test.h5ad")
X = a.X.tocsr(); labels = a.obs[PERT_COL].astype(str).to_numpy(); genes = np.asarray(a.var_names)
Pc = np.asarray(X[labels == CONTROL].sum(0)).ravel().astype(np.float64)
b_ctrl = np.log1p(TS * Pc / Pc.sum())
perts, deltas = [], []
for g in sorted(set(labels) - {CONTROL}):
    P = np.asarray(X[labels == g].sum(0)).ravel().astype(np.float64)
    deltas.append(np.log1p(TS * P / P.sum()) - b_ctrl); perts.append(g)
perts = np.asarray(perts); true = np.vstack(deltas)

donor = Donor.load("data/donors/k562_gwps")
model = TransferModel(donor, genes)
cov = model.covered(perts)
perts_c, true_c = perts[cov], true[cov]
true_cen = true_c - true_c.mean(0, keepdims=True)

m, s = donor.decompose()            # s: (9787, 8248) centered donor effects
print(f"donor s: {s.shape}   target: {true_c.shape}   covered {cov.sum()}/{len(perts)}")

print("\ncomputing SVD of the centered donor effect matrix ...", flush=True)
U, S, Vt = np.linalg.svd(s, full_matrices=False)
var = S**2 / (S**2).sum()
print(f"  top-10 variance explained: {var[:10].sum():.3f}   top-50: {var[:50].sum():.3f}"
      f"   top-200: {var[:200].sum():.3f}")

idx = {p: i for i, p in enumerate(donor.perts)}
rows = [idx[p] for p in perts_c]

def transfer_cos(s_mat):
    pred = np.zeros((len(rows), len(genes)), dtype=np.float64)
    pred[:, model.tcols] = s_mat[np.ix_(rows, model.dcols)]
    pc = pred - pred.mean(0, keepdims=True)
    na = np.linalg.norm(pc, axis=1); nb = np.linalg.norm(true_cen, axis=1)
    ok = (na > 0) & (nb > 0)
    return float(np.mean(np.sum(pc[ok]*true_cen[ok],axis=1)/(na[ok]*nb[ok])))

base = transfer_cos(s)
print(f"\n  rank    var explained   centered transfer cos   vs full")
print(f"  {'full':>6}  {1.0:13.3f}   {base:21.4f}   {'--':>7}")
best = (None, base)
res = [{"rank": "full", "cos": base}]
for k in (5, 10, 20, 30, 50, 75, 100, 150, 200, 400, 800):
    if k >= min(s.shape): break
    sk = (U[:, :k] * S[:k]) @ Vt[:k]
    c = transfer_cos(sk)
    res.append({"rank": k, "cos": c, "var": float(var[:k].sum())})
    flag = "  <-- best" if c > best[1] else ""
    if c > best[1]: best = (k, c)
    print(f"  {k:6d}  {var[:k].sum():13.3f}   {c:21.4f}   {c-base:+7.4f}{flag}")
print(f"\nbest rank {best[0]} -> cos {best[1]:.4f}  (full-rank {base:.4f}, "
      f"{100*(best[1]-base)/abs(base):+.1f}%)")
Path("experiments/004-donor-transfer-h1/lowrank_results.json").write_text(
    json.dumps({"results": res, "best_rank": best[0], "best_cos": best[1], "full_cos": base}, indent=2))
np.save("data/donors/k562_gwps/svd_S.npy", S)
np.save("data/donors/k562_gwps/svd_Vt.npy", Vt[:800].astype(np.float32))
np.save("data/donors/k562_gwps/svd_U.npy", U[:, :800].astype(np.float32))
print("saved SVD factors (top 800) to data/donors/k562_gwps/")
