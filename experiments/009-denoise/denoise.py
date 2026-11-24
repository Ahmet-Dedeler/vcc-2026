"""Experiment 009 — denoise the donor delta without reweighting genes.

Nadig 2025 measured that within K562, the raw log-fold-change correlation between two
replicate experiments is 0.16 while the noise-corrected true correlation is 0.90. Donor
deltas are overwhelmingly sampling noise, so denoising should be the largest available
lever on transfer cosine -- and cosine drives pds, which is 0.017 of our 0.034 gap to #1.

Experiment 004 already killed two approaches: SVD truncation (discards spectrum) and
supervised per-gene weights (overfits, distorts direction). Both failed for reasons that
do NOT apply here:

  A. Empirical-Bayes shrinkage. Per gene, split the observed across-perturbation variance
     into true-effect variance and Poisson sampling variance, and shrink by the ratio.
     Unlike exp 004's weights this is estimated from the DONOR alone -- no target data, so
     nothing to overfit.
  B. kNN smoothing over similar perturbations. Replace each delta with a blend of itself
     and the mean of its nearest neighbours in donor effect space. Averages independent
     noise while preserving the full gene spectrum. This is the "retrieve-and-aggregate"
     shape that TransPert (3rd/338 at VCC 2025) and Genentech's LGR both landed on.
"""
import sys, json
from pathlib import Path
import numpy as np, anndata as ad
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TS

a = ad.read_h5ad("data/harness/h1-test.h5ad"); X = a.X.tocsr()
lab = a.obs["target_gene"].astype(str).to_numpy(); genes = np.asarray(a.var_names)
Pc = np.asarray(X[lab == "non-targeting"].sum(0)).ravel().astype(np.float64)
b = np.log1p(TS * Pc / Pc.sum())
ps, ds = [], []
for g in sorted(set(lab) - {"non-targeting"}):
    P = np.asarray(X[lab == g].sum(0)).ravel().astype(np.float64)
    ds.append(np.log1p(TS * P / P.sum()) - b); ps.append(g)
true = np.vstack(ds); ps = np.asarray(ps)

d = Donor.load("data/donors/k562_gwps")
gi = {g: i for i, g in enumerate(genes)}
tc, dc = [], []
for j, g in enumerate(d.genes):
    t = gi.get(str(g))
    if t is not None: tc.append(t); dc.append(j)
tc, dc = np.asarray(tc), np.asarray(dc)
didx = {p: i for i, p in enumerate(d.perts)}
sel = np.array([p in didx for p in ps])
rows = [didx[p] for p in ps[sel]]
T = true[sel][:, tc]
m, s = d.decompose()
S_all = s[:, dc]
S = S_all[rows]
print(f"{len(rows)} perts x {len(tc)} genes\n")

def cosm(P, Q):
    Pc_ = P - P.mean(0, keepdims=True); Qc = Q - Q.mean(0, keepdims=True)
    na, nb = np.linalg.norm(Pc_, axis=1), np.linalg.norm(Qc, axis=1)
    ok = (na > 0) & (nb > 0)
    return float(np.mean(np.sum(Pc_[ok]*Qc[ok], axis=1)/(na[ok]*nb[ok])))

res = {"baseline": cosm(S, T)}
print(f"  {'baseline':44} {res['baseline']:+.4f}")

# --- A. empirical-Bayes shrinkage -------------------------------------------------
# Poisson noise on a log-pseudobulk delta: var ~ 1/N where N is the gene's count total.
# N ~ n_cells * depth * p_gene, and p_gene comes from the donor's own control profile.
e = np.expm1(d.b_ctrl[dc]); p_gene = e / e.sum()
meta = json.loads(Path("data/donors/k562_gwps/meta.json").read_text())
n_cells = meta["cells_per_pert"]["median"]; depth = 20000.0
sigma2 = 1.0 / np.maximum(n_cells * depth * p_gene, 1e-3)
V = S_all.var(0)
for k in (0.5, 1.0, 2.0):
    tau2 = np.maximum(V - k * sigma2, 0.0)
    w = tau2 / (tau2 + k * sigma2 + 1e-12)
    c = cosm(S * w, T); res[f"EB shrink (sigma2 x{k})"] = c
    print(f"  {'EB shrink, noise scale '+str(k):44} {c:+.4f}   {c-res['baseline']:+.4f}")

# --- B. kNN smoothing over similar perturbations ----------------------------------
U = S_all / (np.linalg.norm(S_all, axis=1, keepdims=True) + 1e-12)
print()
for K in (5, 10, 25, 50):
    # cosine similarity of every evaluated perturbation to the whole donor panel
    sim = U[rows] @ U.T
    sim[np.arange(len(rows)), rows] = -np.inf          # exclude self
    nn = np.argpartition(-sim, K, axis=1)[:, :K]
    nbr = S_all[nn].mean(1)
    for lam in (0.25, 0.5, 0.75):
        c = cosm((1 - lam) * S + lam * nbr, T)
        res[f"kNN K={K} lam={lam}"] = c
        flag = "  <--" if c > max(v for kk, v in res.items() if kk != f"kNN K={K} lam={lam}") else ""
        print(f"  {'kNN smooth K='+str(K)+' lambda='+str(lam):44} {c:+.4f}   {c-res['baseline']:+.4f}{flag}")

best = max(res, key=res.get)
print(f"\nbest: {best} -> {res[best]:+.4f}  vs baseline {res['baseline']:+.4f} "
      f"({100*(res[best]-res['baseline'])/res['baseline']:+.1f}%)")
Path("experiments/009-denoise/results.json").write_text(json.dumps(res, indent=2))
