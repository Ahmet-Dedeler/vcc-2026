"""651 steps 2-3 (PLAN.md, committed 47b018f1 before this script ran).
Step 2 (oracle, descriptive): beta_g from h1-test controls x the observed per-target library shift.
Step 3 (DECISIVE): the same with the shift predicted from K562 GWPS (same target, log UMI vs gem-group NTC mean).
s_Delta for the K562 predictor = sd over all GWPS perturbations with >= 20 cells (a K562-only quantity; pinned here
before running). Analysis only, ~2 GB. Usage: uv run python experiments/651-libsize-filler-sign/run651.py"""
import numpy as np, anndata as ad, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr
from pathlib import Path
W = Path(__file__).resolve().parent / "work"
rng = np.random.default_rng(651)
d = np.load("experiments/213-transport-law/work/de_real_h1.npz", allow_pickle=True)
genes = [str(g) for g in d["genes"]]; perts = [str(p) for p in d["perts"]]
padj = d["real_padj"].copy(); lfc = d["real_lfc"].copy()
P, G = lfc.shape; gi = {g: i for i, g in enumerate(genes)}
own = np.zeros((P, G), bool)
for p, t in enumerate(perts):
    if t in gi: own[p, gi[t]] = True
adj = np.isfinite(lfc) & (lfc != 0) & ~own
null = adj & ~(padj < 0.05)
# ---- beta_g from controls
a = ad.read_h5ad("data/harness/h1-test.h5ad", backed="r")
tt = a.obs["target_gene"].astype(str).values
ci = np.where(tt == "non-targeting")[0]
X = sp.csr_matrix(a[ci].to_memory().X)
L = np.asarray(X.sum(1)).ravel()
vn = {g: i for i, g in enumerate(a.var_names)}
col = np.array([vn[g] for g in genes])
C = (sp.diags(1e6 / L) @ X).tocsc()[:, col]
mu = np.asarray(C.mean(0)).ravel()
x = np.log(L) - np.log(L).mean()
beta = (C.T @ x) / len(x) / np.mean(x * x) / mu
b = (beta - beta.mean()) / beta.std()
np.savez(W / "beta651.npz", genes=np.array(genes), beta=beta, mu=mu)
print("beta: mean %.3f sd %.3f; corr(beta, log mu) %.3f" % (beta.mean(), beta.std(), np.corrcoef(beta, np.log(mu))[0, 1]))
order = np.argsort(-mu)
blk = np.zeros(G, bool); blk[order[:2200]] = True
M = null & blk[None, :]
print("block null genes per target: median %d" % np.median(M.sum(1)))
# ---- CLEAN4 baseline (391)
c4 = np.load("experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
c4g = {str(g): i for i, g in enumerate(c4["genes"])}
pdn = np.array([c4["P_down"][c4g[g]] if g in c4g else np.nan for g in genes])
cov = np.isfinite(pdn)
s_c4 = np.where(cov, np.where(pdn > 0.5, -1, 1), -1)
pz = np.where(cov, pdn, 0.5 + 1e-6) - 0.5
zc4 = pz / pz[cov].std()
truth = np.sign(np.nan_to_num(lfc))
def pnull(S, rows):
    return np.array([(S[p][M[p]] == truth[p][M[p]]).mean() for p in rows])
def boot(dv):
    idx = rng.integers(0, len(dv), (4000, len(dv)))
    m = dv[idx].mean(1); return np.percentile(m, [2.5, 97.5])
def combined(Dhat, sD, lam):
    S = np.empty((P, G), int)
    for p in range(P):
        u = -zc4 + lam * b * Dhat[p] / sD
        S[p] = np.where(u > 0, 1, -1)
    return S
S0 = np.tile(s_c4, (P, 1))
dl = pd.read_csv(W / "dlogL_h1.csv").set_index("target").dlogL
D = np.array([dl[t] for t in perts])
allrows = np.arange(P)
base = pnull(S0, allrows)
print("\n## step 2 (ORACLE: observed h1 shift)")
print("CLEAN4 p_null %.4f" % base.mean())
# per-target regression slope of null-block lfc on beta (expect slope*ln2 ~ Delta_t if the cross-sectional model holds)
sl = np.array([np.polyfit(beta[M[p]], lfc[p, M[p]], 1)[0] for p in range(P)])
print("per-target slope(lfc on beta)*ln2 vs Delta_t: spearman %.3f, OLS slope %.3f" % (spearmanr(sl * np.log(2), D)[0], np.polyfit(D, sl * np.log(2), 1)[0]))
So = np.array([np.where(b * D[p] > 0, 1, -1) for p in range(P)])
r = pnull(So, allrows); print("sign(b*Delta) alone p_null %.4f  (vs CLEAN4 %+.4f)" % (r.mean(), (r - base).mean()))
sDo = D.std()
for lam in (0.5, 1, 2):
    r = pnull(combined(D, sDo, lam), allrows); dv = r - base
    print("combined lam %.1f p_null %.4f  d %+.4f  CI [%+.4f, %+.4f]" % ((lam, r.mean(), dv.mean()) + tuple(boot(dv))))
big = np.abs(D) > 3 * 0.0173
r = pnull(combined(D, sDo, 1), allrows)
print("  lam 1 on the %d targets with |Delta| > 3 se: d %+.4f; others %+.4f" % (big.sum(), (r - base)[big].mean(), (r - base)[~big].mean()))
# ---- step 3: K562 GWPS predictor
print("\n## step 3 (DECISIVE: K562 GWPS predictor)")
g = ad.read_h5ad("data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad", backed="r")
o = g.obs[["gene", "batch", "UMI_count"]].copy()
o["gene"] = o.gene.astype(str); o["lu"] = np.log(o.UMI_count.astype(float))
ntc = o[o.gene == "non-targeting"].groupby("batch", observed=True).lu.mean()
o["rel"] = o.lu - o.batch.map(ntc).astype(float)
k = o[o.gene != "non-targeting"].groupby("gene").rel.agg(["mean", "size"])
k = k[k["size"] >= 20]
sD = k["mean"].std()
print("GWPS perturbations with >= 20 cells: %d; s_Delta %.4f" % (len(k), sD))
rows = np.array([p for p, t in enumerate(perts) if t in k.index])
Dh = np.zeros(P); Dh[rows] = [k.loc[perts[p], "mean"] for p in rows]
rho = spearmanr(Dh[rows], D[rows])[0]
print("covered h1 targets %d / %d; Spearman(K562 Delta-hat, h1 Delta) %.3f" % (len(rows), P, rho))
b0 = pnull(S0, rows)
res = {}
for lam in (0.5, 1, 2):
    r = pnull(combined(Dh, sD, lam), rows); dv = r - b0; ci_ = boot(dv)
    res[lam] = (dv.mean(), ci_)
    print("combined lam %.1f p_null %.4f  d %+.4f  CI [%+.4f, %+.4f]  (CLEAN4 on covered %.4f)" % ((lam, r.mean(), dv.mean()) + tuple(ci_) + (b0.mean(),)))
d1, ci1 = res[1]
go = (rho >= 0.30) and (d1 >= 0.005) and (ci1[0] > 0)
print("\nVERDICT (PLAN bar at lam 1: rho >= 0.30 [%s], d >= +0.005 [%s], CI lower > 0 [%s]): %s" % (
    rho >= 0.30, d1 >= 0.005, ci1[0] > 0, "GO" if go else "CLOSED"))
pd.DataFrame({"target": perts, "dlogL_h1": D, "dhat_k562": np.where(np.isin(allrows, rows), Dh, np.nan), "slope_ln2": sl * np.log(2)}).to_csv(W / "targets651.csv", index=False)
