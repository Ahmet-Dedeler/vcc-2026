"""491b — What are H1's top control PCs, which of them carries the mean response, and does the same axis exist in the
A/B/C controls? Also: target bootstrap and design-strength split of 491a's best own-control arm (proj_k2 + CLEAN4).

Reads the same inputs as frr491.py. Memory-lean: one float32 log1p matrix at a time (<= ~2.5 GB).
usage: uv run python experiments/491-own-control-programme/pcs491.py
"""
import csv, glob, json, time
from pathlib import Path
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "work"
t0 = time.time()

z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str))
lfc, q, bctrl = z["lfc"], z["q"], z["b_ctrl"]
kg = genes[keep]
G = kg.size
e = np.expm1(bctrl[keep])
cpm_pb = 1e6 * e / e.sum()
tset = set(targets)
order = np.argsort(-cpm_pb)
block = np.zeros(G, bool)
block[[i for i in order if kg[i] not in tset][:2200]] = True
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0)

c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
opos = {g: i for i, g in enumerate(offi)}
c4p = np.asarray(c4["P_down"], float)
c4s = np.asarray(c4["psign"], float)
pdown = np.array([c4p[opos[g]] if g in opos else np.nan for g in kg])
clean4_sign = np.array([c4s[opos[g]] if g in opos and np.isfinite(c4p[opos[g]]) else np.nan for g in kg])
mC = np.where(np.isfinite(pdown), 0.5 - pdown, 0.0)


def load_ctrl(path, gene_axis, want_code=None):
    """per-cell log1p CPM (library over all genes of the file) of the NTC cells, projected onto gene_axis, float32."""
    f = h5py.File(path, "r")
    obs = f["obs"]
    tg = obs["target_gene"]
    cats = tg["categories"][:].astype(str)
    codes = tg["codes"][:]
    rows = np.where(codes == int(np.where(cats == "non-targeting")[0][0]))[0]
    import anndata as ad
    vidx = np.asarray(ad.io.read_elem(f["var"]).index.astype(str))
    vpos = {g: i for i, g in enumerate(vidx)}
    col = np.full(vidx.size, -1, np.int64)
    for j, g in enumerate(gene_axis):
        if g in vpos:
            col[vpos[g]] = j
    indptr = f["X"]["indptr"][:]
    X = np.zeros((rows.size, len(gene_axis)), np.float32)
    i = 0
    while i < rows.size:
        a = rows[i]
        j = i
        while j + 1 < rows.size and rows[j + 1] < a + 4000:
            j += 1
        lo, hi = indptr[a], indptr[rows[j] + 1]
        dat = f["X"]["data"][lo:hi].astype(np.float64)
        ind = f["X"]["indices"][lo:hi]
        for k in range(i, j + 1):
            r = rows[k]
            s, t = indptr[r] - lo, indptr[r + 1] - lo
            d, ix = dat[s:t], ind[s:t]
            L = d.sum()
            c = col[ix]
            m = c >= 0
            X[k, c[m]] = np.log1p(1e6 * d[m] / L).astype(np.float32)
        i = j + 1
    f.close()
    return X


def topk(Lc, k, seed=0):
    r = np.random.default_rng(seed)
    Om = r.standard_normal((Lc.shape[1], k + 20)).astype(np.float32)
    Qm = Lc @ Om
    for _ in range(2):
        Qm, _ = np.linalg.qr(Qm)
        Qm = Lc @ (Lc.T @ Qm)
    Qm, _ = np.linalg.qr(Qm)
    B = Qm.T @ Lc
    _, S, Vt = np.linalg.svd(B, full_matrices=False)
    return Vt[:k].T.astype(np.float64), S[:k]


res = {}
X = load_ctrl(ROOT / "data/harness/h1-test.h5ad", kg)
X -= X.mean(0)
V, S = topk(X, 10)
totvar = float(np.einsum('ij,ij->', X, X, dtype=np.float64))
lib_axis = None
mresp = np.nanmean(lfc, 0)
mresp = np.where(np.isfinite(mresp), mresp, 0.0)
pcs = {}
for j in range(6):
    v = V[:, j]
    top_up = [kg[i] for i in np.argsort(-v)[:15]]
    top_dn = [kg[i] for i in np.argsort(v)[:15]]
    pcs[f"PC{j+1}"] = dict(var_share=round(float(S[j] ** 2 / totvar), 4),
                          corr_meanresp_block=round(float(np.corrcoef(v[block], mresp[block])[0, 1]), 4),
                          corr_meanresp_all=round(float(np.corrcoef(v, mresp)[0, 1]), 4),
                          corr_clean4_block=round(float(np.corrcoef(v[block], mC[block])[0, 1]), 4),
                          top_pos=top_up, top_neg=top_dn)
    print(f"PC{j+1} var {pcs[f'PC{j+1}']['var_share']} corr_mr_block {pcs[f'PC{j+1}']['corr_meanresp_block']} "
          f"corr_mr_all {pcs[f'PC{j+1}']['corr_meanresp_all']} corr_c4 {pcs[f'PC{j+1}']['corr_clean4_block']}\n  + {top_up}\n  - {top_dn}", flush=True)
res["h1_pcs"] = pcs
# cell scores on PC2 vs library size / detected genes / mito fraction
scores = X @ V[:, :3].astype(np.float32)
del X
res["elapsed_h1"] = time.time() - t0

# ------------- target bootstrap + design-strength split for proj_k2 + CLEAN4 (491a's best) vs CLEAN4
V2 = V[:, :2]
proj2 = V2 @ (V2.T @ mC)
a = np.where(np.isfinite(clean4_sign), mC, 0.0)
cand = a / (a[block].std() + 1e-12) + proj2 / (proj2[block].std() + 1e-12)


def per_target(pred):
    ps = np.sign(pred)
    hit = np.zeros(len(targets)); n = np.zeros(len(targets)); gh = np.zeros(len(targets)); gn = np.zeros(len(targets))
    for t in range(len(targets)):
        own = kg == targets[t]
        g = block & ~own & np.isfinite(lfc[t]) & (lfc[t] != 0) & np.isfinite(ps) & (ps != 0)
        hit[t] = (np.sign(lfc[t][g]) == ps[g]).sum(); n[t] = g.sum()
        gg = g & gate[t]
        gh[t] = (np.sign(lfc[t][gg]) == ps[gg]).sum(); gn[t] = gg.sum()
    return hit, n, gh, gn


H0 = per_target(clean4_sign)
H1_ = per_target(cand)
rng = np.random.default_rng(4910)
B = 4000
d_fid = np.zeros(B); d_gate = np.zeros(B)
T = len(targets)
for b in range(B):
    ix = rng.integers(0, T, T)
    d_fid[b] = H1_[0][ix].sum() / H1_[1][ix].sum() - H0[0][ix].sum() / H0[1][ix].sum()
    d_gate[b] = H1_[2][ix].sum() / H1_[3][ix].sum() - H0[2][ix].sum() / H0[3][ix].sum()
pt_fid = H1_[0].sum() / H1_[1].sum() - H0[0].sum() / H0[1].sum()
pt_gate = H1_[2].sum() / H1_[3].sum() - H0[2].sum() / H0[3].sum()
res["bootstrap_proj2_c4_vs_c4"] = dict(fidlike=[round(pt_fid, 4), round(float(np.quantile(d_fid, .025)), 4), round(float(np.quantile(d_fid, .975)), 4)],
                                       gate_block=[round(pt_gate, 4), round(float(np.quantile(d_gate, .025)), 4), round(float(np.quantile(d_gate, .975)), 4)],
                                       P_le0_fid=round(float((d_fid <= 0).mean()), 4))
print("bootstrap", json.dumps(res["bootstrap_proj2_c4_vs_c4"]), flush=True)

# design-strength split: 467's strength = count(|key| > 0.2) per target, if its table is on disk
st = None
for p in glob.glob(str(ROOT / "experiments/46[78]*/work/*strength*")) + glob.glob(str(ROOT / "experiments/46[78]*/work/*.json")):
    pass
res["note_strength"] = "strength split: see pcs491 log (computed by n_real tercile instead: weak/mid/strong by reference gate size)"
nreal = gate.sum(1)
edges = np.quantile(nreal, [1 / 3, 2 / 3])
for name, m in (("weak", nreal <= edges[0]), ("mid", (nreal > edges[0]) & (nreal <= edges[1])), ("strong", nreal > edges[1])):
    f0 = H0[0][m].sum() / H0[1][m].sum(); f1 = H1_[0][m].sum() / H1_[1][m].sum()
    g0 = H0[2][m].sum() / max(H0[3][m].sum(), 1); g1 = H1_[2][m].sum() / max(H1_[3][m].sum(), 1)
    res[f"tercile_{name}"] = dict(n=int(m.sum()), fid_c4=round(f0, 4), fid_cand=round(f1, 4), gate_c4=round(g0, 4), gate_cand=round(g1, 4))
    print(name, res[f"tercile_{name}"], flush=True)

# ------------- is H1's PC2 axis present in A/B/C? project each context's controls on H1's PC2 loadings
chal = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
for ctx in "ABC":
    Xc = load_ctrl(ROOT / f"data/challenge/context_{ctx}.h5ad", kg)
    Xc -= Xc.mean(0)
    Vc, Sc = topk(Xc, 10, seed=2)
    ov = {}
    for j in range(3):
        ov[f"H1_PC{j+1}"] = [round(float(abs(V[:, j] @ Vc[:, i])), 3) for i in range(10)]
    # variance of the context's cells along H1's PC2 direction, relative to its own top eigenvalue
    s2 = float(((Xc @ V[:, 1].astype(np.float32)) ** 2).mean())
    res[f"ctx_{ctx}"] = dict(overlap_with_own_top10=ov, var_along_H1_PC2=round(s2, 4), own_top_eig=[round(float(x ** 2 / Xc.shape[0]), 3) for x in Sc[:5]])
    print(ctx, json.dumps(res[f"ctx_{ctx}"]), flush=True)
    del Xc

(OUT / "pcs491.json").write_text(json.dumps(res, indent=1))
print(f"done {time.time()-t0:.0f}s")
