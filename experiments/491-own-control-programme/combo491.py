"""491c — a recipient-side programme for a line with no bank: CLEAN4 + the own-control transport (491a's proj_k2 / k20)
+ 441 ps6's identity axis (I = log mean ctrl_H1 - mean of the A/B/C contexts' log mean ctrl, same chemistry), combined
as z-scores on the block. Read on h1-test with 462f's fid-like instrument and 441's gate-block accuracy, target bootstrap
(4,000 draws), and the control-halves check (programme from half A, truth vs half B) for the own-control part.

usage: uv run python experiments/491-own-control-programme/combo491.py     (<= ~2.5 GB)
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

# ---- H1 controls (log1p CPM on kept genes) and their per-gene arithmetic mean CPM
f = h5py.File(ROOT / "data/harness/h1-test.h5ad", "r")
cats = f["obs"]["target_gene"]["categories"][:].astype(str)
codes = f["obs"]["target_gene"]["codes"][:]
rows = np.where(codes == int(np.where(cats == "non-targeting")[0][0]))[0]
indptr = f["X"]["indptr"][:]
kept_col = np.full(genes.size, -1, np.int64)
kept_col[np.where(keep)[0]] = np.arange(G)
X = np.zeros((rows.size, G), np.float32)
i = 0
while i < rows.size:
    a = rows[i]; j = i
    while j + 1 < rows.size and rows[j + 1] < a + 4000:
        j += 1
    lo, hi = indptr[a], indptr[rows[j] + 1]
    dat = f["X"]["data"][lo:hi].astype(np.float64); ind = f["X"]["indices"][lo:hi]
    for k in range(i, j + 1):
        r = rows[k]; s, t = indptr[r] - lo, indptr[r + 1] - lo
        d, ix = dat[s:t], ind[s:t]; c = kept_col[ix]; m = c >= 0
        X[k, c[m]] = (1e6 * d[m] / d.sum()).astype(np.float32)
    i = j + 1
f.close()
m_all = X.mean(0).astype(np.float64)
rng = np.random.default_rng(491)
perm = rng.permutation(rows.size)
hA, hB = np.sort(perm[: rows.size // 2]), np.sort(perm[rows.size // 2:])
m_B = X[hB].mean(0).astype(np.float64)
eps = 1e-9
with np.errstate(over="ignore", invalid="ignore"):
    m_pert = np.exp2(lfc) * (m_all + eps) - eps
    lfc_B = np.where(np.isfinite(lfc), np.log2((np.maximum(m_pert, 0) + eps) / (m_B + eps)), np.nan)
np.log1p(X, out=X)


def proj(idx, k):
    L = X if idx is None else X[idx]
    mu = L.mean(0)
    Lc = L - mu
    r = np.random.default_rng(0)
    Om = r.standard_normal((G, k + 20)).astype(np.float32)
    Qm = Lc @ Om
    for _ in range(2):
        Qm, _ = np.linalg.qr(Qm)
        Qm = Lc @ (Lc.T @ Qm)
    Qm, _ = np.linalg.qr(Qm)
    _, _, Vt = np.linalg.svd(Qm.T @ Lc, full_matrices=False)
    V = Vt[:k].T.astype(np.float64)
    del Lc
    return V @ (V.T @ mC)


# ---- identity axis: H1 controls vs the A/B/C controls (log of arithmetic mean CPM), 441 ps6's construction
abc = []
for ctx in "ABC":
    mc = np.load(ROOT / f"experiments/491-own-control-programme/work/ctrl_meancpm_{ctx}.npy")
    abc.append(np.array([mc[opos[g]] if g in opos else np.nan for g in kg]))
abc = np.nanmean(np.vstack([np.log(np.maximum(v, 1e-3)) for v in abc]), 0)
ident = np.log(np.maximum(m_all, 1e-3)) - abc
ident = np.where(np.isfinite(ident), ident, 0.0)


def zb(v):
    v = np.where(np.isfinite(v), v, 0.0)
    return v / (v[block].std() + 1e-12)


def per_target(pred, L):
    ps = np.sign(pred)
    hit = np.zeros(len(targets)); n = np.zeros(len(targets)); gh = np.zeros(len(targets)); gn = np.zeros(len(targets))
    for t in range(len(targets)):
        own = kg == targets[t]
        g = block & ~own & np.isfinite(L[t]) & (L[t] != 0) & np.isfinite(ps) & (ps != 0)
        hit[t] = (np.sign(L[t][g]) == ps[g]).sum(); n[t] = g.sum()
        gg = g & gate[t]
        gh[t] = (np.sign(L[t][gg]) == ps[gg]).sum(); gn[t] = gg.sum()
    return hit, n, gh, gn


def boot(H0, H1, B=4000, seed=4911):
    r = np.random.default_rng(seed)
    T = len(targets)
    df, dg = np.empty(B), np.empty(B)
    for b in range(B):
        ix = r.integers(0, T, T)
        df[b] = H1[0][ix].sum() / H1[1][ix].sum() - H0[0][ix].sum() / H0[1][ix].sum()
        dg[b] = H1[2][ix].sum() / H1[3][ix].sum() - H0[2][ix].sum() / H0[3][ix].sum()
    q = lambda a: [round(float(np.quantile(a, x)), 4) for x in (0.025, 0.5, 0.975)]
    return q(df), q(dg), round(float((df <= 0).mean()), 4)


res = {}
c4z = zb(np.where(np.isfinite(clean4_sign), mC, 0.0))
for tag, idx, L in (("all", None, lfc), ("halfA_vs_B", hA, lfc_B)):
    P2, P20 = proj(idx, 2), proj(idx, 20)
    arms = {"clean4": clean4_sign,
            "c4+proj2": c4z + zb(P2), "c4+proj20": c4z + zb(P20),
            "c4+ident(-I) l1": c4z + zb(-ident), "c4+ident(-I) l2": c4z + 2 * zb(-ident),
            "c4+proj2+ident l1": c4z + zb(P2) + zb(-ident), "c4+proj20+ident l1": c4z + zb(P20) + zb(-ident),
            "c4+proj2+ident l2": c4z + zb(P2) + 2 * zb(-ident)}
    H0 = per_target(clean4_sign, L)
    out = {}
    for name, v in arms.items():
        H = per_target(v, L)
        fl, gb = H[0].sum() / H[1].sum(), H[2].sum() / H[3].sum()
        rec = dict(fidlike=round(float(fl), 4), gate_block=round(float(gb), 4))
        if name != "clean4":
            bf, bg, p0 = boot(H0, H)
            rec.update(d_fidlike_ci=bf, d_gate_ci=bg, P_fid_le0=p0)
        out[name] = rec
        print(f"{tag:10s} {name:22s} fid-like {fl:.4f}  gate-block {gb:.4f}  " +
              ("" if name == "clean4" else f"d fid {rec['d_fidlike_ci']} d gate {rec['d_gate_ci']} P(d fid<=0) {rec['P_fid_le0']}"), flush=True)
    res[tag] = out
res["identity_corr_meanresp_block"] = float(np.corrcoef(-ident[block], np.nan_to_num(np.nanmean(lfc, 0))[block])[0, 1])
(OUT / "combo491.json").write_text(json.dumps(res, indent=1))
print(f"identity corr (block) {res['identity_corr_meanresp_block']:.3f}; done {time.time()-t0:.0f}s")
