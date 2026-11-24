"""491a — Does the recipient's own control population carry its per-gene programme once a donor consensus orients it?

441 found single H1 control PCs weakly aligned with H1's mean response (<= 0.18), and 441 ps6 found the control
identity axis worth ~+0.001 on fid-like. Untested there: (1) how much of the mean response lies in the top-k control
SUBSPACE (160 measured 0.71-0.89 at k=50 on Orion, where a technical blob may drive it); (2) transporting the CLEAN4
consensus through the recipient's own covariance: projection P_k m_C, and the fluctuation-response (reweighting) form
Cov_ctrl(y_g / mean_g, h), h = cell score on m_C, which is what a population shift along m_C does to arithmetic means.

Instrument: 462f's fid-like reading (sign agreement with every h1-test target's reference lfc on the 2,200 top
control-CPM block genes, own gene out; KA: CLEAN4 0.5244, LOTO oracle 0.5544) and 441's gate-block accuracy
(significant block genes; CLEAN4 0.580, LOTO oracle 0.645). Control halves: every derived programme is also rebuilt
from control half A and read against lfcs recomputed vs half B, so a shared-control artefact cannot pass.

usage: uv run python experiments/491-own-control-programme/frr491.py   (peak RSS ~2.5 GB)
"""
import csv, glob, json, time
from pathlib import Path
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "work"
OUT.mkdir(exist_ok=True)
t0 = time.time()

# ---------------------------------------------------------------- truth (323's reference table, as 462f reads it)
z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=True)
genes, keep, targets = z["genes"].astype(str), z["keep"], list(z["targets"].astype(str))
lfc, q, bctrl = z["lfc"], z["q"], z["b_ctrl"]
kg = genes[keep]
G = kg.size
e = np.expm1(bctrl[keep])
cpm_ctrl_pb = 1e6 * e / e.sum()
tset = set(targets)
order = np.argsort(-cpm_ctrl_pb)
block = np.zeros(G, bool)
block[[i for i in order if kg[i] not in tset][:2200]] = True

# ---------------------------------------------------------------- CLEAN4 (the incumbent block prior)
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
offi = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
opos = {g: i for i, g in enumerate(offi)}
c4p = np.asarray(c4["P_down"], float)
c4s = np.asarray(c4["psign"], float)
pdown = np.array([c4p[opos[g]] if g in opos else np.nan for g in kg])
clean4_sign = np.array([c4s[opos[g]] if g in opos and np.isfinite(c4p[opos[g]]) else np.nan for g in kg])
mC = np.where(np.isfinite(pdown), 0.5 - pdown, 0.0)  # positive = up; 0 where CLEAN4 is silent

# ---------------------------------------------------------------- h1-test cells (streamed CSR)
f = h5py.File(ROOT / "data/harness/h1-test.h5ad", "r")
cats = f["obs"]["target_gene"]["categories"][:].astype(str)
codes = f["obs"]["target_gene"]["codes"][:]
ntc_code = int(np.where(cats == "non-targeting")[0][0])
indptr = f["X"]["indptr"][:]
kept_col = np.full(genes.size, -1, np.int64)
kept_col[np.where(keep)[0]] = np.arange(G)


def load_rows(rows):
    """per-cell CPM (library over ALL genes, as the scorer normalises), kept genes only, float32."""
    rows = np.sort(rows)
    Yr = np.zeros((rows.size, G), np.float32)
    libs = np.zeros(rows.size)
    # group rows into contiguous file spans of <= 4000 rows to keep reads cheap
    i = 0
    while i < rows.size:
        a = rows[i]
        j = i
        while j + 1 < rows.size and rows[j + 1] < a + 4000:
            j += 1
        b = rows[j] + 1
        lo, hi = indptr[a], indptr[b]
        dat = f["X"]["data"][lo:hi].astype(np.float64)
        ind = f["X"]["indices"][lo:hi]
        for k in range(i, j + 1):
            r = rows[k]
            s, t = indptr[r] - lo, indptr[r + 1] - lo
            d, ix = dat[s:t], ind[s:t]
            L = d.sum()
            libs[k] = L
            c = kept_col[ix]
            m = c >= 0
            Yr[k, c[m]] = (1e6 * d[m] / L).astype(np.float32)
        i = j + 1
    return Yr, libs


ctrl_rows = np.where(codes == ntc_code)[0]
Y, lib = load_rows(ctrl_rows)
print(f"controls {Y.shape}, median lib {np.median(lib):.0f}, {time.time()-t0:.0f}s", flush=True)
m_all = Y.mean(0).astype(np.float64)

# real KA: recompute three targets' lfc from their own cells and compare with the reference table
eps = 1e-9
ka = {}
for t in (0, 40, 84):
    tr = np.where(codes == int(np.where(cats == targets[t])[0][0]))[0]
    Yt, _ = load_rows(tr)
    l = np.log2((Yt.mean(0).astype(np.float64) + eps) / (m_all + eps))
    ok = np.isfinite(lfc[t]) & np.isfinite(l)
    ka[targets[t]] = float(np.max(np.abs(l[ok] - lfc[t][ok])))
    del Yt
f.close()
print("KA lfc recompute max |d|:", ka, flush=True)

# ---------------------------------------------------------------- truth lfcs against control half B
rng = np.random.default_rng(491)
perm = rng.permutation(ctrl_rows.size)
halfA, halfB = np.sort(perm[: ctrl_rows.size // 2]), np.sort(perm[ctrl_rows.size // 2:])
m_B = Y[halfB].mean(0).astype(np.float64)
with np.errstate(over="ignore", invalid="ignore"):
    m_pert = np.exp2(lfc) * (m_all + eps) - eps
    lfc_B = np.log2((np.maximum(m_pert, 0) + eps) / (m_B + eps))
lfc_B = np.where(np.isfinite(lfc), lfc_B, np.nan)
gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0)


def read(pred, L=lfc):
    """fid-like (all block genes vs ref lfc sign) and gate-block accuracy for a per-gene programme (its sign)."""
    ps = np.sign(pred)
    hit = n = gh = gn = 0
    for t in range(len(targets)):
        own = kg == targets[t]
        g = block & ~own & np.isfinite(L[t]) & (L[t] != 0) & np.isfinite(ps) & (ps != 0)
        hit += (np.sign(L[t][g]) == ps[g]).sum()
        n += g.sum()
        gg = g & gate[t]
        gh += (np.sign(L[t][gg]) == ps[gg]).sum()
        gn += gg.sum()
    return dict(fidlike=round(hit / max(n, 1), 4), pairs=int(n), gate_block=round(gh / max(gn, 1), 4), gate_pairs=int(gn))


def loto_oracle(L=lfc):
    tot = np.nansum(L, 0)
    hit = n = gh = gn = 0
    for t in range(len(targets)):
        pred = np.sign((tot - np.nan_to_num(L[t])) / (len(targets) - 1))
        own = kg == targets[t]
        g = block & ~own & np.isfinite(L[t]) & (L[t] != 0) & (pred != 0)
        hit += (np.sign(L[t][g]) == pred[g]).sum()
        n += g.sum()
        gg = g & gate[t]
        gh += (np.sign(L[t][gg]) == pred[gg]).sum()
        gn += gg.sum()
    return dict(fidlike=round(hit / n, 4), pairs=int(n), gate_block=round(gh / gn, 4), gate_pairs=int(gn))


res = {"ka": {"clean4": read(clean4_sign), "oracle_loto": loto_oracle(), "oracle_loto_halfB": loto_oracle(lfc_B),
              "lfc_recompute_max_absdiff": ka}}
print("KA", json.dumps(res["ka"]), flush=True)


def topk_basis(Lc, k, colscale=None, seed=0):
    """randomized SVD (2 power iterations) of Lc * diag(1/colscale); returns an orthonormal gene-space basis."""
    r = np.random.default_rng(seed)
    Om = r.standard_normal((Lc.shape[1], k + 20)).astype(np.float32)
    cs = np.ones(Lc.shape[1], np.float32) if colscale is None else (1.0 / colscale).astype(np.float32)
    Qm = Lc @ (Om * cs[:, None])
    for _ in range(2):
        Qm, _ = np.linalg.qr(Qm)
        Qm = Lc @ (((Lc.T @ Qm) * cs[:, None]) * cs[:, None])
    Qm, _ = np.linalg.qr(Qm)
    B = (Qm.T @ Lc) * cs[None, :]
    _, S, Vt = np.linalg.svd(B, full_matrices=False)
    return Vt[:k].T.astype(np.float64), S[:k]


def programmes(idx):
    Yh = Y if idx is None else Y[idx]
    n = Yh.shape[0]
    mu = Yh.mean(0).astype(np.float64)
    Lc = np.log1p(Yh)
    Lc -= Lc.mean(0)
    sd = (Lc.std(0) + 1e-6).astype(np.float64)
    out = {}
    V, _ = topk_basis(Lc, 200)
    for k in (1, 2, 5, 10, 20, 50, 100, 200):
        Vk = V[:, :k]
        out[f"proj_k{k}"] = Vk @ (Vk.T @ mC)
    Vz, _ = topk_basis(Lc, 100, colscale=sd, seed=1)
    for k in (5, 20, 50, 100):
        Vk = Vz[:, :k]
        out[f"projz_k{k}"] = Vk @ (Vk.T @ mC)
    # fluctuation-response / reweighting: h = cell score on m_C in the standardized frame
    hr = (Lc @ (mC / sd).astype(np.float32)).astype(np.float64)
    hc = (hr - hr.mean()) / (hr.std() + 1e-12)
    h32 = hc.astype(np.float32)
    out["frr_lin"] = (Yh.T @ h32).astype(np.float64) / n / (mu + 1e-3)  # Cov(y_g, h) / mean_g ~ d lfc
    zc = (Lc.T @ h32).astype(np.float64) / n / sd                        # Cov(z_g, h)
    out["frr_z"] = zc
    # remove the diagonal echo Var(z_g) m_C(g) / sd(h): what is left is pure off-diagonal transport
    out["frr_z_offdiag"] = zc - mC / (hr.std() + 1e-12)
    for j in range(3):
        v = V[:, j]
        out[f"pc{j+1}_oriented"] = v * np.sign(v @ mC)
    mresp = np.nanmean(lfc, 0)
    mresp = np.where(np.isfinite(mresp), mresp, 0.0)
    cap = {}
    for k in (1, 5, 10, 20, 50, 100, 200):
        Vk = V[:, :k]
        cap[k] = dict(mean_resp=round(float(np.sum((Vk.T @ mresp) ** 2) / np.sum(mresp ** 2)), 4),
                      mean_resp_block=round(float(np.sum((Vk[block].T @ mresp[block]) ** 2) / np.sum(mresp[block] ** 2)), 4),
                      clean4=round(float(np.sum((Vk.T @ mC) ** 2) / np.sum(mC ** 2)), 4),
                      random=round(k / G, 4))
    corr = {name: round(float(np.corrcoef(v[block], mresp[block])[0, 1]), 4) for name, v in out.items()}
    corr["clean4"] = round(float(np.corrcoef(mC[block], mresp[block])[0, 1]), 4)
    del Lc
    return out, cap, corr


for tag, idx, L in (("all", None, lfc), ("halfA_vs_B", halfA, lfc_B)):
    P, cap, corr = programmes(idx)
    rows = {"clean4": read(clean4_sign, L)}
    a = np.where(np.isfinite(clean4_sign), mC, 0.0)
    za = a / (a[block].std() + 1e-12)
    for name, v in P.items():
        rows[name] = read(v, L)
        zb = v / (v[block].std() + 1e-12)
        for lam in (0.5, 1.0, 2.0):
            rows[f"{name}+c4_l{lam}"] = read(za + lam * zb, L)
    res[tag] = dict(read=rows, capture=cap, corr_with_mean_response_block=corr)
    print(tag, "capture", json.dumps(cap), flush=True)
    print(tag, "corr", json.dumps(corr), flush=True)
    for k, v in sorted(rows.items(), key=lambda kv: -kv[1]["fidlike"])[:10]:
        print(f"  {tag:10s} {k:28s} fidlike {v['fidlike']:.4f}  gate_block {v['gate_block']:.4f}", flush=True)
    for k, v in sorted(rows.items(), key=lambda kv: -kv[1]["gate_block"])[:5]:
        print(f"  {tag:10s} {k:28s} fidlike {v['fidlike']:.4f}  gate_block {v['gate_block']:.4f}  (gate-sorted)", flush=True)
    print(f"  {tag:10s} {'clean4':28s} fidlike {rows['clean4']['fidlike']:.4f}  gate_block {rows['clean4']['gate_block']:.4f}", flush=True)

(OUT / "frr491.json").write_text(json.dumps(res, indent=1))
print(f"done {time.time()-t0:.0f}s")
