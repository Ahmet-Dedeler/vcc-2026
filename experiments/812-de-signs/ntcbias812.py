"""812 H-NTCB — do NTC-guide-specific effects bias the scorer's control mean, and does that bias sign the truth?

Mechanism under test. The scorer's control mean M_g is the average of the NTC guides' cells. A guide with its own effect
on gene g (a seed off-target, a guide-level artefact) shifts M_g by ~(its effect) x (its cell share). Every target's real
lfc on g is measured against that same M_g, so for genes the targets do not move, the truth's sign leans AGAINST the
bias: real sign ~ -sign(M_g - pop_g). pop_g is unknown, but a robust centre of the guide means (a weighted median, a
Huber M-estimate) ignores a single outlying guide, so e_g = M_g - centre_g estimates the guide-driven part of the bias.
Under pure cell sampling (no guide effects) e_g carries no information about M_g - pop_g (Gaussian: Cov(mean - median,
mean) = 0), so accuracy of -sign(e_g) must sit at 0.5 in every |z| bin: that is the falsifier.

Surfaces: h1-test (85 targets) and h1-trainval (126 targets, 618), which share h1-test's 18,400 controls; NTC guide ids
from the 2025 VCC files (`work/stage/h1ctrl_guides.tsv`, 31 guides, 243-1,331 cells each). Truth: the surfaces'
reference DE tables (cell-eval2's own). Our calls: the hook-free zs arms' stats (as surf812).
Known answer (abort on failure): three h1-test targets' lfc recomputed from their cells against M equal the reference
table to <= 1e-4 (float32 cells; 491 matched to 1.1e-5).

usage: uv run python experiments/812-de-signs/ntcbias812.py   -> work/ntcbias812.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ST = HERE / "work/stage"
t0 = time.time()
genes = np.loadtxt(ST / "genes_keep_h1test.txt", dtype=str)
G = genes.size
gpos = {g: i for i, g in enumerate(genes)}

# ------------------------------------------------------------------ control cells (h5py row reads, per-cell CPM over ALL genes)
f = h5py.File(ROOT / "data/harness/h1-test.h5ad", "r")
def _index(grp):
    x = grp[grp.attrs["_index"]]
    if isinstance(x, h5py.Group) and "values" in x:          # nullable-string-array
        return x["values"][:].astype(str)
    if isinstance(x, h5py.Group):                            # categorical
        return x["categories"][:].astype(str)[x["codes"][:]]
    return x[:].astype(str)


var_ix = _index(f["var"])
kept_col = np.array([gpos.get(g, -1) for g in var_ix])
obs = f["obs"]
tg = obs["target_gene"]
cats = tg["categories"][:].astype(str)
codes = tg["codes"][:]
barcodes = _index(obs)
indptr = f["X"]["indptr"][:]


def load_rows(rows):
    rows = np.sort(rows)
    Y = np.zeros((rows.size, G), np.float32)
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
            c = kept_col[ix]
            m = c >= 0
            Y[k, c[m]] = (1e6 * d[m] / d.sum()).astype(np.float32)
        i = j + 1
    return rows, Y


ntc = int(np.where(cats == "non-targeting")[0][0])
rows, Y = load_rows(np.where(codes == ntc)[0])
gmap = pd.read_csv(ST / "h1ctrl_guides.tsv", sep="\t", index_col=0)
guide = gmap.loc[barcodes[rows], "guide_id"].astype(str).to_numpy()
batch = gmap.loc[barcodes[rows], "batch"].astype(str).to_numpy()
N = Y.shape[0]
M = Y.mean(0).astype(np.float64)
print(f"controls {Y.shape} in {time.time()-t0:.0f}s; guides {len(set(guide))}, batches {len(set(batch))}", flush=True)

# KA: recompute three targets' lfc against M and compare with the reference DE table
ref = pq.read_table(ROOT / "experiments/323-current-family-validation/work/reference_de.parquet").to_pandas()
eps = 1e-9
ka = {}
for tname in ("AASS", "TSKU", "GULP1"):
    _, Yt = load_rows(np.where(codes == int(np.where(cats == tname)[0][0]))[0])
    l = np.log2((Yt.mean(0).astype(np.float64) + eps) / (M + eps))
    r = ref[ref.target == tname].set_index("feature")["log2_fold_change"].reindex(genes).to_numpy()
    ok = np.isfinite(r) & np.isfinite(l)
    ka[tname] = float(np.max(np.abs(l[ok] - r[ok])))
print("KA lfc recompute max |d|:", ka, flush=True)
assert max(ka.values()) < 1e-4, f"KA FAIL {ka}"
f.close()


# ------------------------------------------------------------------ group structure and robust centres
def groups(lab):
    u, inv = np.unique(lab, return_inverse=True)
    n = np.bincount(inv)
    S = np.zeros((u.size, G))
    S2 = np.zeros((u.size, G))
    for k in range(u.size):
        Yk = Y[inv == k].astype(np.float64)
        S[k] = Yk.sum(0)
        S2[k] = (Yk ** 2).sum(0)
    m = S / n[:, None]
    v = np.maximum(S2 / n[:, None] - m ** 2, 0) * n[:, None] / np.maximum(n[:, None] - 1, 1)
    return u, n, m, v


def wmedian(m, w):
    """weighted median over axis 0 of m (k x G), weights w (k,)."""
    o = np.argsort(m, 0)
    ms = np.take_along_axis(m, o, 0)
    ws = w[o]
    cw = np.cumsum(ws, 0) / ws.sum(0)
    k = (cw >= 0.5).argmax(0)
    return ms[k, np.arange(m.shape[1])]


def huber(m, se, n, k=2.0, it=30):
    c = wmedian(m, n.astype(float))
    for _ in range(it):
        r = (m - c[None, :]) / np.maximum(se, 1e-12)
        w = np.where(np.abs(r) <= k, 1.0, k / np.maximum(np.abs(r), 1e-12)) * n[:, None]
        c = (w * m).sum(0) / w.sum(0)
    return c


out = dict(ka=ka, n_controls=int(N))
cellvar = Y.astype(np.float64).var(0, ddof=1)
seM = np.sqrt(cellvar / N)
E = {}
for nm, lab in (("guide", guide), ("batch", batch)):
    u, n, m, v = groups(lab)
    se = np.sqrt(v / n[:, None])
    between = (n[:, None] * (m - M[None, :]) ** 2).sum(0) / (u.size - 1)       # n-weighted between-group mean square
    R = between / np.maximum(cellvar, 1e-30)                                  # ~1 under pure cell sampling
    keep5 = M > 5
    out[f"{nm}_R_quantiles"] = {str(q): round(float(np.quantile(R[keep5], q)), 3) for q in (0.1, 0.5, 0.9, 0.99)}
    out[f"{nm}_share_R_gt_2"] = round(float((R[keep5] > 2).mean()), 4)
    E[f"{nm}_wmed"] = (M - wmedian(m, n.astype(float))) / seM
    E[f"{nm}_huber"] = (M - huber(m, se, n)) / seM
    print(nm, u.size, "groups; R q50/q90/q99", out[f"{nm}_R_quantiles"], flush=True)

# ------------------------------------------------------------------ the truth on both surfaces
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
cp = {g: i for i, g in enumerate(map(str, c4["genes"]))}
c4p = np.array([c4["P_down"][cp[g]] if g in cp else np.nan for g in genes])
c4s = np.where(np.isfinite(c4p), np.where(c4p > 0.5, -1.0, 1.0), 0.0)
SURF = {"h1test": (ROOT / "experiments/323-current-family-validation/work/reference_de.parquet",
                   ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz"),
        "h1tv": (ROOT / "experiments/618-h1-trainval-surface/work/ref_tv/reference_de.parquet",
                 ROOT / "experiments/618-h1-trainval-surface/work/g1/stats-base_s0_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz")}
for surf, (refp, stp) in SURF.items():
    targets = np.loadtxt(ST / f"targets_{surf}.txt", dtype=str)
    T = targets.size
    rd = pq.read_table(refp).to_pandas()
    ti = {t: i for i, t in enumerate(targets)}
    rl = np.full((T, G), np.nan)
    rq = np.full((T, G), np.nan)
    a_ = rd.target.map(ti).to_numpy()
    b_ = rd.feature.map(gpos).to_numpy()
    ok = ~(pd.isna(a_) | pd.isna(b_))
    a_ = a_[ok].astype(int)
    b_ = b_[ok].astype(int)
    rl[a_, b_] = rd.log2_fold_change.to_numpy()[ok]
    rq[a_, b_] = rd.p_adj.to_numpy()[ok]
    st = np.load(stp)
    assert list(st["targets"].astype(str)) == list(targets)
    oq, ol = st["q"], st["lfc"]
    own = np.equal.outer(targets, genes)
    adj = np.isfinite(rl) & (rl != 0) & ~own
    sig = (rq < 0.05) & adj
    call = (oq < 0.05) & ~own
    sets = {"nonsig_all": adj & ~sig, "M0_called_nonsig": call & adj & ~sig, "sig_all": sig}
    R = {"real_down_share_nonsig": round(float((rl[adj & ~sig] < 0).mean()), 4)}
    for en, z in E.items():
        pred = -np.sign(z)
        rr = {}
        for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 1e9)):
            gm = (np.abs(z) >= lo) & (np.abs(z) < hi) & (M > 5)
            row = dict(genes=int(gm.sum()))
            for sn, sm in sets.items():
                mm = sm & gm[None, :]
                n = int(mm.sum())
                row[sn] = dict(n=n, acc=round(float((pred[None, :].repeat(T, 0)[mm] == np.sign(rl[mm])).mean()), 4) if n else None,
                               acc_c4=round(float((c4s[None, :].repeat(T, 0)[mm & (c4s != 0)[None, :]] == np.sign(rl[mm & (c4s != 0)[None, :]])).mean()), 4)
                               if (mm & (c4s != 0)[None, :]).sum() else None)
            rr[f"|z|{lo}-{hi:g}"] = row
        # per-gene consistency across targets (the bias is shared): mean real sign of non-sig cells vs -sign(z)
        gm = (np.abs(z) >= 3) & (M > 5)
        rr["corr_z_vs_mean_real_lfc_nonsig"] = round(float(np.corrcoef(z[M > 5], np.nanmean(np.where(adj & ~sig, rl, np.nan), 0)[M > 5])[0, 1]), 4)
        R[en] = rr
    # analytic fid change: override the emitted sign of called non-sig+sig cells by -sign(z) where |z| >= zc (guide_huber)
    from_fid = {}
    for zc in (2.0, 3.0, 5.0):
        z = E["guide_huber"]
        gm = (np.abs(z) >= zc) & (M > 5)
        sel = call & adj
        def fid(sg):
            k = ((np.sign(sg) == np.sign(rl)) & sel).sum(1)
            den = np.maximum(sel.sum(1), sig.sum(1))
            return float(np.nanmean(np.where(den > 0, k / np.maximum(den, 1), np.nan)))
        f0 = fid(ol)
        f1 = fid(np.where(gm[None, :], -np.sign(z)[None, :], ol))
        from_fid[f"zc{zc}"] = dict(genes=int(gm.sum()), called_cells_changed=int((sel & gm[None, :] & (np.sign(ol) != -np.sign(z)[None, :])).sum()),
                                   dfid=round(f1 - f0, 6), d_overall=round((f1 - f0) * 3.32482 / 6, 6))
    R["fid_override_guide_huber"] = from_fid
    out[surf] = R
    print(surf, json.dumps(R)[:3000], flush=True)
json.dump(out, open(HERE / "work/ntcbias812.json", "w"), indent=1)
print(f"done {time.time()-t0:.0f}s")
