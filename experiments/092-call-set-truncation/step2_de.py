"""Step 2: per-perturbation DE on both sides, cell_eval2-exact.
Writes de_tables.npz with per-pert (G,) arrays: p_adj/lfc/p_value for real and pred."""
import os, sys, time, resource, numpy as np, h5py
from scipy import stats
from scipy.stats import false_discovery_control

S = "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/76a42756-c5fb-4398-aea0-8513243571e9/scratchpad"
sys.path.insert(0, S)
from de_kernel import ctrl_tiesums, wilcoxon_z

REAL = "~/Code/virtual-cell-challenge-2026/data/harness/orion-hct116.h5ad"
PRED = "~/Code/virtual-cell-challenge-2026/experiments/031-somatic/out_som_gamma/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
CTRL, PERT_COL, EPS = "non-targeting", "target_gene", 1e-9

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
def log(m): print(f"[{time.strftime('%H:%M:%S')} rss {rss():5.2f}G] {m}", flush=True)

gate = np.load(os.path.join(S, "gate.npz"), allow_pickle=True)
keep = gate["keep"]; var = gate["var"].astype(str); perts = gate["perts"].astype(str)
ctrl_mean_cpm_full = gate["ctrl_mean_cpm"]
kept_genes = var[keep]; G = kept_genes.size
ref_mean = ctrl_mean_cpm_full[keep]                     # control mean CPM, kept genes
keep_idx = np.flatnonzero(keep).astype(np.int64)
log(f"G={G} genes, {len(perts)} perturbations")

ctrl_sorted = np.load(os.path.join(S, "ctrl_log_sorted.npy"))
log(f"ctrl block {ctrl_sorted.shape} loaded")
ctrl_tie = ctrl_tiesums(ctrl_sorted)
log("ctrl tie sums done")

def labels_of(path):
    with h5py.File(path, "r") as f:
        o = f["obs"][PERT_COL]
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in o["categories"][:]])
        return cats[o["codes"][:]]

def read_rows(path, rows):
    """CSR rows -> dense (len(rows), n_genes) float64 raw counts. rows must be sorted."""
    with h5py.File(path, "r") as f:
        Xg = f["X"]; indptr = Xg["indptr"]; nvar = int(Xg.attrs["shape"][1])
        lo, hi = int(rows[0]), int(rows[-1]) + 1
        ip = indptr[lo:hi + 1]
        s, e = int(ip[0]), int(ip[-1])
        data = Xg["data"][s:e].astype(np.float64)
        idx = Xg["indices"][s:e]
        ip = ip - s
        out = np.zeros((len(rows), nvar), dtype=np.float64)
        rowmap = {r: i for i, r in enumerate(rows.tolist())}
        for k in range(hi - lo):
            r = lo + k
            i = rowmap.get(r)
            if i is None: continue
            a, b = ip[k], ip[k + 1]
            out[i, idx[a:b]] = data[a:b]
    return out

def side_stats(dense_counts):
    """counts -> (log1p CPM for kept genes as (G,n) float64, mean CPM over kept genes)."""
    L = dense_counts.sum(axis=1)                       # float64 row sums, ALL genes
    scale = L / 1e6
    scale = scale + (scale == 0)
    cpm = dense_counts / scale[:, None]
    sub = cpm[:, keep_idx]                             # (n, G)
    mean_cpm = sub.mean(axis=0)
    lg = np.ascontiguousarray(np.log1p(sub).T)         # (G, n)
    return lg, mean_cpm

rlab = labels_of(REAL); plab = labels_of(PRED)
out = {}
t0 = time.time()
for pi, pert in enumerate(perts):
    rows_r = np.flatnonzero(rlab == pert); rows_p = np.flatnonzero(plab == pert)
    res = {}
    for side, path, rows in (("real", REAL, rows_r), ("pred", PRED, rows_p)):
        dc = read_rows(path, rows)
        lg, mean_cpm = side_stats(dc)
        del dc
        z = wilcoxon_z(ctrl_sorted, ctrl_tie, lg)
        del lg
        pval = 2 * stats.distributions.norm.sf(np.abs(z))
        lfc = np.log2((mean_cpm + EPS) / (ref_mean + EPS))
        padj = false_discovery_control(pval, method="bh")   # per-target BH over gated genes
        padj = np.where(np.isnan(lfc), 1.0, padj)           # nan_lfc_policy='mask'
        res[side] = (pval, padj, lfc)
    out[pert] = res
    if pi % 10 == 0 or pi == len(perts) - 1:
        log(f"{pi+1}/{len(perts)} {pert}  ({(time.time()-t0)/(pi+1):.2f}s/pert)")

np.savez_compressed(
    os.path.join(S, "de_tables.npz"),
    genes=kept_genes, perts=perts,
    real_p=np.stack([out[p]["real"][0] for p in perts]),
    real_padj=np.stack([out[p]["real"][1] for p in perts]),
    real_lfc=np.stack([out[p]["real"][2] for p in perts]),
    pred_p=np.stack([out[p]["pred"][0] for p in perts]),
    pred_padj=np.stack([out[p]["pred"][1] for p in perts]),
    pred_lfc=np.stack([out[p]["pred"][2] for p in perts]),
)
log(f"saved. peak RSS {rss():.2f} GB, total {time.time()-t0:.0f}s")
