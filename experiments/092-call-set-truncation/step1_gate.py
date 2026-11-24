"""Step 1: validate pairing, build the reference-only CPM gate, cache control block."""
import os, resource, numpy as np, anndata as ad, scipy.sparse as sp

SCRATCH = "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/76a42756-c5fb-4398-aea0-8513243571e9/scratchpad"
REAL = "~/Code/virtual-cell-challenge-2026/data/harness/orion-hct116.h5ad"
PRED = "~/Code/virtual-cell-challenge-2026/experiments/031-somatic/out_som_gamma/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
CTRL = "non-targeting"
PERT_COL = "target_gene"

def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9  # macOS: bytes

def log(msg):
    print(f"[rss {rss():6.2f} GB] {msg}", flush=True)

r = ad.read_h5ad(REAL, backed="r")
p = ad.read_h5ad(PRED, backed="r")
log(f"real {r.shape} pred {p.shape}")
rv = np.asarray(r.var_names, dtype=str); pv = np.asarray(p.var_names, dtype=str)
print("var_names identical & same order:", bool(np.array_equal(rv, pv)))
print("var sets equal:", set(rv) == set(pv))
rl = r.obs[PERT_COL].astype(str).to_numpy(); pl_ = p.obs[PERT_COL].astype(str).to_numpy()
ru, rc = np.unique(rl, return_counts=True); pu, pc = np.unique(pl_, return_counts=True)
print("real groups:", len(ru), "pred groups:", len(pu), "sets equal:", set(ru) == set(pu))
print("control in both:", CTRL in set(ru), CTRL in set(pu))
print("real ctrl n:", rc[ru == CTRL], "pred ctrl n:", pc[pu == CTRL])
print("real pert cells/target: min %d max %d" % (rc[ru != CTRL].min(), rc[ru != CTRL].max()))
print("pred pert cells/target: min %d max %d" % (pc[pu != CTRL].min(), pc[pu != CTRL].max()))
perts = sorted([g for g in ru if g != CTRL])
print("n perturbations:", len(perts))
print("targets that are measured genes:", sum(1 for t in perts if t in set(rv)), "/", len(perts))

# --- control block: the REAL control cells (used as the reference on BOTH sides) ---
cidx = np.flatnonzero(rl == CTRL)
log(f"loading {cidx.size} real control cells")
Xc = r.X[cidx]                      # csr float32
Xc = sp.csr_matrix(Xc)
log(f"control csr nnz={Xc.nnz}")
Xc = Xc.astype(np.float64)          # _to_linear: astype(float)
L = np.zeros(Xc.shape[0], dtype=np.float64)
np.add.reduceat(Xc.data, Xc.indptr[:-1], out=L) if False else None
# exact scanpy _normalize_csr: per-row serial float64 accumulation
for i in range(Xc.shape[0]):
    L[i] = Xc.data[Xc.indptr[i]:Xc.indptr[i+1]].sum()
scale = L / 1e6
scale = scale + (scale == 0)
Xc.data /= np.repeat(scale, np.diff(Xc.indptr))
ctrl_mean_cpm = np.asarray(Xc.mean(axis=0)).ravel()   # float64, matches _grouped_means
keep = ctrl_mean_cpm > 5.0
G = int(keep.sum())
log(f"CPM gate (control-only mean CPM > 5): {G} / {len(rv)} genes kept")
np.savez_compressed(os.path.join(SCRATCH, "gate.npz"),
                    keep=keep, ctrl_mean_cpm=ctrl_mean_cpm, var=rv,
                    ctrl_idx=cidx, perts=np.array(perts, dtype=object).astype(str))
# dense log1p CPM for kept genes, gene-major
sub = Xc[:, keep]              # csr 18400 x G
dense = np.asarray(sub.todense(), dtype=np.float64)   # (18400, G)
np.log1p(dense, out=dense)
log("dense ctrl log1p built")
dense = np.ascontiguousarray(dense.T)                 # (G, 18400)
log("transposed")
np.save(os.path.join(SCRATCH, "ctrl_log_sorted.npy"), np.sort(dense, axis=1))
log("saved sorted ctrl block")
print("peak RSS GB:", rss())
