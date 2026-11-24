"""Cache the REFERENCE DE table for a harness: gate mask, control mean-CPM, and per
perturbation the reference significant set + reference lfc over the gated genes.

This is the object `de_wilcoxon_lfc_nmae` reads: S_p is chosen by the reference ALONE
(BH q<0.05, finite lfc, target gene removed, |S_p| >= 10), so it is identical for every
submission.  Everything downstream in 098 reuses this cache.

    uv run python ref_de.py <h1|orion>
"""
from __future__ import annotations
import os, sys, time
import numpy as np
from scipy.stats import false_discovery_control, norm
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
import anndata as ad
from de_kernel import ctrl_tiesums, wilcoxon_z

WHICH = sys.argv[1]
REAL = {"h1": f"{ROOT}/data/harness/h1-test.h5ad",
        "orion": f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"}[WHICH]
OUT = f"{ROOT}/experiments/098-nmae/work/refde_{WHICH}.npz"
CTRL, MIN_CPM, ALPHA_DE = "non-targeting", 5.0, 0.05

A = ad.read_h5ad(REAL)
lab = np.asarray(A.obs["target_gene"].astype(str).values)
genes = np.asarray([str(g) for g in A.var_names])
X = A.X
X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float64)
print(f"[{WHICH}] cells {X.shape[0]} genes {X.shape[1]} perts {np.unique(lab).size}", flush=True)

def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))

cidx = np.flatnonzero(lab == CTRL)
Cn = norm_log(X[cidx])
cpm = np.expm1(Cn).mean(0)
gate = cpm > MIN_CPM
cs = np.ascontiguousarray(np.sort(Cn[:, gate].T, axis=1)); tie = ctrl_tiesums(cs)
ref_lin = np.expm1(Cn[:, gate]).mean(0)
del Cn
Gk = int(gate.sum()); gk = genes[gate]
gpos = {g: i for i, g in enumerate(gk)}
print(f"[{WHICH}] {cidx.size} controls, {Gk}/{genes.size} genes past {MIN_CPM} CPM", flush=True)

perts = sorted(set(lab) - {CTRL})
names, sig_idx, sig_ptr, LFC, NCELL, QV = [], [], [0], [], [], []
t0 = time.time()
for i, p in enumerate(perts):
    r = np.flatnonzero(lab == p)
    B = norm_log(X[r])
    z = wilcoxon_z(cs, tie, np.ascontiguousarray(B[:, gate].T))
    q = false_discovery_control(2 * norm.sf(np.abs(z)), method="bh")
    lf = np.log2((np.expm1(B[:, gate]).mean(0) + 1e-9) / (ref_lin + 1e-9)); del B
    sig = (q < ALPHA_DE) & np.isfinite(lf) & (lf != 0)
    j = gpos.get(p)
    if j is not None: sig[j] = False
    gs = np.flatnonzero(sig)
    names.append(p); NCELL.append(r.size); LFC.append(lf.astype(np.float32))
    QV.append(q.astype(np.float32))
    sig_idx.append(gs.astype(np.int32)); sig_ptr.append(sig_ptr[-1] + gs.size)
    if (i + 1) % 20 == 0: print(f"  {i+1}/{len(perts)} [{time.time()-t0:.0f}s]", flush=True)

np.savez_compressed(OUT, genes=genes, gate=gate, ref_lin=ref_lin, gk=gk,
                    perts=np.array(names), ncell=np.array(NCELL),
                    lfc=np.vstack(LFC), sig_idx=np.concatenate(sig_idx),
                    sig_ptr=np.array(sig_ptr), n_ctrl=cidx.size, q=np.vstack(QV))
sz = np.diff(sig_ptr)
print(f"[{WHICH}] gate sizes: median {np.median(sz):.0f}  mean {sz.mean():.0f} "
      f"min {sz.min()} max {sz.max()};  rows with >=10: {(sz>=10).sum()}/{len(sz)}")
print(f"[{WHICH}] wrote {OUT}")
