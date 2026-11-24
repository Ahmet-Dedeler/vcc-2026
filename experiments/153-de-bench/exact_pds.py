"""Exact scorer-space PDS: effects are DIFFERENCES OF log1p-PSEUDOBULK, not log2 fold changes.

BRIEF section 1:  b_{p,g} = log(1 + TS * P_{p,g} / sum_g' P_{p,g'}),  TS = 5e4,
P_{p,g} summed over the group's cells.  delta = b_p - b_ctrl, cosine, all 300 panel target
genes removed, PDS_p = 1 - k_p/(n-1) with ties sharing the average position.

This matters here specifically: log1p COMPRESSES low-count genes, while a log2 ratio does not.
The lfc-space proxy therefore over-weights exactly the low-expression coordinates that the
whitening ridge acts on, so it is biased on the axis under test.  Hence this file.
"""
import numpy as np, anndata as ad, scipy.sparse as sp, csv

TS = 5e4

def real_effects(h5ad, ctrl="non-targeting", col="target_gene"):
    A = ad.read_h5ad(h5ad)
    lab = np.asarray(A.obs[col].astype(str).values)
    X = A.X.tocsr()
    genes = np.asarray([str(g) for g in A.var_names])
    def prof(rows):
        P = np.asarray(X[rows].sum(0)).ravel().astype(np.float64)
        return np.log1p(TS * P / max(P.sum(), 1.0))
    b_ctrl = prof(np.flatnonzero(lab == ctrl))
    perts = np.array(sorted(set(lab) - {ctrl}))
    D = np.stack([prof(np.flatnonzero(lab == t)) - b_ctrl for t in perts])
    return genes, perts, D, b_ctrl

def pds(pred, real):
    """pred, real: (n, G) effect matrices on the same axis, panel genes already removed."""
    A = pred / np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), 1e-12)
    B = real / np.maximum(np.linalg.norm(real, axis=1, keepdims=True), 1e-12)
    S = A @ B.T                      # cosine similarity; distance = 1 - this
    own = np.diag(S).copy()
    closer = (S > own[:, None]).sum(1)                    # smaller distance = larger cosine
    tied = (S == own[:, None]).sum(1) - 1
    n = S.shape[0]
    return float(np.mean(1.0 - (closer + 0.5 * tied) / (n - 1)))

if __name__ == "__main__":
    import sys
    genes, perts, D, _ = real_effects(sys.argv[1])
    np.savez_compressed(sys.argv[2], genes=genes, perts=perts, delta=D.astype(np.float32))
    print(f"{sys.argv[1].split('/')[-1]}: {perts.size} perturbations x {genes.size} genes -> {sys.argv[2]}")
