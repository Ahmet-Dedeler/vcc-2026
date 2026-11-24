"""Scanpy-normalization-order DE with production's fast rank-sum kernel.

The order counts / (library / 1e6), not counts * (1e6 / library), matters
for ties. Keep this evaluator separate from the production emitter's control
block: changing the latter would change the prediction being evaluated.
"""
import numpy as np
from scipy.stats import false_discovery_control, norm

from vcc2026 import rc_emit


def normalized_counts(block):
    x = np.asarray(block.toarray() if hasattr(block, "toarray") else block,
                   dtype=np.float64)
    scale = x.sum(1) / 1e6
    scale[scale == 0] = 1.
    return x / scale[:, None]


def control_block(matrix, min_cpm=5., chunk=1000):
    n, g = matrix.shape
    mean = np.zeros(g, dtype=np.float64)
    for lo in range(0, n, chunk):
        mean += normalized_counts(matrix[lo:lo+chunk]).sum(0)
    mean /= n
    keep = mean > min_cpm
    cs = np.empty((int(keep.sum()), n), dtype=np.float64)
    for lo in range(0, n, chunk):
        hi = min(lo+chunk, n)
        cs[:, lo:hi] = np.log1p(normalized_counts(matrix[lo:hi])[:, keep]).T
    cs.sort(axis=1)
    return rc_emit.RcControlBlock(
        keep=keep, cs=cs, tie=rc_emit.ctrl_tiesums(cs), ref_lin=mean[keep],
        n_controls=n, min_cpm=min_cpm)


def de_from_counts(block, control):
    linear = normalized_counts(block)[:, control.keep]
    logged = np.ascontiguousarray(np.log1p(linear).T)
    z = rc_emit.wilcoxon_z(control.cs, control.tie, logged)
    p = 2 * norm.sf(np.abs(z))
    q = false_discovery_control(p, method="bh")
    lfc = np.log2((linear.mean(0)+1e-9)/(control.ref_lin+1e-9))
    return z, p, q, lfc
