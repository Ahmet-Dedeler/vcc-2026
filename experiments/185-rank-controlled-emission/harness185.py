"""185 — shared harness loader: the h1-test control block in the scorer's DE space.

Identical conventions to `experiments/171-reach-prefix/realize171.py` (which inherits
164 -> 157, the board-validated `reach` calibration): control mean-per-cell CPM gate at
5 CPM, sorted log1p(CPM) control block, scanpy-identical tie sums.  Cached to `work/` so the
later scripts do not each pay the 18,400-cell pass.
"""
from __future__ import annotations

import os
import sys
import time

import anndata as ad
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from de_kernel import ctrl_tiesums, wilcoxon_z  # noqa: E402,F401

REAL = os.environ.get("REAL", f"{ROOT}/data/harness/h1-test.h5ad")
CTRL, COL, MIN_CPM = "non-targeting", "target_gene", 5.0
N_EMIT, EMIT_L = 400, 100_000


def dense(A, rows):
    m = A.X[rows]
    return np.asarray(m.todense() if hasattr(m, "todense") else m, dtype=np.float64)


def load_controls(force=False):
    """Returns dict(cs, tie, cpm, keep, genes, ref_lin, b_ctrl, ctrl_lib)."""
    meta_p = os.path.join(WORK, "ctrl_meta.npz")
    cs_p = os.path.join(WORK, "ctrl_cs.npy")
    if not force and os.path.exists(meta_p) and os.path.exists(cs_p):
        z = np.load(meta_p, allow_pickle=True)
        cs = np.load(cs_p, mmap_mode="r")
        cs = np.ascontiguousarray(cs)
        return dict(cs=cs, tie=z["tie"], cpm=z["cpm"], keep=z["keep"], genes=z["genes"],
                    ref_lin=z["ref_lin"], b_ctrl=z["b_ctrl"], ctrl_lib=z["ctrl_lib"])
    t0 = time.time()
    R = ad.read_h5ad(REAL)
    rlab = np.asarray(R.obs[COL].astype(str).values, dtype="U64")
    genes = np.asarray([str(g) for g in R.var_names])
    G_full = genes.size
    ci = np.flatnonzero(rlab == CTRL)
    cpm = np.zeros(G_full)
    libs = []
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        libs.append(s.copy())
        s[s == 0] = 1.0
        cpm += (m * (1e6 / s[:, None])).sum(0)
    cpm /= ci.size
    keep = cpm > MIN_CPM
    G = int(keep.sum())
    cs = np.empty((G, ci.size))
    for i in range(0, ci.size, 2000):
        m = dense(R, ci[i:i + 2000])
        s = m.sum(1)
        s[s == 0] = 1.0
        cs[:, i:i + m.shape[0]] = np.log1p(m[:, keep] * (1e6 / s[:, None])).T
    cs.sort(axis=1)
    cs = np.ascontiguousarray(cs)
    tie = ctrl_tiesums(cs)
    out = dict(cs=cs, tie=tie, cpm=cpm, keep=keep, genes=genes, ref_lin=cpm[keep],
               b_ctrl=np.log1p(cpm), ctrl_lib=np.concatenate(libs))
    os.makedirs(WORK, exist_ok=True)
    np.save(cs_p, cs)
    np.savez_compressed(meta_p, tie=tie, cpm=cpm, keep=keep, genes=genes, ref_lin=cpm[keep],
                        b_ctrl=np.log1p(cpm), ctrl_lib=out["ctrl_lib"])
    print(f"control block built and cached: {ci.size} controls, {G}/{G_full} gated "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return out


def load_totals():
    z = np.load(os.path.join(WORK, "a15_totals.npz"), allow_pickle=True)
    return dict(targets=np.asarray(z["targets"]).astype(str), genes=np.asarray(z["genes"]).astype(str),
                totals=z["totals"], cell_stats=z["cell_stats"])


def load_H():
    """169's recovered design-space effect for the h1-test artifact (|H| is the key)."""
    E = np.load(f"{ROOT}/experiments/169-sharpening-six-member/work/h1test-effect.npz",
                allow_pickle=True)
    return dict(H=E["H"].astype(np.float64), targets=np.asarray(E["targets"]).astype(str),
                genes=np.asarray(E["genes"]).astype(str), gate=E["gate"], b=E["b"])


def de_from_block(B, keep, cs, tie, ref_lin):
    """DE table for a (n, G_full) integer block against the cached control block."""
    from scipy.stats import false_discovery_control, norm
    s = B.sum(1).astype(np.float64)
    s[s == 0] = 1.0
    nb = np.log1p(B[:, keep] * (1e6 / s[:, None]))
    z = wilcoxon_z(cs, tie, np.ascontiguousarray(nb.T))
    p = 2 * norm.sf(np.abs(z))
    q = false_discovery_control(p, method="bh")
    lfc = np.log2((np.expm1(nb).mean(0) + 1e-9) / (ref_lin + 1e-9))
    return z, p, q, lfc
