"""190 — the network library: every prior reduced to one signed weighted matrix N.

    N[t, g]   sign = direction the NETWORK expects gene g to move when target t is KNOCKED DOWN
              |.|  = evidence weight
              0    = the network is silent about this pair

Four families:

  ct-direct     CollecTRI, DIRECT edges only.  Requires the knocked-down gene to itself be a
                CollecTRI TF.  Expected sign = -mor.
  ct-prop       CollecTRI, PROPAGATED through predicted TF expression change.  This is the
                shape `tydymy` published ("signed TF to target correction"): read our own
                predicted change in each TF's transcript, push its regulon by mor times that
                change.  Covers every target, TF or not, and every gene in any regulon.
                    Prop[t, g] = sum_f  a[t, f] * mor(f, g) * w(f, g)
                with a[t, f] the per-target z-scored predicted effect on TF f, optionally
                truncated to the top-K TFs by |a|.
  coexpr        control-derived co-expression GRN, D2R2-shaped.  N = -r(t, g).
                174's known-answer rows: sign accuracy on the pool 0.5158 for -sign(r),
                0.4842 for +sign(r), |r| median 0.0256 at pc0.
  expr          the cheap comparator 189 established: no pair structure at all, just the
                gene's control expression.  Not signed, so it enters only as a key, never
                as an agreement score.

Every builder returns (N, info).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)


# --------------------------------------------------------------------------------------
def _regulon_arrays(targets, genes, weighted=True):
    """CollecTRI as (tf_index_on_gene_axis, gene_index, signed_weight) plus the tf->row map."""
    import prior190
    d = prior190.load_collectri(weighted=weighted)
    gi = {g: j for j, g in enumerate(genes)}
    ti = {t: i for i, t in enumerate(targets)}
    d = d[d["gene"].isin(gi)]
    return d, gi, ti


def ct_direct(targets, genes, H, weighted=True, rng=None, shuffle=None):
    """CollecTRI direct edges.  N[t, g] = -mor * w for t a TF.

    shuffle: None | 'targets' (permute each TF's target list within the gene axis, preserving
    that TF's out-degree and its sign multiset) | 'signs' (permute the mor values globally,
    preserving topology) | 'tfs' (reassign whole regulons to other covered targets).
    """
    d, gi, ti = _regulon_arrays(targets, genes, weighted)
    hit = d[d["tf"].isin(ti)].copy()
    if shuffle == "signs":
        hit["mor"] = rng.permutation(hit["mor"].to_numpy())
    N = np.zeros((len(targets), len(genes)), dtype=np.float64)
    if shuffle == "targets":
        # each TF keeps its out-degree and its (mor, w) multiset; the genes are redrawn
        for tf, grp in hit.groupby("tf"):
            cols = rng.choice(len(genes), size=len(grp), replace=False)
            N[ti[tf], cols] = -grp["mor"].to_numpy() * grp["w"].to_numpy()
    else:
        rows = np.array([ti[t] for t in hit["tf"]])
        cols = np.array([gi[g] for g in hit["gene"]])
        N[rows, cols] = -hit["mor"].to_numpy() * hit["w"].to_numpy()
    if shuffle == "tfs":
        cov = np.flatnonzero((N != 0).any(1))
        N[cov] = N[rng.permutation(cov)]
    return N, {"kind": "ct-direct", "shuffle": shuffle}


def ct_prop(targets, genes, H, topk=None, weighted=True, rng=None, shuffle=None,
            self_only=False):
    """CollecTRI propagated through the predicted TF expression change (tydymy's shape)."""
    d, gi, ti = _regulon_arrays(targets, genes, weighted)
    tfs = sorted(set(d["tf"]) & set(gi))            # TFs measurable on our gene axis
    if not tfs:
        return np.zeros((len(targets), len(genes))), {"kind": "ct-prop", "n_tfs": 0}
    tf_col = np.array([gi[f] for f in tfs])
    fpos = {f: i for i, f in enumerate(tfs)}
    sub = d[d["tf"].isin(fpos)]
    if shuffle == "signs":
        sub = sub.copy()
        sub["mor"] = rng.permutation(sub["mor"].to_numpy())
    R = np.zeros((len(tfs), len(genes)), dtype=np.float64)   # regulon operator
    if shuffle == "targets":
        for f, grp in sub.groupby("tf"):
            cols = rng.choice(len(genes), size=len(grp), replace=False)
            R[fpos[f], cols] = grp["mor"].to_numpy() * grp["w"].to_numpy()
    else:
        R[[fpos[f] for f in sub["tf"]], [gi[g] for g in sub["gene"]]] = (
            sub["mor"].to_numpy() * sub["w"].to_numpy())
    if shuffle == "tfs":
        R = R[rng.permutation(len(tfs))]
    # predicted TF activity change: per-target z-score of H over the TF columns
    A = H[:, tf_col].astype(np.float64)
    if self_only:
        # only the knocked-down gene's own TF entry may act -- isolates the direct channel
        M = np.zeros_like(A)
        for i, t in enumerate(targets):
            if t in fpos:
                M[i, fpos[t]] = A[i, fpos[t]]
        A = M
    sd = A.std(1, keepdims=True)
    A = A / np.where(sd == 0, 1.0, sd)
    if topk is not None and not self_only:
        keep = np.zeros_like(A, dtype=bool)
        idx = np.argsort(-np.abs(A), axis=1)[:, :topk]
        np.put_along_axis(keep, idx, True, axis=1)
        A = np.where(keep, A, 0.0)
    N = A @ R
    return N, {"kind": "ct-prop", "n_tfs": len(tfs), "topk": topk, "shuffle": shuffle}


def coexpr(targets, genes, H, pc=0, topk=None, rng=None, shuffle=None):
    """N = -r(t, g) from the recipient's control cells (174's matrix, sign corrected)."""
    Z = np.load(os.path.join(WORK, "coex190.npz"), allow_pickle=True)
    C = Z[f"C_pc{pc}"].astype(np.float64)
    N = -C                                    # knockdown of an activator lowers its partners
    if topk is not None:
        M = np.zeros(N.shape, bool)
        for p in range(N.shape[0]):
            M[p, np.argsort(-np.abs(N[p]))[:topk]] = True
        N = np.where(M, N, 0.0)
    if shuffle == "targets":
        out = np.zeros_like(N)
        for p in range(N.shape[0]):
            nz = np.flatnonzero(N[p])
            cols = rng.choice(N.shape[1], size=nz.size, replace=False)
            out[p, cols] = N[p, nz]
        N = out
    elif shuffle == "signs":
        s = np.sign(N)
        flat = s.ravel().copy()
        nz = np.flatnonzero(flat)
        flat[nz] = rng.permutation(flat[nz])
        N = np.abs(N) * flat.reshape(N.shape)
    elif shuffle == "tfs":
        N = N[rng.permutation(N.shape[0])]
    return N, {"kind": f"coexpr-pc{pc}", "topk": topk, "shuffle": shuffle}
