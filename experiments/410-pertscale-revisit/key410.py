"""410 — the whitened shared-energy key, clean construction, shared by the PLOLO check and the panel artifact.

    key_p = mean_{i<j} 0.5 * ( <W d_i, d_j> + <W d_j, d_i> )  /  || W m_p ||^2 ,   m_p = mean_i d_i

over the donors carrying target p (>= 2), W = V diag(lambda^(-gamma/2)) V^T on the operator axis (NO on-axis
preserve constant: the key is then exactly homogeneous of degree -1 in W, so a global constant on W rescales every
key identically and the median-relative weight is invariant). d_i are the RAW donor delta rows (bank delta.npy)
mapped to the official axis; inner products run over the operator axis for W-rows and the donor's own coverage for
raw rows (zeros elsewhere). This is the plug-in estimate of the mse-optimal per-target scale s_p = <P_p,T_p>/|P_p|^2
with a held donor standing in for the truth: <W m, T> ~ <W d_i, d_j> for i != j removes each donor's own noise.

weights(key, lo, hi): w = clip(key / median(key over carried), lo, hi); targets carried by < 2 donors get 1.0
(the median, "neutral"). Energy normalisation happens where the effect matrix is known (builder / check).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402
from size410 import donor_rows  # noqa: E402

GAMMA = 0.65


def raw_whiten(op: "plolo.Operator", M: np.ndarray, wcols: np.ndarray, gamma: float = GAMMA) -> np.ndarray:
    """W applied on the operator axis, zeros elsewhere, no preserve constant."""
    op._check_axis(wcols)
    w = op.evals ** (-gamma / 2.0)
    out = np.zeros_like(M, dtype=np.float64)
    out[:, wcols] = ((np.asarray(M, dtype=np.float64)[:, wcols] @ op.evecs) * w) @ op.evecs.T
    return out


def shared_energy_key(targets, donors, op, wcols, cols=None):
    """Return (key, carried, n_donors) for `targets`. `cols` restricts every inner product (None = all)."""
    rows, have = [], []
    for name in donors:
        r, h = donor_rows(name, targets)
        rows.append(r)
        have.append(h)
    H = np.vstack(have)
    Wr = [raw_whiten(op, r, wcols) for r in rows]
    if cols is not None:
        rows = [r[:, cols] for r in rows]
        Wr = [w[:, cols] for w in Wr]
    n = len(targets)
    key = np.full(n, np.nan)
    for p in range(n):
        idx = [k for k in range(len(rows)) if H[k, p]]
        if len(idx) < 2:
            continue
        ips = []
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                ips.append(0.5 * (float(Wr[i][p] @ rows[j][p]) + float(Wr[j][p] @ rows[i][p])))
        m = np.mean([Wr[k][p] for k in idx], axis=0)
        key[p] = float(np.mean(ips)) / max(float(m @ m), 1e-300)
    return key, H.sum(0) >= 2, H.sum(0)


def weights(key, carried, lo, hi):
    med = float(np.median(key[carried]))
    w = np.ones(len(key))
    w[carried] = np.clip(key[carried] / med, lo, hi)
    return w, med
