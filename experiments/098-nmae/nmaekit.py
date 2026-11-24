"""Shared machinery for experiment 098: analytic `de_wilcoxon_lfc_nmae` on a harness.

The metric reads only two things from the submission: the predicted mean expression profile
per perturbation, and nothing else.  Our emitted cells are (quasi) multinomial draws from a
profile with exact totals, so the mean per-cell CPM converges to `1e6 * prof`, which makes the
predicted lfc a DETERMINISTIC function of the profile:

    lfc_hat[p,g] = log2( (1e6 * prof[p,g] + 1e-9) / (ctrl_mean_cpm[g] + 1e-9) )
    prof[p]      = expm1(bctrl + eff[p]) / sum_g expm1(bctrl + eff[p])

`validate.py` checks that against the realised (emitted, rank-summed) lfc.  Everything else in
098 is a transformation of `eff`.
"""
from __future__ import annotations
import os, sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms                                    # noqa: E402
from msekit import CTRL, TS, group_stats                       # noqa: E402
from cell_eval2.prep import bulk_lognorm_means                 # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget       # noqa: E402

REAL = {"h1": f"{ROOT}/data/harness/h1-test.h5ad",
        "orion": f"{ROOT}/data/harness/orion-hct116-thinned.h5ad"}
PRED = {"h1": f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
              "g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        "orion": f"{ROOT}/experiments/095-clean-harness/out/runs/"
                 "g0.65_q1.0_a1.08_spread0.25/pred.h5ad"}
MIN_GATE = 10


class Harness:
    """Everything 098 needs for one harness, loaded once."""

    def __init__(self, which: str):
        self.which = which
        z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{which}.npz", allow_pickle=True)
        self.genes = z["genes"].astype(str)
        self.gate = z["gate"]                       # [G] bool, control CPM > 5
        self.ref_lin = z["ref_lin"]                 # [Gk] control mean CPM
        self.gk = z["gk"].astype(str)
        self.rperts = z["perts"].astype(str)
        self.rlfc = z["lfc"].astype(np.float64)     # [P_real, Gk] reference lfc
        si, sp = z["sig_idx"], z["sig_ptr"]
        self.sig = [si[sp[i]:sp[i + 1]] for i in range(len(self.rperts))]

        rp, rs, rc, rj, rg = group_stats(REAL[which])
        pp, ps, pc, pj, pg = group_stats(PRED[which])
        t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
        for k, c in t["cols"].items():
            t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
        self.dh0 = t["dh"]                          # [P, G] predicted effect, bulk_lognorm space
        self.dreal = t["d"]                         # [P, G] reference effect, same space
        self.bctrl = t["bctrl"]                     # [G] harness control profile, same space
        self.cols = t["cols"]; self.w = t["w"]; self.claim = t["claim"]
        self.tg = np.array([str(x) for x in pp[t["keep"]]])
        self.P, self.G = self.dh0.shape
        assert np.array_equal(rg, self.genes), "gene order mismatch"

        # rows of dh0 that have a reference DE row with a usable gate
        rpos = {p: i for i, p in enumerate(self.rperts)}
        self.rows, self.rrows = [], []
        for i, p in enumerate(self.tg):
            j = rpos.get(p)
            if j is None or self.sig[j].size < MIN_GATE:
                continue
            self.rows.append(i); self.rrows.append(j)
        self.rows = np.asarray(self.rows); self.rrows = np.asarray(self.rrows)
        self.gpos = {g: i for i, g in enumerate(self.gk)}

    # ---- profile -> predicted lfc -----------------------------------------------------
    def lfc_hat(self, eff_rows: np.ndarray, rows: np.ndarray | None = None) -> np.ndarray:
        """`eff_rows` is [n, G] in bulk_lognorm space; returns [n, Gk] predicted lfc."""
        prof = np.expm1(self.bctrl[None, :] + eff_rows)
        np.clip(prof, 0.0, None, out=prof)
        prof /= prof.sum(1, keepdims=True)
        return np.log2((1e6 * prof[:, self.gate] + 1e-9) / (self.ref_lin[None, :] + 1e-9))

    def nmae(self, eff: np.ndarray, per_row: bool = False):
        """`eff` is [P, G] (all builder rows).  Returns mean NMAE over scored rows."""
        L = self.lfc_hat(eff[self.rows])
        out = np.empty(self.rows.size)
        for k, (i, j) in enumerate(zip(self.rows, self.rrows)):
            gs = self.sig[j]
            r = self.rlfc[j][gs]
            out[k] = np.abs(L[k][gs] - r).sum() / np.abs(r).sum()
        return (out, self.rrows) if per_row else float(out.mean())


def a12_eff(H: Harness) -> np.ndarray:
    """a12's effect matrix: no shrinkage, at its own solved cliff amplitude."""
    W = np.ones_like(H.dh0)
    B = float(_across_pert_budget(W * H.dh0, H.cols, H.w))
    s = float(np.sqrt(H.claim / B))
    return s * W * H.dh0, s


# ======================================================================================
#  Exact evaluators on an arbitrary predicted-effect matrix `dh` [P, G]
# ======================================================================================
def budget(H: Harness, dh: np.ndarray) -> float:
    return float(_across_pert_budget(dh, H.cols, H.w))


def AXR(H: Harness, dh: np.ndarray):
    """`(A, X)` of 096's identity, with each row's own target coordinate dropped."""
    d = H.dreal
    a = np.einsum("ij,ij->i", dh, dh)
    x = np.einsum("ij,ij->i", dh, d)
    for k, c in H.cols.items():
        a[k] -= dh[k, c] ** 2
        x[k] -= dh[k, c] * d[k, c]
    return float((a * H.w).sum()), float((x * H.w).sum())


def pds(H: Harness, dh: np.ndarray) -> float:
    """`pds_cosine` of the spec: rank of the correct match by cosine distance, with EVERY
    panel target gene's coordinate removed from both sides."""
    tgt = np.asarray([H.cols[k] for k in sorted(H.cols)], int)
    keep = np.ones(H.G, bool); keep[tgt] = False
    P_ = dh[:, keep]; R_ = H.dreal[:, keep]
    pn = np.linalg.norm(P_, axis=1); rn = np.linalg.norm(R_, axis=1)
    ok = (pn > 0) & (rn > 0)
    C = (P_ @ R_.T) / np.maximum(pn[:, None] * rn[None, :], 1e-300)
    C[~ok, :] = -1.0; C[:, ~ok] = -1.0          # zero norm -> distance 1 -> cosine -1
    n = C.shape[0]
    own = np.diag(C).copy()
    lt = (C > own[:, None]).sum(1)
    eq = (C == own[:, None]).sum(1)
    k = lt + 0.5 * (eq - 1)
    return float(np.mean(1.0 - k / (n - 1)))


def nmae_of(H: Harness, dh: np.ndarray, per_row: bool = False):
    L = H.lfc_hat(dh[H.rows])
    out = np.empty(H.rows.size)
    for k, j in enumerate(H.rrows):
        gs = H.sig[j]; r = H.rlfc[j][gs]
        out[k] = np.abs(L[k][gs] - r).sum() / np.abs(r).sum()
    return out if per_row else float(out.mean())


def mse_cliff(H: Harness, dh: np.ndarray, claim: float, Kc: float, sumD: float):
    """`min_s MSE(s)` at the rho cliff."""
    A, X = AXR(H, dh); B = budget(H, dh)
    s = float(np.sqrt(claim / B))
    return (1.0 + (Kc + claim * (A - B) / B - 2 * s * X) / sumD), A, X, B, s
