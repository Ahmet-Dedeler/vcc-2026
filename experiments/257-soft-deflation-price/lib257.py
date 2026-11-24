"""257 -- a SOFTER CD4 deflation, priced on the transfer factors the 2026-09-05 board point measured.

Everything is 246's construction, imported and not re-derived (`lib246`: 162's PLOLO, 218's
residuals, the shipped stack, `deflate`, `compose`, `score`).  The ONLY new operator is
`deflate_shrink`:

    R_c = R - mean_rows(R)          per-gene centring over the carried rows, support re-zeroed
    U_k = top-k eigenvectors of R_c R_c^T     (246's left-singular directions, `np.linalg.eigh`)
    U_S = the selected columns of U_k         (all of them, or ONE of them for the per-direction cells)
    R_s = R_c - s * U_S U_S^T R_c             support re-zeroed
then RESOLVED sizing: `solve_fraction` on the deflated object, exactly as the builder does it
(249 section 2, `apply_frozen_cd4_residual`: deflate -> size -> merge -> one preserve).

At s = 1 and S = all top-k this is `lib246.deflate(res, "svd", k)` bit for bit (`sweep257.py`
asserts it on every fold).  Also here: 246's energy-matched random control
(`sweep246b.rand_energy_matched`, verbatim), and 188's exact `nmae` identity plus the delivered
amplitude computed exactly as `six246.py` computes them (frame A, no mass fix, wald89 pool).

Read-only everywhere except `experiments/257-soft-deflation-price/work/` (the PLOLO fold caches
are byte copies of 246's, sha-identical, under `work/plolo/`).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E246 = ROOT / "experiments/246-direction-hypotheses"
E249 = ROOT / "experiments/249-cd4-deflate"
E255 = ROOT / "experiments/255-cd4-deflate-k10"
sys.path.insert(0, str(E246))
import lib246 as L  # noqa: E402  (imports plolo / resid218 / compose218; 246's construction)
import plolo  # noqa: E402

# lib246 pins plolo.WORK to 246's directory at import; every cache read here comes from OUR
# byte copies (sha-identical to 246's, verified in RESULTS.md) so nothing outside this
# experiment is ever touched.
plolo.WORK = WORK / "plolo"
L.WORK = WORK

sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402

CANCER = L.CANCER
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CD4_F = L.CD4_F
RAND_SEED = L.RAND_SEED


# ------------------------------------------------------------------------- the operator
def deflate_shrink(res, k: int, s: float = 1.0, dirs=None):
    """`R_c - s * U_S U_S^T R_c` on the carried rows; S = all top-k (dirs None) or the listed
    0-based indices into the top-k ordering (per-direction cells).  Returns (CD4Prepared, audit)."""
    assert k >= 1 and 0.0 <= s <= 1.0
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    e0 = float((R * R).sum())
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    ec = float((Rc * Rc).sum())
    G = Rc @ Rc.T
    w, U = np.linalg.eigh(G)
    order = np.argsort(w)[::-1]
    Uk = U[:, order[:k]]
    wk = w[order[:k]]
    if dirs is not None:
        Uk = Uk[:, list(dirs)]
    if s == 1.0 and dirs is None:
        Rk = Rc - Uk @ (Uk.T @ Rc)              # lib246.deflate's expression, bit for bit
    else:
        Rk = Rc - s * (Uk @ (Uk.T @ Rc))
    Rk[~sup[rows]] = 0.0
    e1 = float((Rk * Rk).sum())
    tot_w = max(float(w.sum()), 1e-300)
    audit = dict(mode="shrink", k=int(k), s=float(s), dirs=(None if dirs is None else [int(d) for d in dirs]),
                 n_carried=int(rows.size), energy_removed_frac=float(1.0 - e1 / e0),
                 top_k_energy_frac_of_centred=float(wk.sum() / tot_w),
                 per_direction_energy_frac_of_centred=[float(x / tot_w) for x in wk],
                 centred_energy_frac_of_raw=float(ec / e0))
    eff[rows] = Rk
    tag = f"shrink{k}s{s:g}" + ("" if dirs is None else "d" + "".join(str(d) for d in dirs))
    return L._as_res(res, eff, tag, k), audit


def rand_energy_matched(res, target_frac: float, seed: int, kmax: int = 260):
    """sweep246b.rand_energy_matched, verbatim: the smallest random subspace dimension whose
    removed energy fraction is >= target_frac (bisection), 246's `rand` mode and seeding."""
    lo, hi = 1, kmax
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        d, aud = L.deflate(res, "rand", mid, seed=seed)
        if aud["energy_removed_frac"] >= target_frac:
            best = (d, aud, mid)
            hi = mid - 1
        else:
            lo = mid + 1
    if best is None:
        d, aud = L.deflate(res, "rand", kmax, seed=seed)
        best = (d, aud, kmax)
    return best


# ------------------------------------------------------------------ construction + score
def resolved_parts(st: dict, d):
    """The package's object: the deflated CD4 residual re-sized by `solve_fraction` (production's
    all-rows rule at f = 0.34), merged with the shipped union at its shipped fraction."""
    fr, per_unit = L.solve_fraction(st["whitened"], d, CD4_F)
    return [(d, fr), (st["union"], st["frac_x"])], float(fr), float(per_unit)


class Identity:
    """188's exact `nmae` identity on the wald89 pool, `key_precision`, delivered amplitude --
    six246.py's construction on frame A, no mass fix -- for any composed effect matrix."""

    def __init__(self, st: dict):
        fold = st["fold"]
        self.wcols = np.asarray(fold["wcols"], dtype=np.int64)
        self.sc = np.asarray(fold["scored"], dtype=np.int64)
        truth_h = np.asarray(fold["truth"], dtype=np.float64)
        self.bA = plolo.challenge_bctrl("A")
        eA = np.expm1(self.bA)
        self.gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (self.bA > 0)
        cpm = 1e6 * np.expm1(fold["b_ctrl"][self.sc]) / float(np.expm1(fold["b_ctrl"]).sum())
        self.LR = to_lfc(truth_h, fold["b_ctrl"], self.gA)[:, self.sc]
        self.pool = topn(np.abs(self.LR) * np.sqrt(np.maximum(cpm, 0))[None, :], POOL_N)
        self.Ts = truth_h[:, self.sc]

    def measure(self, final: np.ndarray) -> dict:
        H = plolo.deliver(final, self.bA, self.gA, None)
        LP = to_lfc(H, self.bA, self.gA)[:, self.sc]
        U, M, D = umd(LP, self.LR, self.pool)
        prec = key_precision(np.abs(H[:, self.sc]), H[:, self.sc], self.Ts, n=POOL_N)
        dn = float(np.median(np.linalg.norm(H[:, self.wcols], axis=1)))
        return dict(U=U, M=M, D=D, nmae_wald89=float(np.mean(1 + (U - 2 * M) / D)),
                    top89_precision=float(prec), delivered_norm=dn)
