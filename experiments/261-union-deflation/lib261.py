"""261 -- deflate the H1-UNION residual's target-shared directions BEFORE its carrying-row sizing, on the
`cd4svd5` base.  MEASUREMENT ONLY.

Everything is 246/257's construction, imported and not re-derived: `lib257` -> `lib246` (162's PLOLO, 218's
residuals and `union_rows`, the shipped stack, `deflate`, `compose`, `score`), `lib257.Identity` (188's exact
`nmae` identity, six246's construction), `lib257.deflate_shrink` (the 257 known-answer cell) and
`lib257.rand_energy_matched` (246's energy-matched random control).  Nothing mathematical is new here:

  * the base is 257's `k5_s1` == 246's `B2_cd4_svd5_resolved` == the shipped `cd4svd5`
    (CD4 deflated at k = 5 by `lib246.deflate(cd4, "svd", 5)`, re-sized by `solve_fraction` at f = 0.34;
    the H1-pair union at its carrying-row fraction for g = 0.28);
  * the arm applies the SAME `lib246.deflate(res, "svd", k)` to the UNION object (centred over ITS carrying
    rows, support re-zeroed) and then re-sizes the deflated union by 218's carrying-row rule
    (`compose218.solve_fraction_carrying`, the builder's `solve_carrying_row_fraction` verbatim) at g --
    deflate BEFORE sizing, exactly the order 249 used for CD4 (the RESOLVED variant);
  * the two projector controls are 246's verbatim (`sweep246b.shuffled_projector`, `sweep246c.split_half_deflate`),
    copied here because those modules run a sweep when imported.

Read-only everywhere except `experiments/261-union-deflation/work/`; the PLOLO fold caches under `work/plolo/`
are byte copies of 257's (= 249's = 246's), sha-listed in `work/plolo/SHA256SUMS.txt`.
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
E257 = ROOT / "experiments/257-soft-deflation-price"
E258 = ROOT / "experiments/258-union-dose-on-cd4svd5"
E259 = ROOT / "experiments/259-union-without-igvf"
sys.path.insert(0, str(E257))
import lib257 as X  # noqa: E402  (imports lib246 as X.L, which imports plolo / resid218 / compose218)
import plolo  # noqa: E402
from compose218 import solve_fraction_carrying  # noqa: E402  (218 section 7.1; == the builder's rule)

L = X.L
# lib246 pins plolo.WORK to 246's directory and lib257 re-pins it to 257's; pin it to OURS so every cache
# read comes from this experiment's byte copies and nothing outside it is ever touched.
plolo.WORK = WORK / "plolo"
L.WORK = WORK
X.WORK = WORK

CANCER = L.CANCER
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
CD4_F = L.CD4_F           # 0.34
H1_G = L.H1_G             # 0.28
RAND_SEED = L.RAND_SEED   # 20260904
SHUF_SEED = L.SHUF_SEED   # 20726


# ------------------------------------------------------------------------------- the base
def cd4svd5_parts(st: dict):
    """The incumbent's object on this fold: 257's `k5_s1` (CD4 deflated k = 5, RESOLVED sizing at f = 0.34)
    merged with the shipped union at its shipped carrying-row fraction (g = 0.28)."""
    d, aud = L.deflate(st["cd4"], "svd", 5)
    fr, pu = L.solve_fraction(st["whitened"], d, CD4_F)
    return d, float(fr), float(pu), aud


def union_parts(st: dict, d_cd4, fr_cd4: float, u, g: float, *, sizing: str = "resolved"):
    """The arm's object: the (deflated) union `u` sized by the carrying-row rule at g against the SAME whitened
    base (the builder's order: deflate -> size -> merge -> one preserve), merged with the deflated CD4 at ITS
    resolved fraction.  `sizing="fixed"` keeps the shipped union fraction (the projector alone, 246's 'fixed')."""
    if sizing == "resolved":
        fu, pu = solve_fraction_carrying(st["whitened"], u, g)
    elif sizing == "fixed":
        fu, pu = st["frac_x"], float("nan")
    else:
        raise ValueError(sizing)
    assert fu > 0.0, "vacuous union"
    return [(d_cd4, fr_cd4), (u, float(fu))], float(fu), float(pu)


# --------------------------------------------------------------------------- the controls
def shuffled_projector(res, k: int, seed: int):
    """sweep246b.shuffled_projector, verbatim: fit the top-k projector on a row-DERANGED copy of the centred
    carried-row matrix, apply it to the REAL row order."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    rng = np.random.default_rng(seed)
    perm = np.arange(rows.size)
    while True:
        rng.shuffle(perm)
        if not np.any(perm == np.arange(rows.size)):
            break
    G = Rc[perm] @ Rc[perm].T
    w, U = np.linalg.eigh(G)
    Uk = U[:, np.argsort(w)[::-1][:k]]          # fitted on shuffled row order
    Rk = Rc - Uk @ (Uk.T @ Rc)                  # applied to the REAL row order
    Rk[~sup[rows]] = 0.0
    e0 = float((R * R).sum())
    eff[rows] = Rk
    audit = dict(mode="svd_shufdef", k=int(k), n_carried=int(rows.size), seed=int(seed),
                 energy_removed_frac=float(1.0 - (Rk * Rk).sum() / e0))
    return L._as_res(res, eff, "svd_shufdef", k), audit


def split_half_deflate(res, k: int, seed: int):
    """sweep246c.split_half_deflate, verbatim: the top-k RIGHT-singular basis fitted on a random half of the
    carried rows, applied to the other half (both ways), so no row is deflated by a basis fitted on itself."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(rows.size)
    halves = (perm[: rows.size // 2], perm[rows.size // 2:])
    e0 = float((R * R).sum())
    out = R.copy()
    for fit, apply in ((halves[0], halves[1]), (halves[1], halves[0])):
        Rf = R[fit]
        mean = Rf.mean(axis=0, keepdims=True)
        Rfc = Rf - mean
        Rfc[~sup[rows][fit]] = 0.0
        G = Rfc @ Rfc.T
        w, U = np.linalg.eigh(G)
        idx = np.argsort(w)[::-1][:k]
        Uk, sk = U[:, idx], np.sqrt(np.maximum(w[idx], 1e-300))
        Vk = (Rfc.T @ Uk) / sk[None, :]                     # (G, k), orthonormal columns
        Ra = R[apply] - mean
        Ra[~sup[rows][apply]] = 0.0
        Ra = Ra - (Ra @ Vk) @ Vk.T
        Ra[~sup[rows][apply]] = 0.0
        out[apply] = Ra
    eff[rows] = out
    audit = dict(mode="svd_splithalf", k=int(k), n_carried=int(rows.size), seed=int(seed),
                 half_sizes=[int(halves[0].size), int(halves[1].size)],
                 energy_removed_frac=float(1.0 - (out * out).sum() / e0))
    return L._as_res(res, eff, "svd_splithalf", k), audit


def spectrum(res, n: int = 12) -> dict:
    """The centred carried-row Gram's eigenvalue shares (what `deflate("svd", k)` removes, direction by
    direction), plus the rank-1 mean's share of the raw energy."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    sup = np.asarray(res.gene_votes) > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    e0 = float((R * R).sum())
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rows]] = 0.0
    ec = float((Rc * Rc).sum())
    w = np.sort(np.linalg.eigvalsh(Rc @ Rc.T))[::-1]
    tot = max(float(w.sum()), 1e-300)
    return dict(n_carried=int(rows.size), raw_energy=e0, centred_energy_frac_of_raw=ec / e0,
                top_shares_of_centred=[float(x / tot) for x in w[:n]],
                cum_shares_of_centred=[float(x) for x in np.cumsum(w[:n]) / tot],
                row_norm_median=float(np.median(np.linalg.norm(np.where(sup[rows], R, 0.0), axis=1))))
