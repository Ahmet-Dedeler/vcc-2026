"""250 — shared construction: 233's shipped stack on any of the five panels, 246's deflation
VERBATIM, 246's controls verbatim, the pinned A/B/C gene-space basis, CD4's partial rows per
panel, and 246's scorer (per-target arrays kept for leave-top-3-out).

Nothing here is re-derived: `fold233.consensus` builds the consensus, `plolo.Operator` applies
the panel's own re-fitted operator (axis-asserted by `fold233.verify_op`), `score233.cd4_prepared`
builds the panel's CD4 residual from the panel's own re-extracted banks, `compose218` sizes and
merges, `E062.add_post_whiten_residual` adds with one global preserve, and
`resid218.score_effects` renders (a23 minnorm + a15) and reads back on the three challenge
frames.  The deflation functions are byte-copies of `lib246.deflate` / `sweep246b` / `sweep246c`
(asserted against 246's json rows by `bridge250.py` and against 249's shipped object by
`basis250.py`).

Sizing, the October distinction (247 §3):
  allrows   production A/B/C: `solve_cd4_fraction` = all-rows median at f = 0.34
  carrying  D/E/F (def_pipeline DEF_ENV): carrying-row median at g_cd4 = 0.3443655278649641
`resolved` = the fraction re-solved on the DEFLATED object (249's production construction:
deflate BEFORE `solve_cd4_dose`); `fixed` = the un-deflated stack's fraction (246's B2 rows).
"""
from __future__ import annotations

import hashlib
import json
import os
import resource
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

import bind250 as B  # noqa: E402  MUST precede score233 / fold233
import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, Operator, solve_fraction  # noqa: E402

resid218, compose218 = S._resid_module()        # pins plolo.WORK to 250's shared cache dir
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

HERE = B.HERE
ROOT = B.ROOT
WORK = B.WORK
OFFICIAL = plolo.OFFICIAL
ALPHA, GAMMA, LAM = S.ALPHA, S.GAMMA, S.LAM
CD4_F, UNION_G, UNION = S.CD4_F, S.UNION_G, S.UNION
G_CD4 = 0.3443655278649641          # 247 §3.2: 0.34 x median_carrying/median_all on A/B/C
RAND_SEED = 20260904                # 246's
SHUF_SEED = 20726                   # 218's / 246's derangement seed
NULL_SEED = 233233                  # 233's deranged-null seed
PROD_CD4_ABC = ROOT / "experiments/249-cd4-deflate/work/stage/cd4-off.npz"   # the shipped object
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")

assert plolo.WORK == S.SHARED and S.SHARED == WORK / "plolo" / "shared", \
    f"cache dir escaped 250: {plolo.WORK}"


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_arr(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


# ------------------------------------------------------------------------------ folds
def load_fold(panel: str, held: str) -> dict:
    """`fold233.build`'s cache-hit branch, reading 233's / 238's cache read-only."""
    cache = B.fold_cache(panel, held)
    assert cache.exists(), f"no fold cache {cache}"
    sources = B.FOLD_SOURCES[held]
    panel_t = B.read_panel(panel)
    panel_set = set(panel_t)
    banks = {n: E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in sources}
    held_bank = E062.Bank.load(held, plolo.DONOR_PATHS[held])
    targets = tuple(g for g in panel_t if g in held_bank.pert_pos)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in panel_t if g in official_pos],
                                 dtype=np.int64)
    donor_union: set[str] = set()
    for n in sources:
        donor_union |= set(banks[n].genes)
    wcols = np.asarray([i for i, g in enumerate(OFFICIAL)
                        if g in donor_union and g not in panel_set], dtype=np.int64)
    z = np.load(cache, allow_pickle=True)
    assert tuple(str(t) for t in z["targets"]) == targets, f"{cache}: target list drifted"
    out = {k: z[k] for k in z.files if k != "audit"}
    out["audit"] = json.loads(str(z["audit"]))
    out.update(panel=panel, held=held, targets=targets, sources=sources,
               panel_positions=panel_positions, wcols=wcols, op_path=B.op_path(panel, held),
               cache=str(cache.relative_to(ROOT)))
    return out


def stage(panel: str, held: str) -> dict:
    """Everything an arm needs; the 2.6 GB operator is released after the whitening."""
    audit = S.bind_panel(panel, resid218)
    fold = load_fold(panel, held)
    targets = fold["targets"]
    wcols = np.asarray(fold["wcols"])
    sc = np.asarray(fold["scored"])
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
    opinfo = fold233.verify_op(panel, held)                     # the axis assertion (233 §2)
    pre, votes = fold233.consensus(fold, pm=True, only=None)
    pre = np.asarray(pre, dtype=np.float64)
    cd4 = S.cd4_prepared(panel, targets)
    h1 = union_rows([resid218.source_residual(n, targets, "gene_norm") for n in UNION])
    op = Operator(fold["op_path"])
    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    # KNOWN-ANSWER: 233's own `build_final` (verbatim) on this operator = the reference stack
    aud233: dict = {}
    final233 = S.build_final(fold, op, pre, cd4, h1, wcols, CD4_F, UNION_G,
                             (resid218, compose218), cd4_sizing="all_rows", audit=aud233)
    del op
    frac_u, pu_u = solve_fraction_carrying(whitened, h1, UNION_G)
    return dict(panel=panel, held=held, fold=fold, targets=targets, truth=truth, pre=pre,
                whitened=whitened, cd4=cd4, h1=h1, frac_u=float(frac_u), pu_u=float(pu_u),
                wcols=wcols, panel_positions=fold["panel_positions"], opinfo=opinfo,
                bind=audit, final233=final233, aud233=aud233)


# ------------------------------------------------------------------------------ sizing
def solve_cd4(whitened, res, sizing: str):
    if sizing == "allrows":
        return solve_fraction(whitened, res, CD4_F)                 # production A/B/C
    if sizing == "carrying":
        return solve_fraction_carrying(whitened, res, G_CD4)        # October, 247 DEF_ENV
    raise ValueError(sizing)


def residual_norms(res) -> dict:
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    n = np.linalg.norm(added, axis=1)
    live = n[n > 0]
    return dict(median_all=float(np.median(n)),
                median_carrying=float(np.median(live)) if live.size else 0.0,
                carrying=int(live.size), counted=int((np.asarray(res.target_condition_count) > 0).sum()))


# --------------------------------------------------------------------------- composition
def merged_many(parts):
    """246 `lib246.merged_many`, verbatim: ONE add, ONE global preserve; == `merged` for two."""
    parts = [(r, f) for r, f in parts if r is not None and f != 0.0]
    assert parts, "nothing to merge"
    base_r, base_f = parts[0]
    out = merged(base_r, base_f, None, 0.0)
    eff = np.asarray(out.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(out.gene_votes) > 0
    counts = np.asarray(out.target_condition_count).astype(np.int64)
    for r, f in parts[1:]:
        v = np.asarray(r.gene_votes) > 0
        eff += np.where(v, np.asarray(r.mean_effect, dtype=np.float64), 0.0) * f
        votes |= v
        counts = np.maximum(counts, np.asarray(r.target_condition_count).astype(np.int64))
    return E062.CD4Prepared(mode="merged", mean_effect=eff,
                            gene_votes=votes.astype(np.uint8),
                            target_condition_count=np.minimum(counts, 255).astype(np.uint8),
                            condition_details=[])


def compose(st: dict, cd4_res, frac_cd4: float, extra=()):
    """whitened + ONE merged residual (CD4 @ frac + H1 union @ its fraction [+ extra parts]),
    ONE preserve, panel columns zeroed -- 218 §7.2 / 222's builder / 233's `build_final` tail."""
    parts = [(cd4_res, frac_cd4), (st["h1"], st["frac_u"])] + list(extra)
    final, preserve = E062.add_post_whiten_residual(st["whitened"], merged_many(parts), 1.0,
                                                    st["wcols"])
    final[:, st["panel_positions"]] = 0.0
    return final, float(preserve["post_add_global_preserve_scale"])


def derange(final: np.ndarray, seed: int = NULL_SEED) -> np.ndarray:
    """233's `null` arm: the predicted rows deranged (score233.main, verbatim)."""
    rng = np.random.default_rng(seed)
    n = final.shape[0]
    perm = rng.permutation(n)
    for _ in range(50):
        if not (perm == np.arange(n)).any():
            break
        perm = rng.permutation(n)
    return final[perm]


def score(st: dict, final: np.ndarray) -> dict:
    """246's scorer verbatim (`resid218.score_effects`): a23 minnorm (pds / cos / cos_agg) and
    a15 (q_w = `a15_sgn_w`), mean over challenge frames A/B/C, per-target arrays kept."""
    return resid218.score_effects(final, st["fold"], st["truth"])


# ----------------------------------------------------------------------------- deflation
def _carried(res):
    return np.asarray(res.target_condition_count) > 0


def _as_res(res, eff, mode, k):
    return E062.CD4Prepared(mode=f"{res.mode}|{mode}{k if k else ''}", mean_effect=eff,
                            gene_votes=np.asarray(res.gene_votes).copy(),
                            target_condition_count=np.asarray(
                                res.target_condition_count).copy(),
                            condition_details=res.condition_details)


def deflate(res, mode: str, k: int = 0, *, seed: int = RAND_SEED):
    """`lib246.deflate`, verbatim (modes none / centre / svd / rand; the `op` mode is not used).

    'svd': per-gene centring over the CARRIED rows (count > 0: complete rows plus the zero
    partial rows, as production's `freeze_complete_mean3` hands them), re-zero outside the
    gene support, remove the top-k left-singular directions of the centred matrix, re-zero.
    == the builder's `deflate_cd4_residual` (249 §1; asserted in basis250.py on the shipped
    object)."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(_carried(res))
    audit = dict(mode=mode, k=int(k), n_carried=int(rows.size))
    if mode == "none" or rows.size == 0:
        audit["energy_removed_frac"] = 0.0
        return _as_res(res, eff, mode, k), audit
    R = eff[rows]
    e0 = float((R * R).sum())
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    if mode == "centre":
        R = Rc
    elif mode in ("svd", "rand"):
        assert k > 0
        if mode == "svd":
            G = Rc @ Rc.T
            w, U = np.linalg.eigh(G)
            Uk = U[:, np.argsort(w)[::-1][:k]]
            audit["top_k_energy_frac_of_centred"] = float(
                np.sort(w)[::-1][:k].sum() / max(w.sum(), 1e-300))
            audit["basis_u_k_sha256"] = sha256_arr(Uk)
        else:
            rng = np.random.default_rng(seed + 1000 * k)
            Uk, _ = np.linalg.qr(rng.standard_normal((rows.size, k)))
        R = Rc - Uk @ (Uk.T @ Rc)
        R[~sup[rows]] = 0.0
    else:
        raise ValueError(mode)
    e1 = float((R * R).sum())
    audit["energy_removed_frac"] = float(1.0 - e1 / e0) if e0 > 0 else 0.0
    eff[rows] = R
    return _as_res(res, eff, mode, k), audit


def rand_energy_matched(res, target_frac: float, seed: int, kmax: int = 260):
    """246 `sweep246b.rand_energy_matched`, verbatim: the smallest random subspace dimension
    whose removed energy >= what svd-k removes (bisection)."""
    lo, hi = 1, kmax
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        d, aud = deflate(res, "rand", mid, seed=seed)
        if aud["energy_removed_frac"] >= target_frac:
            best = (d, aud, mid)
            hi = mid - 1
        else:
            lo = mid + 1
    if best is None:
        d, aud = deflate(res, "rand", kmax, seed=seed)
        best = (d, aud, kmax)
    return best


def shuffled_projector(res, k: int, seed: int):
    """246 `sweep246b.shuffled_projector`, verbatim: top-k fitted on a row-deranged copy,
    its left projector applied to the REAL row order."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(_carried(res))
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
    Uk = U[:, np.argsort(w)[::-1][:k]]
    Rk = Rc - Uk @ (Uk.T @ Rc)
    Rk[~sup[rows]] = 0.0
    e0 = float((R * R).sum())
    eff[rows] = Rk
    audit = dict(mode="svd_shufdef", k=int(k), n_carried=int(rows.size),
                 energy_removed_frac=float(1.0 - (Rk * Rk).sum() / e0))
    return _as_res(res, eff, "svd_shufdef", k), audit


def split_half_deflate(res, k: int, seed: int):
    """246 `sweep246c.split_half_deflate`, verbatim: basis fitted on a random half of the
    carried rows (gene-space right-singular vectors), applied to the OTHER half, both ways."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(_carried(res))
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
        Vk = (Rfc.T @ Uk) / sk[None, :]
        Ra = R[apply] - mean
        Ra[~sup[rows][apply]] = 0.0
        Ra = Ra - (Ra @ Vk) @ Vk.T
        Ra[~sup[rows][apply]] = 0.0
        out[apply] = Ra
    eff[rows] = out
    audit = dict(mode="svd_splithalf", k=int(k), n_carried=int(rows.size),
                 energy_removed_frac=float(1.0 - (out * out).sum() / e0))
    return _as_res(res, eff, "svd_splithalf", k), audit


# ------------------------------------------------------------- the gene-space basis (Phase 3)
def centred_rows(res):
    """(rows, R, Rc, sup_rows): the carried rows, raw and centred on their own mean."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    sup = np.asarray(res.gene_votes) > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(_carried(res))
    R = eff[rows]
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rows]] = 0.0
    return rows, R, Rc, sup[rows]


def gene_basis(res, k: int):
    """V_k (genes x k, orthonormal): the top-k RIGHT-singular directions of the centred carried
    rows -- the object that is comparable ACROSS panels (the left vectors U_k live in row space
    and change with the rows).  V_k = R_c^T U_k / s_k, exactly 246's split-half construction."""
    rows, R, Rc, _ = centred_rows(res)
    G = Rc @ Rc.T
    w, U = np.linalg.eigh(G)
    order = np.argsort(w)[::-1]
    idx = order[:k]
    Uk, sk = U[:, idx], np.sqrt(np.maximum(w[idx], 1e-300))
    Vk = (Rc.T @ Uk) / sk[None, :]
    tot = float(max(w.sum(), 1e-300))
    audit = dict(k=int(k), n_carried=int(rows.size),
                 top_k_energy_frac_of_centred=float(w[idx].sum() / tot),
                 spectrum_frac=[float(x / tot) for x in w[order[:15]]],
                 centred_energy=tot, raw_energy=float((R * R).sum()),
                 basis_u_k_sha256=sha256_arr(Uk),
                 orthonormality_err=float(np.abs(Vk.T @ Vk - np.eye(k)).max()))
    return Vk, Uk, audit


def deflate_gene_basis(res, V, tag: str):
    """Deflate `res` with a gene-space basis fitted ELSEWHERE (the pinned A/B/C basis):
    centre on the residual's OWN carried rows, R_k = R_c (I - V V^T), re-zero outside support."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    sup = np.asarray(res.gene_votes) > 0
    eff[~sup] = 0.0
    rows, R, Rc, sup_rows = centred_rows(res)
    e0 = float((R * R).sum())
    ec = float((Rc * Rc).sum())
    proj = Rc @ V
    captured = float((proj * proj).sum() / max(ec, 1e-300))
    Rk = Rc - proj @ V.T
    Rk[~sup_rows] = 0.0
    e1 = float((Rk * Rk).sum())
    eff[rows] = Rk
    audit = dict(mode=f"pinned_{tag}", k=int(V.shape[1]), n_carried=int(rows.size),
                 energy_removed_frac=float(1.0 - e1 / e0) if e0 > 0 else 0.0,
                 captured_frac_of_centred=captured)
    return _as_res(res, eff, f"pin{tag}", V.shape[1]), audit


def principal_angles_deg(Va, Vb) -> list[float]:
    s = np.linalg.svd(Va.T @ Vb, compute_uv=False)
    return [float(np.degrees(np.arccos(min(1.0, max(-1.0, x))))) for x in s]


def production_cd4_abc():
    """The shipped A/B/C CD4 residual as the builder held it before deflation (249 `cd4-off.npz`:
    300 rows, 291 carried, 276 complete)."""
    z = np.load(PROD_CD4_ABC)
    return E062.CD4Prepared(mode="production_abc_off", mean_effect=np.asarray(z["mean_effect"], np.float64),
                            gene_votes=np.asarray(z["gene_votes"], np.uint8),
                            target_condition_count=np.asarray(z["target_condition_count"], np.uint8),
                            condition_details=[]), tuple(str(t) for t in z["targets"])


# ---------------------------------------------------------------- CD4 partial rows (Phase 4)
def cd4_partial(panel: str, targets, cd4_root: Path | None = None):
    """246 `lib246.cd4_partial(equalise=False)`, parameterised by panel: CD4 rows carried by ONE
    or TWO conditions -- the rows production's complete-mean3 rule zeroes -- as a residual that
    is zero on the complete rows."""
    root = cd4_root or B.cd4_root(panel)
    banks = [E062.Bank.load(n, root / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, OFFICIAL, set(B.read_panel(panel)), "gene_norm",
                                 chunk_rows=256, max_rss_gib=10.0)
    count = np.asarray(prep.target_condition_count, dtype=np.int64)
    complete = count == len(CD4_NAMES)
    partial = (count > 0) & ~complete
    eff = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    norms_c = np.linalg.norm(np.where(votes[complete] > 0, eff[complete], 0.0), axis=1)
    eff[~partial] = 0.0
    votes[~partial] = 0
    cnt = np.where(partial, count, 0).astype(np.uint8)
    return E062.CD4Prepared(mode="cd4_partial", mean_effect=eff, gene_votes=votes,
                            target_condition_count=cnt,
                            condition_details=[dict(n_partial=int(partial.sum()),
                                                    n_one=int((count == 1).sum()),
                                                    n_two=int((count == 2).sum()),
                                                    median_complete_norm=float(
                                                        np.median(norms_c)) if norms_c.size else 0.0)])
