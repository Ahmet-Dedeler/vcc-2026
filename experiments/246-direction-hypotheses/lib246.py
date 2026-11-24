"""246 — shared construction: the shipped stack on 162's PLOLO, plus the deflation arms.

Everything here is 218's / 222's construction, imported, not re-derived:
  * `plolo`            162's instrument, imported verbatim (fold caches copied from 222's
                       `work/plolo/`, sha-identical to 218's).
  * `resid218`         `score_effects`, `source_residual`, `self_residual`, `shuffled`,
                       ALPHA / GAMMA / LAM.
  * `compose218`       `union_rows`, `solve_fraction_carrying`, `merged`, `restrict_rows`.

The ONLY new objects are the deflation operators in `deflate()` (the specific-subspace
relaxation the mandate names: per-gene centring across the carried rows, then removal of the
top-k right-singular directions of the residual matrix; plus the operator's own top-k
eigenvectors, and a random-subspace control), `merged_many()` (218's `merged()` for three
residuals) and `cd4_partial()` (the CD4 rows production discards because fewer than three
conditions carry the target).

Known-answer rows are asserted by the caller (`sweep246.py`).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E218 = ROOT / "experiments/218-independent-donors"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

# resid218 / compose218 set `plolo.WORK` to 218's directory at import; pin it back to OURS
# so every cache read and (never expected) cache write stays inside this experiment.
from resid218 import (ALPHA, GAMMA, LAM, CD4_NAMES, CD4_ROOT,  # noqa: E402
                      score_effects, self_residual, shuffled, source_residual)
from compose218 import (merged, restrict_rows, solve_fraction_carrying,  # noqa: E402
                        union_rows)

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, OFFICIAL, PANEL_SET, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, solve_fraction)

CD4_F = 0.34            # board incumbent `cd4f034`, kept by `h1union`
H1_G = 0.28             # shipped `h1union` dose
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
SHUF_SEED = 20726       # 218's derangement seed
RAND_SEED = 20260904


# ----------------------------------------------------------------------------- staging
def op_path_for(held: str) -> Path:
    return (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
            else PROD_OPS / "whitener_union_pruned.npz")


def stage(held: str, kmax_op: int = 20) -> dict:
    """Everything an arm needs, with the 2.5 GB operator released after use."""
    fold = build_fold(held)
    targets = fold["targets"]
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    cd4 = cd4_for(targets)
    op = Operator(op_path_for(held))
    whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)
    order = np.argsort(op.evals)[::-1][:kmax_op]
    op_top = np.ascontiguousarray(op.evecs[:, order])          # (n_wcols, kmax) float64
    op_top_evals = np.asarray(op.evals[order], dtype=np.float64)
    del op
    union = union_rows([source_residual(n, targets, "gene_norm")
                        for n in ("gse295214_h1", "h1_2025")])
    frac_cd4, _ = solve_fraction(whitened, cd4, CD4_F)
    frac_x, _ = solve_fraction_carrying(whitened, union, H1_G)
    return dict(held=held, fold=fold, targets=targets, truth=truth, pre=pre, cd4=cd4,
                whitened=whitened, union=union, frac_cd4=float(frac_cd4),
                frac_x=float(frac_x), op_top=op_top, op_top_evals=op_top_evals)


def cd4_partial(targets, *, equalise: bool):
    """CD4 rows carried by ONE or TWO conditions -- the rows `cd4_for`'s complete-mean3 rule
    zeroes.  Production's own `prepare_cd4_mean` averages over whichever conditions carry the
    target; only the caller (`plolo.cd4_for`) discards incomplete rows.  Returned as a
    residual that is zero on the complete rows, so it can ride on top of the shipped stack.

    `equalise=True` rescales every partial row to the median norm of the complete rows, so a
    one-condition row (not shrunk by averaging) does not carry more weight than a two-condition
    row.  With carrying-row sizing a global scalar cancels; this per-row-class one does not.
    """
    banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, OFFICIAL, PANEL_SET, "gene_norm",
                                 chunk_rows=256, max_rss_gib=40.0)
    count = np.asarray(prep.target_condition_count, dtype=np.int64)
    complete = count == len(CD4_NAMES)
    partial = (count > 0) & ~complete
    eff = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    norms_c = np.linalg.norm(np.where(votes[complete] > 0, eff[complete], 0.0), axis=1)
    if equalise:
        ref = float(np.median(norms_c))
        for i in np.flatnonzero(partial):
            n = float(np.linalg.norm(np.where(votes[i] > 0, eff[i], 0.0)))
            if n > 0:
                eff[i] *= ref / n
    eff[~partial] = 0.0
    votes[~partial] = 0
    cnt = np.where(partial, count, 0).astype(np.uint8)
    return E062.CD4Prepared(mode=f"cd4_partial{'_eq' if equalise else ''}",
                            mean_effect=eff, gene_votes=votes, target_condition_count=cnt,
                            condition_details=[dict(n_partial=int(partial.sum()),
                                                    n_one=int((count == 1).sum()),
                                                    n_two=int((count == 2).sum()),
                                                    median_complete_norm=float(
                                                        np.median(norms_c)))])


# ---------------------------------------------------------------------------- deflation
def _carried(res):
    return np.asarray(res.target_condition_count) > 0


def deflate(res, mode: str, k: int = 0, *, seed: int = RAND_SEED, op_top=None,
            wcols=None):
    """Remove the component SHARED across the carried rows from an un-whitened residual.

    mode:
      'none'    identity (k ignored) -- the known-answer construction.
      'centre'  subtract the per-gene mean over the carried rows (rank-1, the across-target
                mean), on the residual's own support.
      'svd'     'centre', then remove the top-k left-singular directions of the centred
                carried-row matrix R_c:  R_k = (I - U_k U_k^T) R_c, which equals removing the
                top-k RIGHT-singular directions (the rows' shared gene-space axes).
      'rand'    'centre', then remove a RANDOM k-dimensional subspace of R_c's row space --
                209's permutation-style control: energy removal without direction selection.
      'op'      remove the whitening operator's own top-k eigenvectors (the donors' dominant
                PCs, what gamma = 0.65 suppresses hardest) from the residual on the operator
                columns.  No centring.
    Returns a new CD4Prepared and an audit dict (energy removed, subspace dims).
    """
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
    if mode == "op":
        assert op_top is not None and wcols is not None and k > 0
        V = op_top[:, :k]                                     # (n_wcols, k), orthonormal
        Rw = R[:, wcols]
        R[:, wcols] = Rw - (Rw @ V) @ V.T
        R[~sup[rows]] = 0.0
    else:
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


def _as_res(res, eff, mode, k):
    return E062.CD4Prepared(mode=f"{res.mode}|{mode}{k if k else ''}", mean_effect=eff,
                            gene_votes=np.asarray(res.gene_votes).copy(),
                            target_condition_count=np.asarray(
                                res.target_condition_count).copy(),
                            condition_details=res.condition_details)


# --------------------------------------------------------------------------- composition
def merged_many(parts):
    """218 §7.2's `merged()` for any number of (residual, fraction) pairs: ONE add, ONE
    global preserve.  Reduces to `compose218.merged` for two parts."""
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


def compose(st: dict, parts) -> np.ndarray:
    """Whitened base + ONE merged residual add + ONE preserve; panel columns zeroed."""
    fold = st["fold"]
    if not [p for p in parts if p[0] is not None and p[1] != 0.0]:
        final = np.asarray(st["whitened"], dtype=np.float64).copy()
    else:
        final, _ = E062.add_post_whiten_residual(
            st["whitened"], merged_many(parts), 1.0, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    return final


def score(st: dict, final: np.ndarray) -> dict:
    return score_effects(final, st["fold"], st["truth"])


def stack_parts(st: dict):
    return [(st["cd4"], st["frac_cd4"]), (st["union"], st["frac_x"])]


def residual_norms(res):
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    n = np.linalg.norm(added, axis=1)
    live = n[n > 0]
    return dict(median_all=float(np.median(n)), median_carrying=float(np.median(live))
                if live.size else 0.0, carrying=int(live.size))
