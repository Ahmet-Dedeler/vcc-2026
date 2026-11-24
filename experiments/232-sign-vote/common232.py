"""232 — shared machinery: the fold/source loader, the production stack, and the VOTE.

THE OBJECT UNDER TEST
---------------------
`E062.build_baseline` is an equal-weight MEAN over the sources that carry a target:

    pre[i, j] = ( sum_s E_s[i, j] ) / |{s : target i available in s}|

The denominator is a positive scalar, so

    sign(pre[i, j]) == sign( sum_s E_s[i, j] )

i.e. production's per-gene sign is the sign of the SUM over covering donors -- a
MAGNITUDE-weighted vote in which one large donor outvotes two agreeing ones.  204 swept 29
donor weightings; every one of them is a weighted mean, so all 29 share this sign rule.
This module builds the alternatives: the unweighted majority, and confidence-weighted votes
whose weight is a per-(donor, gene) quantity rather than a per-gene one.

WHY NOT 189's `theta`.  `theta_g` is a per-GENE reliability from a hierarchical logit fit.
In a weighted sign vote the weight enters as `sign( sum_s w_s(i,j) * sign(e_s[i,j]) )`; a
weight that does not depend on `s` factors straight out of the sum and the vote is
unchanged.  **A per-gene weight cancels exactly, at every arity.**  So `theta` cannot serve
here and the confidence weights below are all per-(donor, gene) or per-(donor, target, gene).

COVERAGE.  On the official 18,533-gene axis a donor that does not carry gene `j` writes a
literal `0.0` into its effect row.  That zero is NO VOTE, not "no response".  Every
statistic here uses `cov[s, i, j] = target_available[s, i] & gene_support[s, j]`.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   deliver, endpoints, solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from resid218 import ALPHA, GAMMA, LAM, score_effects, source_residual  # noqa: E402
from compose218 import merged, solve_fraction_carrying, union_rows  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments/188-magnitude-calibration"))
sys.path.insert(0, str(ROOT / "experiments/203-reopen-125"))
from mapfit import to_lfc  # noqa: E402
from phase23 import POOL_N, key_precision, topn, umd  # noqa: E402

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALL_FOLDS = CANCER + ("kolf2",)
CD4_F, UNION_G = 0.34, 0.28
UNION = ("gse295214_h1", "h1_2025")
PROD204 = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10)   # 204's geometry, for KA-C

# 223 `work/focus223.json` via 228 `arms228.KA` -- the shipped stack's own (cos, pds)
KA_STACK = {
    ("k562_gwps", "cd4"): (0.03635247771108033, 0.7958811888512827),
    ("k562_gwps", "base222"): (0.03943565390646537, 0.8053968077826482),
    ("orion_hct116", "cd4"): (0.02968920878381447, 0.790140396929563),
    ("orion_hct116", "base222"): (0.032080093772601896, 0.7922125442469077),
    ("orion_hek293t", "cd4"): (0.024558952922866917, 0.7605504090145728),
    ("orion_hek293t", "base222"): (0.026156872210626163, 0.7636228552854588),
    ("kolf2", "cd4"): (0.009016327626527735, 0.6145272211448684),
    ("kolf2", "base222"): (0.012520242539233312, 0.6361551030668678),
}
# 204 RESULTS.md §4, "Absolute equal-weight baseline `pds`", 204's geometry
KA_204_PDS = dict(k562_gwps=0.78912, orion_hct116=0.79182,
                  orion_hek293t=0.77738, kolf2=0.61013)


# ----------------------------------------------------------------------- loading
def load_sources(held: str) -> dict:
    z = np.load(WORK / f"sources-{held}.npz", allow_pickle=True)
    out = dict(
        names=[str(x) for x in z["names"]],
        E=np.asarray(z["effects"], dtype=np.float64),
        avail=np.asarray(z["available"], dtype=bool),
        supp=np.asarray(z["support"], dtype=bool),
        votes=np.asarray(z["votes"], dtype=np.int64),
        targets=tuple(str(t) for t in z["targets"]),
        ka_rel=float(z["ka_rel"]))
    if "extra_names" in z.files:
        out.update(extra_names=[str(x) for x in z["extra_names"]],
                   extra_E=np.asarray(z["extra_effects"], dtype=np.float64),
                   extra_avail=np.asarray(z["extra_available"], dtype=bool),
                   extra_supp=np.asarray(z["extra_support"], dtype=bool))
    return out


def coverage(src: dict) -> np.ndarray:
    """`cov[s, i, j]` -- does donor `s` have a vote on (target i, gene j)?"""
    return src["avail"][:, :, None] & src["supp"][:, None, :]


def stack_sources(src: dict, use_extra: bool = False):
    """Return (names, E, avail, supp) optionally including the non-production voters."""
    if not use_extra or "extra_E" not in src:
        return src["names"], src["E"], src["avail"], src["supp"]
    return (src["names"] + src["extra_names"],
            np.concatenate([src["E"], src["extra_E"]], 0),
            np.concatenate([src["avail"], src["extra_avail"]], 0),
            np.concatenate([src["supp"], src["extra_supp"]], 0))


def weighted_baseline(E, avail, w, fallback=None) -> np.ndarray:
    """204's `weighted_baseline`; `w = 1` reproduces `E062.build_baseline` exactly.

    A target no source carries takes 062's fallback row
    (`EFFECT_SCALE * mean_response / |mean_response|`), which 204 deliberately leaves
    UNWEIGHTED.  It is therefore identical in every arm, so it is copied from the
    incumbent `pre` rather than recomputed -- exact, and it cannot drift between arms.
    """
    num = np.zeros(E.shape[1:], dtype=np.float64)
    den = np.zeros(E.shape[1], dtype=np.float64)
    for s in range(E.shape[0]):
        rows = avail[s]
        num[rows] += float(w[s]) * E[s][rows]
        den += np.where(rows, float(w[s]), 0.0)
    covered = den > 0
    num[covered] /= den[covered, None]
    if (~covered).any():
        assert fallback is not None, "uncovered target and no fallback row supplied"
        num[~covered] = np.asarray(fallback, dtype=np.float64)[~covered]
    return num


# --------------------------------------------------------------- the vote itself
def vote_signs(E, avail, supp, weights=None) -> tuple[np.ndarray, np.ndarray]:
    """(vote sign, covering-donor count) on every (target, gene).

    `weights` is `None` for the unweighted majority, or an array broadcastable to
    (S, T, G) of NON-NEGATIVE per-donor confidences.  Ties return 0.
    """
    cov = avail[:, :, None] & supp[:, None, :]
    s = np.sign(E) * cov
    if weights is None:
        tot = s.sum(0)
    else:
        tot = (s * weights).sum(0)
    return np.sign(tot), cov.sum(0).astype(np.int16)


def z_weights(E, avail, supp) -> np.ndarray:
    """`|e_s(i,j)| / sd_s(j)` -- the donor's own response standardised by its own
    per-gene spread across targets.  Per (donor, target, gene), so it does NOT cancel.
    It differs from `sign(mean)` because `sign(mean)` compares RAW magnitudes across
    donors whose per-gene scales differ; this compares each donor to its own null."""
    out = np.zeros_like(E)
    for s in range(E.shape[0]):
        rows = avail[s]
        sd = E[s][rows].std(0)
        sd = np.where(sd > 0, sd, np.inf)
        out[s] = np.abs(E[s]) / sd[None, :]
    return out * (avail[:, :, None] & supp[:, None, :])


def rho_weights(E, avail, supp) -> np.ndarray:
    """Per-(donor, gene) leave-one-in agreement: corr across targets between donor `s`'s
    column `j` and the MEAN of the other covering donors' column `j`.

    At arity 2 this is symmetric (both donors get the same number) and therefore cancels
    exactly, like `theta`.  It is only informative at arity >= 3.  Reported, not hidden."""
    S = E.shape[0]
    out = np.zeros((S, E.shape[2]), dtype=np.float64)
    for s in range(S):
        others = [t for t in range(S) if t != s]
        rows = avail[s] & np.logical_and.reduce([avail[t] for t in others])
        if rows.sum() < 5:
            continue
        a = E[s][rows]
        b = np.mean([E[t][rows] for t in others], axis=0)
        a = a - a.mean(0)
        b = b - b.mean(0)
        na, nb = np.linalg.norm(a, axis=0), np.linalg.norm(b, axis=0)
        good = (na > 0) & (nb > 0)
        r = np.zeros(E.shape[2])
        r[good] = (a[:, good] * b[:, good]).sum(0) / (na[good] * nb[good])
        out[s] = np.clip(r, 0.0, None)
    return out[:, None, :]


def vote_bundle(E, avail, supp) -> dict:
    """Every sign rule and its CONFIDENCE, on the full official axis.

    `margin` is the dose knob: the vote's normalised decisiveness in [0, 1].
      majority   |pos - neg| / (pos + neg)
      z-vote     |sum_s z_s sign(e_s)| / sum_s z_s
    At arity 2 the majority margin is identically 0 on every contested pair -- a tie is
    not a rule -- which is why the dose is swept on the z-vote.
    """
    cov = avail[:, :, None] & supp[:, None, :]
    s = (np.sign(E) * cov).astype(np.int8)
    pos = (s > 0).sum(0).astype(np.int16)
    neg = (s < 0).sum(0).astype(np.int16)
    ncov = cov.sum(0).astype(np.int16)
    zw = z_weights(E, avail, supp)
    zsum = zw.sum(0)
    zscore = (s * zw).sum(0)
    del zw
    return dict(
        ncov=ncov, pos=pos, neg=neg,
        s_mean=np.sign(E.sum(0)).astype(np.int8),
        s_maj=np.sign(pos.astype(np.int32) - neg).astype(np.int8),
        s_z=np.sign(zscore).astype(np.int8),
        margin_maj=np.abs(pos.astype(np.float64) - neg) / np.maximum(pos + neg, 1),
        margin_z=np.abs(zscore) / np.maximum(zsum, 1e-300),
        disagree=(pos > 0) & (neg > 0))


def flip_mask(vb: dict, rule: str, tau: float) -> tuple[np.ndarray, np.ndarray]:
    """The coordinates the arm is allowed to touch, and the sign it imposes there.

    A coordinate qualifies only if (a) at least TWO donors have a vote on it -- a gene
    one donor carries has no vote and is LEFT ALONE -- (b) the rule returns a non-zero
    sign, (c) that sign DISAGREES with `sign(mean)`, the incumbent, and (d) the vote is
    decisive at dose `tau`.
    """
    s_r = vb[f"s_{rule}"]
    margin = vb[f"margin_{rule}"]
    m = (vb["ncov"] >= 2) & (s_r != 0) & (s_r != vb["s_mean"]) & (margin >= tau)
    return m, s_r


def impose_sign(final: np.ndarray, mask: np.ndarray, s_target: np.ndarray) -> np.ndarray:
    """`|effect|` is bit-identical to the incumbent's; only signs move.  Asserted.

    A coordinate whose target sign is 0 is LEFT ALONE -- imposing a zero sign would zero
    the magnitude, which is a magnitude arm and §8.1d says those are dead.  This is the
    only place the invariant could have been broken and it is checked here, elementwise.
    """
    eff = mask & (s_target != 0)
    out = final.copy()
    out[eff] = np.abs(final[eff]) * s_target[eff]
    assert np.array_equal(np.abs(out), np.abs(final)), "MAGNITUDE INVARIANCE VIOLATED"
    return out


# -------------------------------------------------------------- production stack
class FoldBundle:
    """Everything the production stack needs for one fold, built once."""

    def __init__(self, held: str, limit_gib: float = 7.5):
        self.held = held
        self.fold = build_fold(held, limit_gib=limit_gib)
        f = self.fold
        self.wcols = np.asarray(f["wcols"])
        self.sc = np.asarray(f["scored"])
        self.truth_h = np.asarray(f["truth"], dtype=np.float64)
        self.truth = self.truth_h[:, self.sc]
        self.pre = np.asarray(f["pre"], dtype=np.float64)
        self.targets = f["targets"]
        self.cd4 = cd4_for(self.targets)
        self.h1 = union_rows([source_residual(n, self.targets, "gene_norm")
                              for n in UNION])
        self.op = Operator(PROD_OPS / "whitener_union_pruned.npz"
                           if FOLDS[held]["op_kind"] == "production"
                           else Path(FOLDS[held]["op"]))

    def whitened(self, pre=None, alpha_scale: float = 1.0, gamma=GAMMA, lam=LAM,
                 alpha=ALPHA):
        p = self.pre if pre is None else pre
        return self.op.apply(p * (alpha * alpha_scale), self.wcols, gamma, lam)

    def final(self, pre=None, alpha_scale: float = 1.0, stack: str = "base222"):
        """The shipped stack: whiten, then CD4 @0.34 + H1 union @0.28, panel zeroed."""
        wh = self.whitened(pre, alpha_scale)
        frac_cd4, _ = solve_fraction(wh, self.cd4, CD4_F)
        if stack == "base222":
            frac_u, _ = solve_fraction_carrying(wh, self.h1, UNION_G)
            res = merged(self.cd4, frac_cd4, self.h1, frac_u)
        elif stack == "cd4":
            res = merged(self.cd4, frac_cd4, None, 0.0)
        else:
            raise ValueError(stack)
        fin, _ = E062.add_post_whiten_residual(wh, res, 1.0, self.wcols)
        fin[:, self.fold["panel_positions"]] = 0.0
        return fin

    def final_gamma(self, gamma: float, pre=None):
        """The shipped stack with ONE knob moved: the whitening exponent.  Every other
        input -- alpha, lam, the operator, CD4 @ 0.34, H1 union @ 0.28, the panel zeroing,
        the render -- is the incumbent's."""
        wh = self.whitened(pre, 1.0, gamma=gamma)
        frac_cd4, _ = solve_fraction(wh, self.cd4, CD4_F)
        frac_u, _ = solve_fraction_carrying(wh, self.h1, UNION_G)
        fin, _ = E062.add_post_whiten_residual(
            wh, merged(self.cd4, frac_cd4, self.h1, frac_u), 1.0, self.wcols)
        fin[:, self.fold["panel_positions"]] = 0.0
        return fin

    def final_204(self, pre=None):
        """204's geometry: alpha 0.982, gamma 0.65, lam 0, CD4 @0.10, NO H1 union."""
        wh = self.op.apply((self.pre if pre is None else pre) * PROD204["alpha"],
                           self.wcols, PROD204["gamma"], PROD204["lam"])
        frac, _ = solve_fraction(wh, self.cd4, PROD204["cd4_f"])
        fin, _ = E062.add_post_whiten_residual(wh, self.cd4, frac, self.wcols)
        fin[:, self.fold["panel_positions"]] = 0.0
        return fin

    def score(self, final) -> dict:
        e = score_effects(final, self.fold, self.truth)
        return {k: v for k, v in e.items() if not k.startswith("_")}


# ------------------------------------------------------------------- the DE pool
_GATE_CACHE: dict[str, np.ndarray] = {}


def ctx_gate(ctx: str) -> tuple[np.ndarray, np.ndarray]:
    if ctx not in _GATE_CACHE:
        b = plolo.challenge_bctrl(ctx)
        e = np.expm1(b)
        _GATE_CACHE[ctx] = (b, (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (b > 0))
    return _GATE_CACHE[ctx]


def de_pools(bundle: FoldBundle) -> dict[str, np.ndarray]:
    """188's three pools at the reference's BH depth (89 genes/row), six228's exact
    construction.  `wald89` is the one 188 validated on h1-test."""
    _, gA = ctx_gate("A")
    sc = bundle.sc
    b = bundle.fold["b_ctrl"]
    cpm_sc = 1e6 * np.expm1(b[sc]) / float(np.expm1(b).sum())
    LR = to_lfc(bundle.truth_h, b, gA)[:, sc]
    return {"wald89": topn(np.abs(LR) * np.sqrt(np.maximum(cpm_sc, 0)), POOL_N),
            "lfc89": topn(np.abs(LR), POOL_N),
            "h89": topn(np.abs(bundle.truth), POOL_N)}


def weighted_sign_accuracy(pred_sign, truth, pool, restrict=None) -> tuple[float, float]:
    """|t|-weighted sign accuracy over `pool` (optionally further restricted).

    Returns (accuracy, total weight).  `key_precision`'s convention: a zero predicted
    sign never counts as correct.  Targets whose restricted pool weight is 0 are dropped
    from the mean, and their count is reported by the caller.
    """
    vals, wsum = [], 0.0
    for i in range(truth.shape[0]):
        m = pool[i] if restrict is None else (pool[i] & restrict[i])
        if not m.any():
            continue
        w = np.abs(truth[i, m])
        s = float(w.sum())
        if s <= 0:
            continue
        ok = (pred_sign[i, m] != 0) & (pred_sign[i, m] == np.sign(truth[i, m]))
        vals.append(float(np.where(ok, w, 0.0).sum() / s))
        wsum += s
    return (float(np.mean(vals)) if vals else float("nan")), wsum


# ------------------------------------------- the shipped call set (228 §5, `rc_callset`)
RC_HEAD, RC_CPM = 400, 2200
MEDIAN_LIBRARY, CELLS = 20109.0, 400
N_REALS = (89, 200, 400)
SCORER_MIN_CPM = 5.0


def control_block(context: str) -> tuple[np.ndarray, np.ndarray]:
    """(keep, ref_lin) for one challenge context -- 228 `ctrl228.control_block`, cached
    inside 232's own work directory.  Both are properties of the CONTROL cells alone, so
    both are invariant to every arm: that is why 2,200 of the 2,600 called genes cannot move."""
    import h5py
    cache = WORK / f"ctrl-{context}.npz"
    if cache.exists():
        z = np.load(cache)
        return np.asarray(z["keep"], bool), np.asarray(z["ref_lin"], np.float64)
    path = ROOT / f"data/challenge/context_{context}.h5ad"
    with h5py.File(path, "r") as h:
        indptr = np.asarray(h["X/indptr"][:])
        n_cells = len(indptr) - 1
        n_genes = int(h["X"].attrs["shape"][1])
        acc = np.zeros(n_genes, dtype=np.float64)
        for s in range(0, n_cells, 2000):
            e = min(s + 2000, n_cells)
            lo, hi = int(indptr[s]), int(indptr[e])
            idx = np.asarray(h["X/indices"][lo:hi])
            dat = np.asarray(h["X/data"][lo:hi], dtype=np.float64)
            rows = np.repeat(np.arange(s, e), np.diff(indptr[s:e + 1]))
            tot = np.bincount(rows - s, weights=dat, minlength=e - s)
            acc += np.bincount(idx, weights=1e6 * dat / tot[rows - s], minlength=n_genes)
    ref_lin = acc / n_cells
    keep = ref_lin > SCORER_MIN_CPM
    np.savez_compressed(cache, keep=keep, ref_lin=ref_lin)
    return keep, ref_lin


def call_set(effect_row, gate10, keep, ref_lin):
    """`rc_emit.CallSet(400, 2200, 1.90)` exactly: |effect| head, then control-CPM block."""
    score = np.where(gate10 & keep, np.abs(np.nan_to_num(effect_row)), -np.inf)
    head = np.argsort(-score, kind="stable")[:RC_HEAD]
    rest = np.where(keep, ref_lin, -np.inf).copy()
    rest[head] = -np.inf
    cpm = np.argsort(-rest, kind="stable")[:RC_CPM]
    call = np.zeros(len(score), dtype=bool)
    call[head] = True
    call[cpm] = True
    return call, head


def block_digest(profile: np.ndarray) -> str:
    import hashlib
    tot = np.rint(CELLS * MEDIAN_LIBRARY * profile).astype(np.int64)
    return hashlib.sha256(tot.tobytes()).hexdigest()[:16]


def pool_weight_share(truth, pool, restrict) -> float:
    """Per-target |t|-weighted share of the pool that falls inside `restrict`, macro-averaged."""
    vals = []
    for i in range(truth.shape[0]):
        w = np.abs(truth[i, pool[i]])
        s = float(w.sum())
        if s <= 0:
            continue
        vals.append(float(np.abs(truth[i, pool[i] & restrict[i]]).sum() / s))
    return float(np.mean(vals)) if vals else float("nan")
