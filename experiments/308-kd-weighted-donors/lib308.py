"""308 -- |kd|-WEIGHTED DONORS: the weight constructions, frozen here before any fold was scored.

The feature (kd308.py): |kd|(t, d) = |delta_d[row t, column of gene t]|, donor screen d's own measured on-target
knockdown of target t, in d's own log1p(50,000 p) units.  Read from d's bank ONLY.  The production prepared rows
zero every panel column, so this coordinate is absent from the incumbent's base and cannot leak from it.

Cuts (both reported, both frozen before scoring):
  * ABS  = 293's absolute cut |kd| >= 1.0 (the challenge library's own call);
  * REL  = 296's library-relative cut C2 = the donor library's OWN median |kd| over every perturbation whose own
    gene is on its axis (k562_gwps 0.6571, orion_hct116 0.1563, orion_hek293t 0.1476) -- 296 measured that the
    absolute cut is NOT comparable across libraries and that the library-relative one is the transferable number.

The arms (all on top of cd4svd5; one variable each; W is (n, J), NaN rows keep the cached mean row):
  1. gate:<cut>          w_j = 1[|kd|_j >= cut_j], renormalised over the survivors; if EVERY carrying donor fails,
                         the row falls back to the cached equal-weight mean (the fallback rate is reported).
  2. soft:<cut>:p<P>     w_j ∝ (max(|kd|_j, EPS) / cut_j)^P, one parameter P, swept coarsely.  P = 0 is the
                         incumbent's equal weight (a known answer).  The standardisation by the donor's own cut
                         makes the weight scale-free across libraries, which the raw |kd| is not.
  3. <arm>:dis           the same restricted to the DISAGREEING set (>= 1 carrying donor at or above its cut and
                         >= 1 below); every other target keeps the incumbent's row verbatim.

Controls, frozen with the arms:
  * cd4svd5 / mean_explicit          the incumbent, and it rebuilt through the selection path (0.0e+00)
  * relw                             271's reliability-weighted mean, w_j ∝ R_j/(1-R_j) with R = SB(split-half
                                     cosine of the donor's own batch-parity halves); est271's `weights()` verbatim
                                     on 271's half rows -- THE DECISIVE CONTROL
  * relw_glob                        271's donor-level scalar reliability weights (the 019 / 142 class)
  * normrel1.5                       097 route (a) a6: w_j ∝ (||row_j(t)|| / median_t ||row_j(t)||)^(-1.5)
  * <arm>:shuf:s<k>                  |kd| permuted ACROSS TARGETS within each donor, 3 seeds
  * <arm>:rand:s<k>                  random weights matched to the arm's realised weight ENTROPY, 3 seeds
  * <arm>:rev                        reversed |kd|: the FAILING construct is weighted more (must lose)
Every weighted row is NORM-MATCHED to the cached mean row's norm (263's `:nm`, shape-only) unless `:lit`.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
sys.path.insert(0, str(E263))
sys.path.insert(0, str(E271))
import lib263 as Y  # noqa: E402  (-> lib257 -> lib246 -> plolo / resid218 / compose218)
import est271 as E  # noqa: E402  (read-only import; 271's frozen estimators)

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK

EPS = 1e-3                       # the floor under |kd| in the soft weight (a zero knockdown gets a tiny weight, not 0)
POWERS = (0.5, 1.0, 2.0, 4.0)    # the coarse sweep of the soft weight's one parameter
SEEDS = (0, 1, 2)
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def load_feature():
    z = np.load(WORK / "kd308_feature.npz", allow_pickle=True)
    return dict(panel=[str(t) for t in z["panel"]], donors=[str(d) for d in z["donors"]],
                kd=np.asarray(z["kd"], float), abs_kd=np.asarray(z["abs_kd"], float),
                measured=np.asarray(z["measured"], bool), cells=np.asarray(z["cells"], float),
                cut_rel=np.asarray(z["cut_rel"], float), cut_abs=float(z["cut_abs"]))


def kd_for_targets(targets: list[str], donors: list[str], feat: dict) -> dict:
    """The SAME feature for an arbitrary target list (h1-test's 85 2025-library targets are not the 2026 panel):
    |kd|(t, d) read straight from donor d's bank.  Identical construction to kd308.py; the cuts are the donors'
    own library-wide constants, already computed there."""
    import sys as _sys
    _sys.path.insert(0, str(HERE))
    import kd308 as K
    plolo = Y.plolo
    E062 = plolo.E062
    di = {d: j for j, d in enumerate(feat["donors"])}
    A = np.full((len(donors), len(targets)), np.nan)
    M = np.zeros((len(donors), len(targets)), bool)
    for j, name in enumerate(donors):
        bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
        gpos = K.own_column_map(bank)
        pairs, idx = [], []
        for i, t in enumerate(targets):
            r, c = bank.pert_pos.get(t), gpos.get(t)
            if r is not None and c is not None:
                pairs.append((int(r), int(c)))
                idx.append(i)
        if pairs:
            A[j, np.asarray(idx)] = np.abs(K.read_diagonal(bank, pairs))
            M[j, np.asarray(idx)] = True
        del bank
    cols = np.asarray([di[d] for d in donors], dtype=np.int64)
    return dict(A=A, M=M, cells=np.full_like(A, np.nan),
                cut_abs=np.full(len(donors), feat["cut_abs"]), cut_rel=feat["cut_rel"][cols])


def fold_feature(feat: dict, targets: list[str], donors: list[str]) -> dict:
    """Project the panel-wide feature onto ONE fold's target list and donor order.  (J, n) arrays."""
    ti = {t: i for i, t in enumerate(feat["panel"])}
    di = {d: j for j, d in enumerate(feat["donors"])}
    rows = np.asarray([ti[t] for t in targets], dtype=np.int64)
    cols = np.asarray([di[d] for d in donors], dtype=np.int64)
    A = feat["abs_kd"][np.ix_(rows, cols)].T                    # (J, n)
    M = feat["measured"][np.ix_(rows, cols)].T
    C = feat["cells"][np.ix_(rows, cols)].T
    return dict(A=A, M=M, cells=C, cut_abs=np.full(len(cols), feat["cut_abs"]),
                cut_rel=feat["cut_rel"][cols])


# ------------------------------------------------------------------------------- weight constructions
def _fill_missing(A: np.ndarray, M: np.ndarray, avail: np.ndarray) -> np.ndarray:
    """A carrying donor with NO measured |kd| (its bank does not measure the target's own gene, or does not carry
    the perturbation) gets that donor's MEDIAN |kd| over its measured targets -- the neutral fill, so missingness
    never becomes a covert weight."""
    out = A.copy()
    for j in range(A.shape[0]):
        bad = avail[j] & ~M[j]
        if bad.any():
            out[j, bad] = float(np.median(A[j, M[j]])) if M[j].any() else 1.0
    return out


def gate_weights(A: np.ndarray, avail: np.ndarray, cut: np.ndarray, *, gate: np.ndarray | None = None,
                 reverse: bool = False) -> tuple[np.ndarray, dict]:
    """Arm 1.  w_j = 1[|kd|_j >= cut_j] over the carrying donors, renormalised; all-fail -> keep the mean row."""
    J, n = avail.shape
    W = np.full((n, J), np.nan)
    pass_ = A >= cut[:, None]
    if reverse:
        pass_ = ~pass_
    n_fallback = n_acted = 0
    act = np.ones(n, bool) if gate is None else gate
    for t in np.flatnonzero((avail.sum(0) >= 2) & act):
        js = np.flatnonzero(avail[:, t])
        keep = js[pass_[js, t]]
        if keep.size == 0:
            n_fallback += 1                      # every carrying donor is below the cut: the incumbent's row stands
            continue
        if keep.size == js.size:
            continue                              # no donor dropped: the incumbent's row, left as the cached row
        W[t] = 0.0
        W[t, keep] = 1.0 / keep.size
        n_acted += 1
    return W, dict(n_fallback=n_fallback, n_acted=n_acted)


def soft_weights(A: np.ndarray, avail: np.ndarray, cut: np.ndarray, power: float, *,
                 gate: np.ndarray | None = None, reverse: bool = False) -> tuple[np.ndarray, dict]:
    """Arm 2.  w_j ∝ (max(|kd|_j, EPS) / cut_j)^power (reverse: ^-power)."""
    J, n = avail.shape
    W = np.full((n, J), np.nan)
    p = -power if reverse else power
    S = np.maximum(A, EPS) / cut[:, None]
    act = np.ones(n, bool) if gate is None else gate
    idx = np.flatnonzero((avail.sum(0) >= 2) & act)
    for t in idx:
        js = np.flatnonzero(avail[:, t])
        w = S[js, t] ** p
        if not np.isfinite(w).all() or w.sum() <= 0:
            continue
        W[t] = 0.0
        W[t, js] = w / w.sum()
    return W, dict(n_acted=int(len(idx)))


def entropy_of(W: np.ndarray) -> np.ndarray:
    """Per-target Shannon entropy (nats) of the weight row; NaN where the row is not set."""
    out = np.full(W.shape[0], np.nan)
    live = np.isfinite(W).all(1)
    P = np.clip(W[live], 1e-12, 1.0)
    out[live] = -(P * np.log(P)).sum(1)
    return out


def entropy_matched_random(W: np.ndarray, avail: np.ndarray, seed: int) -> np.ndarray:
    """The random control: for every target the arm touches, draw a weight row with the SAME entropy as the arm's
    row but a random assignment of which donor gets which weight (J = 2: the same {w, 1-w} pair, randomly ordered;
    J > 2: the arm's own sorted weight vector randomly permuted).  Entropy is matched exactly by construction."""
    rng = np.random.default_rng(seed)
    out = np.full_like(W, np.nan)
    for t in np.flatnonzero(np.isfinite(W).all(1)):
        js = np.flatnonzero(avail[:, t])
        w = np.asarray(W[t, js], float)
        out[t] = 0.0
        out[t, js] = w[rng.permutation(w.size)]
    return out


def shuffle_kd(A: np.ndarray, M: np.ndarray, seed: int) -> np.ndarray:
    """The `:shuf` control: each donor's |kd| vector permuted across ITS measured targets (within donor)."""
    rng = np.random.default_rng(seed)
    out = A.copy()
    for j in range(A.shape[0]):
        idx = np.flatnonzero(M[j])
        out[j, idx] = A[j, idx][rng.permutation(idx.size)]
    return out


def normrel_weights(rows: list[np.ndarray], avail: np.ndarray, wcols: np.ndarray, power: float = -1.5) -> np.ndarray:
    """097 route (a) a6: w_j ∝ (||row_j(t)||_wcols / median_t ||row_j(t)||)^power."""
    J, n = avail.shape
    K = np.zeros((J, n))
    for j in range(J):
        nj = np.linalg.norm(rows[j][:, wcols], axis=1)
        med = float(np.median(nj[avail[j]])) or 1.0
        K[j] = nj / med
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(avail.sum(0) >= 2):
        js = np.flatnonzero(avail[:, t])
        k = np.maximum(K[js, t], 1e-9) ** power
        W[t] = 0.0
        W[t, js] = k / k.sum()
    return W
