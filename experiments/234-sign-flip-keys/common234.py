"""234 — shared machinery for the SIGN-FLIP KEY screen.

THE QUESTION, AS ARITHMETIC
---------------------------
232 measured the exchange rate: **+0.01 of |t|-weighted sign accuracy `q_w` is worth
+0.00399 overall** (k_cos 1.384; +0.00379 at k_cos 1.000).  Rule 13's gate is +0.0019, so a
sign arm must deliver **Dq_w = +0.0048 ... +0.0050**.

A SIGN-FLIP arm picks a set `S` of (target, gene) coordinates by a CONFIDENCE KEY and
replaces `effect -> |effect| * (-sign(effect))` there.  Magnitudes are bit-identical, so
§8.1d classifies it as a sign arm and 232 §11(a) pins `jac` at exactly zero.  Its effect on
`q_w` is EXACT and needs no scoring run:

    Dq_w = mean_i [ sum_{j in pool_i & S_i} w_ij (1 - 2 c_ij) / sum_{j in pool_i} w_ij ]

with `w = |truth|` and `c = 1[our delivered sign is right]`.  So a flip pays **iff the key
selects a set whose sign accuracy is below 0.5**, and the whole screen is a question about
one number per key: how far below 0.5 can it push a set that still carries mass?

This module extracts, per PLOLO fold, the truth / pool / delivered-sign / key tensors that
the screen needs, and nothing else.  All heavy objects (the 17,942^2 whitener) live only
inside `extract234.py`, one fold per process.

KEYS.  Everything already on disk that is per-(target, gene) or broadcasts to it:
  absH        |delivered effect| -- 220's best key, corr +0.0743 with correctness
  absfin      |final effect| (pre-render)
  abspre      |consensus effect| (pre-whitening)
  key189      |H| * exp(0.5 * theta) -- 189's shipped ordering key
  theta       189's per-GENE reliability prior (broadcast)
  expr        control log-expression (broadcast) / cpm
  zpred       |final| / sd_across_targets(final) -- our own prediction standardised
  absH_rank   within-target rank fraction of |H| (scale-free version of absH)
  ncov        number of covering donors
  margin_maj  |pos-neg|/(pos+neg) -- at arity 2 this IS "unanimous vs contested"
  margin_z    232's z-vote normalised margin
  coh         |sum_s e_s| / sum_s |e_s| -- donor coherence
  tdonor      |mean_s e_s| / sem_s(e_s) -- across-donor t
  minz/sumz   min / sum of |e_s|/sd_s(gene) over covering donors
  whitenkept  1[sign(pre) == sign(final)] -- did the operator rewrite this sign?
  rownorm     per-TARGET ||H|| (a per-target key still moves a macro-averaged q_w)
  nl_targ     per-TARGET number of live genes
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
WORK.mkdir(parents=True, exist_ok=True)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALL_FOLDS = CANCER + ("kolf2",)

# ---- 232's exchange rate (RESULTS.md §5).  The ONLY PLOLO -> overall conversion used here.
RATE_PER_QW = {"k1.384": 0.003990 / 0.01, "k1.000": 0.003785 / 0.01}
GATE = 0.0019
QW_NEEDED = {k: GATE / v for k, v in RATE_PER_QW.items()}      # +0.00476 / +0.00502

# ---- 232's own measured rows, used as known-answer targets
KA232_EXCHANGE = dict(composite_k1384=0.003990, composite_k1000=0.003785,
                      d_nmae=-0.002573, r_nmae=-0.995, r_composite=0.975,
                      o_pds=0.002434, o_nmae=0.000719)
KA232_CEILING = dict(disagree_mass=0.3627, q_contested=0.52292, ceiling=0.17325)
# 232's BASE rows in `work/arm232.json` -- the delivered |t|-weighted sign accuracy on
# `wald89` that the whole exchange rate is regressed against.  (NOT §2b's table, which is
# `sign232.py`'s consensus-level statistic and carries the uncovered-target fallback rows;
# the two differ by ~0.002 and only THIS one is the object 232 priced.)
KA232_QW_DELIVERED = dict(k562_gwps=0.5762383750822039, orion_hct116=0.56494474812357,
                          orion_hek293t=0.5529074867709532, kolf2=0.516146835399458)
KA232_NMAE_BASE = dict(k562_gwps=0.954307, orion_hct116=0.995679, orion_hek293t=0.996245)
KA232_PDS_BASE = dict(k562_gwps=0.805397, orion_hct116=0.792213, orion_hek293t=0.763623)
KA228_E = 0.072
KA220_CORR = dict(absH=0.0743, theta=0.002, expr=0.0080)
KA220_AUC = dict(absH=0.5350, expr=0.5049)

SEED = 20260904


def qw_to_overall(dqw: float, k: str = "k1.384") -> float:
    """232 §5's exchange rate.  PLOLO units in, overall points out."""
    return float(dqw) * RATE_PER_QW[k]


# ------------------------------------------------------------------ the exact Dq_w algebra
def prep_weights(truth: np.ndarray, pool: np.ndarray):
    """Per-target |t| weights normalised inside the pool, and the correctness bit.

    Returns (`w`, `den`) where `w[i, j] = |truth[i, j]| * pool[i, j]` and `den[i]` is the
    row's pool weight.  Rows with zero pool weight are dropped by `qw` / `dqw`.
    """
    w = np.where(pool, np.abs(truth), 0.0)
    den = w.sum(1)
    return w, den


def qw(sgn: np.ndarray, truth: np.ndarray, w: np.ndarray, den: np.ndarray) -> float:
    """|t|-weighted sign accuracy, macro-averaged over targets.

    232 `weighted_sign_accuracy`'s convention exactly: a ZERO predicted sign is never
    correct, and a target with no pool weight is dropped from the mean.
    """
    ok = (sgn != 0) & (sgn == np.sign(truth))
    live = den > 0
    return float(np.mean((w * ok).sum(1)[live] / den[live]))


def dqw(flip: np.ndarray, correct: np.ndarray, w: np.ndarray, den: np.ndarray) -> float:
    """EXACT change in `qw` when the coordinates in `flip` have their sign reversed.

    A flipped correct coordinate loses its weight; a flipped incorrect one gains it.  A
    coordinate whose delivered sign is ZERO gains weight when flipped only if the flip
    gives it the truth's sign -- but `-0 == 0`, so a zero stays zero and gains nothing.
    That is why `correct` must be the strict bit and `zero_sign` is excluded upstream.
    """
    live = den > 0
    num = (w * flip * (1.0 - 2.0 * correct)).sum(1)
    return float(np.mean(num[live] / den[live]))


def mass_share(sel: np.ndarray, w: np.ndarray, den: np.ndarray) -> float:
    """Macro-averaged |t|-weighted share of the pool that `sel` covers."""
    live = den > 0
    return float(np.mean((w * sel).sum(1)[live] / den[live]))


def acc_in(sel: np.ndarray, correct: np.ndarray, w: np.ndarray) -> float:
    """|t|-weighted sign accuracy INSIDE `sel`, pooled (not macro) -- the diagnostic that
    says whether a flip can pay at all: it pays iff this is below 0.5."""
    num = float((w * sel * correct).sum())
    den = float((w * sel).sum())
    return num / den if den > 0 else float("nan")


# ------------------------------------------------------------------------- dose machinery
def row_lowest(key: np.ndarray, frac: float, allowed: np.ndarray | None = None):
    """Select the `frac` lowest-key coordinates in each row (optionally within `allowed`).

    The selection uses only the KEY, never the pool or the truth -- the pool is the
    reference's DE set and is not knowable at prediction time.
    """
    T, G = key.shape
    out = np.zeros((T, G), dtype=bool)
    for i in range(T):
        idx = np.flatnonzero(allowed[i]) if allowed is not None else np.arange(G)
        if idx.size == 0:
            continue
        k = int(round(frac * idx.size))
        if k <= 0:
            continue
        k = min(k, idx.size)
        sub = idx[np.argsort(key[i, idx], kind="stable")[:k]]
        out[i, sub] = True
    return out


def row_rank_frac(key: np.ndarray) -> np.ndarray:
    """Within-row rank of `key` mapped to [0, 1) -- a scale-free version of any key."""
    order = np.argsort(key, axis=1, kind="stable")
    r = np.empty_like(order)
    np.put_along_axis(r, order, np.arange(key.shape[1])[None, :], axis=1)
    return r.astype(np.float32) / key.shape[1]
