"""277 -- shared construction: the low-conservation CLASS and the class-conditional renderings.

Nothing mathematical is new.  Everything is imported:
  * `lib257` -> `lib246` -> `plolo` / `resid218` / `compose218`   (162's PLOLO, 218's residuals, the shipped
    stack, `deflate`, `compose`, `score`, `solve_fraction`, 188's exact `nmae` identity via `lib257.Identity`).
  * `op272.DampedOperator`                                        (272's A' eigenvalue floor max(lambda, c n)).
  * `boot265.boot` / `boot_folds`                                 (265's paired target bootstrap).
  * `experiments/276-conservation-map/work/conservation276.csv`   (276's per-target `C_lolo_<fold>`).

The ONLY new objects here are (a) the CLASS selector, (b) `mix_rows` -- a CD4 residual whose rows come from
two prepared residuals by a boolean mask, and (c) `mixed_whiten` -- one whitening in which the class rows
carry the damped eigenvalue weights and every other row the incumbent weights, with production's SINGLE
on-axis preserve.  Frozen in `work/DEFINITION-277.md` (sha `05a89391...`) before anything was scored.

Read-only everywhere except `experiments/277-conservation-class-rendering/work/` (the PLOLO fold caches
under `work/plolo/` are byte copies of 262's, `shasum -c` 7/7 OK).
"""
from __future__ import annotations

import csv
import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E257 = ROOT / "experiments/257-soft-deflation-price"
E262 = ROOT / "experiments/262-where-pds-is-lost"
E272 = ROOT / "experiments/272-noise-corrected-operator"
E273 = ROOT / "experiments/273-flex-reference-instrument"
E276 = ROOT / "experiments/276-conservation-map"
E265 = ROOT / "experiments/265-majority-permutation-test"

for p in (E257, E272, E265):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import lib257 as X  # noqa: E402  (imports lib246 as X.L -> plolo)
import plolo  # noqa: E402
import op272 as O  # noqa: E402
from boot265 import boot, boot_folds  # noqa: E402

L = X.L
# every cache read comes from THIS experiment's byte copies; nothing outside it is touched
plolo.WORK = WORK / "plolo"
L.WORK = WORK
X.WORK = WORK
E062, E064 = plolo.E062, plolo.E064

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
CD4_F = L.CD4_F                     # 0.34
GAMMA, LAM = L.GAMMA, L.LAM
SEED = 20260906
B_BOOT = 2000
CSV276 = E276 / "work/conservation276.csv"
KEYS = ("cos_agg", "a23_cos", "pds")


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def row_cos(A, B):
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


# ------------------------------------------------------------------------------- the CLASS
def c_lolo(fold: str, targets) -> np.ndarray:
    """276's leave-the-held-line-out conservation score for this fold's targets; NaN where 276 has none."""
    col = f"C_lolo_{fold}"
    tab: dict[str, float] = {}
    with open(CSV276) as fh:
        for r in csv.DictReader(fh):
            v = r[col]
            if v not in ("", "nan"):
                tab[r["target"]] = float(v)
    return np.asarray([tab.get(str(t), np.nan) for t in targets], dtype=np.float64)


def bottom_class(score: np.ndarray, targets, frac: float, *, top: bool = False) -> np.ndarray:
    """DEFINITION-277 section 1: the floor(frac * n_eligible) eligible rows with the SMALLEST (or, with
    top=True, the LARGEST) score; ties broken by ascending target symbol.  Ineligible rows (NaN) are never
    in the class."""
    elig = np.flatnonzero(np.isfinite(score))
    n_class = int(np.floor(frac * elig.size))
    key = [(-score[i] if top else score[i], str(targets[i])) for i in elig]
    order = sorted(range(elig.size), key=lambda j: key[j])
    sel = elig[np.asarray(order[:n_class], dtype=np.int64)]
    m = np.zeros(score.size, dtype=bool)
    m[sel] = True
    return m


def random_class(score: np.ndarray, n_class: int, seed: int) -> np.ndarray:
    elig = np.flatnonzero(np.isfinite(score))
    rng = np.random.default_rng(seed)
    sel = rng.choice(elig, size=n_class, replace=False)
    m = np.zeros(score.size, dtype=bool)
    m[sel] = True
    return m


# ------------------------------------------------------------------------- class-conditional renderings
def mix_rows(res_class, res_rest, mask: np.ndarray, tag: str):
    """A CD4 residual whose rows come from `res_class` where `mask` and from `res_rest` elsewhere.  The two
    inputs must share the support (`lib246._as_res` copies votes / counts through `deflate`, so `cd4` and
    its deflations do)."""
    a = np.asarray(res_class.mean_effect, dtype=np.float64)
    b = np.asarray(res_rest.mean_effect, dtype=np.float64)
    va = np.asarray(res_class.gene_votes)
    vb = np.asarray(res_rest.gene_votes)
    assert a.shape == b.shape and np.array_equal(va, vb), "mix_rows: inputs do not share the support"
    eff = np.where(mask[:, None], a, b)
    return E062.CD4Prepared(mode=f"{res_rest.mode}|mix:{tag}", mean_effect=eff,
                            gene_votes=vb.copy(),
                            target_condition_count=np.asarray(res_rest.target_condition_count).copy(),
                            condition_details=res_rest.condition_details)


def zero_rows(res, mask: np.ndarray, tag: str):
    """The residual with the masked rows removed entirely (effect, votes AND count zeroed)."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes).copy()
    cnt = np.asarray(res.target_condition_count).copy()
    eff[mask] = 0.0
    votes[mask] = 0
    cnt[mask] = 0
    return E062.CD4Prepared(mode=f"{res.mode}|zero:{tag}", mean_effect=eff, gene_votes=votes,
                            target_condition_count=cnt, condition_details=res.condition_details)


def mixed_whiten(M, wcols, evecs, w_class: np.ndarray, w_rest: np.ndarray, mask: np.ndarray):
    """`plolo.Operator.apply`'s expression with a PER-ROW weight choice and production's SINGLE on-axis
    preserve: rows in `mask` are transformed with `w_class`, every other row with `w_rest`, then one scalar
    restores the median on-axis row norm.  `mask.sum() == 0` (or `w_class is w_rest`) is `Operator.apply`
    bit for bit."""
    out = np.asarray(M, dtype=np.float64).copy()
    w = np.asarray(wcols, dtype=np.int64)
    norms_in = np.linalg.norm(out[:, w], axis=1)
    before = float(np.median(norms_in))
    assert before > 0.0
    C = out[:, w] @ evecs                                   # (n_rows, n_cols) in the eigenbasis
    Cw = np.where(mask[:, None], C * w_class[None, :], C * w_rest[None, :])
    t = Cw @ evecs.T
    after = float(np.median(np.linalg.norm(t, axis=1)))
    assert after > 0.0
    out[:, w] = t * (before / after)
    return out


# ------------------------------------------------------------------------------------- groups + retention
def groups_262(base_cos: np.ndarray) -> dict[str, np.ndarray]:
    return {"e (base.truth <= 0)": base_cos <= 0,
            "marginal (0, 0.05]": (base_cos > 0) & (base_cos <= 0.05),
            "saturated > 0.05": base_cos > 0.05}


def lt3o(diff: np.ndarray) -> dict:
    """Leave-top-3-out: the fold-mean per-target difference recomputed after deleting the three targets
    with the largest positive difference.  Retention = the LT3O mean / the full mean."""
    d = np.asarray(diff, dtype=np.float64)
    full = float(d.mean())
    if d.size <= 3:
        return dict(full=full, lt3o=float("nan"), retention=float("nan"))
    drop = np.argsort(d)[::-1][:3]
    keep = np.ones(d.size, dtype=bool)
    keep[drop] = False
    cut = float(d[keep].mean())
    return dict(full=full, lt3o=cut, retention=float(cut / full) if full != 0 else float("nan"),
                dropped_idx=[int(i) for i in drop], dropped_delta=[float(d[i]) for i in drop])


def fmt(b: dict) -> str:
    return f"{b['point']:+.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}]"
