"""285 -- the shared response object, the two readouts, the nested ladder design and the bootstrap.

Nothing mathematical is new.  The profile / residual arithmetic is `truth273.py`'s and `lib280.py`'s; the
per-target rank statistic is 218's `_per_target`, imported verbatim (273 KA-2: == cell-eval2's `pds_cosine`
to 5.6e-17).
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/285-reliability-benchmark"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
# The scorer's per-target rank statistic, taken VERBATIM out of 218's `resid218.py` (273's KA-2 proved that
# function equals cell-eval2 0.16.0's `pds_cosine` on a cloud-scored artifact to 5.6e-17).  The function text is
# extracted and exec'd on its own so that none of 218's module-level code (plolo, fold caches) runs here.
_R218 = ROOT / "experiments/218-independent-donors/resid218.py"
_SRC = _R218.read_text()
_I = _SRC.index("def _per_target(")
_J = _SRC.index("\ndef ", _I + 1)
PER_TARGET_SRC = _SRC[_I:_J]
PER_TARGET_SHA = hashlib.sha256(PER_TARGET_SRC.encode()).hexdigest()
_NS: dict = {"np": np}
exec(PER_TARGET_SRC, _NS)                  # noqa: S102 -- verbatim archive code, hash recorded
per_target = _NS["_per_target"]

TS = 5e4
SEED = 20260906
NBOOT = 20000
RUNGS = (25, 50, 100, 200, 400)
NDRAW = 8


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------ the response object
def prof(S: np.ndarray) -> np.ndarray:
    """log1p(TS * S / sum S), row-wise (2-D) or on a vector."""
    S = np.asarray(S, dtype=np.float64)
    if S.ndim == 1:
        t = S.sum()
        return np.log1p(TS * S / t) if t > 0 else np.zeros_like(S)
    t = S.sum(1, keepdims=True)
    return np.log1p(TS * S / np.where(t > 0, t, 1.0))


def delta(S: np.ndarray, ctrl_S: np.ndarray | None) -> np.ndarray:
    """The view's effect rows: profile minus the view's own control profile (or zero when no control exists;
    the residual below is algebraically independent of that choice)."""
    P = prof(S)
    return P - (prof(ctrl_S)[None, :] if ctrl_S is not None else 0.0)


def resid(D: np.ndarray) -> np.ndarray:
    """159's panel-centred residual: minus the mean over the set's targets, on the SAME side."""
    return D - D.mean(0, keepdims=True)


def row_cos(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


def cos_agg(P: np.ndarray, T: np.ndarray) -> float:
    a = float((P * P).sum()); x = float((P * T).sum()); d = float((T * T).sum())
    return float(np.sqrt(max(0.0, x * x / (a * d)))) if a > 0 and d > 0 else 0.0


def views(SA, SB, cA=None, cB=None, scored=None, mA=None, mB=None):
    """Both readouts for one pair of views.  `scored` = the column mask (target genes already zeroed/removed).

    `mA` / `mB` (optional, in delta space on the scored columns) replace the side's own mean as the centring
    vector.  Use them whenever the rows passed are a SUBSET of the source's targets: centring a 4-row matrix on
    its own mean leaves almost nothing, and centring it on the full panel's mean (273's `halves_sub`
    convention) keeps the residual space the same object across units.

    Returns dict(cos_resid, pds, cos_own, cos_agg): `cos_resid` per target on the panel-centred residual;
    `pds` = mean of pds(A->B) and pds(B->A) with every target PASSED as a distractor; `cos_own` the
    per-target cosine on the UN-centred effects (what the scorer reads)."""
    DA, DB = delta(SA, cA), delta(SB, cB)
    if scored is not None:
        DA, DB = DA[:, scored], DB[:, scored]
    RA = resid(DA) if mA is None else DA - np.asarray(mA)[None, :]
    RB = resid(DB) if mB is None else DB - np.asarray(mB)[None, :]
    cr = row_cos(RA, RB)
    ok = np.isfinite(cr) & (np.linalg.norm(DA, axis=1) > 0) & (np.linalg.norm(DB, axis=1) > 0)
    pab = np.full(len(DA), np.nan); pba = np.full(len(DA), np.nan); cown = np.full(len(DA), np.nan)
    if ok.sum() >= 3:
        p1, c1 = per_target(DA[ok], DB[ok])
        p2, _ = per_target(DB[ok], DA[ok])
        pab[ok], pba[ok], cown[ok] = p1, p2, c1
    return dict(cos_resid=cr, pds=(pab + pba) / 2.0, pds_ab=pab, pds_ba=pba, cos_own=cown,
                cos_agg=cos_agg(DA[ok], DB[ok]) if ok.sum() >= 3 else np.nan, n_ok=int(ok.sum()))


def thin(S: np.ndarray, frac: float, rng) -> np.ndarray:
    """Binomial thinning of integer count sums (lib280.thin, vectorised over rows)."""
    if frac >= 1.0:
        return np.asarray(S, dtype=np.float64)
    return rng.binomial(np.rint(np.asarray(S)).astype(np.int64), frac).astype(np.float64)


# ------------------------------------------------------------------------------------ the nested ladder design
def ladder_blocks(rungs=RUNGS) -> list[tuple[int, int]]:
    """Nested rank blocks: side A of rung n = ranks [0, n), side B = [n, 2n).  One random ranking per draw
    delivers every rung, so the rungs of one draw are nested and share cells across rungs by construction."""
    edges = [0]
    for n in rungs:
        edges += [n, 2 * n]
    edges = sorted(set(edges))
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def sides_from_blocks(block_sums: dict[tuple[int, int], np.ndarray], n: int) -> tuple[np.ndarray, np.ndarray]:
    """Sum the nested blocks into side A = [0, n) and side B = [n, 2n)."""
    A = sum(v for (lo, hi), v in block_sums.items() if hi <= n)
    B = sum(v for (lo, hi), v in block_sums.items() if lo >= n and hi <= 2 * n)
    return A, B


def rank_assign(counts: np.ndarray, rng) -> list[np.ndarray]:
    """Per group, a random permutation of its member positions (returns the rank of each member)."""
    return [rng.permutation(c) for c in counts]


# ------------------------------------------------------------------------------------ statistics
def stats(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=np.float64)
    v = x[np.isfinite(x)]
    if v.size == 0:
        return dict(n=0)
    return dict(n=int(v.size), median=float(np.median(v)), mean=float(v.mean()),
                p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                share_le0=float((v <= 0).mean()), share_ge005=float((v >= 0.05).mean()))


def boot_median(x: np.ndarray, n: int = NBOOT, seed: int = SEED):
    v = np.asarray(x, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return dict(point=np.nan, lo=np.nan, hi=np.nan, n=0)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, v.size, size=(n, v.size))
    d = np.median(v[idx], axis=1)
    return dict(point=float(np.median(v)), lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)), n=int(v.size))


def boot_median_diff(x: np.ndarray, y: np.ndarray, n: int = NBOOT, seed: int = SEED):
    """PAIRED target bootstrap of median(x) - median(y) on the targets where both are finite."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size == 0:
        return dict(point=np.nan, lo=np.nan, hi=np.nan, n=0, p_le0=np.nan)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n, x.size))
    d = np.median(x[idx], axis=1) - np.median(y[idx], axis=1)
    return dict(point=float(np.median(x) - np.median(y)), lo=float(np.percentile(d, 2.5)),
                hi=float(np.percentile(d, 97.5)), n=int(x.size), p_le0=float((d <= 0).mean()))


def spearman_brown(r: float, k: float) -> float:
    """Reliability of a k-times-longer measurement."""
    return float(k * r / (1.0 + (k - 1.0) * r)) if np.isfinite(r) else np.nan


def fmt(b: dict) -> str:
    return f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]"


class Logger:
    def __init__(self, path: Path):
        self.path, self.lines, self.t0 = path, [], time.time()

    def __call__(self, s: str = "") -> None:
        print(s, flush=True)
        self.lines.append(s)

    def save(self) -> None:
        self.path.write_text("\n".join(self.lines) + "\n")

    def dump(self, obj, path: Path) -> None:
        path.write_text(json.dumps(obj, indent=1, default=float))
