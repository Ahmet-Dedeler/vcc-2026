"""303 — shared, frozen definitions for claim C7 (second overexpression library).

PRE-ACCESS SAFETY CONTRACT
--------------------------
Nothing in this module reads `X` or `layers/control` for any reserve (R5 =
the 510 `Nakatake2020_hESC_TF_screen` rows).  The only row-value reads it
performs are

  * `layers/measure_mask` (any row) — explicitly free per the R5 declaration;
  * `X` / `layers/control` rows whose `obs.dataset == 'Pan2025_hESC'` — the
    already-used bank source, explicitly free per the R5 declaration.

`read_source_rows()` enforces that with an allow-list argument.  The single
reserve read lives in `access303.py` and is gated by a coordinator release
JSON; this module cannot perform it.

The arithmetic here is 115's `extract_compact_donor.py` transform, restated so
the Nakatake rows and the `gse295214_h1` bank rows sit in ONE convention:

    source X          = log1p(CPM10k perturbed) - log1p(CPM10k control)   (per row)
    layers/control    = log1p(CPM10k control)                             (per row)
    perturbed_log     = control + X
    *_cpm             = expm1(max(*_log, 0))                              -> sums to 10_000
    collapse          = sum duplicate source symbols in LINEAR cpm space
                        (a canonical gene is measured iff EVERY source column
                         mapping onto it is measured)
    delta             = log1p(5 * collapsed_perturbed) - log1p(5 * collapsed_control)
    b_ctrl            = log1p(5 * collapsed_control)

`delta` is the bank's `delta.npy` semantics ("TS50000 log1p perturbation
effect"); coordinates outside the row's measured support are written as exact
zero by 115 and are NOT expression values.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/303-second-overexpression-library"
WORK = HERE / "work"

# ---------------------------------------------------------------- inputs ----
DB = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
DB_SHA = "0754d6206fe045744dfd9b21eb37c8e3435152c76b158b9810a84c64a029c0ac"  # 115 meta.json
DB_BYTES = 4_425_260_319

H1BANK = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
CANON = ROOT / "data/challenge/gene_names.csv"
HELD279 = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
ZNF279 = ROOT / "experiments/279-znf-state-dependence/work/znf279.csv"
PERSRC290 = ROOT / "experiments/290-remainder-response-features/work/persrc290.npz"
SYM302 = ROOT / "experiments/302-orf-vs-crispri-symmetry/work/symmetry302.csv"
SUMS302 = ROOT / "experiments/302-orf-vs-crispri-symmetry/work/sums302.npz"

RESERVE_DATASET = "Nakatake2020_hESC_TF_screen"
BANK_DATASET = "Pan2025_hESC"
JOUNG_DATASET = "Joung2023_hESC_TF_screen"
FREE_DATASETS = (BANK_DATASET,)          # value rows this module may read

SOURCE_TARGET_SUM = 10_000.0
PROJECT_TARGET_SUM = 50_000.0
TARGET_SUM_ATOL = 1.0

SEED = 20260907
N_BOOT = 20_000
N_NULL = 2_000       # target-shuffled pairing draws (>= the draft's 200)
BLOCK = 64           # rows per h5py read block


class Refuse(RuntimeError):
    """A frozen precondition failed.  Never caught inside a measurement path."""


def sha256_file(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_strings(values) -> str:
    """Order-sensitive hash of a list of strings (newline-joined, utf-8)."""
    return sha256_bytes(("\n".join(str(v) for v in values) + "\n").encode())


def sha256_array(a: np.ndarray) -> str:
    a = np.ascontiguousarray(a)
    return sha256_bytes(a.tobytes() + f"|{a.dtype.str}|{a.shape}".encode())


# --------------------------------------------------------------- metadata ---
@dataclass(frozen=True)
class DbMeta:
    n_obs: int
    n_var: int
    dataset: np.ndarray        # (n_obs,) str
    obs_id: np.ndarray         # (n_obs,) str
    pert_symbol: np.ndarray    # (n_obs,) str  = obs.perturbed_gene_name
    pert_gene: np.ndarray      # (n_obs,) str  = obs.perturbed_gene (ENSG)
    condition: np.ndarray      # (n_obs,) str
    sign: np.ndarray           # (n_obs,) int64
    var_id: np.ndarray         # (n_var,) str  (var index, 'human_id')
    var_symbol: np.ndarray     # (n_var,) str  (var.gene_name)


def read_db_meta(path: Path = DB) -> DbMeta:
    """obs/var only.  No X, no control, no mask."""
    with h5py.File(path, "r") as f:
        def cat(name: str) -> np.ndarray:
            c = np.asarray(f[f"obs/{name}/categories"][:]).astype(str)
            return c[f[f"obs/{name}/codes"][:]]
        n_obs, n_var = [int(x) for x in f["X"].shape]
        m = DbMeta(
            n_obs=n_obs, n_var=n_var,
            dataset=cat("dataset"),
            obs_id=np.asarray(f["obs/id"][:]).astype(str),
            pert_symbol=cat("perturbed_gene_name"),
            pert_gene=cat("perturbed_gene"),
            condition=cat("condition"),
            sign=f["obs/perturbation_sign"][:].astype(np.int64),
            var_id=np.asarray(f["var/human_id"][:]).astype(str),
            var_symbol=np.asarray(
                np.asarray(f["var/gene_name/categories"][:]).astype(str)[
                    f["var/gene_name/codes"][:]]),
        )
    return m


def read_bank() -> dict:
    perts = [x.strip() for x in (H1BANK / "perts.txt").read_text().splitlines() if x.strip()]
    genes = [x.strip() for x in (H1BANK / "genes.txt").read_text().splitlines() if x.strip()]
    return dict(perts=perts, genes=genes, row_of={p: i for i, p in enumerate(perts)},
                col_of={g: i for i, g in enumerate(genes)})


def canonical_genes() -> list[str]:
    lines = CANON.read_text().splitlines()
    return [x.strip() for x in lines[1:] if x.strip()]


# ------------------------------------------------------------- projection ---
@dataclass(frozen=True)
class Projection:
    """115's source->canonical gene projection, rebuilt from var metadata."""
    canonical: tuple[str, ...]
    single_src: np.ndarray
    single_dst: np.ndarray
    dup_groups: tuple[tuple[int, tuple[int, ...]], ...]
    covered: np.ndarray

    @property
    def n_canonical(self) -> int:
        return len(self.canonical)


def build_projection(var_symbol: np.ndarray, canonical: list[str]) -> Projection:
    pos = {g: i for i, g in enumerate(canonical)}
    dst = np.asarray([pos.get(s, -1) for s in var_symbol], np.int32)
    grouped: dict[int, list[int]] = {}
    for i, d in enumerate(dst):
        if d >= 0:
            grouped.setdefault(int(d), []).append(i)
    singles, dups = [], []
    for d, srcs in sorted(grouped.items()):
        (singles if len(srcs) == 1 else dups).append((srcs[0], d) if len(srcs) == 1 else (d, tuple(srcs)))
    covered = np.zeros(len(canonical), np.bool_)
    covered[np.fromiter(grouped, np.int64)] = True
    return Projection(
        canonical=tuple(canonical),
        single_src=np.asarray([s for s, _ in singles], np.int64),
        single_dst=np.asarray([d for _, d in singles], np.int64),
        dup_groups=tuple(dups),
        covered=covered,
    )


def collapse_linear(values: np.ndarray, measured: np.ndarray, proj: Projection):
    """115's `_collapse_linear`: sum duplicates in linear space; AND their masks."""
    if values.shape != measured.shape or values.ndim != 2:
        raise Refuse("linear values and measurement mask disagree")
    rows = values.shape[0]
    out = np.zeros((rows, proj.n_canonical), np.float64)
    msk = np.zeros((rows, proj.n_canonical), np.bool_)
    if proj.single_src.size:
        out[:, proj.single_dst] = values[:, proj.single_src]
        msk[:, proj.single_dst] = measured[:, proj.single_src]
    for d, srcs in proj.dup_groups:
        s = np.asarray(srcs, np.int64)
        out[:, d] = values[:, s].sum(axis=1, dtype=np.float64)
        msk[:, d] = measured[:, s].all(axis=1)
    return out, msk


def collapse_mask_only(measured: np.ndarray, proj: Projection) -> np.ndarray:
    """Mask-only projection — the FREE preflight path (no X, no control)."""
    if measured.ndim != 2:
        raise Refuse("mask must be 2-D")
    msk = np.zeros((measured.shape[0], proj.n_canonical), np.bool_)
    if proj.single_src.size:
        msk[:, proj.single_dst] = measured[:, proj.single_src]
    for d, srcs in proj.dup_groups:
        msk[:, d] = measured[:, np.asarray(srcs, np.int64)].all(axis=1)
    return msk


# ------------------------------------------------------------ the transform -
def bank_transform(effect: np.ndarray, control_log: np.ndarray, mask: np.ndarray,
                   proj: Projection, labels, strict_sums: bool = True):
    """115's arithmetic, verbatim, for an arbitrary row block.

    Returns (delta_canonical float32, measured_canonical bool, diagnostics).
    `strict_sums=False` records the CPM10k sum drift instead of raising — used
    only to REPORT a precondition failure for a dataset 115 never validated.
    """
    if not (effect.shape == control_log.shape == mask.shape):
        raise Refuse("X / control / measure_mask block shapes disagree")
    for name, a in (("X", effect), ("control", control_log), ("measure_mask", mask)):
        if not np.all(np.isfinite(a)):
            raise Refuse(f"{name} contains nonfinite values near {labels[0]}")
    if np.any((mask != 0.0) & (mask != 1.0)):
        raise Refuse(f"measure_mask is not exactly binary near {labels[0]}")
    measured_src = mask.astype(np.bool_, copy=False)
    perturbed_log = control_log + effect
    min_c, min_p = float(np.min(control_log)), float(np.min(perturbed_log))
    if min_c < -1e-6 or min_p < -1e-6:
        raise Refuse(f"log1p values imply negative CPM near {labels[0]}: "
                     f"control={min_c}, perturbed={min_p}")
    c_cpm = np.expm1(np.maximum(control_log, 0.0), dtype=np.float64)
    p_cpm = np.expm1(np.maximum(perturbed_log, 0.0), dtype=np.float64)
    c_sum = c_cpm.sum(axis=1, dtype=np.float64)
    p_sum = p_cpm.sum(axis=1, dtype=np.float64)
    drift = dict(control_min=float(np.min(c_sum)), control_max=float(np.max(c_sum)),
                 perturbed_min=float(np.min(p_sum)), perturbed_max=float(np.max(p_sum)))
    bad = (np.abs(c_sum - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL) | \
          (np.abs(p_sum - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL)
    if strict_sums and bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise Refuse(f"CPM10k sum drift at {labels[i]}: control={c_sum[i]:.6f} "
                     f"perturbed={p_sum[i]:.6f}")
    cp, m1 = collapse_linear(p_cpm, measured_src, proj)
    cc, m2 = collapse_linear(c_cpm, measured_src, proj)
    if not np.array_equal(m1, m2):
        raise Refuse("perturbed and control masks disagree after collapse")
    k = PROJECT_TARGET_SUM / SOURCE_TARGET_SUM
    delta = (np.log1p(k * cp) - np.log1p(k * cc)).astype(np.float32)
    delta[~m1] = 0.0
    return delta, m1, dict(drift=drift, n_bad_rows=int(bad.sum()))


def control_profile(control_log_row: np.ndarray, proj: Projection) -> np.ndarray:
    """115's `b_ctrl`: log1p(5 * collapsed control cpm) from ONE control row."""
    c = np.expm1(np.maximum(np.asarray(control_log_row, np.float64), 0.0), dtype=np.float64)
    s = float(c.sum())
    if abs(s - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL:
        raise Refuse(f"reference control CPM10k sum drift: {s:.9f}")
    coll = np.zeros(proj.n_canonical, np.float64)
    if proj.single_src.size:
        coll[proj.single_dst] = c[proj.single_src]
    for d, srcs in proj.dup_groups:
        coll[d] = c[np.asarray(srcs, np.int64)].sum(dtype=np.float64)
    return np.log1p((PROJECT_TARGET_SUM / SOURCE_TARGET_SUM) * coll).astype(np.float32)


# ------------------------------------------------------------- safe reader ---
def read_source_rows(path: Path, rows: np.ndarray, datasets_allowed: tuple[str, ...],
                     dataset_of: np.ndarray, what: tuple[str, ...] = ("X", "layers/control",
                                                                      "layers/measure_mask")):
    """Row-blocked read with a dataset allow-list.

    Refuses outright if any requested row belongs to a dataset not on the
    allow-list.  `layers/measure_mask` is free for every dataset (the R5
    declaration), so it is exempted from the allow-list check.
    """
    rows = np.asarray(rows, np.int64)
    if rows.ndim != 1 or rows.size == 0:
        raise Refuse("row index list must be a non-empty 1-D array")
    if np.any(np.diff(rows) <= 0):
        raise Refuse("row index list must be strictly increasing")
    value_reads = [w for w in what if w != "layers/measure_mask"]
    if value_reads:
        bad = sorted(set(dataset_of[rows]) - set(datasets_allowed))
        if bad:
            raise Refuse(f"read of {value_reads} refused for dataset(s) {bad}; "
                         f"allowed: {list(datasets_allowed)}")
    out = {w: [] for w in what}
    with h5py.File(path, "r") as f:
        for lo in range(0, rows.size, BLOCK):
            sel = rows[lo:lo + BLOCK]
            for w in what:
                out[w].append(np.asarray(f[w][sel, :], np.float64))
    return {w: np.concatenate(v, axis=0) for w, v in out.items()}


# ------------------------------------------------------------- statistics ---
def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na == 0.0 or nb == 0.0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    from scipy.stats import rankdata
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan")
    a, b = rankdata(x[ok]), rankdata(y[ok])
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else float("nan")


def partial_spearman(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> float:
    from scipy.stats import rankdata
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    if ok.sum() < Z.shape[1] + 3:
        return float("nan")
    R = np.column_stack([rankdata(Z[ok, j]) for j in range(Z.shape[1])] + [np.ones(int(ok.sum()))])
    a, b = rankdata(x[ok]), rankdata(y[ok])
    ra = a - R @ np.linalg.lstsq(R, a, rcond=None)[0]
    rb = b - R @ np.linalg.lstsq(R, b, rcond=None)[0]
    d = np.linalg.norm(ra) * np.linalg.norm(rb)
    return float(ra @ rb / d) if d > 0 else float("nan")


def boot_ci(fn, n: int, draws: int, seed: int):
    """Paired unit (TF) bootstrap.  Returns (lo, hi, P(stat <= 0), point-free)."""
    rng = np.random.default_rng(seed)
    vals = np.empty(draws)
    for d in range(draws):
        vals[d] = fn(rng.integers(0, n, n))
    v = vals[np.isfinite(vals)]
    if v.size == 0:
        return float("nan"), float("nan"), float("nan")
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)),
            float((v <= 0).mean()))


def holm(pvals: dict[str, float], alpha: float = 0.05) -> dict[str, dict]:
    """Holm step-down at familywise `alpha`.  Returns per-cell adjusted p and reject."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    out, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        adj = min(1.0, max(running, (m - i) * p))
        running = adj
        out[k] = dict(p_raw=float(p), p_holm=float(adj), reject=bool(adj <= alpha))
    return out


def perm_p(observed: float, null: np.ndarray, side: str) -> float:
    """One-sided permutation p with the +1 correction.  side in {'less','greater'}."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if null.size == 0 or not np.isfinite(observed):
        return float("nan")
    if side == "less":
        k = int((null <= observed).sum())
    elif side == "greater":
        k = int((null >= observed).sum())
    else:
        raise Refuse(f"bad side {side!r}")
    return float((k + 1) / (null.size + 1))


def dump_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, default=_default) + "\n")


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))
