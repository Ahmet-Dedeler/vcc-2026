"""052 — source-only rank-64 geometry versus full-rank union whitening.

The implementation is deliberately isolated and memory-bounded. Large donor
arrays, the temporary basis, and the covariance are memory-mapped. No
single-cell counts, cell-eval2 process, submission, or production artifact is
created.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import shutil
import signal
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CACHE_044 = ROOT / "experiments/044-csls-distill/cache"
SCRATCH = HERE / ".scratch"

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
OUTER_DONORS = tuple(DONOR_PATHS)
PRIMARY_VARIANTS = (
    "identity",
    "union_only",
    "subspace_only",
    "union_then_subspace",
    "subspace_then_union",
)
ELIGIBLE_OVER_UNION = (
    "subspace_only",
    "union_then_subspace",
    "subspace_then_union",
)
NULL_VARIANTS = (
    "random_subspace_only",
    "union_then_random",
    "random_then_union",
)
MATCHED_NULL_VARIANT = {
    "subspace_only": "random_subspace_only",
    "union_then_subspace": "union_then_random",
    "subspace_then_union": "random_then_union",
}


def rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / 1024**3
    return value * 1024 / 1024**3


def enforce_soft_limit(stage: str, limit_gib: float) -> None:
    peak = rss_gib()
    if peak >= limit_gib:
        raise MemoryError(
            f"soft stop at {stage}: peak RSS {peak:.2f} GiB >= {limit_gib:.2f} GiB"
        )


def read_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


def stable_seed(parts: Iterable[str], extra: int = 0) -> int:
    payload = "|".join(parts) + f"|{extra}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:4], "little")


def row_normalize(x: np.ndarray) -> np.ndarray:
    x64 = np.asarray(x, dtype=np.float64)
    return x64 / (np.linalg.norm(x64, axis=1, keepdims=True) + 1e-12)


def pds_per_query(
    prediction: np.ndarray,
    truth: np.ndarray,
) -> np.ndarray:
    if prediction.ndim != 2 or truth.ndim != 2:
        raise ValueError("PDS inputs must be two-dimensional")
    if prediction.shape != truth.shape:
        raise ValueError(
            "PDS requires aligned prediction/truth matrices, got "
            f"{prediction.shape} and {truth.shape}"
        )
    if len(truth) < 2:
        raise ValueError("PDS requires at least two perturbations")
    q = row_normalize(prediction)
    t = row_normalize(truth)
    with np.errstate(all="ignore"):
        scores = q @ t.T
    if not np.isfinite(scores).all():
        raise FloatingPointError("non-finite PDS scores")
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    rank = (scores > own).sum(1) + 0.5 * ((scores == own).sum(1) - 1)
    return 1.0 - rank / (n - 1)


def pds(prediction: np.ndarray, truth: np.ndarray) -> float:
    return float(np.mean(pds_per_query(prediction, truth)))


def split_folds(n: int, folds: int, seed: int) -> list[np.ndarray]:
    order = np.random.default_rng(seed).permutation(n)
    return [np.sort(part) for part in np.array_split(order, folds)]


def close_memmap(array: np.ndarray | None) -> None:
    if isinstance(array, np.memmap) and getattr(array, "_mmap", None) is not None:
        array.flush()
        array._mmap.close()


def clean_scratch() -> None:
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)


def install_signal_cleanup() -> None:
    def handler(signum: int, _frame: object) -> None:
        clean_scratch()
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, handler)
    signal.signal(signal.SIGINT, handler)


@dataclass
class Donor:
    name: str
    path: Path
    genes: list[str]
    perts: list[str]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    supported_official: np.ndarray
    supported_donor: np.ndarray
    background_rows: np.ndarray
    center: np.ndarray


@dataclass
class UnionOperator:
    positions: np.ndarray
    evals: np.ndarray
    evecs: np.ndarray
    gamma: float

    def apply(self, matrix: np.ndarray) -> np.ndarray:
        out = np.asarray(matrix, dtype=np.float32).copy()
        block = out[:, self.positions]
        before = float(np.median(np.linalg.norm(block, axis=1))) or 1.0
        weights = np.power(self.evals, -self.gamma / 2.0).astype(np.float32)
        with np.errstate(all="ignore"):
            transformed = ((block @ self.evecs) * weights) @ self.evecs.T
        if not np.isfinite(transformed).all():
            raise FloatingPointError("non-finite union-whitened values")
        after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
        out[:, self.positions] = transformed * (before / after)
        return out


@dataclass
class SubspaceOperator:
    positions: np.ndarray
    basis: np.ndarray
    weights: np.ndarray

    def apply(self, matrix: np.ndarray) -> np.ndarray:
        out = np.asarray(matrix, dtype=np.float32).copy()
        block = out[:, self.positions]
        before = float(np.median(np.linalg.norm(block, axis=1))) or 1.0
        with np.errstate(all="ignore"):
            scores = block @ self.basis
            transformed = block + (scores * (self.weights - 1.0)) @ self.basis.T
        if not np.isfinite(transformed).all():
            raise FloatingPointError("non-finite subspace values")
        after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
        out[:, self.positions] = transformed * (before / after)
        return out


@dataclass(frozen=True)
class FoldAxes:
    score_genes: list[str]
    union_genes: list[str]
    subspace_genes: list[str]
    score_official: np.ndarray
    union_official: np.ndarray
    union_positions: np.ndarray
    subspace_positions: np.ndarray
    subspace_columns: np.ndarray
    inner_score_positions: np.ndarray


def validate_design() -> tuple[list[str], list[str], list[str], list[str]]:
    official = list(
        pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
    )
    vcc = list(
        pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
    )
    train = read_lines(CACHE_044 / "train_perts.txt")
    outer = read_lines(CACHE_044 / "outer_perts.txt")
    master = read_lines(CACHE_044 / "master_genes.txt")
    train_mask = np.load(CACHE_044 / "train_mask.npy")
    subspace_genes = [
        gene for gene, keep in zip(master, train_mask, strict=True) if keep
    ]
    if len(official) != 18_533:
        raise RuntimeError(f"expected 18,533 official genes, got {len(official)}")
    if len(vcc) != 300 or len(train) != 600 or len(outer) != 259:
        raise RuntimeError("unexpected panel sizes")
    if set(vcc) & set(train):
        raise RuntimeError("VCC and inner-selection perturbations overlap")
    if not set(outer) <= set(vcc):
        raise RuntimeError("outer panel is not a VCC subset")
    if len(subspace_genes) != 5_582:
        raise RuntimeError(
            f"expected the experiment-048 5,582-gene subspace, got "
            f"{len(subspace_genes)}"
        )
    return official, vcc, train, subspace_genes


def build_fold_axes(
    official: list[str],
    vcc: list[str],
    train: list[str],
    subspace_genes: list[str],
    held_supported_official: np.ndarray,
    source_supported_official: Iterable[np.ndarray],
) -> FoldAxes:
    """Construct score and transform axes with target coordinates isolated."""
    official_array = np.asarray(official)
    official_pos = {gene: index for index, gene in enumerate(official)}
    if len(official_pos) != len(official):
        raise RuntimeError("official gene axis contains duplicate symbols")

    vcc_set = set(vcc)
    train_set = set(train)
    excluded_coordinates = vcc_set | train_set
    held_support = set(official_array[held_supported_official])
    source_union = set().union(
        *(set(official_array[supported]) for supported in source_supported_official)
    )

    # The official outer score removes only the 300 VCC targets. The 600
    # inner-panel target coordinates remain in this matrix but neither
    # operator may alter them.
    score_genes = [
        gene for gene in official if gene in held_support and gene not in vcc_set
    ]
    score_pos = {gene: index for index, gene in enumerate(score_genes)}
    union_genes = [
        gene
        for gene in score_genes
        if gene in source_union and gene not in excluded_coordinates
    ]
    safe_subspace = [
        gene
        for gene in subspace_genes
        if (
            gene in score_pos
            and gene in source_union
            and gene not in excluded_coordinates
        )
    ]
    union_local_pos = {gene: index for index, gene in enumerate(union_genes)}
    missing_from_union = set(safe_subspace) - set(union_local_pos)
    if missing_from_union:
        raise RuntimeError(
            "subspace coordinates are not contained in the union axis: "
            f"{sorted(missing_from_union)[:5]}"
        )
    if set(union_genes) & excluded_coordinates:
        raise RuntimeError("target-gene coordinate entered the union operator")
    if set(safe_subspace) & excluded_coordinates:
        raise RuntimeError("target-gene coordinate entered the subspace operator")

    return FoldAxes(
        score_genes=score_genes,
        union_genes=union_genes,
        subspace_genes=safe_subspace,
        score_official=np.fromiter(
            (official_pos[gene] for gene in score_genes),
            dtype=np.int64,
            count=len(score_genes),
        ),
        union_official=np.fromiter(
            (official_pos[gene] for gene in union_genes),
            dtype=np.int64,
            count=len(union_genes),
        ),
        union_positions=np.fromiter(
            (score_pos[gene] for gene in union_genes),
            dtype=np.int64,
            count=len(union_genes),
        ),
        subspace_positions=np.fromiter(
            (score_pos[gene] for gene in safe_subspace),
            dtype=np.int64,
            count=len(safe_subspace),
        ),
        subspace_columns=np.fromiter(
            (union_local_pos[gene] for gene in safe_subspace),
            dtype=np.int64,
            count=len(safe_subspace),
        ),
        inner_score_positions=np.fromiter(
            (
                index
                for index, gene in enumerate(score_genes)
                if gene not in excluded_coordinates
            ),
            dtype=np.int64,
        ),
    )


def stream_center(
    delta: np.ndarray,
    rows: np.ndarray,
    donor_cols: np.ndarray,
    chunk_rows: int,
) -> np.ndarray:
    total = np.zeros(len(donor_cols), dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        rr = rows[lo : lo + chunk_rows]
        block = np.asarray(delta[np.ix_(rr, donor_cols)], dtype=np.float32)
        total += block.sum(axis=0, dtype=np.float64)
    return (total / len(rows)).astype(np.float32)


def prepare_donors(
    official: list[str],
    excluded_perts: set[str],
    chunk_rows: int,
    soft_limit: float,
) -> dict[str, Donor]:
    official_pos = {gene: index for index, gene in enumerate(official)}
    donors = {}
    for name, path in DONOR_PATHS.items():
        genes = read_lines(path / "genes.txt")
        perts = read_lines(path / "perts.txt")
        gene_pos = {gene: index for index, gene in enumerate(genes)}
        pert_pos = {pert: index for index, pert in enumerate(perts)}
        supported = [gene for gene in genes if gene in official_pos]
        supported_official = np.fromiter(
            (official_pos[gene] for gene in supported),
            dtype=np.int64,
            count=len(supported),
        )
        supported_donor = np.fromiter(
            (gene_pos[gene] for gene in supported),
            dtype=np.int64,
            count=len(supported),
        )
        background_rows = np.fromiter(
            (index for index, pert in enumerate(perts) if pert not in excluded_perts),
            dtype=np.int64,
        )
        delta = np.load(path / "delta.npy", mmap_mode="r")
        center = stream_center(delta, background_rows, supported_donor, chunk_rows)
        close_memmap(delta)
        del delta
        gc.collect()
        donors[name] = Donor(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos=pert_pos,
            supported_official=supported_official,
            supported_donor=supported_donor,
            background_rows=background_rows,
            center=center,
        )
        print(
            f"  prepared {name:16} official={len(supported):5d} "
            f"background={len(background_rows):5d} peak={rss_gib():.2f} GiB",
            flush=True,
        )
        enforce_soft_limit(f"prepare {name}", soft_limit)
    return donors


def donor_axis_map(
    donor: Donor,
    axis_official: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    support_map = {
        int(official_index): (int(donor_col), donor.center[index])
        for index, (official_index, donor_col) in enumerate(
            zip(
                donor.supported_official,
                donor.supported_donor,
                strict=True,
            )
        )
    }
    dst = []
    src = []
    center = []
    for axis_index, official_index in enumerate(axis_official):
        found = support_map.get(int(official_index))
        if found is None:
            continue
        donor_col, value = found
        dst.append(axis_index)
        src.append(donor_col)
        center.append(value)
    return (
        np.asarray(dst, dtype=np.int64),
        np.asarray(src, dtype=np.int64),
        np.asarray(center, dtype=np.float32),
    )


def choose_background_rows(
    donor: Donor,
    maximum: int,
    seed: int,
) -> np.ndarray:
    rows = donor.background_rows
    if len(rows) <= maximum:
        return rows.copy()
    chosen = np.random.default_rng(seed).choice(rows, maximum, replace=False)
    return np.sort(chosen.astype(np.int64))


def fill_basis(
    path: Path,
    sources: tuple[str, ...],
    donors: dict[str, Donor],
    union_official: np.ndarray,
    basis_per_donor: int,
    operator_seed: int,
    chunk_rows: int,
    soft_limit: float,
) -> tuple[np.memmap, dict[str, float], float, dict[str, int]]:
    selections = {
        name: choose_background_rows(
            donors[name],
            basis_per_donor,
            stable_seed(("union-basis", name), operator_seed),
        )
        for name in sources
    }
    total_rows = sum(len(rows) for rows in selections.values())
    basis = np.lib.format.open_memmap(
        path,
        mode="w+",
        dtype=np.float32,
        shape=(total_rows, len(union_official)),
    )
    basis[:] = 0.0
    donor_ranges = {}
    donor_scales = {}
    cursor = 0
    for name in sources:
        donor = donors[name]
        rows = selections[name]
        dst, src, center = donor_axis_map(donor, union_official)
        norms = np.empty(len(rows), dtype=np.float64)
        delta = np.load(donor.path / "delta.npy", mmap_mode="r")
        start = cursor
        for lo in range(0, len(rows), chunk_rows):
            hi = min(len(rows), lo + chunk_rows)
            rr = rows[lo:hi]
            block = np.asarray(delta[np.ix_(rr, src)], dtype=np.float32)
            block -= center
            view = basis[cursor + lo : cursor + hi]
            view[:, dst] = block
            norms[lo:hi] = np.linalg.norm(view, axis=1)
        close_memmap(delta)
        del delta
        donor_scales[name] = float(np.median(norms)) or 1.0
        cursor += len(rows)
        donor_ranges[name] = (start, cursor)
        enforce_soft_limit(f"fill basis {name}", soft_limit)

    reference = float(np.median(list(donor_scales.values())))
    for name, (lo, hi) in donor_ranges.items():
        factor = reference / donor_scales[name]
        for start in range(lo, hi, chunk_rows):
            stop = min(hi, start + chunk_rows)
            basis[start:stop] *= factor
    basis.flush()
    return (
        basis,
        donor_scales,
        reference,
        {name: len(rows) for name, rows in selections.items()},
    )


def fit_rank64_subspace(
    basis: np.ndarray,
    subspace_columns: np.ndarray,
    seed: int,
    rank: int = 64,
    gamma: float = 2.0,
    oversample: int = 16,
    power_iterations: int = 2,
    chunk_rows: int = 256,
) -> tuple[np.ndarray, np.ndarray]:
    """Randomized top-PC fit from the background basis only."""
    width = rank + oversample
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((len(subspace_columns), width), dtype=np.float32)
    y = np.empty((len(basis), width), dtype=np.float32)
    for lo in range(0, len(basis), chunk_rows):
        hi = min(len(basis), lo + chunk_rows)
        with np.errstate(all="ignore"):
            y[lo:hi] = basis[lo:hi, subspace_columns] @ omega
    for _ in range(power_iterations):
        q, _ = np.linalg.qr(y, mode="reduced")
        z = np.zeros((len(subspace_columns), width), dtype=np.float32)
        for lo in range(0, len(basis), chunk_rows):
            hi = min(len(basis), lo + chunk_rows)
            with np.errstate(all="ignore"):
                z += basis[lo:hi, subspace_columns].T @ q[lo:hi]
        for lo in range(0, len(basis), chunk_rows):
            hi = min(len(basis), lo + chunk_rows)
            with np.errstate(all="ignore"):
                y[lo:hi] = basis[lo:hi, subspace_columns] @ z
    q, _ = np.linalg.qr(y, mode="reduced")
    small = np.zeros((width, len(subspace_columns)), dtype=np.float32)
    for lo in range(0, len(basis), chunk_rows):
        hi = min(len(basis), lo + chunk_rows)
        with np.errstate(all="ignore"):
            small += q[lo:hi].T @ basis[lo:hi, subspace_columns]
    _, singular, vt = np.linalg.svd(small, full_matrices=False)
    vectors = vt[:rank].T.astype(np.float32, copy=False)
    evals = np.square(singular[:rank].astype(np.float64)) / len(basis)
    tail = max(float(evals[-1]), float(evals[0]) * 1e-12, 1e-12)
    weights = np.power(evals / tail, -gamma / 2.0).astype(np.float32)

    if not all(np.isfinite(value).all() for value in (vectors, evals, weights)):
        raise FloatingPointError("non-finite subspace fit")
    return vectors, weights


def build_covariance(
    basis: np.memmap,
    covariance_path: Path,
    soft_limit: float,
) -> np.memmap:
    """Build the union covariance while the source basis is memory-mapped."""
    covariance = np.lib.format.open_memmap(
        covariance_path,
        mode="w+",
        dtype=np.float32,
        shape=(basis.shape[1], basis.shape[1]),
        fortran_order=True,
    )
    print(
        f"    covariance {basis.shape[1]:,}x{basis.shape[1]:,} "
        f"from {basis.shape[0]:,} rows",
        flush=True,
    )
    with np.errstate(all="ignore"):
        np.matmul(basis.T, basis, out=covariance)
    covariance /= float(len(basis))
    covariance.flush()
    if not np.isfinite(covariance[::257, ::257]).all():
        raise FloatingPointError("non-finite covariance sample")
    enforce_soft_limit("union covariance", soft_limit)
    return covariance


def diagonalize_union(
    covariance: np.memmap,
    covariance_path: Path,
    soft_limit: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Diagonalize after the much larger source basis has been unmapped."""
    print("    full float32 symmetric eigendecomposition (evr)", flush=True)
    evals, evecs = eigh(
        covariance,
        overwrite_a=True,
        check_finite=False,
        driver="evr",
    )
    close_memmap(covariance)
    del covariance
    covariance_path.unlink(missing_ok=True)
    gc.collect()
    evals = np.maximum(
        evals.astype(np.float64),
        float(evals[-1]) * 1e-12,
    )
    evecs = np.asarray(evecs, dtype=np.float32)
    if not np.isfinite(evecs[::257, ::257]).all():
        raise FloatingPointError("non-finite union eigenvectors")
    enforce_soft_limit("union eigendecomposition", soft_limit)
    return evals, evecs


def extract_panel(
    donor: Donor,
    perturbations: list[str],
    score_official: np.ndarray,
    chunk_rows: int,
    *,
    center_on_background: bool = True,
) -> np.ndarray:
    dst, src, center = donor_axis_map(donor, score_official)
    rows = np.fromiter(
        (donor.pert_pos[pert] for pert in perturbations),
        dtype=np.int64,
        count=len(perturbations),
    )
    out = np.zeros((len(rows), len(score_official)), dtype=np.float32)
    delta = np.load(donor.path / "delta.npy", mmap_mode="r")
    for lo in range(0, len(rows), chunk_rows):
        hi = min(len(rows), lo + chunk_rows)
        block = np.asarray(delta[np.ix_(rows[lo:hi], src)], dtype=np.float32)
        if center_on_background:
            block -= center
        out[lo:hi, dst] = block
    close_memmap(delta)
    del delta
    return out


def transform_primary_variants(
    query: np.ndarray,
    union: UnionOperator,
    subspace: SubspaceOperator,
) -> dict[str, np.ndarray]:
    union_only = union.apply(query)
    subspace_only = subspace.apply(query)
    return {
        "identity": query,
        "union_only": union_only,
        "subspace_only": subspace_only,
        "union_then_subspace": subspace.apply(union_only),
        "subspace_then_union": union.apply(subspace_only),
    }


def transform_null_variants(
    query: np.ndarray,
    union: UnionOperator,
    random_subspace: SubspaceOperator,
    *,
    union_only: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    if union_only is None:
        union_only = union.apply(query)
    random_only = random_subspace.apply(query)
    return {
        "random_subspace_only": random_only,
        "union_then_random": random_subspace.apply(union_only),
        "random_then_union": union.apply(random_only),
    }


def score_inner_seed(
    transformed_by_direction: dict[str, dict[str, np.ndarray]],
    truth_by_direction: dict[str, np.ndarray],
    inner_score_positions: np.ndarray,
    seed: int,
    folds: int,
) -> dict:
    rows = {name: [] for name in (*PRIMARY_VARIANTS, *NULL_VARIANTS)}
    direction_means = {name: {} for name in (*PRIMARY_VARIANTS, *NULL_VARIANTS)}
    for direction, transformed in transformed_by_direction.items():
        truth = truth_by_direction[direction]
        score_truth = truth[:, inner_score_positions]
        fold_rows = split_folds(
            len(truth),
            folds,
            stable_seed(("inner-folds", direction), seed),
        )
        for name, prediction in transformed.items():
            # Fold membership controls only which query rows are aggregated.
            # Every query is ranked against the complete 600-perturbation
            # candidate catalog, matching the competition retrieval task.
            query_scores = pds_per_query(
                prediction[:, inner_score_positions],
                score_truth,
            )
            direction_means[name][direction] = float(np.mean(query_scores))
            for fold_index, selected_rows in enumerate(fold_rows):
                value = float(np.mean(query_scores[selected_rows]))
                rows[name].append(
                    {
                        "direction": direction,
                        "fold": fold_index,
                        "score": value,
                    }
                )

    summary = {}
    union_scores = np.asarray([row["score"] for row in rows["union_only"]])
    for name, values in rows.items():
        scores = np.asarray([row["score"] for row in values])
        deltas = scores - union_scores
        summary[name] = {
            "mean": float(np.mean(scores)),
            "delta_vs_union": float(np.mean(deltas)),
            "positive_folds_vs_union": int(np.sum(deltas > 0)),
            "direction_delta_vs_union": {
                direction: (
                    direction_means[name][direction]
                    - direction_means["union_only"][direction]
                )
                for direction in direction_means[name]
            },
            "fold_scores": [float(value) for value in scores],
            "fold_deltas_vs_union": [float(value) for value in deltas],
        }
    return summary


def select_over_union(
    summary: dict[str, dict],
    min_gain: float,
    min_positive_folds: int,
) -> str:
    eligible = []
    for name in ELIGIBLE_OVER_UNION:
        row = summary[name]
        if row["delta_vs_union"] < min_gain:
            continue
        if row["positive_folds_vs_union"] < min_positive_folds:
            continue
        if not all(delta > 0 for delta in row["direction_delta_vs_union"].values()):
            continue
        eligible.append(name)
    if not eligible:
        return "union_only"
    complexity = {
        "subspace_only": 0,
        "union_then_subspace": 1,
        "subspace_then_union": 2,
    }
    best = max(summary[name]["delta_vs_union"] for name in eligible)
    near = [name for name in eligible if summary[name]["delta_vs_union"] >= best - 5e-4]
    return min(near, key=lambda name: complexity[name])


def run_outer(
    outer: str,
    donors: dict[str, Donor],
    official: list[str],
    vcc: list[str],
    train: list[str],
    subspace_genes: list[str],
    args: argparse.Namespace,
) -> dict:
    sources = tuple(name for name in OUTER_DONORS if name != outer)
    excluded_coordinates = set(vcc) | set(train)
    axes = build_fold_axes(
        official,
        vcc,
        train,
        subspace_genes,
        donors[outer].supported_official,
        (donors[name].supported_official for name in sources),
    )

    fold_scratch = SCRATCH / outer
    fold_scratch.mkdir(parents=True, exist_ok=True)
    basis_path = fold_scratch / "basis.npy"
    covariance_path = fold_scratch / "covariance.npy"
    basis = None
    covariance = None
    try:
        print(
            f"\n=== outer {outer}: score={len(axes.score_genes):,} "
            f"union={len(axes.union_genes):,} "
            f"subspace={len(axes.subspace_genes):,} ===",
            flush=True,
        )
        basis, donor_scales, reference, basis_rows = fill_basis(
            basis_path,
            sources,
            donors,
            axes.union_official,
            args.basis_per_donor,
            args.operator_seed,
            args.chunk_rows,
            args.soft_rss_gib,
        )
        print(
            f"    basis={basis.shape} rows={basis_rows} peak={rss_gib():.2f} GiB",
            flush=True,
        )
        sub_basis, sub_weights = fit_rank64_subspace(
            basis,
            axes.subspace_columns,
            stable_seed(("subspace", outer), args.operator_seed),
            chunk_rows=args.chunk_rows,
        )
        covariance = build_covariance(
            basis,
            covariance_path,
            args.soft_rss_gib,
        )
        close_memmap(basis)
        del basis
        basis = None
        basis_path.unlink(missing_ok=True)
        gc.collect()
        evals, evecs = diagonalize_union(
            covariance,
            covariance_path,
            args.soft_rss_gib,
        )
        covariance = None

        union_operator = UnionOperator(
            positions=axes.union_positions,
            evals=evals,
            evecs=evecs,
            gamma=args.union_gamma,
        )
        subspace_operator = SubspaceOperator(
            positions=axes.subspace_positions,
            basis=sub_basis,
            weights=sub_weights,
        )

        # Freeze every selection decision before loading the outer panel.
        transformed_deterministic = {}
        truth_by_direction = {}
        train_panels = {
            name: extract_panel(
                donors[name],
                train,
                axes.score_official,
                args.chunk_rows,
            )
            for name in sources
        }
        for valid in sources:
            predictor = next(name for name in sources if name != valid)
            query = (reference / donor_scales[predictor]) * train_panels[predictor]
            truth = train_panels[valid]
            variants = transform_primary_variants(
                query,
                union_operator,
                subspace_operator,
            )
            direction = f"{predictor}->{valid}"
            transformed_deterministic[direction] = variants
            truth_by_direction[direction] = truth

        seed_results = {}
        random_bases = {}
        for seed in args.seeds:
            rng = np.random.default_rng(stable_seed(("random-null", outer), seed))
            random_matrix = rng.standard_normal(sub_basis.shape, dtype=np.float32)
            random_basis, _ = np.linalg.qr(random_matrix, mode="reduced")
            random_operator = SubspaceOperator(
                axes.subspace_positions,
                random_basis.astype(np.float32, copy=False),
                sub_weights,
            )
            random_bases[str(seed)] = random_operator.basis

            transformed = {}
            for direction, deterministic in transformed_deterministic.items():
                predictor = direction.split("->", 1)[0]
                query = (reference / donor_scales[predictor]) * train_panels[predictor]
                null_variants = transform_null_variants(
                    query,
                    union_operator,
                    random_operator,
                    union_only=deterministic["union_only"],
                )
                transformed[direction] = {
                    **deterministic,
                    **null_variants,
                }

            inner = score_inner_seed(
                transformed,
                truth_by_direction,
                axes.inner_score_positions,
                seed,
                args.folds,
            )
            selected = select_over_union(
                inner,
                args.min_inner_gain,
                args.min_positive_folds,
            )
            seed_results[str(seed)] = {
                "inner": inner,
                "selected": selected,
            }
            del transformed, null_variants, random_operator
            gc.collect()

        # The held panel is touched only after all inner selections are frozen.
        outer_perts = read_lines(CACHE_044 / "outer_perts.txt")
        outer_source_panels = {
            name: extract_panel(
                donors[name],
                outer_perts,
                axes.score_official,
                args.chunk_rows,
            )
            for name in sources
        }
        outer_query = sum(
            (reference / donor_scales[name]) * outer_source_panels[name]
            for name in sources
        ) / len(sources)
        outer_truth = extract_panel(
            donors[outer],
            outer_perts,
            axes.score_official,
            args.chunk_rows,
            center_on_background=False,
        )
        outer_base_variants = transform_primary_variants(
            outer_query,
            union_operator,
            subspace_operator,
        )
        outer_primary_scores = {
            name: pds(matrix, outer_truth)
            for name, matrix in outer_base_variants.items()
        }

        for seed in args.seeds:
            seed_key = str(seed)
            random_operator = SubspaceOperator(
                axes.subspace_positions,
                random_bases[seed_key],
                sub_weights,
            )
            outer_null_variants = transform_null_variants(
                outer_query,
                union_operator,
                random_operator,
                union_only=outer_base_variants["union_only"],
            )
            outer_null_scores = {
                name: pds(outer_null_variants[name], outer_truth)
                for name in NULL_VARIANTS
            }
            selected = seed_results[seed_key]["selected"]
            matched_null = MATCHED_NULL_VARIANT.get(selected)
            seed_results[seed_key].update(
                {
                    "outer_selected_score": outer_primary_scores[selected],
                    "outer_selected_delta_vs_union": (
                        outer_primary_scores[selected]
                        - outer_primary_scores["union_only"]
                    ),
                    "outer_null_scores": outer_null_scores,
                    "matched_null_variant": matched_null,
                    "outer_matched_null_delta_vs_union": (
                        0.0
                        if matched_null is None
                        else (
                            outer_null_scores[matched_null]
                            - outer_primary_scores["union_only"]
                        )
                    ),
                }
            )
            del outer_null_variants, random_operator
            gc.collect()

        result = {
            "outer_donor": outer,
            "source_donors": list(sources),
            "axes": {
                "outer_score_genes": len(axes.score_genes),
                "inner_score_genes": len(axes.inner_score_positions),
                "union_transform_genes": len(axes.union_genes),
                "subspace_transform_genes": len(axes.subspace_genes),
                "excluded_target_coordinates": len(excluded_coordinates),
            },
            "basis": {
                "rows_by_donor": basis_rows,
                "rows_total": int(sum(basis_rows.values())),
                "donor_scales": donor_scales,
                "reference_scale": reference,
            },
            "outer_scores": outer_primary_scores,
            "seeds": seed_results,
        }
        print(
            "    outer "
            + " ".join(
                f"{name}={outer_primary_scores[name]:.5f}" for name in PRIMARY_VARIANTS
            ),
            flush=True,
        )
        for seed in args.seeds:
            row = seed_results[str(seed)]
            print(
                f"    seed={seed} selected={row['selected']} "
                f"delta_vs_union={row['outer_selected_delta_vs_union']:+.5f}",
                flush=True,
            )
        del (
            evecs,
            evals,
            union_operator,
            subspace_operator,
            train_panels,
            outer_source_panels,
            outer_query,
            outer_truth,
            transformed_deterministic,
            truth_by_direction,
            outer_base_variants,
            random_bases,
        )
        gc.collect()
        enforce_soft_limit(f"finish outer {outer}", args.soft_rss_gib)
        return result
    finally:
        if basis is not None:
            close_memmap(basis)
        if covariance is not None:
            close_memmap(covariance)
        clean_scratch()


def summarize_folds(
    folds: dict[str, dict],
    outer_names: list[str],
    seeds: list[int],
) -> tuple[dict[str, float], dict[str, dict], dict]:
    """Summarize frozen outer evaluations without selecting on outer scores."""
    fixed_macro = {
        variant: float(
            np.mean([folds[outer]["outer_scores"][variant] for outer in outer_names])
        )
        for variant in PRIMARY_VARIANTS
    }
    fixed_deltas = {
        variant: score - fixed_macro["union_only"]
        for variant, score in fixed_macro.items()
    }
    seed_macro = {}
    for seed in seeds:
        seed_key = str(seed)
        selected_scores = [
            folds[outer]["seeds"][seed_key]["outer_selected_score"]
            for outer in outer_names
        ]
        selected_deltas = [
            folds[outer]["seeds"][seed_key]["outer_selected_delta_vs_union"]
            for outer in outer_names
        ]
        matched_null_deltas = [
            folds[outer]["seeds"][seed_key]["outer_matched_null_delta_vs_union"]
            for outer in outer_names
        ]
        matched_advantages = np.asarray(selected_deltas) - np.asarray(
            matched_null_deltas
        )
        selected_variants = {
            outer: folds[outer]["seeds"][seed_key]["selected"] for outer in outer_names
        }
        seed_macro[seed_key] = {
            "selected_score": float(np.mean(selected_scores)),
            "selected_delta_vs_union": float(np.mean(selected_deltas)),
            "positive_outer_folds_vs_union": int(
                np.sum(np.asarray(selected_deltas) > 0)
            ),
            "selected_variants_by_outer": selected_variants,
            "all_selected_non_union": all(
                name != "union_only" for name in selected_variants.values()
            ),
            "matched_random_delta_vs_union": float(np.mean(matched_null_deltas)),
            "positive_outer_folds_vs_matched_random": int(
                np.sum(matched_advantages > 0)
            ),
        }

    nested_gate = all(
        seed_macro[str(seed)]["all_selected_non_union"]
        and seed_macro[str(seed)]["selected_delta_vs_union"] > 0
        and seed_macro[str(seed)]["positive_outer_folds_vs_union"] == len(outer_names)
        and seed_macro[str(seed)]["positive_outer_folds_vs_matched_random"]
        == len(outer_names)
        for seed in seeds
    )
    head_to_head = {
        "fixed_outer_scores_are_diagnostic_only": True,
        "fixed_pipeline_deltas_vs_union": fixed_deltas,
        "nested_selection_gate": bool(nested_gate),
    }
    return fixed_macro, seed_macro, head_to_head


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, nargs="+", default=[52, 152])
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--basis-per-donor", type=int, default=12_000)
    parser.add_argument("--operator-seed", type=int, default=0)
    parser.add_argument("--union-gamma", type=float, default=0.65)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--soft-rss-gib", type=float, default=3.5)
    parser.add_argument("--min-inner-gain", type=float, default=0.001)
    parser.add_argument("--min-positive-folds", type=int, default=6)
    parser.add_argument(
        "--outer",
        choices=OUTER_DONORS,
        nargs="+",
        default=list(OUTER_DONORS),
    )
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    if args.soft_rss_gib <= 0:
        parser.error("--soft-rss-gib must be positive")
    if len(args.seeds) != 2:
        parser.error("exactly two seeds are required")

    install_signal_cleanup()
    clean_scratch()
    official, vcc, train, subspace_genes = validate_design()
    excluded_perts = set(vcc) | set(train)
    print(
        f"design: outer={args.outer} seeds={args.seeds} "
        f"excluded_perts={len(excluded_perts)}",
        flush=True,
    )
    donors = prepare_donors(
        official,
        excluded_perts,
        args.chunk_rows,
        args.soft_rss_gib,
    )

    folds = {}
    try:
        for outer in args.outer:
            folds[outer] = run_outer(
                outer,
                donors,
                official,
                vcc,
                train,
                subspace_genes,
                args,
            )
    finally:
        clean_scratch()

    fixed_macro, seed_macro, head_to_head = summarize_folds(
        folds,
        args.outer,
        args.seeds,
    )
    result = {
        "design": {
            "outer_donors": list(args.outer),
            "seeds": args.seeds,
            "outer_vcc_perturbations": 259,
            "inner_non_vcc_perturbations": 600,
            "excluded_target_coordinates": 900,
            "union_gamma": args.union_gamma,
            "subspace_rank": 64,
            "subspace_gamma": 2.0,
            "basis_per_donor_cap": args.basis_per_donor,
            "serial": True,
            "memory_mapped_large_arrays": True,
            "full_gene_eigendecompositions": len(args.outer),
            "held_donor_used_in_fit_or_selection": False,
            "cell_eval2_runs": 0,
            "submission_files_touched": 0,
        },
        "folds": folds,
        "fixed_macro": fixed_macro,
        "seed_macro": seed_macro,
        "head_to_head": head_to_head,
        "peak_rss_gib": rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("\nfixed macro PDS", flush=True)
    for name in PRIMARY_VARIANTS:
        print(f"  {name:24} {fixed_macro[name]:.6f}", flush=True)
    print(
        f"nested selection gate={result['head_to_head']['nested_selection_gate']}",
        flush=True,
    )
    print(
        f"wrote {args.output}; peak RSS={rss_gib():.2f} GiB",
        flush=True,
    )


if __name__ == "__main__":
    main()
