"""048 — strict source-only low-rank geometry under nested LODO.

The target-effect PCA basis from experiment 046 cannot be used: permutation
controls proved that its gain came from target covariance geometry rather than
paired context adaptation. This experiment tests low-rank transforms whose
basis, spectrum, and hyperparameters are learned using source donors and
non-VCC perturbations only.

Inputs are memory-mapped. No single-cell counts, cell-eval2 invocation,
submission file, or production model artifact is created.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Prevent a BLAS thread storm before NumPy is imported.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CACHE_044 = ROOT / "experiments/044-csls-distill/cache"

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
OUTER_DONORS = tuple(DONOR_PATHS)
CURRENT_COMPARATOR_KEY = "source:spectral:r=256:gamma=0.65:ridge=0"


def rss_gib() -> float:
    """Peak process RSS in GiB on macOS and Linux."""
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value / 1024**3
    return value * 1024 / 1024**3


def enforce_soft_memory_limit(stage: str, limit_gib: float) -> None:
    peak = rss_gib()
    if peak >= limit_gib:
        raise MemoryError(
            f"soft memory stop at {stage}: peak RSS {peak:.2f} GiB "
            f">= {limit_gib:.2f} GiB"
        )


def read_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


def stable_seed(parts: Iterable[str], extra: int = 0) -> int:
    payload = "|".join(parts) + f"|{extra}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:4], "little")


def row_normalize(x: np.ndarray) -> np.ndarray:
    x64 = np.asarray(x, dtype=np.float64)
    return x64 / (np.linalg.norm(x64, axis=1, keepdims=True) + 1e-12)


def pds_from_scores(scores: np.ndarray) -> float:
    """Exact VCC rank statistic for a square cosine-similarity matrix."""
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    rank = (scores > own).sum(1) + 0.5 * ((scores == own).sum(1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def split_folds(n: int, folds: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    return [np.sort(part) for part in np.array_split(order, folds)]


@dataclass(frozen=True)
class Candidate:
    basis_kind: str
    family: str
    rank: int
    gamma: float = 0.0
    ridge: float = 0.0
    eta: float = 0.0

    @property
    def key(self) -> str:
        if self.family == "spectral":
            return (
                f"{self.basis_kind}:spectral:r={self.rank}:"
                f"gamma={self.gamma:g}:ridge={self.ridge:g}"
            )
        if self.family == "uniform":
            return (
                f"{self.basis_kind}:uniform:r={self.rank}:eta={self.eta:g}"
            )
        raise ValueError(self.family)

    @property
    def complexity(self) -> tuple:
        family_order = {"uniform": 0, "spectral": 1}
        return (
            self.rank,
            family_order[self.family],
            self.gamma,
            self.eta,
            self.ridge,
        )


@dataclass
class DonorBank:
    name: str
    train_raw: np.ndarray
    train_src: np.ndarray
    outer_raw: np.ndarray
    outer_src: np.ndarray
    normalization_count: int
    normalization_floor: float


@dataclass
class Subspace:
    source_basis: np.ndarray
    source_evals: np.ndarray
    random_basis: np.ndarray


@dataclass
class Pair:
    full_query: np.ndarray
    full_truth: np.ndarray
    local_query: np.ndarray
    local_truth: np.ndarray


def validate_design_files() -> tuple[list[str], list[str], list[str], np.ndarray]:
    required = [
        CACHE_044 / "master_genes.txt",
        CACHE_044 / "train_perts.txt",
        CACHE_044 / "outer_perts.txt",
        CACHE_044 / "train_mask.npy",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"experiment 044 design cache incomplete: {missing}")

    genes = read_lines(CACHE_044 / "master_genes.txt")
    train_perts = read_lines(CACHE_044 / "train_perts.txt")
    outer_perts = read_lines(CACHE_044 / "outer_perts.txt")
    train_mask = np.load(CACHE_044 / "train_mask.npy")
    vcc_panel = set(
        pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
    )

    if len(train_perts) != 600:
        raise RuntimeError(f"expected 600 non-VCC train perturbations, got {len(train_perts)}")
    if len(outer_perts) != 259:
        raise RuntimeError(f"expected 259 outer VCC perturbations, got {len(outer_perts)}")
    if set(train_perts) & vcc_panel:
        raise RuntimeError("training perturbations overlap the VCC panel")
    if not set(outer_perts) <= vcc_panel:
        raise RuntimeError("outer perturbations are not a subset of the VCC panel")
    if set(train_perts) & set(outer_perts):
        raise RuntimeError("training and outer perturbations overlap")
    if train_mask.shape != (len(genes),):
        raise RuntimeError("train feature mask shape mismatch")

    gpos = {gene: i for i, gene in enumerate(genes)}
    if any(gene in gpos for gene in outer_perts):
        raise RuntimeError("VCC target coordinate leaked into the score axis")
    if any(train_mask[gpos[p]] for p in train_perts if p in gpos):
        raise RuntimeError("training target coordinate leaked into shaping axis")
    return genes, train_perts, outer_perts, train_mask.astype(bool, copy=False)


def _stream_background_stats(
    delta: np.ndarray,
    donor_perts: list[str],
    cols: np.ndarray,
    excluded_perts: set[str],
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray, int, float]:
    """Center and scale away from both selection and VCC perturbations."""
    total = np.zeros(len(cols), dtype=np.float64)
    count = 0
    for lo in range(0, delta.shape[0], chunk_rows):
        hi = min(delta.shape[0], lo + chunk_rows)
        keep = np.fromiter(
            (donor_perts[i] not in excluded_perts for i in range(lo, hi)),
            dtype=bool,
            count=hi - lo,
        )
        if not keep.any():
            continue
        block = np.asarray(delta[lo:hi, cols], dtype=np.float32)[keep]
        total += block.sum(axis=0, dtype=np.float64)
        count += int(keep.sum())
    if count == 0:
        raise RuntimeError("donor has no background perturbations")
    mean = total / count

    absolute = np.zeros(len(cols), dtype=np.float64)
    for lo in range(0, delta.shape[0], chunk_rows):
        hi = min(delta.shape[0], lo + chunk_rows)
        keep = np.fromiter(
            (donor_perts[i] not in excluded_perts for i in range(lo, hi)),
            dtype=bool,
            count=hi - lo,
        )
        if not keep.any():
            continue
        block = np.asarray(delta[lo:hi, cols], dtype=np.float32)[keep]
        absolute += np.abs(block - mean).sum(axis=0, dtype=np.float64)
    mean_abs = absolute / count
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25.0)) if live.size else 1e-9
    scale = np.maximum(mean_abs, floor)
    return mean.astype(np.float32), scale.astype(np.float32), count, floor


def _extract_rows(
    delta: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray]:
    raw = np.empty((len(rows), len(cols)), dtype=np.float32)
    src = np.empty_like(raw)
    for lo in range(0, len(rows), chunk_rows):
        hi = min(len(rows), lo + chunk_rows)
        rr = rows[lo:hi]
        block = np.asarray(delta[np.ix_(rr, cols)], dtype=np.float32)
        centered = block - mean
        raw[lo:hi] = centered
        src[lo:hi] = centered / scale
    return raw, src


def prepare_donor(
    name: str,
    path: Path,
    master_genes: list[str],
    train_perts: list[str],
    outer_perts: list[str],
    vcc_panel: set[str],
    chunk_rows: int,
) -> DonorBank:
    donor_genes = read_lines(path / "genes.txt")
    donor_perts = read_lines(path / "perts.txt")
    delta = np.load(path / "delta.npy", mmap_mode="r")
    if delta.shape != (len(donor_perts), len(donor_genes)):
        raise RuntimeError(
            f"{name} delta shape {delta.shape} != "
            f"({len(donor_perts)}, {len(donor_genes)})"
        )

    gpos = {gene: index for index, gene in enumerate(donor_genes)}
    ppos = {pert: index for index, pert in enumerate(donor_perts)}
    missing_genes = [gene for gene in master_genes if gene not in gpos]
    missing_train = [pert for pert in train_perts if pert not in ppos]
    missing_outer = [pert for pert in outer_perts if pert not in ppos]
    if missing_genes or missing_train or missing_outer:
        raise RuntimeError(
            f"{name} missing genes={len(missing_genes)} "
            f"train={len(missing_train)} outer={len(missing_outer)}"
        )

    cols = np.fromiter(
        (gpos[gene] for gene in master_genes),
        dtype=np.int64,
        count=len(master_genes),
    )
    excluded_perts = set(train_perts) | vcc_panel
    mean, scale, normalization_count, floor = _stream_background_stats(
        delta, donor_perts, cols, excluded_perts, chunk_rows
    )
    train_rows = np.fromiter(
        (ppos[pert] for pert in train_perts),
        dtype=np.int64,
        count=len(train_perts),
    )
    outer_rows = np.fromiter(
        (ppos[pert] for pert in outer_perts),
        dtype=np.int64,
        count=len(outer_perts),
    )
    train_raw, train_src = _extract_rows(
        delta, train_rows, cols, mean, scale, chunk_rows
    )
    outer_raw, outer_src = _extract_rows(
        delta, outer_rows, cols, mean, scale, chunk_rows
    )
    del delta
    return DonorBank(
        name=name,
        train_raw=train_raw,
        train_src=train_src,
        outer_raw=outer_raw,
        outer_src=outer_src,
        normalization_count=normalization_count,
        normalization_floor=floor,
    )


def fit_subspace(
    x: np.ndarray,
    max_rank: int,
    seed: int,
) -> Subspace:
    """Fit exact source PCs through the small sample-space Gram matrix."""
    x32 = np.asarray(x, dtype=np.float32)
    with np.errstate(all="ignore"):
        gram = x32 @ x32.T
        gram_evals, gram_evecs = np.linalg.eigh(gram)
    order = np.argsort(gram_evals)[::-1]
    gram_evals = np.maximum(gram_evals[order], 0.0)
    gram_evecs = gram_evecs[:, order]
    live = gram_evals > max(float(gram_evals[0]) * 1e-12, 1e-12)
    rank = min(max_rank, int(live.sum()))
    if rank < max_rank:
        raise RuntimeError(
            f"source matrix rank {rank} is below requested rank {max_rank}"
        )
    gram_evals = gram_evals[:rank]
    gram_evecs = gram_evecs[:, :rank]
    with np.errstate(all="ignore"):
        source_basis = (
            x32.T @ gram_evecs
        ) / np.sqrt(gram_evals)[None, :]
    source_basis = source_basis.astype(np.float32, copy=False)
    source_evals = gram_evals.astype(np.float64) / max(len(x32), 1)

    rng = np.random.default_rng(seed)
    random_matrix = rng.standard_normal(
        (x32.shape[1], rank), dtype=np.float32
    )
    random_basis, _ = np.linalg.qr(random_matrix, mode="reduced")
    if not (
        np.isfinite(source_basis).all()
        and np.isfinite(source_evals).all()
        and np.isfinite(random_basis).all()
    ):
        raise FloatingPointError("non-finite fitted subspace")
    return Subspace(
        source_basis=source_basis,
        source_evals=source_evals,
        random_basis=random_basis.astype(np.float32, copy=False),
    )


def candidate_weights(candidate: Candidate, evals: np.ndarray) -> np.ndarray:
    if candidate.family == "uniform":
        return np.full(candidate.rank, 1.0 - candidate.eta, dtype=np.float64)
    if candidate.family != "spectral":
        raise ValueError(candidate.family)
    local = np.asarray(evals[: candidate.rank], dtype=np.float64)
    tail = max(float(local[-1]), float(local[0]) * 1e-12, 1e-12)
    ratio = (local + candidate.ridge * tail) / (
        tail * (1.0 + candidate.ridge)
    )
    return ratio ** (-candidate.gamma / 2.0)


def score_candidates(
    pair: Pair,
    subspace: Subspace,
    candidates: list[Candidate],
    matrix_keys: set[str] | None = None,
) -> tuple[float, dict[str, float], dict[str, np.ndarray]]:
    matrix_keys = matrix_keys or set()
    with np.errstate(all="ignore"):
        base_scores = pair.full_query @ pair.full_truth.T
        zt_source = pair.local_truth @ subspace.source_basis
        zq_source = pair.local_query @ subspace.source_basis
        zt_random = pair.local_truth @ subspace.random_basis
        zq_random = pair.local_query @ subspace.random_basis
    if not all(
        np.isfinite(value).all()
        for value in (
            base_scores,
            zt_source,
            zq_source,
            zt_random,
            zq_random,
        )
    ):
        raise FloatingPointError("non-finite score projection")
    baseline = pds_from_scores(base_scores)
    scores: dict[str, float] = {}
    matrices: dict[str, np.ndarray] = {}
    for candidate in candidates:
        if candidate.basis_kind == "source":
            zq, zt = zq_source, zt_source
        elif candidate.basis_kind == "random":
            zq, zt = zq_random, zt_random
        else:
            raise ValueError(candidate.basis_kind)
        weights = candidate_weights(candidate, subspace.source_evals)
        rank = candidate.rank
        with np.errstate(all="ignore"):
            transformed_scores = base_scores + (
                zq[:, :rank] * (weights - 1.0)
            ) @ zt[:, :rank].T
        if not np.isfinite(transformed_scores).all():
            raise FloatingPointError(f"non-finite scores for {candidate.key}")
        scores[candidate.key] = pds_from_scores(transformed_scores)
        if candidate.key in matrix_keys:
            matrices[candidate.key] = transformed_scores
    return baseline, scores, matrices


def make_inner_pair(
    predictor: DonorBank,
    truth: DonorBank,
    rows: np.ndarray,
    local_feat: np.ndarray,
) -> Pair:
    q = row_normalize(predictor.train_src[np.ix_(rows, local_feat)])
    t = row_normalize(truth.train_raw[np.ix_(rows, local_feat)])
    return Pair(q, t, q, t)


def source_scales(
    banks: dict[str, DonorBank],
    names: tuple[str, ...],
    local_feat: np.ndarray,
) -> dict[str, float]:
    return {
        name: float(
            np.median(
                np.linalg.norm(
                    banks[name].train_src[:, local_feat],
                    axis=1,
                )
            )
        )
        or 1.0
        for name in names
    }


def scaled_fit_matrix(
    banks: dict[str, DonorBank],
    names: tuple[str, ...],
    rows: np.ndarray,
    local_feat: np.ndarray,
) -> np.ndarray:
    scales = source_scales(banks, names, local_feat)
    reference = float(np.median(list(scales.values())))
    parts = [
        (reference / scales[name])
        * banks[name].train_src[np.ix_(rows, local_feat)]
        for name in names
    ]
    return np.vstack(parts).astype(np.float32, copy=False)


def make_outer_pair(
    banks: dict[str, DonorBank],
    target: str,
    sources: tuple[str, ...],
    local_feat: np.ndarray,
) -> Pair:
    scales = source_scales(banks, sources, local_feat)
    reference = float(np.median(list(scales.values())))
    q = sum(
        (reference / scales[name]) * banks[name].outer_src
        for name in sources
    ) / len(sources)
    t = banks[target].outer_raw
    qn = row_normalize(q)
    tn = row_normalize(t)
    return Pair(
        full_query=qn,
        full_truth=tn,
        local_query=qn[:, local_feat],
        local_truth=tn[:, local_feat],
    )


def build_candidates() -> tuple[list[Candidate], list[Candidate]]:
    ranks = (8, 16, 32, 64, 128, 256, 384)
    source = [
        Candidate("source", "spectral", rank, gamma=gamma, ridge=ridge)
        for rank in ranks
        for gamma in (
            0.25,
            0.5,
            0.65,
            0.85,
            1.0,
            1.25,
            1.5,
            2.0,
            3.0,
            4.0,
            8.0,
        )
        for ridge in (0.0, 0.1)
    ]
    source.extend(
        Candidate("source", "uniform", rank, eta=eta)
        for rank in ranks
        for eta in (0.25, 0.5, 0.75, 1.0)
    )
    null = [
        Candidate("random", "spectral", rank, gamma=0.65, ridge=0.0)
        for rank in (32, 128, 256, 384)
    ]
    null.extend(
        Candidate("random", "uniform", rank, eta=0.5)
        for rank in (32, 128, 256, 384)
    )
    keys = [candidate.key for candidate in [*source, *null]]
    if len(keys) != len(set(keys)):
        raise RuntimeError("candidate keys are not unique")
    if CURRENT_COMPARATOR_KEY not in set(keys):
        raise RuntimeError("current whitening comparator missing from grid")
    return source, null


def summarize_inner(
    records: dict[str, list[dict]],
    baselines: list[dict],
    source_names: tuple[str, ...],
) -> dict[str, dict]:
    baseline_values = np.asarray([row["score"] for row in baselines])
    out = {}
    for key, rows in records.items():
        scores = np.asarray([row["score"] for row in rows])
        deltas = scores - baseline_values
        donor_delta = {}
        for name in source_names:
            idx = [i for i, row in enumerate(rows) if row["valid_donor"] == name]
            donor_delta[name] = float(np.mean(deltas[idx]))
        out[key] = {
            "mean": float(np.mean(scores)),
            "delta": float(np.mean(deltas)),
            "positive_folds": int(np.sum(deltas > 0)),
            "donor_delta": donor_delta,
            "fold_scores": [float(value) for value in scores],
            "fold_deltas": [float(value) for value in deltas],
        }
    return out


def select_candidate(
    summary: dict[str, dict],
    candidates: list[Candidate],
    min_gain: float,
    min_positive_folds: int,
    tie_tolerance: float,
    require_each_donor: bool,
) -> str:
    by_key = {candidate.key: candidate for candidate in candidates}
    eligible = []
    for candidate in candidates:
        row = summary[candidate.key]
        if row["delta"] < min_gain:
            continue
        if row["positive_folds"] < min_positive_folds:
            continue
        if require_each_donor and not all(
            value > 0 for value in row["donor_delta"].values()
        ):
            continue
        eligible.append((candidate.key, row["delta"]))
    if not eligible:
        return "identity"
    best = max(delta for _, delta in eligible)
    near = [
        key for key, delta in eligible if delta >= best - tie_tolerance
    ]
    return min(near, key=lambda key: by_key[key].complexity)


def pooled_source_only_summary(
    folds: dict[str, dict],
    candidates: list[Candidate],
) -> dict[str, dict]:
    """Pool the six ordered source-to-source non-VCC validation directions."""
    pooled = {}
    for candidate in candidates:
        fold_deltas = []
        donor_deltas = {name: [] for name in OUTER_DONORS}
        for outer_row in folds.values():
            row = outer_row["inner"]["all_candidates"][candidate.key]
            fold_deltas.extend(row["fold_deltas"])
            for donor, delta in row["donor_delta"].items():
                donor_deltas[donor].append(delta)
        pooled[candidate.key] = {
            "delta": float(np.mean(fold_deltas)),
            "positive_folds": int(np.sum(np.asarray(fold_deltas) > 0)),
            "donor_delta": {
                donor: float(np.mean(values))
                for donor, values in donor_deltas.items()
            },
            "fold_deltas": [float(value) for value in fold_deltas],
        }
    return pooled


def run_outer_fold(
    outer: str,
    banks: dict[str, DonorBank],
    local_feat: np.ndarray,
    source_candidates: list[Candidate],
    null_candidates: list[Candidate],
    args: argparse.Namespace,
) -> tuple[dict, dict]:
    sources = tuple(name for name in OUTER_DONORS if name != outer)
    all_candidates = [*source_candidates, *null_candidates]
    max_rank = max(candidate.rank for candidate in all_candidates)
    records = {candidate.key: [] for candidate in all_candidates}
    baselines: list[dict] = []

    for valid_donor in sources:
        fit_donor = next(name for name in sources if name != valid_donor)
        fold_indices = split_folds(
            len(banks[fit_donor].train_src),
            args.folds,
            stable_seed(("inner", outer, valid_donor), args.seed),
        )
        all_rows = np.arange(len(banks[fit_donor].train_src))
        for fold_number, valid_rows in enumerate(fold_indices):
            fit_rows = np.setdiff1d(all_rows, valid_rows, assume_unique=True)
            xfit = scaled_fit_matrix(
                banks, (fit_donor,), fit_rows, local_feat
            )
            subspace = fit_subspace(
                xfit,
                max_rank,
                stable_seed(
                    ("basis", outer, valid_donor, fit_donor),
                    args.seed + fold_number,
                ),
            )
            pair = make_inner_pair(
                banks[fit_donor],
                banks[valid_donor],
                valid_rows,
                local_feat,
            )
            baseline, scores, _ = score_candidates(
                pair, subspace, all_candidates
            )
            baselines.append(
                {
                    "valid_donor": valid_donor,
                    "fit_donor": fit_donor,
                    "fold": fold_number,
                    "score": baseline,
                }
            )
            for candidate in all_candidates:
                records[candidate.key].append(
                    {
                        "valid_donor": valid_donor,
                        "fit_donor": fit_donor,
                        "fold": fold_number,
                        "score": scores[candidate.key],
                    }
                )
            print(
                f"    inner outer={outer:16} valid={valid_donor:16} "
                f"fold={fold_number + 1}/{args.folds} "
                f"identity={baseline:.5f} peak={rss_gib():.2f} GiB",
                flush=True,
            )
            del xfit, subspace, pair
            enforce_soft_memory_limit("inner fold", args.soft_rss_gib)

    summary = summarize_inner(records, baselines, sources)
    selected = select_candidate(
        summary,
        source_candidates,
        args.min_inner_gain,
        args.min_positive_folds,
        args.tie_tolerance,
        require_each_donor=True,
    )
    selected_null = select_candidate(
        summary,
        null_candidates,
        0.0,
        args.min_positive_folds,
        args.tie_tolerance,
        require_each_donor=False,
    )

    all_rows = np.arange(len(next(iter(banks.values())).train_src))
    xfit = scaled_fit_matrix(banks, sources, all_rows, local_feat)
    final_seed = stable_seed(("final", outer, *sources), args.seed)
    final_subspace = fit_subspace(xfit, max_rank, final_seed)
    outer_pair = make_outer_pair(banks, outer, sources, local_feat)
    audit_key = selected if selected != "identity" else CURRENT_COMPARATOR_KEY
    outer_identity, outer_scores, score_matrices = score_candidates(
        outer_pair,
        final_subspace,
        all_candidates,
        matrix_keys={audit_key},
    )
    selected_score = (
        outer_identity if selected == "identity" else outer_scores[selected]
    )
    null_score = (
        outer_identity
        if selected_null == "identity"
        else outer_scores[selected_null]
    )

    # Perturbation labels are absent from the covariance transform. Refit after
    # row permutation and verify that the induced scores are numerically stable.
    rng = np.random.default_rng(stable_seed(("permutation", outer), args.seed))
    permuted = fit_subspace(
        xfit[rng.permutation(len(xfit))],
        max_rank,
        final_seed,
    )
    _, permuted_scores, permuted_matrices = score_candidates(
        outer_pair,
        permuted,
        [
            next(
                candidate
                for candidate in all_candidates
                if candidate.key == audit_key
            )
        ],
        matrix_keys={audit_key},
    )
    permutation_audit = {
        "candidate": audit_key,
        "original_pds": outer_scores[audit_key],
        "permuted_pds": permuted_scores[audit_key],
        "absolute_pds_difference": abs(
            outer_scores[audit_key] - permuted_scores[audit_key]
        ),
        "max_score_matrix_difference": float(
            np.max(
                np.abs(
                    score_matrices[audit_key]
                    - permuted_matrices[audit_key]
                )
            )
        ),
    }

    result = {
        "outer_donor": outer,
        "source_donors": list(sources),
        "inner": {
            "identity_mean": float(np.mean([row["score"] for row in baselines])),
            "selected": selected,
            "selected_summary": (
                None if selected == "identity" else summary[selected]
            ),
            "current_comparator": {
                "candidate": CURRENT_COMPARATOR_KEY,
                **summary[CURRENT_COMPARATOR_KEY],
            },
            "selected_null": selected_null,
            "selected_null_summary": (
                None if selected_null == "identity" else summary[selected_null]
            ),
            "all_candidates": summary,
        },
        "outer": {
            "identity": outer_identity,
            "selected": selected_score,
            "selected_delta": selected_score - outer_identity,
            "current_comparator": outer_scores[CURRENT_COMPARATOR_KEY],
            "current_comparator_delta": (
                outer_scores[CURRENT_COMPARATOR_KEY] - outer_identity
            ),
            "selected_null": null_score,
            "selected_null_delta": null_score - outer_identity,
            "all_scores": outer_scores,
        },
        "permutation_audit": permutation_audit,
    }
    print(
        f"  selected={selected}\n"
        f"  outer identity={outer_identity:.5f} "
        f"selected={selected_score:.5f} ({selected_score-outer_identity:+.5f}) "
        f"current={outer_scores[CURRENT_COMPARATOR_KEY]:.5f} "
        f"null={null_score:.5f}",
        flush=True,
    )
    del xfit, final_subspace, permuted, outer_pair
    enforce_soft_memory_limit("outer fold", args.soft_rss_gib)
    return result, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=48)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--soft-rss-gib", type=float, default=3.5)
    parser.add_argument("--min-inner-gain", type=float, default=0.001)
    parser.add_argument("--min-positive-folds", type=int, default=6)
    parser.add_argument("--tie-tolerance", type=float, default=5e-4)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    if args.folds < 2:
        parser.error("--folds must be at least 2")
    if args.soft_rss_gib >= 4.0:
        parser.error("--soft-rss-gib must stay below the external 4 GiB cap")

    master_genes, train_perts, outer_perts, train_mask = validate_design_files()
    local_feat = np.flatnonzero(train_mask).astype(np.int64)
    vcc_panel = set(
        pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
    )
    print(
        f"design: donors={len(OUTER_DONORS)} train={len(train_perts)} "
        f"outer={len(outer_perts)} score_genes={len(master_genes)} "
        f"shape_genes={len(local_feat)}",
        flush=True,
    )

    banks = {}
    for name, path in DONOR_PATHS.items():
        banks[name] = prepare_donor(
            name,
            path,
            master_genes,
            train_perts,
            outer_perts,
            vcc_panel,
            args.chunk_rows,
        )
        print(
            f"  prepared {name:16} "
            f"background={banks[name].normalization_count:5d} "
            f"floor={banks[name].normalization_floor:.6g} "
            f"peak={rss_gib():.2f} GiB",
            flush=True,
        )
        enforce_soft_memory_limit(f"prepare {name}", args.soft_rss_gib)

    source_candidates, null_candidates = build_candidates()
    print(
        f"grid: source={len(source_candidates)} null={len(null_candidates)} "
        f"max_rank={max(c.rank for c in [*source_candidates, *null_candidates])}",
        flush=True,
    )

    folds = {}
    for outer in OUTER_DONORS:
        print(f"\n=== outer holdout {outer} ===", flush=True)
        folds[outer], _ = run_outer_fold(
            outer,
            banks,
            local_feat,
            source_candidates,
            null_candidates,
            args,
        )

    identity = float(np.mean([
        folds[name]["outer"]["identity"] for name in OUTER_DONORS
    ]))
    selected = float(np.mean([
        folds[name]["outer"]["selected"] for name in OUTER_DONORS
    ]))
    current = float(np.mean([
        folds[name]["outer"]["current_comparator"] for name in OUTER_DONORS
    ]))
    null = float(np.mean([
        folds[name]["outer"]["selected_null"] for name in OUTER_DONORS
    ]))
    deployment_summary = pooled_source_only_summary(folds, source_candidates)
    deployment_selected = select_candidate(
        deployment_summary,
        source_candidates,
        args.min_inner_gain,
        args.min_positive_folds * len(OUTER_DONORS),
        args.tie_tolerance,
        require_each_donor=True,
    )
    deployment_outer_scores = {
        name: (
            folds[name]["outer"]["identity"]
            if deployment_selected == "identity"
            else folds[name]["outer"]["all_scores"][deployment_selected]
        )
        for name in OUTER_DONORS
    }
    deployment_outer_macro = float(np.mean(list(deployment_outer_scores.values())))
    result = {
        "design": {
            "outer_donors": list(OUTER_DONORS),
            "training_perturbations": len(train_perts),
            "outer_vcc_perturbations": len(outer_perts),
            "train_outer_overlap": len(set(train_perts) & set(outer_perts)),
            "score_genes": len(master_genes),
            "shaping_genes": len(local_feat),
            "inner_donor_directions_per_outer": 2,
            "perturbation_folds_per_direction": args.folds,
            "source_candidate_count": len(source_candidates),
            "null_candidate_count": len(null_candidates),
            "preprocessing": (
                "center and gene-normalize each donor on background "
                "perturbations excluding both the 600 selection perturbations "
                "and all VCC targets, streamed from memory-mapped delta.npy"
            ),
            "fit_and_selection_rule": (
                "outer donor excluded; source PCA fit on one source donor's "
                "non-VCC training folds; candidate selected by held-out non-VCC "
                "PDS against the other source donor"
            ),
            "outer_test_rule": (
                "selected transform refit on both source donors' 600 non-VCC "
                "perturbations, then scored once on 259 outer VCC perturbations"
            ),
            "full_gene_eigendecompositions": 0,
            "cell_eval2_runs": 0,
            "submission_files_touched": 0,
        },
        "candidate_grid": {
            "source": [candidate.key for candidate in source_candidates],
            "null": [candidate.key for candidate in null_candidates],
            "current_comparator": CURRENT_COMPARATOR_KEY,
        },
        "folds": folds,
        "macro": {
            "identity": identity,
            "nested_selected": selected,
            "nested_selected_delta": selected - identity,
            "current_comparator": current,
            "current_comparator_delta": current - identity,
            "selected_random_null": null,
            "selected_random_null_delta": null - identity,
            "selected_positive_outer_folds": int(sum(
                folds[name]["outer"]["selected_delta"] > 0
                for name in OUTER_DONORS
            )),
        },
        "deployment_source_only_selection": {
            "selected": deployment_selected,
            "selection_summary": (
                None
                if deployment_selected == "identity"
                else deployment_summary[deployment_selected]
            ),
            "selection_uses": (
                "only the 30 held-out non-VCC source-to-source folds; no VCC "
                "panel score enters selection"
            ),
            "vcc_panel_diagnostic_warning": (
                "These VCC-panel scores are not the strict nested estimate "
                "because all three donors contributed non-VCC data to the "
                "global hyperparameter selection."
            ),
            "vcc_panel_scores": deployment_outer_scores,
            "vcc_panel_macro": deployment_outer_macro,
            "vcc_panel_delta": deployment_outer_macro - identity,
        },
        "verdict_gate": {
            "criterion": (
                "nested-selected source-only transform must improve macro PDS, "
                "be positive in every outer donor, and beat the selected "
                "random-subspace control"
            ),
            "passed": bool(
                selected > identity
                and all(
                    folds[name]["outer"]["selected_delta"] > 0
                    for name in OUTER_DONORS
                )
                and selected > null
            ),
        },
        "peak_rss_gib": rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"\nmacro identity={identity:.5f} selected={selected:.5f} "
        f"delta={selected-identity:+.5f} current={current:.5f} "
        f"null={null:.5f} gate={result['verdict_gate']['passed']}",
        flush=True,
    )
    print(
        f"wrote {args.output}; peak RSS={rss_gib():.2f} GiB",
        flush=True,
    )


if __name__ == "__main__":
    main()
