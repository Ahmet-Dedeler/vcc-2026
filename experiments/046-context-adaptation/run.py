"""046 — supervised target-context adaptation on public donor effects.

For each outer target donor, learn a low-capacity map

    consensus(other public donors) -> target donor

using 600 common, non-VCC perturbations. Hyperparameters are selected only on
held-out folds of those training perturbations. The selected map is then
evaluated once on the 259 VCC-panel perturbations shared by K562-gwps,
Orion-HCT116, and Orion-HEK293T.

The experiment reuses experiment 044's compact memory-mapped cache. It does not
load single-cell data, run cell-eval2, fit a full eigensystem, or submit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Set before importing NumPy so one laptop process cannot create a BLAS storm.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CACHE = ROOT / "experiments/044-csls-distill/cache"

OUTER_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
JURKAT = "nadig_jurkat"
ALL_CACHED_DONORS = (*OUTER_DONORS, "rpe1", JURKAT)


def rss_gib() -> float:
    """Peak process RSS in GiB on macOS."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def read_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


def stable_seed(parts: Iterable[str], extra: int = 0) -> int:
    text = "|".join(parts) + f"|{extra}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], "little")


def row_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-12)


def pds_from_scores(scores: np.ndarray) -> float:
    """VCC PDS rank statistic for a square similarity matrix."""
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    rank = (scores > own).sum(1) + 0.5 * ((scores == own).sum(1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def randomized_right_basis(
    x: np.ndarray,
    rank: int,
    seed: int,
    oversample: int = 8,
    power_iters: int = 1,
) -> np.ndarray:
    """Approximate top right singular vectors without a gene-gene covariance."""
    x = np.asarray(x, dtype=np.float32)
    width = min(rank + oversample, min(x.shape))
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((x.shape[1], width), dtype=np.float32)
    y = x @ omega
    for _ in range(power_iters):
        q, _ = np.linalg.qr(y, mode="reduced")
        y = x @ (x.T @ q)
    q, _ = np.linalg.qr(y, mode="reduced")
    small = q.T @ x
    _, _, vt = np.linalg.svd(small, full_matrices=False)
    return vt[:rank].T.astype(np.float64)


def ridge_solve(
    x: np.ndarray,
    y: np.ndarray,
    ridge: float,
) -> np.ndarray:
    """Ridge regression with a scale relative to the latent design energy."""
    gram = x.T @ x
    scale = float(np.trace(gram) / max(gram.shape[0], 1)) or 1.0
    return np.linalg.solve(
        gram + ridge * scale * np.eye(gram.shape[0]),
        x.T @ y,
    )


@dataclass(frozen=True)
class Candidate:
    family: str
    ridge: float = 0.0
    rank: int = 0
    residual_rank: int = 0
    eta: float = 1.0
    clip: str = ""

    @property
    def key(self) -> str:
        bits = [self.family]
        if self.rank:
            bits.append(f"r={self.rank}")
        if self.residual_rank:
            bits.append(f"rr={self.residual_rank}")
        if self.ridge:
            bits.append(f"ridge={self.ridge:g}")
        if self.eta != 1.0:
            bits.append(f"eta={self.eta:g}")
        if self.clip:
            bits.append(f"clip={self.clip}")
        return ":".join(bits)

    @property
    def complexity(self) -> tuple:
        family_order = {
            "scalar": 0,
            "residual_mean": 1,
            "diagonal": 2,
            "pca_ridge": 3,
            "residual_program": 4,
        }
        return (
            family_order[self.family],
            self.rank + self.residual_rank,
            -self.ridge,
            abs(self.eta),
            self.clip,
        )


@dataclass
class Fitted:
    candidate: Candidate
    scalar: float | None = None
    diagonal: np.ndarray | None = None
    basis: np.ndarray | None = None
    mapping: np.ndarray | None = None
    x_mean: np.ndarray | None = None
    residual_mean: np.ndarray | None = None
    source_basis: np.ndarray | None = None
    residual_basis: np.ndarray | None = None


@dataclass
class Pair:
    full_query: np.ndarray
    full_truth: np.ndarray
    local_query: np.ndarray
    local_truth: np.ndarray
    local_norm: np.ndarray


class Banks:
    """Memory-mapped matrices from experiment 044."""

    def __init__(self) -> None:
        self.arrays: dict[tuple[str, str, str], np.ndarray] = {}
        for name in ALL_CACHED_DONORS:
            for kind in ("raw", "src"):
                self.arrays[(name, kind, "train")] = np.load(
                    CACHE / f"{name}_{kind}_train.npy", mmap_mode="r"
                )
                if name in OUTER_DONORS:
                    self.arrays[(name, kind, "outer")] = np.load(
                        CACHE / f"{name}_{kind}_outer.npy", mmap_mode="r"
                    )

    def get(self, name: str, kind: str, panel: str) -> np.ndarray:
        return self.arrays[(name, kind, panel)]


def validate_cache() -> tuple[list[str], list[str], list[str], np.ndarray]:
    required = [
        CACHE / "meta.json",
        CACHE / "master_genes.txt",
        CACHE / "train_perts.txt",
        CACHE / "outer_perts.txt",
        CACHE / "train_mask.npy",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"experiment 044 cache is incomplete: {missing}")

    meta = json.loads((CACHE / "meta.json").read_text())
    if meta.get("train_count") != 600:
        raise RuntimeError(f"expected 600 cached training perturbations, got {meta}")
    genes = read_lines(CACHE / "master_genes.txt")
    train_perts = read_lines(CACHE / "train_perts.txt")
    outer_perts = read_lines(CACHE / "outer_perts.txt")
    train_mask = np.load(CACHE / "train_mask.npy")
    if len(outer_perts) != 259:
        raise RuntimeError(f"expected 259 outer VCC perturbations, got {len(outer_perts)}")
    if set(train_perts) & set(outer_perts):
        raise RuntimeError("training perturbations overlap the outer VCC panel")
    if train_mask.shape != (len(genes),):
        raise RuntimeError("cached training feature mask has the wrong shape")
    gene_set = set(genes)
    if any(p in gene_set for p in outer_perts):
        raise RuntimeError("VCC target coordinate leaked into the cached score axis")
    gene_pos = {gene: index for index, gene in enumerate(genes)}
    if any(train_mask[gene_pos[p]] for p in train_perts if p in gene_set):
        raise RuntimeError("training target coordinate leaked into the adaptation axis")
    return genes, train_perts, outer_perts, train_mask


def source_scales(
    banks: Banks,
    sources: tuple[str, ...],
    kind: str,
    local_feat: np.ndarray,
) -> dict[str, float]:
    scales = {}
    for name in sources:
        x = np.asarray(
            banks.get(name, kind, "train")[:, local_feat],
            dtype=np.float64,
        )
        scales[name] = float(np.median(np.linalg.norm(x, axis=1))) or 1.0
    return scales


def consensus(
    banks: Banks,
    sources: tuple[str, ...],
    kind: str,
    panel: str,
    feat: np.ndarray,
    norm_feat: np.ndarray,
) -> np.ndarray:
    scales = source_scales(banks, sources, kind, norm_feat)
    reference = float(np.median(list(scales.values())))
    parts = [
        (reference / scales[name])
        * np.asarray(banks.get(name, kind, panel)[:, feat], dtype=np.float64)
        for name in sources
    ]
    return sum(parts) / len(parts)


def make_train_pair(
    banks: Banks,
    target: str,
    sources: tuple[str, ...],
    kind: str,
    local_feat: np.ndarray,
) -> Pair:
    q = consensus(banks, sources, kind, "train", local_feat, local_feat)
    t = np.asarray(
        banks.get(target, "raw", "train")[:, local_feat],
        dtype=np.float64,
    )
    q = row_normalize(q)
    t = row_normalize(t)
    return Pair(q, t, q, t, np.ones(len(q), dtype=np.float64))


def make_outer_pair(
    banks: Banks,
    target: str,
    sources: tuple[str, ...],
    kind: str,
    local_feat: np.ndarray,
    all_feat: np.ndarray,
) -> Pair:
    q = consensus(banks, sources, kind, "outer", all_feat, local_feat)
    t = np.asarray(
        banks.get(target, "raw", "outer")[:, all_feat],
        dtype=np.float64,
    )
    q = row_normalize(q)
    t = row_normalize(t)
    q_block = q[:, local_feat]
    block_norm = np.linalg.norm(q_block, axis=1)
    q_local = q_block / (block_norm[:, None] + 1e-12)
    t_local = t[:, local_feat]
    return Pair(q, t, q_local, t_local, block_norm)


def candidate_grid() -> list[Candidate]:
    ridges = (0.3, 3.0, 30.0)
    etas = (0.25, 0.5, 1.0)
    candidates = [Candidate("scalar", ridge=3.0)]
    candidates.extend(Candidate("residual_mean", eta=eta) for eta in etas)
    candidates.extend(
        Candidate("diagonal", ridge=ridge, eta=eta, clip=clip)
        for ridge in ridges
        for eta in etas
        for clip in ("positive2", "signed3")
    )
    candidates.extend(
        Candidate("pca_ridge", ridge=ridge, rank=rank, eta=eta)
        for rank in (4, 8, 16, 32)
        for ridge in ridges
        for eta in etas
    )
    candidates.extend(
        Candidate(
            "residual_program",
            ridge=ridge,
            rank=rank,
            residual_rank=residual_rank,
            eta=eta,
        )
        for rank in (8, 16, 32)
        for residual_rank in (4, 8, 16)
        for ridge in ridges
        for eta in etas
    )
    return candidates


def fit_models(
    x: np.ndarray,
    y: np.ndarray,
    candidates: list[Candidate],
    seed: int,
) -> dict[str, Fitted]:
    """Fit all candidates while sharing randomized bases across the grid."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    residual = y - x
    models: dict[str, Fitted] = {}

    # Scalar ridge-to-identity is a sanity check. It cannot change PDS when the
    # transformed block is the entire score axis.
    scalar_candidates = [c for c in candidates if c.family == "scalar"]
    if scalar_candidates:
        scalar_candidate = scalar_candidates[0]
        den = float(np.sum(x * x))
        cross = float(np.sum(x * y))
        scalar = (cross + scalar_candidate.ridge * den) / (
            den + scalar_candidate.ridge * den
        )
        models[scalar_candidate.key] = Fitted(
            scalar_candidate,
            scalar=scalar,
        )

    # Gene-wise shrinkage towards the identity map.
    den_g = np.sum(x * x, axis=0)
    cross_g = np.sum(x * y, axis=0)
    live = den_g[den_g > 1e-12]
    scale_g = float(np.median(live)) if live.size else 1.0
    for candidate in (c for c in candidates if c.family == "diagonal"):
        coef = (cross_g + candidate.ridge * scale_g) / (
            den_g + candidate.ridge * scale_g
        )
        if candidate.clip == "positive2":
            coef = np.clip(coef, 0.0, 2.0)
        elif candidate.clip == "signed3":
            coef = np.clip(coef, -1.0, 3.0)
        else:
            raise ValueError(candidate.clip)
        models[candidate.key] = Fitted(candidate, diagonal=coef)

    # A target-specific mean residual is the lowest-capacity additive map.
    residual_mean = residual.mean(axis=0)
    for candidate in (c for c in candidates if c.family == "residual_mean"):
        models[candidate.key] = Fitted(
            candidate,
            residual_mean=residual_mean,
        )

    pca_candidates = [c for c in candidates if c.family == "pca_ridge"]
    if pca_candidates:
        max_pca = max(c.rank for c in pca_candidates)
        joint = np.vstack([x, y]).astype(np.float32, copy=False)
        pca_basis = randomized_right_basis(joint, max_pca, seed + 11)
        del joint
        for rank in sorted({c.rank for c in pca_candidates}):
            basis = pca_basis[:, :rank]
            zx = x @ basis
            zr = residual @ basis
            for ridge in sorted({
                c.ridge for c in pca_candidates if c.rank == rank
            }):
                mapping = ridge_solve(zx, zr, ridge)
                for candidate in (
                    c for c in pca_candidates
                    if c.rank == rank and c.ridge == ridge
                ):
                    models[candidate.key] = Fitted(
                        candidate,
                        basis=basis,
                        mapping=mapping,
                    )

    residual_candidates = [
        c for c in candidates if c.family == "residual_program"
    ]
    if residual_candidates:
        x_mean = x.mean(axis=0)
        residual_center = residual - residual_mean
        x_center = x - x_mean
        max_source_rank = max(c.rank for c in residual_candidates)
        max_residual_rank = max(c.residual_rank for c in residual_candidates)
        source_basis = randomized_right_basis(
            x_center, max_source_rank, seed + 23
        )
        residual_basis = randomized_right_basis(
            residual_center, max_residual_rank, seed + 37
        )
        for rank in sorted({c.rank for c in residual_candidates}):
            vx = source_basis[:, :rank]
            zx = x_center @ vx
            for residual_rank in sorted({
                c.residual_rank for c in residual_candidates if c.rank == rank
            }):
                vr = residual_basis[:, :residual_rank]
                zr = residual_center @ vr
                for ridge in sorted({
                    c.ridge for c in residual_candidates
                    if c.rank == rank and c.residual_rank == residual_rank
                }):
                    mapping = ridge_solve(zx, zr, ridge)
                    for candidate in (
                        c for c in residual_candidates
                        if c.rank == rank
                        and c.residual_rank == residual_rank
                        and c.ridge == ridge
                    ):
                        models[candidate.key] = Fitted(
                            candidate,
                            x_mean=x_mean,
                            residual_mean=residual_mean,
                            source_basis=vx,
                            residual_basis=vr,
                            mapping=mapping,
                        )

    if len(models) != len(candidates):
        raise RuntimeError(f"fit {len(models)} models for {len(candidates)} candidates")
    return models


def model_scores(pair: Pair, model: Fitted) -> np.ndarray:
    """Cosine ranking scores without materializing normalized adapted vectors.

    Query normalization divides every entry in one row by the same positive
    number, so it cannot change that row's rank. Truth vectors are already
    normalized on the complete score axis.
    """
    base = pair.full_query @ pair.full_truth.T
    x = pair.local_query
    y = pair.local_truth
    scale = pair.local_norm[:, None]
    c = model.candidate

    if c.family == "scalar":
        delta = (float(model.scalar) - 1.0) * x
        return base + (scale * delta) @ y.T
    if c.family == "diagonal":
        factor = 1.0 + c.eta * (model.diagonal - 1.0)
        delta = x * (factor - 1.0)
        return base + (scale * delta) @ y.T
    if c.family == "residual_mean":
        correction = c.eta * model.residual_mean
        return base + scale * (correction @ y.T)[None, :]
    if c.family == "pca_ridge":
        zx = x @ model.basis
        zy = y @ model.basis
        latent_delta = c.eta * (zx @ model.mapping)
        return base + (scale * latent_delta) @ zy.T
    if c.family == "residual_program":
        zx = (x - model.x_mean) @ model.source_basis
        zy = y @ model.residual_basis
        latent_delta = zx @ model.mapping
        mean_scores = model.residual_mean @ y.T
        return base + c.eta * (
            scale * mean_scores[None, :]
            + (scale * latent_delta) @ zy.T
        )
    raise ValueError(c.family)


def apply_local(x: np.ndarray, model: Fitted) -> np.ndarray:
    """Materialize the adapted local vectors for diagnostics or serialization."""
    c = model.candidate
    if c.family == "scalar":
        return float(model.scalar) * x
    if c.family == "diagonal":
        factor = 1.0 + c.eta * (model.diagonal - 1.0)
        return x * factor
    if c.family == "residual_mean":
        return x + c.eta * model.residual_mean
    if c.family == "pca_ridge":
        return x + c.eta * (((x @ model.basis) @ model.mapping) @ model.basis.T)
    if c.family == "residual_program":
        predicted = (
            (x - model.x_mean)
            @ model.source_basis
            @ model.mapping
            @ model.residual_basis.T
        )
        return x + c.eta * (model.residual_mean + predicted)
    raise ValueError(c.family)


def split_folds(n: int, folds: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    return [np.sort(part) for part in np.array_split(order, folds)]


def cv_target(
    pair: Pair,
    target: str,
    candidates: list[Candidate],
    folds: int,
    seed: int,
) -> dict:
    fold_indices = split_folds(len(pair.full_query), folds, seed)
    candidate_scores = {candidate.key: [] for candidate in candidates}
    baseline_scores = []
    all_rows = np.arange(len(pair.full_query))

    for fold_index, valid in enumerate(fold_indices):
        fit = np.setdiff1d(all_rows, valid, assume_unique=True)
        models = fit_models(
            pair.local_query[fit],
            pair.local_truth[fit],
            candidates,
            stable_seed((target, str(seed)), extra=fold_index),
        )
        valid_pair = Pair(
            pair.full_query[valid],
            pair.full_truth[valid],
            pair.local_query[valid],
            pair.local_truth[valid],
            pair.local_norm[valid],
        )
        baseline = pds_from_scores(
            valid_pair.full_query @ valid_pair.full_truth.T
        )
        baseline_scores.append(baseline)
        for key, model in models.items():
            candidate_scores[key].append(pds_from_scores(model_scores(valid_pair, model)))
        print(
            f"    {target:16} fold {fold_index + 1}/{folds} "
            f"baseline={baseline:.5f} peak_rss={rss_gib():.2f} GiB",
            flush=True,
        )

    base_mean = float(np.mean(baseline_scores))
    summary = {}
    for key, scores in candidate_scores.items():
        deltas = np.asarray(scores) - np.asarray(baseline_scores)
        summary[key] = {
            "mean": float(np.mean(scores)),
            "delta": float(np.mean(deltas)),
            "folds": [float(value) for value in scores],
            "fold_deltas": [float(value) for value in deltas],
            "positive_folds": int(np.sum(deltas > 0)),
        }
    return {
        "baseline_mean": base_mean,
        "baseline_folds": baseline_scores,
        "candidates": summary,
    }


def select_candidate(
    cv: dict[str, dict],
    targets: tuple[str, ...],
    candidates: list[Candidate],
    min_gain: float,
    tie_tolerance: float,
    require_each_target: bool,
) -> str:
    by_key = {candidate.key: candidate for candidate in candidates}
    eligible = []
    for candidate in candidates:
        gains = [cv[target]["candidates"][candidate.key]["delta"] for target in targets]
        gain = float(np.mean(gains))
        if gain < min_gain:
            continue
        if require_each_target and not all(value > 0 for value in gains):
            continue
        eligible.append((candidate.key, gain))
    if not eligible:
        return "baseline"
    best_gain = max(gain for _, gain in eligible)
    near = [
        key for key, gain in eligible
        if gain >= best_gain - tie_tolerance
    ]
    return min(near, key=lambda key: by_key[key].complexity)


def select_for_one_target(
    result: dict,
    candidates: list[Candidate],
    min_gain: float,
    tie_tolerance: float,
    min_positive_folds: int,
) -> str:
    by_key = {candidate.key: candidate for candidate in candidates}
    eligible = []
    for candidate in candidates:
        row = result["candidates"][candidate.key]
        if row["delta"] < min_gain:
            continue
        if row["positive_folds"] < min_positive_folds:
            continue
        eligible.append((candidate.key, row["delta"]))
    if not eligible:
        return "baseline"
    best_gain = max(gain for _, gain in eligible)
    near = [
        key for key, gain in eligible
        if gain >= best_gain - tie_tolerance
    ]
    return min(near, key=lambda key: by_key[key].complexity)


def outer_scores(
    pair: Pair,
    models: dict[str, Fitted],
) -> tuple[float, dict[str, float]]:
    baseline = pds_from_scores(pair.full_query @ pair.full_truth.T)
    scores = {
        key: pds_from_scores(model_scores(pair, model))
        for key, model in models.items()
    }
    return baseline, scores


def serialize_model(
    model: Fitted,
    path: Path,
    metadata: dict,
) -> None:
    arrays = {
        "scalar": np.asarray(
            [] if model.scalar is None else [model.scalar],
            dtype=np.float64,
        ),
        "diagonal": np.asarray(
            [] if model.diagonal is None else model.diagonal,
            dtype=np.float32,
        ),
        "basis": np.asarray(
            [] if model.basis is None else model.basis,
            dtype=np.float32,
        ),
        "mapping": np.asarray(
            [] if model.mapping is None else model.mapping,
            dtype=np.float64,
        ),
        "x_mean": np.asarray(
            [] if model.x_mean is None else model.x_mean,
            dtype=np.float32,
        ),
        "residual_mean": np.asarray(
            [] if model.residual_mean is None else model.residual_mean,
            dtype=np.float32,
        ),
        "source_basis": np.asarray(
            [] if model.source_basis is None else model.source_basis,
            dtype=np.float32,
        ),
        "residual_basis": np.asarray(
            [] if model.residual_basis is None else model.residual_basis,
            dtype=np.float32,
        ),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **arrays)
    path.with_suffix(".json").write_text(json.dumps(metadata, indent=2))


def run_representation(
    banks: Banks,
    kind: str,
    candidates: list[Candidate],
    local_feat: np.ndarray,
    all_feat: np.ndarray,
    train_perts: list[str],
    outer_perts: list[str],
    args: argparse.Namespace,
) -> tuple[dict, dict[str, Pair], dict[str, dict[str, Fitted]]]:
    print(f"\n=== source representation: {kind} ===", flush=True)
    train_pairs = {}
    outer_pairs = {}
    cv = {}

    for target in OUTER_DONORS:
        sources = tuple(name for name in OUTER_DONORS if name != target)
        train_pairs[target] = make_train_pair(
            banks, target, sources, kind, local_feat
        )
        outer_pairs[target] = make_outer_pair(
            banks, target, sources, kind, local_feat, all_feat
        )
        print(
            f"  CV target={target}, sources={sources}, "
            f"train={len(train_perts)}x{len(local_feat)}",
            flush=True,
        )
        cv[target] = cv_target(
            train_pairs[target],
            target,
            candidates,
            args.folds,
            args.seed + stable_seed((kind, target)),
        )

    target_selected = {
        target: select_for_one_target(
            cv[target],
            candidates,
            args.min_cv_gain,
            args.tie_tolerance,
            args.min_positive_folds,
        )
        for target in OUTER_DONORS
    }
    donor_nested_selected = {
        target: select_candidate(
            cv,
            tuple(name for name in OUTER_DONORS if name != target),
            candidates,
            args.min_cv_gain,
            args.tie_tolerance,
            require_each_target=True,
        )
        for target in OUTER_DONORS
    }
    global_selected = select_candidate(
        cv,
        OUTER_DONORS,
        candidates,
        args.min_cv_gain,
        args.tie_tolerance,
        require_each_target=True,
    )

    fitted_full = {}
    outer = {}
    fixed_scores = {candidate.key: {} for candidate in candidates}
    for target in OUTER_DONORS:
        fitted_full[target] = fit_models(
            train_pairs[target].local_query,
            train_pairs[target].local_truth,
            candidates,
            stable_seed((kind, target, "full"), extra=args.seed),
        )
        baseline, scores = outer_scores(outer_pairs[target], fitted_full[target])
        for key, score in scores.items():
            fixed_scores[key][target] = score

        def selected_score(key: str) -> float:
            return baseline if key == "baseline" else scores[key]

        target_key = target_selected[target]
        nested_key = donor_nested_selected[target]
        outer[target] = {
            "baseline": baseline,
            "target_cv_selected": target_key,
            "target_cv_score": selected_score(target_key),
            "target_cv_delta": selected_score(target_key) - baseline,
            "donor_nested_selected": nested_key,
            "donor_nested_score": selected_score(nested_key),
            "donor_nested_delta": selected_score(nested_key) - baseline,
            "global_selected": global_selected,
            "global_score": selected_score(global_selected),
            "global_delta": selected_score(global_selected) - baseline,
        }
        print(
            f"  OUTER {target:16} base={baseline:.5f} "
            f"target-CV={target_key} {outer[target]['target_cv_delta']:+.5f} "
            f"donor-nested={nested_key} {outer[target]['donor_nested_delta']:+.5f} "
            f"global={global_selected} {outer[target]['global_delta']:+.5f}",
            flush=True,
        )

    def macro(field: str) -> float:
        return float(np.mean([outer[target][field] for target in OUTER_DONORS]))

    fixed_summary = {}
    baseline_macro = macro("baseline")
    for key, folds in fixed_scores.items():
        value = float(np.mean(list(folds.values())))
        fixed_summary[key] = {
            **folds,
            "macro": value,
            "delta": value - baseline_macro,
        }
    optimistic = sorted(
        fixed_summary.items(),
        key=lambda item: item[1]["macro"],
        reverse=True,
    )

    result = {
        "source_kind": kind,
        "design": {
            "training_perturbations": len(train_perts),
            "outer_vcc_perturbations": len(outer_perts),
            "overlap": len(set(train_perts) & set(outer_perts)),
            "score_genes": len(all_feat),
            "adaptation_genes": len(local_feat),
            "folds": args.folds,
            "candidate_count": len(candidates),
            "full_eigendecompositions": 0,
        },
        "cv": cv,
        "selection": {
            "target_specific": target_selected,
            "donor_nested": donor_nested_selected,
            "global": global_selected,
        },
        "outer": outer,
        "macro": {
            "baseline": baseline_macro,
            "target_cv": macro("target_cv_score"),
            "target_cv_delta": macro("target_cv_delta"),
            "donor_nested": macro("donor_nested_score"),
            "donor_nested_delta": macro("donor_nested_delta"),
            "global": macro("global_score"),
            "global_delta": macro("global_delta"),
        },
        "strict_transfer_gate": {
            "criterion": (
                "target-CV and donor-nested deltas must both be positive in all "
                f"three outer donors, and both macro gains must be >= {args.ship_gain:g}"
            ),
            "target_cv_all_positive": all(
                outer[target]["target_cv_delta"] > 0 for target in OUTER_DONORS
            ),
            "donor_nested_all_positive": all(
                outer[target]["donor_nested_delta"] > 0 for target in OUTER_DONORS
            ),
            "passed": (
                all(outer[target]["target_cv_delta"] > 0 for target in OUTER_DONORS)
                and all(
                    outer[target]["donor_nested_delta"] > 0
                    for target in OUTER_DONORS
                )
                and macro("target_cv_delta") >= args.ship_gain
                and macro("donor_nested_delta") >= args.ship_gain
            ),
        },
        "fixed_grid_outer_optimistic": {
            "warning": (
                "Exploratory only. These candidates are ranked on the same VCC "
                "outer folds and cannot be used as an unbiased selection result."
            ),
            "top": [
                {"candidate": key, **row}
                for key, row in optimistic[:20]
            ],
        },
    }
    return result, train_pairs, fitted_full


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-kinds",
        nargs="+",
        choices=["raw", "src"],
        default=["raw", "src"],
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=46)
    parser.add_argument("--min-cv-gain", type=float, default=0.001)
    parser.add_argument("--tie-tolerance", type=float, default=0.0005)
    parser.add_argument("--min-positive-folds", type=int, default=3)
    parser.add_argument("--ship-gain", type=float, default=0.002)
    parser.add_argument(
        "--fit-jurkat",
        action="store_true",
        help=(
            "write a Jurkat artifact only if the outer transfer gate passes and "
            "results_controls.json independently establishes label dependence"
        ),
    )
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()

    genes, train_perts, outer_perts, train_mask = validate_cache()
    all_feat = np.arange(len(genes), dtype=np.int64)
    local_feat = np.flatnonzero(train_mask).astype(np.int64)
    candidates = candidate_grid()
    banks = Banks()
    print(
        f"cache: train={len(train_perts)}, outer={len(outer_perts)}, "
        f"score_genes={len(all_feat)}, adaptation_genes={len(local_feat)}\n"
        f"candidate grid={len(candidates)}, peak_rss={rss_gib():.2f} GiB",
        flush=True,
    )

    results = {}
    for kind in args.source_kinds:
        result, _, _ = run_representation(
            banks,
            kind,
            candidates,
            local_feat,
            all_feat,
            train_perts,
            outer_perts,
            args,
        )
        results[kind] = result

    passed = [
        kind for kind in args.source_kinds
        if results[kind]["strict_transfer_gate"]["passed"]
    ]
    artifact = None
    control_path = HERE / "results_controls.json"
    control_gate = None
    if control_path.exists():
        control_gate = json.loads(control_path.read_text()).get(
            "context_specificity_gate", {}
        )
    may_fit_jurkat = (
        bool(passed)
        and args.fit_jurkat
        and bool(control_gate)
        and control_gate.get("passed") is True
    )
    if may_fit_jurkat:
        best_kind = max(
            passed,
            key=lambda kind: results[kind]["macro"]["donor_nested_delta"],
        )
        # The adaptation procedure transferred. Select Jurkat's candidate using
        # Jurkat's own non-VCC perturbation CV, exactly as for each outer donor.
        sources = OUTER_DONORS
        jurkat_pair = make_train_pair(
            banks, JURKAT, sources, best_kind, local_feat
        )
        print(
            f"\nstrict gate passed for {best_kind}; fitting analogous "
            f"{sources} -> {JURKAT} map",
            flush=True,
        )
        jurkat_cv = cv_target(
            jurkat_pair,
            JURKAT,
            candidates,
            args.folds,
            args.seed + stable_seed((best_kind, JURKAT)),
        )
        jurkat_selected = select_for_one_target(
            jurkat_cv,
            candidates,
            args.min_cv_gain,
            args.tie_tolerance,
            args.min_positive_folds,
        )
        if jurkat_selected != "baseline":
            models = fit_models(
                jurkat_pair.local_query,
                jurkat_pair.local_truth,
                candidates,
                stable_seed((best_kind, JURKAT, "full"), extra=args.seed),
            )
            model = models[jurkat_selected]
            artifact_path = HERE / "artifacts" / f"{best_kind}_to_jurkat.npz"
            metadata = {
                "experiment": 46,
                "source_kind": best_kind,
                "sources": list(sources),
                "target": JURKAT,
                "candidate": jurkat_selected,
                "training_perturbations": train_perts,
                "adaptation_genes": [genes[index] for index in local_feat],
                "score_axis_sha256": hashlib.sha256(
                    "\n".join(genes).encode()
                ).hexdigest(),
                "training_perturbations_sha256": hashlib.sha256(
                    "\n".join(train_perts).encode()
                ).hexdigest(),
                "source_scales": source_scales(
                    banks, sources, best_kind, local_feat
                ),
                "source_scale_reference": float(np.median(list(
                    source_scales(
                        banks, sources, best_kind, local_feat
                    ).values()
                ))),
                "input_contract": (
                    "For each source, use the centered donor effect from "
                    "Donor.decompose(); multiply by source_scale_reference / "
                    "source_scales[source]; average sources; restrict to "
                    "adaptation_genes; row-normalize; apply the stored map. "
                    "Embed the adapted block back with its original block norm."
                ),
                "target_coordinate_policy": (
                    "all 600 training targets and all 300 VCC targets excluded "
                    "before fitting the dense/low-rank map"
                ),
                "jurkat_cv": jurkat_cv,
            }
            serialize_model(model, artifact_path, metadata)
            artifact = str(artifact_path.relative_to(ROOT))
            print(f"  wrote Jurkat artifact: {artifact}", flush=True)
        else:
            print(
                "  outer gate passed, but Jurkat CV selected baseline; "
                "no Jurkat artifact written",
                flush=True,
            )
            results[best_kind]["jurkat_cv"] = jurkat_cv
    elif passed:
        print(
            f"\nouter transfer gate passed for {passed}, but no Jurkat artifact "
            "was written. --fit-jurkat and a passing label-dependence control "
            "are both required.",
            flush=True,
        )

    payload = {
        "experiment": 46,
        "cache": str(CACHE.relative_to(ROOT)),
        "representations": results,
        "jurkat_artifact": artifact,
        "context_specificity_gate": control_gate,
        "peak_rss_gib": rss_gib(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2))
    print(
        f"\nwrote {args.output}; peak RSS={rss_gib():.2f} GiB; "
        f"jurkat_artifact={artifact}",
        flush=True,
    )


if __name__ == "__main__":
    main()
