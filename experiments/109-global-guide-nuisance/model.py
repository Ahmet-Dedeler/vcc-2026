#!/usr/bin/env python3
"""Pure numerical model for experiment 109.

The only external code loaded here is experiment 083's source-only numerical
helper module. No expression file, held-context truth, network resource, or
challenge artifact is opened by importing this module.
"""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.util
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP083_MODEL_PATH = ROOT / "experiments/083-guide-latent-decontamination/model.py"
SPLIT_SEED = 1_092_026
STABLE_SINGULAR_RATIO = 1e-6
EPS = 1e-12


def sha256_axis(values: Sequence[str]) -> str:
    """Hash one ordered string axis, including its terminal delimiter."""
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


@dataclasses.dataclass(frozen=True)
class AxisContract:
    """Exact count/hash contract for an ordered string axis."""

    count: int
    sha256: str

    @classmethod
    def from_values(cls, values: Sequence[str]) -> AxisContract:
        ordered = tuple(map(str, values))
        return cls(count=len(ordered), sha256=sha256_axis(ordered))


def validate_axis(
    values: Sequence[str],
    contract: AxisContract,
    *,
    label: str,
    require_unique: bool = True,
    require_sorted: bool = False,
) -> tuple[str, ...]:
    """Authenticate one axis before it controls fitting or centering."""
    ordered = tuple(map(str, values))
    if len(ordered) != int(contract.count):
        raise ValueError(f"{label} count changed: {len(ordered)} != {contract.count}")
    if sha256_axis(ordered) != str(contract.sha256):
        raise ValueError(f"{label} hash changed")
    if require_unique and len(set(ordered)) != len(ordered):
        raise ValueError(f"{label} contains duplicates")
    if require_sorted and ordered != tuple(sorted(ordered)):
        raise ValueError(f"{label} must be sorted")
    return ordered


@dataclasses.dataclass(frozen=True)
class AuthenticatedResponseAxis:
    """Response coordinates derived from exact gene and target axes."""

    gene_axis: tuple[str, ...]
    panel_targets: tuple[str, ...]
    development_targets: tuple[str, ...]
    response_mask: np.ndarray
    gene_axis_contract: AxisContract
    panel_targets_contract: AxisContract
    development_targets_contract: AxisContract
    response_axis_contract: AxisContract


def authenticate_response_axis(
    gene_axis: Sequence[str],
    panel_targets: Sequence[str],
    development_targets: Sequence[str],
    *,
    gene_axis_contract: AxisContract,
    panel_targets_contract: AxisContract,
    development_targets_contract: AxisContract,
    response_axis_contract: AxisContract,
) -> AuthenticatedResponseAxis:
    """Construct the fit mask only from pinned full axes.

    The mask is never accepted directly from a caller. Every panel and
    development target coordinate must exist on the gene axis and is removed
    before the learned response is formed.
    """
    genes = validate_axis(gene_axis, gene_axis_contract, label="gene axis")
    panel = validate_axis(
        panel_targets,
        panel_targets_contract,
        label="panel target axis",
        require_sorted=True,
    )
    development = validate_axis(
        development_targets,
        development_targets_contract,
        label="development target axis",
        require_sorted=True,
    )
    overlap = set(panel) & set(development)
    if overlap:
        raise ValueError(
            f"panel/development target axes overlap: {sorted(overlap)[:5]}"
        )
    excluded = set(panel) | set(development)
    mask = np.asarray([gene not in excluded for gene in genes], dtype=bool)
    if not np.any(mask):
        raise ValueError("response axis is empty after exclusions")
    response = tuple(gene for gene, keep in zip(genes, mask, strict=True) if keep)
    validate_axis(
        response,
        response_axis_contract,
        label="retained response axis",
    )
    mask.setflags(write=False)
    return AuthenticatedResponseAxis(
        gene_axis=genes,
        panel_targets=panel,
        development_targets=development,
        response_mask=mask,
        gene_axis_contract=gene_axis_contract,
        panel_targets_contract=panel_targets_contract,
        development_targets_contract=development_targets_contract,
        response_axis_contract=response_axis_contract,
    )


def _load_exp083_model():
    name = "exp109_reused_exp083_model"
    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, EXP083_MODEL_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {EXP083_MODEL_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP083 = _load_exp083_model()


def deterministic_target_split(
    targets: Iterable[str],
    *,
    seed: int = SPLIT_SEED,
) -> dict[str, tuple[str, ...]]:
    """Reuse experiment 083's deterministic target-level split."""
    return EXP083.deterministic_target_split(targets, seed=seed)


def assert_split_contract(
    splits: Mapping[str, Sequence[str]],
    panel_targets: Iterable[str],
) -> None:
    """Reuse experiment 083's disjointness and panel-leakage guard."""
    EXP083.assert_split_contract(splits, panel_targets)


def _center_by_target(values: np.ndarray, labels: Sequence[str]) -> np.ndarray:
    matrix = np.asarray(values, dtype=np.float64)
    target_labels = np.asarray(tuple(map(str, labels)), dtype=object)
    if matrix.ndim != 2 or len(matrix) != len(target_labels):
        raise ValueError("value/target row mismatch")
    centered = np.empty_like(matrix)
    for target in sorted(set(target_labels.tolist())):
        rows = np.flatnonzero(target_labels == target)
        centered[rows] = matrix[rows] - matrix[rows].mean(axis=0, keepdims=True)
    return centered


def _fit_rows(
    guide_targets: Sequence[str],
    identifiable: np.ndarray,
    fit_targets: Iterable[str],
) -> tuple[np.ndarray, tuple[str, ...]]:
    labels = tuple(map(str, guide_targets))
    identifiable = np.asarray(identifiable, dtype=bool)
    if identifiable.shape != (len(labels),):
        raise ValueError("target/identifiable length mismatch")
    fit_set = set(map(str, fit_targets))
    grouped = EXP083.target_rows(labels)
    eligible = {
        target
        for target in fit_set
        if target in grouped and int(np.sum(identifiable[grouped[target]])) >= 2
    }
    rows = np.asarray(
        [
            index
            for index, target in enumerate(labels)
            if identifiable[index] and target in eligible
        ],
        dtype=np.int64,
    )
    if len(rows) < 4 or len(eligible) < 2:
        raise ValueError("too few identifiable multi-guide fit targets")
    return rows, tuple(labels[index] for index in rows)


@dataclasses.dataclass(frozen=True, order=True)
class Candidate:
    response_rank: int
    ridge: float
    eta: float
    weight_power: float
    huber_c: float = 2.0

    @property
    def is_identity(self) -> bool:
        return self.eta == 0.0


def candidate_grid() -> tuple[Candidate, ...]:
    """The complete frozen candidate grid, including one exact identity."""
    values = [Candidate(0, 0.0, 0.0, 0.0, 2.0)]
    for rank in (8, 16, 32):
        for ridge in (0.1, 1.0, 10.0, 100.0):
            for eta in (0.25, 0.5, 0.75, 1.0):
                values.append(Candidate(rank, ridge, eta, 0.0, 2.0))
    return tuple(values)


@dataclasses.dataclass(frozen=True)
class StableSubspaceModel:
    """Low-rank response map restricted to identifiable feature contrasts."""

    candidate: Candidate
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    feature_keep: np.ndarray
    stable_feature_basis: np.ndarray
    stable_singular_values: np.ndarray
    response_mask: np.ndarray
    gene_axis_sha256: str
    panel_targets_sha256: str
    development_targets_sha256: str
    response_axis_sha256: str
    fit_targets_sha256: str
    response_basis: np.ndarray
    coefficients: np.ndarray

    @property
    def stable_rank(self) -> int:
        return int(self.stable_feature_basis.shape[1])

    @property
    def response_width(self) -> int:
        return int(np.sum(self.response_mask))

    def standardized_features(self, features: np.ndarray) -> np.ndarray:
        values = np.asarray(features, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != len(self.feature_keep):
            raise ValueError("feature matrix shape mismatch")
        kept = values[:, self.feature_keep]
        return (kept - self.feature_mean) / self.feature_scale

    def stable_scores(self, features: np.ndarray) -> np.ndarray:
        return self.standardized_features(features) @ self.stable_feature_basis

    def projected_standardized_features(self, features: np.ndarray) -> np.ndarray:
        scores = self.stable_scores(features)
        return scores @ self.stable_feature_basis.T

    def nuisance(self, features: np.ndarray) -> np.ndarray:
        if self.candidate.is_identity:
            return np.zeros((len(features), self.response_width), dtype=np.float64)
        return self.stable_scores(features) @ self.coefficients @ self.response_basis


def fit_stable_subspace_model(
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    identifiable: np.ndarray,
    fit_targets: Iterable[str],
    response_axis: AuthenticatedResponseAxis,
    candidate: Candidate,
    *,
    stable_ratio: float = STABLE_SINGULAR_RATIO,
    seed: int = SPLIT_SEED,
    fixed_response_basis: np.ndarray | None = None,
) -> StableSubspaceModel:
    """Fit only on rows belonging to the supplied fit targets.

    Feature standardization, feature SVD, response SVD, and ridge regression
    are all isolated to fit-target guide rows.
    """
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    if not isinstance(response_axis, AuthenticatedResponseAxis):
        raise TypeError("response_axis must be authenticated")
    verified_response_axis = authenticate_response_axis(
        response_axis.gene_axis,
        response_axis.panel_targets,
        response_axis.development_targets,
        gene_axis_contract=response_axis.gene_axis_contract,
        panel_targets_contract=response_axis.panel_targets_contract,
        development_targets_contract=response_axis.development_targets_contract,
        response_axis_contract=response_axis.response_axis_contract,
    )
    if not np.array_equal(
        response_axis.response_mask,
        verified_response_axis.response_mask,
    ):
        raise ValueError("caller-supplied response mask differs from pinned axes")
    response_axis = verified_response_axis
    response_mask = np.asarray(response_axis.response_mask, dtype=bool)
    if effects.ndim != 2 or features.ndim != 2:
        raise ValueError("effects/features must be matrices")
    if len(effects) != len(features) or len(effects) != len(guide_targets):
        raise ValueError("guide arrays have inconsistent rows")
    if effects.shape[1] != len(response_axis.gene_axis):
        raise ValueError("effect and authenticated gene axes differ")
    if response_mask.shape != (effects.shape[1],):
        raise ValueError("authenticated response mask shape changed")
    if not 0.0 < stable_ratio < 1.0:
        raise ValueError("stable_ratio must be in (0, 1)")

    requested_fit_targets = tuple(sorted(set(map(str, fit_targets))))
    rows, labels = _fit_rows(
        guide_targets,
        identifiable,
        requested_fit_targets,
    )
    observed_fit_targets = tuple(sorted(set(labels)))
    if observed_fit_targets != requested_fit_targets:
        missing = sorted(set(requested_fit_targets) - set(observed_fit_targets))
        raise ValueError(
            "one or more requested fit targets lack two identifiable guides: "
            f"{missing[:5]}"
        )
    fit_features = np.asarray(features[rows], dtype=np.float64)
    feature_mean_all = fit_features.mean(axis=0)
    feature_scale_all = fit_features.std(axis=0)
    feature_keep = np.isfinite(feature_scale_all) & (feature_scale_all > 1e-8)
    if not np.any(feature_keep):
        raise ValueError("no varying fit features")
    feature_mean = feature_mean_all[feature_keep]
    feature_scale = feature_scale_all[feature_keep]
    standardized = (fit_features[:, feature_keep] - feature_mean) / feature_scale
    centered_x = _center_by_target(standardized, labels)

    _, singular_values, vt = np.linalg.svd(centered_x, full_matrices=False)
    if len(singular_values) == 0 or singular_values[0] <= EPS:
        raise ValueError("within-target feature contrasts have zero rank")
    stable = singular_values / singular_values[0] > stable_ratio
    if not np.any(stable):
        raise ValueError("no stable within-target feature directions")
    stable_basis = np.asarray(vt[stable].T, dtype=np.float64)
    stable_singular_values = np.asarray(singular_values[stable], dtype=np.float64)
    stable_x = centered_x @ stable_basis

    response = np.asarray(
        effects[np.ix_(rows, np.flatnonzero(response_mask))],
        dtype=np.float64,
    )
    centered_y = _center_by_target(response, labels)

    if candidate.is_identity:
        response_basis = np.zeros((0, response.shape[1]), dtype=np.float64)
        coefficients = np.zeros((stable_basis.shape[1], 0), dtype=np.float64)
    else:
        max_rank = min(centered_y.shape)
        if not 0 < candidate.response_rank <= max_rank:
            raise ValueError(
                f"response rank {candidate.response_rank} outside 1..{max_rank}"
            )
        if fixed_response_basis is None:
            response_basis = EXP083.randomized_right_basis(
                centered_y,
                candidate.response_rank,
                seed=EXP083.stable_seed(
                    seed, "exp109-response-basis", candidate.response_rank
                ),
            )
        else:
            response_basis = np.asarray(
                fixed_response_basis[: candidate.response_rank],
                dtype=np.float64,
            )
            if response_basis.shape != (
                candidate.response_rank,
                response.shape[1],
            ):
                raise ValueError("fixed response basis shape mismatch")
        scores = centered_y @ response_basis.T
        gram = stable_x.T @ stable_x
        gram.flat[:: len(gram) + 1] += candidate.ridge
        coefficients = np.linalg.solve(gram, stable_x.T @ scores)

    return StableSubspaceModel(
        candidate=candidate,
        feature_mean=feature_mean,
        feature_scale=feature_scale,
        feature_keep=feature_keep,
        stable_feature_basis=stable_basis,
        stable_singular_values=stable_singular_values,
        response_mask=response_mask,
        gene_axis_sha256=response_axis.gene_axis_contract.sha256,
        panel_targets_sha256=response_axis.panel_targets_contract.sha256,
        development_targets_sha256=(response_axis.development_targets_contract.sha256),
        response_axis_sha256=response_axis.response_axis_contract.sha256,
        fit_targets_sha256=sha256_axis(requested_fit_targets),
        response_basis=response_basis,
        coefficients=coefficients,
    )


def centered_singleton_corrections(
    nuisance: np.ndarray,
    target_ids: Sequence[str],
    singleton_supported: np.ndarray,
    *,
    eta: float,
    centering_axis_contract: AxisContract,
) -> np.ndarray:
    """Return globally centered corrections on one authenticated population."""
    values = np.asarray(nuisance, dtype=np.float64)
    targets = tuple(map(str, target_ids))
    support = np.asarray(singleton_supported, dtype=bool)
    if (
        values.ndim != 2
        or support.shape != (len(values),)
        or len(targets) != len(values)
    ):
        raise ValueError("nuisance/target/support shape mismatch")
    if len(set(targets)) != len(targets):
        raise ValueError("target-level centering rows contain duplicates")
    if not np.isfinite(values).all():
        raise ValueError("non-finite nuisance prediction")
    output = np.zeros_like(values)
    selected = np.flatnonzero(support)
    selected_targets = tuple(targets[index] for index in selected)
    validate_axis(
        selected_targets,
        centering_axis_contract,
        label="singleton centering axis",
        require_sorted=True,
    )
    if eta == 0.0:
        return output
    if len(selected) < 2:
        raise ValueError("need at least two singleton targets for global centering")
    block = values[selected]
    # Remove the arbitrary gauge before averaging. Centering direct absolute
    # predictions is mathematically equivalent but loses avoidable precision
    # when a large constant response vector is added to every row.
    relative = block - block[:1]
    output[selected] = -float(eta) * (relative - relative.mean(axis=0, keepdims=True))

    # Reclose after floating-point centering without concentrating roundoff in
    # one biological target. Three residual-distribution passes keep both the
    # zero-mean and gauge-invariance audits below their frozen tolerances at
    # the full 6,689-target centering scale.
    for _ in range(3):
        residual = np.sum(output[selected], axis=0, dtype=np.float64)
        output[selected] -= residual / len(selected)
    output[~support] = 0.0
    return output


def deploy_singleton_corrections(
    fitted: StableSubspaceModel,
    target_features: np.ndarray,
    target_ids: Sequence[str],
    singleton_supported: np.ndarray,
    *,
    centering_axis_contract: AxisContract,
    eta: float | None = None,
) -> np.ndarray:
    """Predict and globally center target-level singleton corrections."""
    chosen_eta = fitted.candidate.eta if eta is None else float(eta)
    nuisance = fitted.nuisance(target_features)
    return centered_singleton_corrections(
        nuisance,
        target_ids,
        singleton_supported,
        eta=chosen_eta,
        centering_axis_contract=centering_axis_contract,
    )


def singleton_nuisance_mean(
    fitted: StableSubspaceModel,
    target_features: np.ndarray,
    target_ids: Sequence[str],
    singleton_supported: np.ndarray,
    *,
    centering_axis_contract: AxisContract,
) -> np.ndarray:
    """Return the nuisance mean over the exact deployment centering axis."""
    features = np.asarray(target_features)
    targets = tuple(map(str, target_ids))
    support = np.asarray(singleton_supported, dtype=bool)
    if (
        features.ndim != 2
        or len(features) != len(targets)
        or support.shape != (len(targets),)
    ):
        raise ValueError("singleton feature/target/support shape mismatch")
    if len(set(targets)) != len(targets):
        raise ValueError("target-level singleton feature rows contain duplicates")
    selected = np.flatnonzero(support)
    validate_axis(
        tuple(targets[index] for index in selected),
        centering_axis_contract,
        label="singleton centering axis",
        require_sorted=True,
    )
    nuisance = fitted.nuisance(features[selected])
    if not np.isfinite(nuisance).all():
        raise ValueError("non-finite singleton nuisance prediction")
    return nuisance.mean(axis=0, dtype=np.float64)


def expand_response_corrections(
    fitted: StableSubspaceModel,
    response_corrections: np.ndarray,
) -> np.ndarray:
    """Insert learned response coordinates; excluded coordinates stay zero."""
    values = np.asarray(response_corrections, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != fitted.response_width:
        raise ValueError("response correction width mismatch")
    output = np.zeros((len(values), len(fitted.response_mask)), dtype=np.float64)
    output[:, fitted.response_mask] = values
    return output


def apply_k562_slot_retention(
    correction: np.ndarray,
    perturbmatch_target_support: np.ndarray,
    perturbmatch_gene_support: np.ndarray,
) -> np.ndarray:
    """Retain 50% only on the target-by-gene PerturbMatch intersection."""
    values = np.asarray(correction, dtype=np.float64)
    target_support = np.asarray(perturbmatch_target_support, dtype=bool)
    gene_support = np.asarray(perturbmatch_gene_support, dtype=bool)
    if values.ndim != 2:
        raise ValueError("correction must be a target-by-gene matrix")
    if target_support.shape != (len(values),):
        raise ValueError("correction/PerturbMatch target-axis mismatch")
    if values.shape[1] != len(gene_support):
        raise ValueError("correction/PerturbMatch gene-axis mismatch")
    output = values.copy()
    shared = target_support[:, None] & gene_support[None, :]
    output[shared] *= 0.5
    return output


def held_guide_predictions(
    fitted: StableSubspaceModel,
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    evaluation_targets: Iterable[str],
    singleton_mean_nuisance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...], np.ndarray]:
    """Score translation-invariant recovery of guide-specific residuals.

    The source-side estimand is the raw held-guide residual relative to its
    raw sibling center. The candidate is the portion of that residual removed
    by deployment cleaning. A target-constant nuisance prediction therefore
    receives exactly zero credit.
    """
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    counts = np.asarray(cell_counts)
    identifiable = np.asarray(identifiable, dtype=bool)
    nuisance_mean = np.asarray(singleton_mean_nuisance, dtype=np.float64)
    if nuisance_mean.shape != (fitted.response_width,):
        raise ValueError("singleton nuisance mean width mismatch")
    grouped = EXP083.target_rows(guide_targets)
    mask_columns = np.flatnonzero(fitted.response_mask)
    truth_rows: list[np.ndarray] = []
    baseline_rows: list[np.ndarray] = []
    candidate_rows: list[np.ndarray] = []
    labels: list[str] = []
    held_rows: list[int] = []
    for target in sorted(set(map(str, evaluation_targets))):
        rows = grouped.get(target)
        if rows is None:
            continue
        rows = rows[identifiable[rows]]
        if len(rows) < 2:
            continue
        local_effects = np.asarray(
            effects[np.ix_(rows, mask_columns)], dtype=np.float64
        )
        local_nuisance = fitted.nuisance(features[rows])
        for held_position, held_row in enumerate(rows):
            sibling_positions = np.arange(len(rows)) != held_position
            sibling_rows = rows[sibling_positions]
            raw_sibling_center = EXP083.robust_center(
                local_effects[sibling_positions],
                counts[sibling_rows],
                weight_power=fitted.candidate.weight_power,
                huber_c=fitted.candidate.huber_c,
            )
            cleaned = local_effects - fitted.candidate.eta * (
                local_nuisance - nuisance_mean[None, :]
            )
            clean_sibling_center = EXP083.robust_center(
                cleaned[sibling_positions],
                counts[sibling_rows],
                weight_power=fitted.candidate.weight_power,
                huber_c=fitted.candidate.huber_c,
            )
            raw_residual = local_effects[held_position] - raw_sibling_center
            clean_residual = cleaned[held_position] - clean_sibling_center
            recovered_nuisance = raw_residual - clean_residual
            truth_rows.append(raw_residual)
            baseline_rows.append(np.zeros_like(raw_residual))
            candidate_rows.append(recovered_nuisance)
            labels.append(target)
            held_rows.append(int(held_row))
    if not labels:
        empty = np.empty((0, fitted.response_width), dtype=np.float64)
        return empty, empty.copy(), empty.copy(), (), np.empty(0, dtype=np.int64)
    return (
        np.stack(truth_rows),
        np.stack(baseline_rows),
        np.stack(candidate_rows),
        tuple(labels),
        np.asarray(held_rows, dtype=np.int64),
    )


def held_guide_reconstruction_diagnostic(
    fitted: StableSubspaceModel,
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    evaluation_targets: Iterable[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    """Legacy raw-guide reconstruction diagnostic; never use for selection."""
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    counts = np.asarray(cell_counts)
    identifiable = np.asarray(identifiable, dtype=bool)
    grouped = EXP083.target_rows(guide_targets)
    mask_columns = np.flatnonzero(fitted.response_mask)
    truth_rows: list[np.ndarray] = []
    baseline_rows: list[np.ndarray] = []
    candidate_rows: list[np.ndarray] = []
    labels: list[str] = []
    for target in sorted(set(map(str, evaluation_targets))):
        rows = grouped.get(target)
        if rows is None:
            continue
        rows = rows[identifiable[rows]]
        if len(rows) < 2:
            continue
        local_effects = np.asarray(
            effects[np.ix_(rows, mask_columns)], dtype=np.float64
        )
        local_nuisance = fitted.nuisance(features[rows])
        for held_position, held_row in enumerate(rows):
            sibling_positions = np.arange(len(rows)) != held_position
            sibling_rows = rows[sibling_positions]
            baseline = EXP083.robust_center(
                local_effects[sibling_positions],
                counts[sibling_rows],
                weight_power=fitted.candidate.weight_power,
                huber_c=fitted.candidate.huber_c,
            )
            corrected_siblings = (
                local_effects[sibling_positions]
                - fitted.candidate.eta * local_nuisance[sibling_positions]
            )
            target_center = EXP083.robust_center(
                corrected_siblings,
                counts[sibling_rows],
                weight_power=fitted.candidate.weight_power,
                huber_c=fitted.candidate.huber_c,
            )
            candidate = (
                target_center + fitted.candidate.eta * local_nuisance[held_position]
            )
            truth_rows.append(local_effects[held_position])
            baseline_rows.append(baseline)
            candidate_rows.append(candidate)
            labels.append(target)
    if not labels:
        empty = np.empty((0, fitted.response_width), dtype=np.float64)
        return empty, empty.copy(), empty.copy(), ()
    return (
        np.stack(truth_rows),
        np.stack(baseline_rows),
        np.stack(candidate_rows),
        tuple(labels),
    )


def held_guide_direction_preserving_shrink(
    fitted: StableSubspaceModel,
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    evaluation_targets: Iterable[str],
    singleton_mean_nuisance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...], np.ndarray]:
    """Match cleaned row norms while preserving each raw guide's direction."""
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    counts = np.asarray(cell_counts)
    identifiable = np.asarray(identifiable, dtype=bool)
    nuisance_mean = np.asarray(singleton_mean_nuisance, dtype=np.float64)
    if nuisance_mean.shape != (fitted.response_width,):
        raise ValueError("singleton nuisance mean width mismatch")
    grouped = EXP083.target_rows(guide_targets)
    mask_columns = np.flatnonzero(fitted.response_mask)
    truth_rows: list[np.ndarray] = []
    baseline_rows: list[np.ndarray] = []
    shrink_rows: list[np.ndarray] = []
    labels: list[str] = []
    held_rows: list[int] = []
    for target in sorted(set(map(str, evaluation_targets))):
        rows = grouped.get(target)
        if rows is None:
            continue
        rows = rows[identifiable[rows]]
        if len(rows) < 2:
            continue
        local_effects = np.asarray(
            effects[np.ix_(rows, mask_columns)],
            dtype=np.float64,
        )
        local_nuisance = fitted.nuisance(features[rows])
        cleaned = local_effects - fitted.candidate.eta * (
            local_nuisance - nuisance_mean[None, :]
        )
        shrunken = np.stack(
            [
                EXP083.shrink_control(original, candidate)
                for original, candidate in zip(
                    local_effects,
                    cleaned,
                    strict=True,
                )
            ]
        )
        for held_position, held_row in enumerate(rows):
            sibling_positions = np.arange(len(rows)) != held_position
            sibling_rows = rows[sibling_positions]
            raw_sibling_center = EXP083.robust_center(
                local_effects[sibling_positions],
                counts[sibling_rows],
                weight_power=0.0,
                huber_c=fitted.candidate.huber_c,
            )
            shrink_sibling_center = EXP083.robust_center(
                shrunken[sibling_positions],
                counts[sibling_rows],
                weight_power=0.0,
                huber_c=fitted.candidate.huber_c,
            )
            raw_residual = local_effects[held_position] - raw_sibling_center
            shrink_residual = shrunken[held_position] - shrink_sibling_center
            truth_rows.append(raw_residual)
            baseline_rows.append(np.zeros_like(raw_residual))
            shrink_rows.append(raw_residual - shrink_residual)
            labels.append(target)
            held_rows.append(int(held_row))
    if not labels:
        empty = np.empty((0, fitted.response_width), dtype=np.float64)
        return empty, empty.copy(), empty.copy(), (), np.empty(0, dtype=np.int64)
    return (
        np.stack(truth_rows),
        np.stack(baseline_rows),
        np.stack(shrink_rows),
        tuple(labels),
        np.asarray(held_rows, dtype=np.int64),
    )


def gauge_invariance_tolerance(*arrays: np.ndarray) -> float:
    """Scale-aware roundoff allowance for a bounded finite gauge shift."""
    scale = max(
        1.0,
        *(
            float(np.max(np.abs(np.asarray(value, dtype=np.float64))))
            for value in arrays
            if np.asarray(value).size
        ),
    )
    return 128.0 * np.finfo(np.float64).eps * scale


def target_macro_metrics(
    truth: np.ndarray,
    baseline: np.ndarray,
    candidate: np.ndarray,
    labels: Sequence[str],
) -> dict[str, float]:
    """Average guide-level improvements within target, then across targets."""
    truth = np.asarray(truth, dtype=np.float64)
    baseline = np.asarray(baseline, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    labels = tuple(map(str, labels))
    if truth.shape != baseline.shape or truth.shape != candidate.shape:
        raise ValueError("metric matrix shape mismatch")
    if len(truth) != len(labels):
        raise ValueError("metric label length mismatch")
    if not labels:
        return {
            "targets": 0,
            "guide_holdouts": 0,
            "cosine_delta": float("nan"),
            "nmse_delta": float("nan"),
        }
    baseline_cos = EXP083.row_cosines(truth, baseline)
    candidate_cos = EXP083.row_cosines(truth, candidate)
    denominator = np.maximum(np.sum(truth * truth, axis=1), EPS)
    baseline_nmse = np.sum((baseline - truth) ** 2, axis=1) / denominator
    candidate_nmse = np.sum((candidate - truth) ** 2, axis=1) / denominator
    target_cos: list[float] = []
    target_nmse: list[float] = []
    label_array = np.asarray(labels, dtype=object)
    for target in sorted(set(labels)):
        rows = label_array == target
        target_cos.append(float(np.mean(candidate_cos[rows] - baseline_cos[rows])))
        target_nmse.append(float(np.mean(candidate_nmse[rows] - baseline_nmse[rows])))
    return {
        "targets": len(target_cos),
        "guide_holdouts": len(labels),
        "cosine_delta": float(np.mean(target_cos)),
        "nmse_delta": float(np.mean(target_nmse)),
    }


def per_guide_metric_deltas(
    truth: np.ndarray,
    candidate: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the exact guide-level deltas underlying source-test metrics."""
    truth = np.asarray(truth, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    if truth.ndim != 2 or candidate.shape != truth.shape:
        raise ValueError("source-test truth/candidate shape mismatch")
    if not np.isfinite(truth).all() or not np.isfinite(candidate).all():
        raise ValueError("source-test arrays contain non-finite values")
    baseline = np.zeros_like(truth)
    cosine_delta = EXP083.row_cosines(truth, candidate) - EXP083.row_cosines(
        truth,
        baseline,
    )
    denominator = np.maximum(np.sum(truth * truth, axis=1), EPS)
    baseline_nmse = np.sum((baseline - truth) ** 2, axis=1) / denominator
    candidate_nmse = np.sum((candidate - truth) ** 2, axis=1) / denominator
    return cosine_delta, candidate_nmse - baseline_nmse


def stable_score_extrapolation_diagnostics(
    fitted: StableSubspaceModel,
    reference_features: np.ndarray,
    query_features: np.ndarray,
) -> dict[str, float | int]:
    """Measure singleton score leverage relative to development-guide scores.

    This is diagnostic only. Projection fraction cannot detect extrapolation
    inside the retained subspace, so the frozen report also records robust
    score-space distances and tail exceedance.
    """
    reference = fitted.stable_scores(reference_features)
    query = fitted.stable_scores(query_features)
    if reference.ndim != 2 or query.ndim != 2:
        raise ValueError("stable-score inputs must be matrices")
    if reference.shape[1] != query.shape[1] or reference.shape[1] == 0:
        raise ValueError("stable-score dimensions differ or are empty")
    if len(reference) < 2 or len(query) < 1:
        raise ValueError("insufficient rows for extrapolation diagnostics")
    if not np.isfinite(reference).all() or not np.isfinite(query).all():
        raise ValueError("non-finite stable scores")

    center = np.median(reference, axis=0)
    scale = np.median(np.abs(reference - center), axis=0) * 1.4826
    fallback = np.std(reference, axis=0)
    scale = np.where(scale > EPS, scale, fallback)
    scale = np.where(scale > EPS, scale, 1.0)
    reference_distance = np.sqrt(np.mean(((reference - center) / scale) ** 2, axis=1))
    query_z = (query - center) / scale
    query_distance = np.sqrt(np.mean(query_z**2, axis=1))
    reference_p99 = float(np.quantile(reference_distance, 0.99))
    return {
        "reference_rows": len(reference),
        "query_rows": len(query),
        "stable_rank": int(reference.shape[1]),
        "reference_distance_p99": reference_p99,
        "query_distance_p50": float(np.quantile(query_distance, 0.50)),
        "query_distance_p95": float(np.quantile(query_distance, 0.95)),
        "query_distance_p99": float(np.quantile(query_distance, 0.99)),
        "query_distance_max": float(np.max(query_distance)),
        "query_fraction_above_reference_p99": float(
            np.mean(query_distance > reference_p99)
        ),
        "query_max_abs_robust_z": float(np.max(np.abs(query_z))),
    }


def clone_with_candidate(
    fitted: StableSubspaceModel,
    candidate: Candidate,
) -> StableSubspaceModel:
    """Reuse a fitted rank/ridge map for another eta/weight candidate."""
    if candidate.response_rank != fitted.candidate.response_rank:
        raise ValueError("cannot clone across response ranks")
    if candidate.ridge != fitted.candidate.ridge:
        raise ValueError("cannot clone across ridge strengths")
    return dataclasses.replace(fitted, candidate=candidate)
