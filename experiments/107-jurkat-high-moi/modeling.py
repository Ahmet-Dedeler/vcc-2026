#!/usr/bin/env python3
"""Leakage-safe estimator and frozen selector for Experiment 107.

Production entry points derive target eligibility from guide-level evidence.
They do not accept caller-provided support masks, and channel-held-out support
is recomputed from training channels only. Private ``*_nonproduction`` helpers
keep the numerical algebra independently testable with synthetic inputs.
"""

from __future__ import annotations

import dataclasses
from typing import Sequence

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg


MAX_RESPONSE_RANK = 64
IDENTIFIABLE_RELATIVE_TOLERANCE = 1e-6
TRANSFERABLE_PROJECTION_DIAGONAL = 0.99
FROZEN_RESPONSE_RANK = 64
FROZEN_LAMBDA_GRID = (1.0, 3.0, 10.0, 30.0, 100.0)
FROZEN_CHANNELS = tuple(range(1, 17))
FROZEN_HELD_CHANNEL_FOLDS = (
    (1, 5, 9, 13),
    (2, 6, 10, 14),
    (3, 7, 11, 15),
    (4, 8, 12, 16),
)
MIN_DISTINCT_CALLED_GUIDES = 3
MIN_GUIDES_WITH_THREE_UNTREATED_CELLS = 2
MIN_UNTREATED_CARRIER_CELLS = 20
MIN_TRAINING_CHANNELS = 4
LAMBDA_RELATIVE_GAIN_TOLERANCE = 0.001
SPARSE_SVD_SEED = 1_072_026
PROJECTION_POSTCONDITION_TOLERANCE = 1e-10


class ModelingError(ValueError):
    """Raised when an estimator input violates the frozen algebra contract."""


def _as_float_matrix(
    values: np.ndarray,
    *,
    name: str,
    rows: int | None = None,
) -> np.ndarray:
    matrix = np.asarray(values, dtype=np.float64)
    if matrix.ndim != 2:
        raise ModelingError(f"{name} must be a two-dimensional matrix")
    if rows is not None and matrix.shape[0] != rows:
        raise ModelingError(f"{name} must have {rows} rows")
    if not np.all(np.isfinite(matrix)):
        raise ModelingError(f"{name} contains non-finite values")
    return matrix


def _as_channels(channels: Sequence[object], rows: int) -> np.ndarray:
    labels = np.asarray(channels)
    if labels.ndim != 1 or labels.shape[0] != rows:
        raise ModelingError("channels must be one-dimensional and match matrix rows")
    if rows == 0:
        raise ModelingError("at least one row is required")
    return labels


def _as_binary_target_design(values: np.ndarray, *, rows: int) -> np.ndarray:
    matrix = _as_float_matrix(values, name="binary_target_design", rows=rows)
    if matrix.shape[1] == 0:
        raise ModelingError("binary_target_design must contain at least one target")
    if not np.all((matrix == 0.0) | (matrix == 1.0)):
        raise ModelingError(
            "binary_target_design must contain only exact zeros and ones"
        )
    return matrix


def _as_nonnegative_integer_matrix(
    values: np.ndarray,
    *,
    name: str,
    rows: int | None = None,
    columns: int | None = None,
) -> np.ndarray:
    matrix = np.asarray(values)
    if matrix.ndim != 2:
        raise ModelingError(f"{name} must be a two-dimensional matrix")
    if rows is not None and matrix.shape[0] != rows:
        raise ModelingError(f"{name} must have {rows} rows")
    if columns is not None and matrix.shape[1] != columns:
        raise ModelingError(f"{name} must have {columns} columns")
    try:
        finite = np.asarray(matrix, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ModelingError(f"{name} must contain finite integer counts") from error
    if not np.all(np.isfinite(finite)):
        raise ModelingError(f"{name} contains non-finite values")
    if np.any(finite < 0.0) or np.any(finite != np.floor(finite)):
        raise ModelingError(f"{name} must contain non-negative integer counts")
    if np.any(finite > np.iinfo(np.int64).max):
        raise ModelingError(f"{name} contains counts too large for int64")
    return finite.astype(np.int64)


def _as_guide_target_indices(
    values: np.ndarray,
    *,
    target_count: int,
) -> np.ndarray:
    indices = np.asarray(values)
    if indices.ndim != 1 or indices.size == 0:
        raise ModelingError("guide_target_indices must be a non-empty vector")
    try:
        finite = np.asarray(indices, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ModelingError("guide_target_indices must contain integers") from error
    if not np.all(np.isfinite(finite)) or np.any(finite != np.floor(finite)):
        raise ModelingError("guide_target_indices must contain finite integers")
    integer = finite.astype(np.int64)
    if np.any(integer < -1) or np.any(integer >= target_count):
        raise ModelingError(
            "guide_target_indices must be -1 or a valid target column index"
        )
    return integer


def _matmul(left: np.ndarray, right: np.ndarray, *, name: str) -> np.ndarray:
    """Multiply finite matrices and reject a genuinely non-finite result.

    Some macOS Accelerate builds leave floating-point status flags set after
    otherwise valid BLAS calls, which can make NumPy emit false overflow
    warnings on later ``matmul`` operations. Suppress those flags locally and
    enforce the condition that matters: the returned matrix must be finite.
    """

    with np.errstate(all="ignore"):
        result = np.matmul(left, right)
    if not np.all(np.isfinite(result)):
        raise ModelingError(f"{name} produced non-finite values")
    return result


def _stable_channel_rows(channels: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return channel row indices in deterministic first-appearance order."""

    rows_by_label: dict[object, list[int]] = {}
    for row, raw_label in enumerate(channels.tolist()):
        try:
            rows_by_label.setdefault(raw_label, []).append(row)
        except TypeError as error:
            raise ModelingError("channel labels must be hashable") from error
    return tuple(np.asarray(rows, dtype=np.int64) for rows in rows_by_label.values())


def _integer_channel_tuple(
    channels: Sequence[object],
    *,
    name: str,
) -> tuple[int, ...]:
    result: list[int] = []
    for raw in np.asarray(channels, dtype=object).tolist():
        if isinstance(raw, bool) or not isinstance(raw, (int, np.integer)):
            raise ModelingError(f"{name} must contain integer channel labels")
        result.append(int(raw))
    if len(set(result)) != len(result):
        raise ModelingError(f"{name} must not contain duplicate channel labels")
    return tuple(result)


@dataclasses.dataclass(frozen=True)
class GuideEvidence:
    """Sparse called-guide presence for the untreated estimator population."""

    guide_presence: sparse.csr_matrix
    guide_target_indices: np.ndarray


@dataclasses.dataclass(frozen=True)
class TargetSupportEvidence:
    """Exact frozen support counts computed from specified training channels."""

    training_channels: np.ndarray
    distinct_called_guides: np.ndarray
    guides_with_three_untreated_cells: np.ndarray
    untreated_carrier_cells: np.ndarray
    training_channel_count: np.ndarray
    supported: np.ndarray


@dataclasses.dataclass(frozen=True)
class DesignGeometryAudit:
    """Response-free identifiability geometry for one training population."""

    training_channels: np.ndarray
    evidence_support: np.ndarray
    zero_variance: np.ndarray
    nuisance_aliased: np.ndarray
    fitted: np.ndarray
    target_scales: np.ndarray
    nuisance_target_slopes: np.ndarray
    nullspace_basis: np.ndarray
    projection_diagonal: np.ndarray
    singular_values: np.ndarray
    retained_singular_directions: np.ndarray
    effective_rank: int
    condition_number: float


@dataclasses.dataclass(frozen=True)
class ProductionDesignAudit:
    """Frozen response-free design audits for all folds and the final fit."""

    held_channel_folds: tuple[tuple[int, ...], ...]
    fold_support: tuple[TargetSupportEvidence, ...]
    fold_geometry: tuple[DesignGeometryAudit, ...]
    final_support: TargetSupportEvidence
    final_geometry: DesignGeometryAudit


@dataclasses.dataclass(frozen=True)
class _ValidatedGuideEvidence:
    guide_presence: sparse.csr_matrix
    guide_target_indices: np.ndarray


def _as_csr_binary_design(
    values: sparse.spmatrix,
    *,
    name: str,
    rows: int | None = None,
    columns: int | None = None,
) -> sparse.csr_matrix:
    if not sparse.isspmatrix_csr(values):
        raise ModelingError(f"{name} must be a scipy CSR matrix")
    matrix = values
    if matrix.ndim != 2:
        raise ModelingError(f"{name} must be two-dimensional")
    if rows is not None and matrix.shape[0] != rows:
        raise ModelingError(f"{name} must have {rows} rows")
    if columns is not None and matrix.shape[1] != columns:
        raise ModelingError(f"{name} must have {columns} columns")
    if not matrix.has_canonical_format:
        raise ModelingError(f"{name} must have canonical sorted CSR indices")
    if matrix.data.size and (
        not np.all(np.isfinite(matrix.data)) or not np.all(matrix.data == 1)
    ):
        raise ModelingError(
            f"{name} must contain only implicit zeros and explicit ones"
        )
    return matrix.astype(np.float64, copy=False)


def _validate_guide_evidence(
    evidence: GuideEvidence | None,
    *,
    rows: int,
    target_count: int,
) -> _ValidatedGuideEvidence:
    if not isinstance(evidence, GuideEvidence):
        raise ModelingError("guide_evidence is required and must be GuideEvidence")
    guide_target_indices = _as_guide_target_indices(
        evidence.guide_target_indices,
        target_count=target_count,
    )
    guide_presence = _as_csr_binary_design(
        evidence.guide_presence,
        name="guide_presence",
        rows=rows,
        columns=guide_target_indices.size,
    )
    return _ValidatedGuideEvidence(
        guide_presence=guide_presence,
        guide_target_indices=guide_target_indices,
    )


def _validate_guide_target_collapse(
    evidence: _ValidatedGuideEvidence,
    binary_target_design: sparse.csr_matrix,
) -> None:
    eligible = evidence.guide_target_indices >= 0
    guide_rows = np.flatnonzero(eligible)
    mapping = sparse.csr_matrix(
        (
            np.ones(guide_rows.size, dtype=np.float64),
            (
                guide_rows,
                evidence.guide_target_indices[eligible],
            ),
        ),
        shape=(
            evidence.guide_target_indices.size,
            binary_target_design.shape[1],
        ),
    )
    collapsed = (evidence.guide_presence @ mapping).tocsr()
    collapsed.sum_duplicates()
    collapsed.eliminate_zeros()
    collapsed.data.fill(1.0)
    if collapsed.shape != binary_target_design.shape:
        raise ModelingError("guide collapse and target design shapes differ")
    if (collapsed != binary_target_design).nnz:
        raise ModelingError(
            "binary_target_design does not equal eligible guide collapse"
        )


def _derive_target_support_validated(
    evidence: _ValidatedGuideEvidence,
    *,
    binary_target_design: sparse.csr_matrix,
    channels: np.ndarray,
    training_channels: Sequence[object],
    target_count: int,
) -> TargetSupportEvidence:
    requested_channels = _integer_channel_tuple(
        training_channels,
        name="training_channels",
    )
    if not requested_channels:
        raise ModelingError("training_channels must not be empty")
    data_channels = set(_stable_integer_data_channels(channels))
    if not set(requested_channels).issubset(data_channels):
        raise ModelingError(
            "training_channels contain a channel absent from model inputs"
        )
    training_rows = np.isin(channels, requested_channels)
    if not np.any(training_rows):
        raise ModelingError("training_channels select no cells")
    guide_totals = np.asarray(
        evidence.guide_presence[training_rows].sum(axis=0)
    ).ravel()
    guide_targets = evidence.guide_target_indices
    eligible = guide_targets >= 0
    distinct_called = np.bincount(
        guide_targets[eligible & (guide_totals > 0)],
        minlength=target_count,
    )
    guides_with_three = np.bincount(
        guide_targets[eligible & (guide_totals >= 3)],
        minlength=target_count,
    )
    untreated_carriers = (
        np.asarray(binary_target_design[training_rows].sum(axis=0))
        .ravel()
        .astype(np.int64, copy=False)
    )
    training_channel_count = np.zeros(target_count, dtype=np.int64)
    for channel in requested_channels:
        channel_carriers = np.asarray(
            binary_target_design[channels == channel].sum(axis=0)
        ).ravel()
        training_channel_count += channel_carriers > 0
    supported = (
        (distinct_called >= MIN_DISTINCT_CALLED_GUIDES)
        & (guides_with_three >= MIN_GUIDES_WITH_THREE_UNTREATED_CELLS)
        & (untreated_carriers >= MIN_UNTREATED_CARRIER_CELLS)
        & (training_channel_count >= MIN_TRAINING_CHANNELS)
    )
    return TargetSupportEvidence(
        training_channels=np.asarray(requested_channels, dtype=np.int64),
        distinct_called_guides=distinct_called.astype(np.int64, copy=False),
        guides_with_three_untreated_cells=guides_with_three.astype(
            np.int64,
            copy=False,
        ),
        untreated_carrier_cells=untreated_carriers,
        training_channel_count=training_channel_count.astype(
            np.int64,
            copy=False,
        ),
        supported=supported,
    )


def derive_target_support(
    guide_evidence: GuideEvidence,
    binary_target_design: sparse.csr_matrix,
    channels: Sequence[object],
    *,
    training_channels: Sequence[object],
) -> TargetSupportEvidence:
    """Derive the four frozen support gates from guide evidence only."""

    labels = _as_channels(channels, binary_target_design.shape[0])
    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
        rows=labels.size,
    )
    validated = _validate_guide_evidence(
        guide_evidence,
        rows=targets.shape[0],
        target_count=targets.shape[1],
    )
    _validate_guide_target_collapse(validated, targets)
    return _derive_target_support_validated(
        validated,
        binary_target_design=targets,
        channels=labels,
        training_channels=training_channels,
        target_count=targets.shape[1],
    )


def within_channel_center(
    values: np.ndarray,
    channels: Sequence[object],
) -> np.ndarray:
    """Subtract each column's own mean within every complete channel.

    The input is converted to float64 and never modified. Means are accumulated
    in original row order, while channels are visited in first-appearance
    order, making repeated calls deterministic.
    """

    matrix = _as_float_matrix(values, name="values")
    labels = _as_channels(channels, matrix.shape[0])
    centered = matrix.copy()
    for rows in _stable_channel_rows(labels):
        centered[rows] -= matrix[rows].mean(axis=0, keepdims=True)
    return centered


@dataclasses.dataclass(frozen=True)
class ProjectionDiagnostics:
    """SVD audit record for one nuisance-residualized target design."""

    projection: np.ndarray
    singular_values: np.ndarray
    retained: np.ndarray
    diagonal: np.ndarray
    effective_rank: int
    condition_number: float
    relative_tolerance: float


def identifiable_row_space(
    residualized_target_design: np.ndarray,
    *,
    relative_tolerance: float = IDENTIFIABLE_RELATIVE_TOLERANCE,
) -> ProjectionDiagnostics:
    """Compute the identifiable target row-space projector in float64.

    A right-singular direction is retained exactly when
    ``sigma / sigma_max > relative_tolerance``.
    """

    design = _as_float_matrix(
        residualized_target_design,
        name="residualized_target_design",
    )
    if design.shape[1] == 0:
        raise ModelingError("residualized_target_design must contain a target")
    if not np.isfinite(relative_tolerance) or relative_tolerance < 0:
        raise ModelingError("relative_tolerance must be finite and non-negative")

    _, singular_values, right = np.linalg.svd(design, full_matrices=False)
    retained = np.zeros(singular_values.shape, dtype=bool)
    if singular_values.size and singular_values[0] > 0.0:
        retained = singular_values / singular_values[0] > relative_tolerance

    identifiable = right[retained]
    projection = _matmul(
        identifiable.T,
        identifiable,
        name="identifiable projection",
    )
    # Symmetrization removes harmless last-bit BLAS asymmetry and gives a
    # deterministic diagonal for the transferability gate.
    projection = 0.5 * (projection + projection.T)
    effective_rank = int(np.count_nonzero(retained))
    condition_number = (
        float(singular_values[0] / singular_values[retained][-1])
        if effective_rank
        else float("inf")
    )
    return ProjectionDiagnostics(
        projection=projection,
        singular_values=singular_values,
        retained=retained,
        diagonal=np.diag(projection).copy(),
        effective_rank=effective_rank,
        condition_number=condition_number,
        relative_tolerance=float(relative_tolerance),
    )


def project_identifiable_coefficients(
    residualized_target_design: np.ndarray,
    coefficients_binary_units: np.ndarray,
    *,
    relative_tolerance: float = IDENTIFIABLE_RELATIVE_TOLERANCE,
) -> tuple[np.ndarray, ProjectionDiagnostics]:
    """Remove every coefficient component in the target-design nullspace."""

    coefficients = _as_float_matrix(
        coefficients_binary_units,
        name="coefficients_binary_units",
    )
    if coefficients.shape[0] != np.asarray(residualized_target_design).shape[1]:
        raise ModelingError(
            "coefficient rows must match residualized target-design columns"
        )
    diagnostics = identifiable_row_space(
        residualized_target_design,
        relative_tolerance=relative_tolerance,
    )
    return (
        _matmul(
            diagnostics.projection,
            coefficients,
            name="coefficient projection",
        ),
        diagnostics,
    )


@dataclasses.dataclass(frozen=True)
class SupportDiagnostics:
    """Per-target support and identifiability decisions from training only."""

    evidence_support: np.ndarray
    zero_variance: np.ndarray
    nuisance_aliased: np.ndarray
    fitted: np.ndarray
    projection_diagonal: np.ndarray
    transferable: np.ndarray
    target_scales: np.ndarray
    singular_values: np.ndarray
    retained_singular_directions: np.ndarray
    effective_rank: int
    condition_number: float
    relative_singular_tolerance: float
    projection_diagonal_threshold: float
    nullspace_annihilation_error: float = 0.0
    projection_idempotence_error: float = 0.0
    gauge_shift_invariance_error: float = 0.0

    @property
    def unsupported(self) -> np.ndarray:
        return ~self.transferable

    @property
    def requested_support(self) -> np.ndarray:
        """Backward-readable alias; support is evidence-derived in public fits."""

        return self.evidence_support


@dataclasses.dataclass(frozen=True)
class HeldChannelResult:
    """Held-channel transforms and diagnostics; none are used during fitting."""

    centered_response: np.ndarray
    centered_target_design: np.ndarray | None
    centered_nuisance: np.ndarray
    residual_response: np.ndarray
    residual_target_design: np.ndarray | None
    predicted_residual_response: np.ndarray
    unmodeled_target_burden: np.ndarray


@dataclasses.dataclass(frozen=True)
class RidgeSolveDiagnostics:
    """Deterministic batched-CG convergence record for one ridge penalty."""

    ridge_lambda: float
    iterations: int
    operator_applications: int
    warm_started: bool
    converged: np.ndarray
    rhs_norms: np.ndarray
    residual_norms: np.ndarray
    relative_residual_norms: np.ndarray
    absolute_tolerance: float
    relative_tolerance: float


@dataclasses.dataclass(frozen=True)
class FittedEstimator:
    """One estimator of identifiable high-MOI conditional target contrasts.

    Exported rows are minimum-norm components in the identifiable row space.
    They are not claimed to be uniquely identified singleton marginal effects.
    """

    ridge_lambda: float
    requested_response_rank: int
    effective_response_rank: int
    nuisance_response_slopes: np.ndarray
    nuisance_target_slopes: np.ndarray
    response_basis: np.ndarray
    coefficients_binary_units: np.ndarray
    identifiable_projection: np.ndarray | None
    nullspace_basis: np.ndarray
    projected_coefficients: np.ndarray
    exported_effects: np.ndarray
    diagnostics: SupportDiagnostics
    ridge_solve_diagnostics: RidgeSolveDiagnostics | None = None

    def project_coefficients(self, coefficients: np.ndarray) -> np.ndarray:
        """Project another coefficient matrix through this fit's row space."""

        matrix = _as_float_matrix(coefficients, name="coefficients")
        if matrix.shape[0] != self.exported_effects.shape[0]:
            raise ModelingError("coefficient rows do not match fitted targets")
        if self.identifiable_projection is not None:
            return _matmul(
                self.identifiable_projection,
                matrix,
                name="fitted coefficient projection",
            )
        if not self.nullspace_basis.size:
            return matrix.copy()
        return matrix - _matmul(
            self.nullspace_basis.T,
            _matmul(
                self.nullspace_basis,
                matrix,
                name="nullspace coefficient coordinates",
            ),
            name="implicit coefficient projection",
        )

    def predict_held(
        self,
        binary_target_design: np.ndarray | sparse.csr_matrix,
        nuisance: np.ndarray,
        channels: Sequence[object],
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        """Predict held residual response using training-fitted quantities only.

        Returns ``(prediction, residualized_target_design,
        unmodeled_target_burden)``.
        """

        if sparse.issparse(binary_target_design):
            targets_sparse = _as_csr_binary_design(
                binary_target_design,
                name="binary_target_design",
            )
            nuisance_matrix = _as_float_matrix(
                nuisance,
                name="nuisance",
                rows=targets_sparse.shape[0],
            )
            labels = _as_channels(channels, targets_sparse.shape[0])
            if targets_sparse.shape[1] != self.nuisance_target_slopes.shape[1]:
                raise ModelingError("held target columns do not match fitted targets")
            centered_nuisance = within_channel_center(nuisance_matrix, labels)
            raw_prediction = np.asarray(
                targets_sparse @ self.exported_effects,
                dtype=np.float64,
            )
            prediction = within_channel_center(raw_prediction, labels)
            if centered_nuisance.shape[1]:
                prediction -= _matmul(
                    centered_nuisance,
                    _matmul(
                        self.nuisance_target_slopes,
                        self.exported_effects,
                        name="held sparse nuisance coefficients",
                    ),
                    name="held sparse nuisance prediction",
                )
            unmodeled = np.asarray(
                targets_sparse[:, self.diagnostics.unsupported].sum(axis=1)
            ).ravel()
            return prediction, None, unmodeled

        targets = _as_binary_target_design(
            binary_target_design,
            rows=binary_target_design.shape[0],
        )
        nuisance_matrix = _as_float_matrix(
            nuisance,
            name="nuisance",
            rows=targets.shape[0],
        )
        labels = _as_channels(channels, targets.shape[0])
        if targets.shape[1] != self.nuisance_target_slopes.shape[1]:
            raise ModelingError("held target columns do not match fitted targets")
        if nuisance_matrix.shape[1] != self.nuisance_target_slopes.shape[0]:
            raise ModelingError("held nuisance columns do not match fitted nuisance")

        centered_targets = within_channel_center(targets, labels)
        centered_nuisance = within_channel_center(nuisance_matrix, labels)
        residual_targets = centered_targets - _matmul(
            centered_nuisance,
            self.nuisance_target_slopes,
            name="held target nuisance prediction",
        )
        prediction = _matmul(
            residual_targets,
            self.exported_effects,
            name="held residual response prediction",
        )
        unmodeled = np.sum(
            np.abs(residual_targets[:, self.diagnostics.unsupported]),
            axis=1,
        )
        return prediction, residual_targets, unmodeled

    def evaluate_held(
        self,
        response: np.ndarray,
        binary_target_design: np.ndarray | sparse.csr_matrix,
        nuisance: np.ndarray,
        channels: Sequence[object],
    ) -> HeldChannelResult:
        """Transform held response without refitting any estimator quantity."""

        response_matrix = _as_float_matrix(response, name="response")
        if sparse.issparse(binary_target_design):
            targets_sparse = _as_csr_binary_design(
                binary_target_design,
                name="binary_target_design",
                rows=response_matrix.shape[0],
            )
            nuisance_matrix = _as_float_matrix(
                nuisance,
                name="nuisance",
                rows=response_matrix.shape[0],
            )
            labels = _as_channels(channels, response_matrix.shape[0])
            return _evaluate_sparse_held(
                self,
                response_matrix,
                targets_sparse,
                nuisance_matrix,
                labels,
            )
        targets = _as_binary_target_design(
            binary_target_design,
            rows=response_matrix.shape[0],
        )
        nuisance_matrix = _as_float_matrix(
            nuisance,
            name="nuisance",
            rows=response_matrix.shape[0],
        )
        labels = _as_channels(channels, response_matrix.shape[0])
        if response_matrix.shape[1] != self.nuisance_response_slopes.shape[1]:
            raise ModelingError("held response columns do not match fitted response")

        centered_response = within_channel_center(response_matrix, labels)
        centered_targets = within_channel_center(targets, labels)
        centered_nuisance = within_channel_center(nuisance_matrix, labels)
        residual_response = centered_response - _matmul(
            centered_nuisance,
            self.nuisance_response_slopes,
            name="held response nuisance prediction",
        )
        residual_targets = centered_targets - _matmul(
            centered_nuisance,
            self.nuisance_target_slopes,
            name="held target nuisance prediction",
        )
        prediction = _matmul(
            residual_targets,
            self.exported_effects,
            name="held residual response prediction",
        )
        unmodeled = np.sum(
            np.abs(residual_targets[:, self.diagnostics.unsupported]),
            axis=1,
        )
        return HeldChannelResult(
            centered_response=centered_response,
            centered_target_design=centered_targets,
            centered_nuisance=centered_nuisance,
            residual_response=residual_response,
            residual_target_design=residual_targets,
            predicted_residual_response=prediction,
            unmodeled_target_burden=unmodeled,
        )

    def apply_exact_fallback(
        self,
        incumbent: np.ndarray,
        *,
        response_supported: np.ndarray | None = None,
    ) -> np.ndarray:
        """Insert transferable source effects and preserve all else bit-exactly."""

        support = self.diagnostics.transferable[:, None]
        if response_supported is not None:
            response_mask = np.asarray(response_supported, dtype=bool)
            if response_mask.shape != (self.exported_effects.shape[1],):
                raise ModelingError("response_supported has the wrong shape")
            support = support & response_mask[None, :]
        return exact_unsupported_fallback(
            incumbent,
            self.exported_effects,
            support,
        )


def _relative_max_abs(error: np.ndarray, reference: np.ndarray) -> float:
    numerator = float(np.max(np.abs(error), initial=0.0))
    denominator = max(
        1.0,
        float(np.max(np.abs(reference), initial=0.0)),
    )
    return numerator / denominator


def validate_production_projection(
    coefficients: np.ndarray,
    projected: np.ndarray,
    nullspace: np.ndarray,
    *,
    tolerance: float = PROJECTION_POSTCONDITION_TOLERANCE,
) -> tuple[float, float, float]:
    """Enforce nullspace, idempotence, and gauge invariance postconditions."""

    raw = _as_float_matrix(coefficients, name="projection coefficients")
    observed = _as_float_matrix(
        projected,
        name="projected coefficients",
        rows=raw.shape[0],
    )
    if observed.shape != raw.shape:
        raise ModelingError("projected coefficients have the wrong shape")
    basis = _as_float_matrix(nullspace, name="projection nullspace")
    if basis.shape[1] != raw.shape[0]:
        raise ModelingError("projection nullspace has the wrong target axis")
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ModelingError("projection tolerance must be finite and non-negative")

    if basis.size:
        null_coordinates = _matmul(
            basis,
            observed,
            name="projected nullspace coordinates",
        )
        annihilation_error = _relative_max_abs(null_coordinates, observed)
        reprojection = observed - _matmul(
            basis.T,
            null_coordinates,
            name="projected coefficient reprojection",
        )
        idempotence_error = _relative_max_abs(reprojection - observed, observed)

        gauge_coordinates = np.where(
            np.arange(basis.shape[0], dtype=np.int64) % 2 == 0,
            1.0,
            -1.0,
        )[:, None]
        gauge_shift = _matmul(
            basis.T,
            gauge_coordinates,
            name="deterministic nullspace gauge shift",
        )
        shifted = raw[:, :1] + gauge_shift
        shifted_projection = shifted - _matmul(
            basis.T,
            _matmul(
                basis,
                shifted,
                name="shifted nullspace coordinates",
            ),
            name="shifted identifiable projection",
        )
        gauge_error = _relative_max_abs(
            shifted_projection - observed[:, :1],
            observed[:, :1],
        )
    else:
        annihilation_error = 0.0
        idempotence_error = _relative_max_abs(observed - raw, observed)
        gauge_error = idempotence_error

    failures = {
        "nullspace annihilation": annihilation_error,
        "projection idempotence": idempotence_error,
        "gauge-shift invariance": gauge_error,
    }
    for name, error in failures.items():
        if error > tolerance:
            raise ModelingError(
                f"production {name} error {error:.3e} exceeds {tolerance:.1e}"
            )
    return annihilation_error, idempotence_error, gauge_error


def _least_squares_pseudoinverse(design: np.ndarray) -> np.ndarray:
    """Build the float64 least-squares operator used for every nuisance fit."""

    matrix = _as_float_matrix(design, name="least_squares_design")
    if matrix.shape[1] == 0:
        return np.zeros((0, matrix.shape[0]), dtype=np.float64)
    left, singular_values, right = np.linalg.svd(matrix, full_matrices=False)
    if not singular_values.size:
        return np.zeros((matrix.shape[1], matrix.shape[0]), dtype=np.float64)
    cutoff = np.finfo(np.float64).eps * max(matrix.shape) * float(singular_values[0])
    retained = singular_values > cutoff
    inverse = np.zeros_like(singular_values)
    inverse[retained] = 1.0 / singular_values[retained]
    return _matmul(
        right.T * inverse[None, :],
        left.T,
        name="least-squares pseudoinverse",
    )


def _nuisance_slopes(nuisance: np.ndarray, outcome: np.ndarray) -> np.ndarray:
    if nuisance.shape[1] == 0:
        return np.zeros((0, outcome.shape[1]), dtype=np.float64)
    return _matmul(
        _least_squares_pseudoinverse(nuisance),
        outcome,
        name="nuisance least-squares slopes",
    )


def _canonical_response_basis(
    residual_response: np.ndarray,
    requested_rank: int,
) -> np.ndarray:
    _, singular_values, right = np.linalg.svd(
        residual_response,
        full_matrices=False,
    )
    if not singular_values.size or singular_values[0] == 0.0:
        return np.zeros((0, residual_response.shape[1]), dtype=np.float64)
    numerical_tolerance = (
        np.finfo(np.float64).eps * max(residual_response.shape) * singular_values[0]
    )
    effective_rank = min(
        requested_rank,
        int(np.count_nonzero(singular_values > numerical_tolerance)),
    )
    basis = right[:effective_rank].copy()
    # SVD signs are arbitrary. Canonical signs keep serialized diagnostics
    # stable while leaving the reconstructed coefficient matrix unchanged.
    for component in basis:
        pivot = int(np.argmax(np.abs(component)))
        if component[pivot] < 0.0:
            component *= -1.0
    return basis


def _canonical_response_basis_from_covariance(
    residual_response: np.ndarray,
    requested_rank: int,
) -> np.ndarray:
    """Fit the response basis without allocating the tall left SVD factor."""

    covariance = _matmul(
        residual_response.T,
        residual_response,
        name="response covariance",
    )
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    singular_values = np.sqrt(eigenvalues)
    if not singular_values.size or singular_values[0] == 0.0:
        return np.zeros((0, residual_response.shape[1]), dtype=np.float64)
    numerical_tolerance = (
        np.finfo(np.float64).eps * max(residual_response.shape) * singular_values[0]
    )
    effective_rank = min(
        requested_rank,
        int(np.count_nonzero(singular_values > numerical_tolerance)),
    )
    basis = eigenvectors[:, order[:effective_rank]].T.copy()
    for component in basis:
        pivot = int(np.argmax(np.abs(component)))
        if component[pivot] < 0.0:
            component *= -1.0
    return basis


def _validate_fit_controls(
    ridge_lambda: float,
    rank: int,
    relative_singular_tolerance: float,
    projection_diagonal_threshold: float,
) -> tuple[float, int]:
    penalty = float(ridge_lambda)
    if not np.isfinite(penalty) or penalty < 0.0:
        raise ModelingError("ridge_lambda must be finite and non-negative")
    if isinstance(rank, bool) or int(rank) != rank or not 1 <= int(rank) <= 64:
        raise ModelingError("rank must be an integer from 1 through 64")
    if (
        not np.isfinite(relative_singular_tolerance)
        or relative_singular_tolerance < 0.0
    ):
        raise ModelingError(
            "relative_singular_tolerance must be finite and non-negative"
        )
    if (
        not np.isfinite(projection_diagonal_threshold)
        or not 0.0 <= projection_diagonal_threshold <= 1.0
    ):
        raise ModelingError(
            "projection_diagonal_threshold must be between zero and one"
        )
    return penalty, int(rank)


def _validate_frozen_fit_controls(
    ridge_lambda: float,
    rank: int,
) -> tuple[float, int]:
    penalty, requested_rank = _validate_fit_controls(
        ridge_lambda,
        rank,
        IDENTIFIABLE_RELATIVE_TOLERANCE,
        TRANSFERABLE_PROJECTION_DIAGONAL,
    )
    if penalty not in FROZEN_LAMBDA_GRID:
        raise ModelingError(
            f"production ridge_lambda must be one of {FROZEN_LAMBDA_GRID}"
        )
    if requested_rank != FROZEN_RESPONSE_RANK:
        raise ModelingError(f"production response rank must be {FROZEN_RESPONSE_RANK}")
    return penalty, requested_rank


def _stable_integer_data_channels(channels: np.ndarray) -> tuple[int, ...]:
    result: list[int] = []
    seen: set[int] = set()
    for raw in channels.tolist():
        if isinstance(raw, bool) or not isinstance(raw, (int, np.integer)):
            raise ModelingError("production channels must be integer labels")
        channel = int(raw)
        if channel not in seen:
            result.append(channel)
            seen.add(channel)
    return tuple(result)


def _validate_evidence_matches_design(
    evidence: _ValidatedGuideEvidence,
    binary_target_design: sparse.csr_matrix,
) -> None:
    _validate_guide_target_collapse(evidence, binary_target_design)


def _validate_frozen_selector_controls(
    lambda_grid: Sequence[float],
    rank: int,
    folds: Sequence[Sequence[object]],
) -> tuple[tuple[float, ...], int, tuple[tuple[int, ...], ...]]:
    try:
        penalties = tuple(float(value) for value in lambda_grid)
    except (TypeError, ValueError) as error:
        raise ModelingError("lambda_grid must contain numeric values") from error
    if penalties != FROZEN_LAMBDA_GRID:
        raise ModelingError(f"production lambda_grid must equal {FROZEN_LAMBDA_GRID}")
    if isinstance(rank, bool) or int(rank) != rank:
        raise ModelingError("production response rank must be an integer")
    requested_rank = int(rank)
    if requested_rank != FROZEN_RESPONSE_RANK:
        raise ModelingError(f"production response rank must be {FROZEN_RESPONSE_RANK}")
    try:
        normalized_folds = tuple(
            _integer_channel_tuple(fold, name="held fold") for fold in folds
        )
    except TypeError as error:
        raise ModelingError("folds must be a sequence of channel sequences") from error
    if normalized_folds != FROZEN_HELD_CHANNEL_FOLDS:
        raise ModelingError(f"production folds must equal {FROZEN_HELD_CHANNEL_FOLDS}")
    return penalties, requested_rank, normalized_folds


def _validate_frozen_production_channels(
    labels: np.ndarray,
) -> None:
    data_channels = _stable_integer_data_channels(labels)
    if set(data_channels) != set(FROZEN_CHANNELS):
        raise ModelingError(
            f"production data must contain exactly channels {FROZEN_CHANNELS}"
        )


def _fit_centered(
    centered_response: np.ndarray,
    centered_targets: np.ndarray,
    centered_nuisance: np.ndarray,
    *,
    ridge_lambda: float,
    rank: int,
    support: np.ndarray,
    relative_singular_tolerance: float,
    projection_diagonal_threshold: float,
) -> FittedEstimator:
    response_slopes = _nuisance_slopes(centered_nuisance, centered_response)
    target_slopes = _nuisance_slopes(centered_nuisance, centered_targets)
    residual_response = centered_response - _matmul(
        centered_nuisance,
        response_slopes,
        name="training response nuisance prediction",
    )
    residual_targets = centered_targets - _matmul(
        centered_nuisance,
        target_slopes,
        name="training target nuisance prediction",
    )

    centered_scales = np.std(centered_targets, axis=0, ddof=0)
    residual_scales = np.std(residual_targets, axis=0, ddof=0)
    scale_reference = max(
        1.0,
        float(np.max(centered_scales, initial=0.0)),
        float(np.max(residual_scales, initial=0.0)),
    )
    scale_tolerance = (
        np.finfo(np.float64).eps * max(centered_targets.shape) * scale_reference
    )
    zero_variance = centered_scales <= scale_tolerance
    nuisance_aliased = (~zero_variance) & (residual_scales <= scale_tolerance)
    fitted = support & ~zero_variance & ~nuisance_aliased

    target_count = centered_targets.shape[1]
    response_count = centered_response.shape[1]
    target_scales = np.full(target_count, np.nan, dtype=np.float64)
    target_scales[fitted] = residual_scales[fitted]
    coefficients = np.zeros((target_count, response_count), dtype=np.float64)
    projected = np.zeros_like(coefficients)
    full_projection = np.zeros((target_count, target_count), dtype=np.float64)
    projection_diagonal = np.zeros(target_count, dtype=np.float64)
    response_basis = np.zeros((0, response_count), dtype=np.float64)
    singular_values = np.zeros(0, dtype=np.float64)
    retained = np.zeros(0, dtype=bool)
    effective_rank = 0
    condition_number = float("inf")

    if np.any(fitted):
        active_targets = residual_targets[:, fitted]
        scaled_targets = active_targets / target_scales[fitted]
        response_basis = _canonical_response_basis(residual_response, rank)
        if response_basis.shape[0]:
            response_scores = _matmul(
                residual_response,
                response_basis.T,
                name="training response scores",
            )
            if ridge_lambda == 0.0:
                scaled_coefficients = np.linalg.lstsq(
                    scaled_targets,
                    response_scores,
                    rcond=None,
                )[0]
            else:
                gram = _matmul(
                    scaled_targets.T,
                    scaled_targets,
                    name="scaled target Gram matrix",
                )
                gram.flat[:: gram.shape[0] + 1] += ridge_lambda
                scaled_coefficients = np.linalg.solve(
                    gram,
                    _matmul(
                        scaled_targets.T,
                        response_scores,
                        name="target-response cross product",
                    ),
                )
            active_coefficients = (
                _matmul(
                    scaled_coefficients,
                    response_basis,
                    name="binary-unit coefficient reconstruction",
                )
            ) / target_scales[fitted, None]
        else:
            active_coefficients = np.zeros(
                (int(np.count_nonzero(fitted)), response_count),
                dtype=np.float64,
            )

        active_projected, projection = project_identifiable_coefficients(
            active_targets,
            active_coefficients,
            relative_tolerance=relative_singular_tolerance,
        )
        active_indices = np.flatnonzero(fitted)
        coefficients[fitted] = active_coefficients
        projected[fitted] = active_projected
        full_projection[np.ix_(active_indices, active_indices)] = projection.projection
        projection_diagonal[fitted] = projection.diagonal
        singular_values = projection.singular_values
        retained = projection.retained
        effective_rank = projection.effective_rank
        condition_number = projection.condition_number

    transferable = fitted & (projection_diagonal >= projection_diagonal_threshold)
    exported = np.zeros_like(projected)
    exported[transferable] = projected[transferable]
    diagnostics = SupportDiagnostics(
        evidence_support=support.copy(),
        zero_variance=zero_variance,
        nuisance_aliased=nuisance_aliased,
        fitted=fitted,
        projection_diagonal=projection_diagonal,
        transferable=transferable,
        target_scales=target_scales,
        singular_values=singular_values,
        retained_singular_directions=retained,
        effective_rank=effective_rank,
        condition_number=condition_number,
        relative_singular_tolerance=float(relative_singular_tolerance),
        projection_diagonal_threshold=float(projection_diagonal_threshold),
    )
    return FittedEstimator(
        ridge_lambda=ridge_lambda,
        requested_response_rank=rank,
        effective_response_rank=response_basis.shape[0],
        nuisance_response_slopes=response_slopes,
        nuisance_target_slopes=target_slopes,
        response_basis=response_basis,
        coefficients_binary_units=coefficients,
        identifiable_projection=full_projection,
        nullspace_basis=np.zeros((0, target_count), dtype=np.float64),
        projected_coefficients=projected,
        exported_effects=exported,
        diagnostics=diagnostics,
    )


def _fit_estimator_nonproduction(
    response: np.ndarray,
    binary_target_design: np.ndarray,
    nuisance: np.ndarray,
    channels: Sequence[object],
    *,
    ridge_lambda: float,
    rank: int = MAX_RESPONSE_RANK,
    support: np.ndarray,
    relative_singular_tolerance: float = IDENTIFIABLE_RELATIVE_TOLERANCE,
    projection_diagonal_threshold: float = TRANSFERABLE_PROJECTION_DIAGONAL,
) -> FittedEstimator:
    """Exercise low-level algebra with explicit controls; not for production."""

    penalty, requested_rank = _validate_fit_controls(
        ridge_lambda,
        rank,
        relative_singular_tolerance,
        projection_diagonal_threshold,
    )
    response_matrix = _as_float_matrix(response, name="response")
    targets = _as_binary_target_design(
        binary_target_design,
        rows=response_matrix.shape[0],
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=response_matrix.shape[0],
    )
    labels = _as_channels(channels, response_matrix.shape[0])
    support_mask = np.asarray(support)
    if support_mask.dtype != np.bool_ or support_mask.shape != (targets.shape[1],):
        raise ModelingError("support must be a Boolean value per target")

    return _fit_centered(
        within_channel_center(response_matrix, labels),
        within_channel_center(targets, labels),
        within_channel_center(nuisance_matrix, labels),
        ridge_lambda=penalty,
        rank=requested_rank,
        support=support_mask,
        relative_singular_tolerance=relative_singular_tolerance,
        projection_diagonal_threshold=projection_diagonal_threshold,
    )


@dataclasses.dataclass(frozen=True)
class _SparseResidualizedDesign:
    """Implicit ``(I-P_nuisance) C_channel X`` target design."""

    target_design: sparse.csr_matrix
    channels: np.ndarray
    centered_nuisance: np.ndarray
    nuisance_pseudoinverse: np.ndarray
    target_scales: np.ndarray

    def _residualize_dense(self, values: np.ndarray) -> np.ndarray:
        centered = within_channel_center(values, self.channels)
        if not self.centered_nuisance.shape[1]:
            return centered
        slopes = _matmul(
            self.nuisance_pseudoinverse,
            centered,
            name="implicit nuisance slopes",
        )
        return centered - _matmul(
            self.centered_nuisance,
            slopes,
            name="implicit nuisance prediction",
        )

    def matmat(self, coefficients: np.ndarray) -> np.ndarray:
        matrix = _as_float_matrix(
            coefficients,
            name="sparse design coefficients",
            rows=self.target_design.shape[1],
        )
        scaled = matrix / self.target_scales[:, None]
        raw = np.asarray(self.target_design @ scaled, dtype=np.float64)
        return self._residualize_dense(raw)

    def rmatmat(self, values: np.ndarray) -> np.ndarray:
        matrix = _as_float_matrix(
            values,
            name="sparse design outcome",
            rows=self.target_design.shape[0],
        )
        residual = self._residualize_dense(matrix)
        result = self.target_design.T @ residual
        return np.asarray(result, dtype=np.float64) / self.target_scales[:, None]

    def linear_operator(self) -> sparse_linalg.LinearOperator:
        rows, columns = self.target_design.shape

        def matvec(vector: np.ndarray) -> np.ndarray:
            return self.matmat(np.asarray(vector)[:, None])[:, 0]

        def rmatvec(vector: np.ndarray) -> np.ndarray:
            return self.rmatmat(np.asarray(vector)[:, None])[:, 0]

        return sparse_linalg.LinearOperator(
            shape=(rows, columns),
            matvec=matvec,
            rmatvec=rmatvec,
            matmat=self.matmat,
            rmatmat=self.rmatmat,
            dtype=np.float64,
        )


def _sparse_target_statistics(
    target_design: sparse.csr_matrix,
    centered_nuisance: np.ndarray,
    channels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows = target_design.shape[0]
    counts = np.asarray(target_design.sum(axis=0)).ravel()
    centered_sum_squares = counts.astype(np.float64, copy=True)
    for channel_rows in _stable_channel_rows(channels):
        channel_counts = np.asarray(target_design[channel_rows].sum(axis=0)).ravel()
        centered_sum_squares -= channel_counts * channel_counts / channel_rows.size
    if centered_nuisance.shape[1]:
        nuisance_pseudoinverse = _least_squares_pseudoinverse(centered_nuisance)
        nuisance_cross = np.asarray(
            centered_nuisance.T @ target_design,
            dtype=np.float64,
        )
        target_slopes = np.asarray(
            target_design.T @ nuisance_pseudoinverse.T,
            dtype=np.float64,
        ).T
        residual_sum_squares = centered_sum_squares - np.sum(
            target_slopes * nuisance_cross,
            axis=0,
        )
    else:
        nuisance_pseudoinverse = np.zeros(
            (0, target_design.shape[0]),
            dtype=np.float64,
        )
        target_slopes = np.zeros(
            (0, target_design.shape[1]),
            dtype=np.float64,
        )
        residual_sum_squares = centered_sum_squares.copy()
    centered_sum_squares = np.maximum(centered_sum_squares, 0.0)
    residual_sum_squares = np.maximum(residual_sum_squares, 0.0)
    return (
        np.sqrt(centered_sum_squares / rows),
        np.sqrt(residual_sum_squares / rows),
        nuisance_pseudoinverse,
        target_slopes,
    )


def _canonicalize_row_signs(matrix: np.ndarray) -> np.ndarray:
    result = np.asarray(matrix, dtype=np.float64).copy()
    for row in result:
        pivot = int(np.argmax(np.abs(row)))
        if row[pivot] < 0.0:
            row *= -1.0
    return result


def _validate_sparse_nullspace(
    operator: _SparseResidualizedDesign,
    nullspace: np.ndarray,
    *,
    largest: float,
    threshold: float,
) -> None:
    if not nullspace.size:
        return
    if not np.all(np.isfinite(nullspace)):
        raise ModelingError("sparse target-design nullspace is non-finite")
    gram = _matmul(
        nullspace,
        nullspace.T,
        name="sparse nullspace orthonormality",
    )
    if not np.allclose(
        gram,
        np.eye(nullspace.shape[0]),
        rtol=0.0,
        atol=1e-10,
    ):
        raise ModelingError("sparse target-design nullspace is not orthonormal")
    residual = operator.matmat(nullspace.T)
    relative = np.linalg.norm(residual, axis=0) / largest
    if np.any(relative > threshold / largest + 1e-10):
        raise ModelingError(
            "sparse target-design nullspace violates the singular-value cutoff"
        )
    gauge_projection = nullspace.T - _matmul(
        nullspace.T,
        gram,
        name="sparse nullspace gauge projection",
    )
    scale = max(1.0, float(np.max(np.abs(nullspace), initial=0.0)))
    if float(np.max(np.abs(gauge_projection), initial=0.0)) > 1e-10 * scale:
        raise ModelingError("sparse target-design gauge projection is unstable")


def _sparse_identifiable_nullspace(
    operator: _SparseResidualizedDesign,
) -> tuple[np.ndarray, np.ndarray, int, float]:
    """Return near-null right vectors without constructing an N-by-T matrix."""

    columns = operator.target_design.shape[1]
    if columns == 0:
        return (
            np.zeros((0, 0), dtype=np.float64),
            np.zeros(0, dtype=np.float64),
            0,
            float("inf"),
        )
    linear = operator.linear_operator()
    if columns <= 256:
        dense = operator.matmat(np.eye(columns, dtype=np.float64))
        _, singular_values, right = np.linalg.svd(dense, full_matrices=False)
        retained = (
            singular_values / singular_values[0] > IDENTIFIABLE_RELATIVE_TOLERANCE
            if singular_values.size and singular_values[0] > 0.0
            else np.zeros(singular_values.shape, dtype=bool)
        )
        nullspace = right[~retained]
        effective_rank = int(np.count_nonzero(retained))
        condition = (
            float(singular_values[0] / singular_values[retained][-1])
            if effective_rank
            else float("inf")
        )
        return nullspace, singular_values, effective_rank, condition

    try:
        largest = float(
            sparse_linalg.svds(
                linear,
                k=1,
                which="LM",
                return_singular_vectors=False,
                rng=np.random.default_rng(SPARSE_SVD_SEED),
            )[0]
        )
    except Exception as error:
        raise ModelingError(
            "failed to establish sparse target-design largest singular value"
        ) from error
    if not np.isfinite(largest) or largest <= 0.0:
        raise ModelingError("sparse target design has no identifiable direction")
    threshold = largest * IDENTIFIABLE_RELATIVE_TOLERANCE
    k = min(8, columns - 1)
    while True:
        try:
            _, smallest, right = sparse_linalg.svds(
                linear,
                k=k,
                which="SM",
                return_singular_vectors=True,
                rng=np.random.default_rng(SPARSE_SVD_SEED),
            )
        except Exception as error:
            raise ModelingError(
                "failed to establish sparse target-design nullspace"
            ) from error
        order = np.argsort(smallest)
        smallest = smallest[order]
        right = right[order]
        null_mask = smallest <= threshold
        if np.any(~null_mask) or k == columns - 1:
            nullspace = _canonicalize_row_signs(right[null_mask])
            _validate_sparse_nullspace(
                operator,
                nullspace,
                largest=largest,
                threshold=threshold,
            )
            retained_minimum = (
                float(smallest[~null_mask][0]) if np.any(~null_mask) else largest
            )
            effective_rank = columns - nullspace.shape[0]
            singular_values = np.r_[largest, smallest]
            return (
                nullspace,
                singular_values,
                effective_rank,
                float(largest / retained_minimum),
            )
        k = min(columns - 1, k * 2)


@dataclasses.dataclass(frozen=True)
class _ValidatedDesignGeometry:
    """One audited design geometry bound to its exact sparse inputs."""

    audit: DesignGeometryAudit
    centered_nuisance: np.ndarray
    nuisance_pseudoinverse: np.ndarray
    active_indices: np.ndarray
    active_operator: _SparseResidualizedDesign | None
    active_nullspace: np.ndarray


@dataclasses.dataclass(frozen=True)
class _ValidatedProductionDesignAudit:
    """Validated reusable design geometry for every production population."""

    audit: ProductionDesignAudit
    fold_geometry: tuple[_ValidatedDesignGeometry, ...]
    final_geometry: _ValidatedDesignGeometry


@dataclasses.dataclass(frozen=True)
class _PreparedSparseFitGeometry:
    """Response and prevalidated design geometry for one training population."""

    requested_response_rank: int
    nuisance_response_slopes: np.ndarray
    nuisance_target_slopes: np.ndarray
    response_basis: np.ndarray
    evidence_support: np.ndarray
    zero_variance: np.ndarray
    nuisance_aliased: np.ndarray
    fitted: np.ndarray
    active_indices: np.ndarray
    target_scales: np.ndarray
    active_operator: _SparseResidualizedDesign | None
    ridge_rhs: np.ndarray
    active_nullspace: np.ndarray
    full_nullspace: np.ndarray
    projection_diagonal: np.ndarray
    singular_values: np.ndarray
    retained_singular_directions: np.ndarray
    effective_rank: int
    condition_number: float
    target_count: int
    response_count: int


def _validated_support_mask(values: np.ndarray, target_count: int) -> np.ndarray:
    support = np.asarray(values)
    if support.dtype != np.bool_ or support.shape != (target_count,):
        raise ModelingError("design support must be a Boolean value per target")
    return support.astype(bool, copy=False)


def _compute_sparse_design_geometry(
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    *,
    support: np.ndarray,
    training_channels: Sequence[object],
) -> _ValidatedDesignGeometry:
    """Compute response-free target identifiability and reusable operators."""

    targets = _as_csr_binary_design(
        target_design,
        name="binary_target_design",
    )
    labels = _as_channels(channels, targets.shape[0])
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=targets.shape[0],
    )
    expected_channels = _integer_channel_tuple(
        training_channels,
        name="design audit training_channels",
    )
    if _stable_integer_data_channels(labels) != expected_channels:
        raise ModelingError("design audit rows do not match its training channels")
    support_mask = _validated_support_mask(support, targets.shape[1])

    centered_nuisance = within_channel_center(nuisance_matrix, labels)
    (
        centered_scales,
        residual_scales,
        nuisance_pseudoinverse,
        target_slopes,
    ) = _sparse_target_statistics(
        targets,
        centered_nuisance,
        labels,
    )
    scale_reference = max(
        1.0,
        float(np.max(centered_scales, initial=0.0)),
        float(np.max(residual_scales, initial=0.0)),
    )
    scale_tolerance = np.finfo(np.float64).eps * max(targets.shape) * scale_reference
    zero_variance = centered_scales <= scale_tolerance
    nuisance_aliased = (~zero_variance) & (residual_scales <= scale_tolerance)
    fitted = support_mask & ~zero_variance & ~nuisance_aliased
    active_indices = np.flatnonzero(fitted)
    target_count = targets.shape[1]
    target_scales = np.full(target_count, np.nan, dtype=np.float64)
    target_scales[fitted] = residual_scales[fitted]
    projection_diagonal = np.zeros(target_count, dtype=np.float64)
    full_nullspace = np.zeros((0, target_count), dtype=np.float64)
    active_nullspace = np.zeros((0, active_indices.size), dtype=np.float64)
    active_operator: _SparseResidualizedDesign | None = None
    singular_values = np.zeros(0, dtype=np.float64)
    effective_rank = 0
    condition_number = float("inf")

    if active_indices.size:
        active_operator = _SparseResidualizedDesign(
            target_design=targets[:, active_indices].tocsr(),
            channels=labels,
            centered_nuisance=centered_nuisance,
            nuisance_pseudoinverse=nuisance_pseudoinverse,
            target_scales=residual_scales[active_indices],
        )
        unscaled_operator = dataclasses.replace(
            active_operator,
            target_scales=np.ones(active_indices.size, dtype=np.float64),
        )
        (
            active_nullspace,
            singular_values,
            effective_rank,
            condition_number,
        ) = _sparse_identifiable_nullspace(unscaled_operator)
        if active_nullspace.size:
            full_nullspace = np.zeros(
                (active_nullspace.shape[0], target_count),
                dtype=np.float64,
            )
            full_nullspace[:, active_indices] = active_nullspace
        projection_diagonal[active_indices] = 1.0 - np.sum(
            active_nullspace * active_nullspace,
            axis=0,
        )

    retained_singular_directions = (
        singular_values / singular_values[0] > IDENTIFIABLE_RELATIVE_TOLERANCE
        if singular_values.size and singular_values[0] > 0.0
        else np.zeros(singular_values.shape, dtype=bool)
    )
    audit = DesignGeometryAudit(
        training_channels=np.asarray(expected_channels, dtype=np.int64),
        evidence_support=support_mask.copy(),
        zero_variance=zero_variance,
        nuisance_aliased=nuisance_aliased,
        fitted=fitted,
        target_scales=target_scales,
        nuisance_target_slopes=target_slopes,
        nullspace_basis=full_nullspace,
        projection_diagonal=projection_diagonal,
        singular_values=singular_values,
        retained_singular_directions=retained_singular_directions,
        effective_rank=effective_rank,
        condition_number=condition_number,
    )
    return _ValidatedDesignGeometry(
        audit=audit,
        centered_nuisance=centered_nuisance,
        nuisance_pseudoinverse=nuisance_pseudoinverse,
        active_indices=active_indices,
        active_operator=active_operator,
        active_nullspace=active_nullspace,
    )


def _validate_sparse_design_geometry(
    audit: DesignGeometryAudit,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    *,
    support: np.ndarray,
    training_channels: Sequence[object],
) -> _ValidatedDesignGeometry:
    """Bind a serialized design audit to exact design inputs without refitting SVD."""

    if not isinstance(audit, DesignGeometryAudit):
        raise ModelingError("design audit geometry has the wrong type")
    targets = _as_csr_binary_design(
        target_design,
        name="binary_target_design",
    )
    labels = _as_channels(channels, targets.shape[0])
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=targets.shape[0],
    )
    expected_channels = np.asarray(
        _integer_channel_tuple(
            training_channels,
            name="design audit training_channels",
        ),
        dtype=np.int64,
    )
    if _stable_integer_data_channels(labels) != tuple(expected_channels.tolist()):
        raise ModelingError("design audit rows do not match its training channels")
    if not np.array_equal(audit.training_channels, expected_channels):
        raise ModelingError("design audit training-channel axis changed")
    support_mask = _validated_support_mask(support, targets.shape[1])
    target_count = targets.shape[1]
    boolean_fields = {
        "evidence_support": audit.evidence_support,
        "zero_variance": audit.zero_variance,
        "nuisance_aliased": audit.nuisance_aliased,
        "fitted": audit.fitted,
        "retained_singular_directions": audit.retained_singular_directions,
    }
    for name, value in boolean_fields.items():
        array = np.asarray(value)
        expected_shape = (
            np.asarray(audit.singular_values).shape
            if name == "retained_singular_directions"
            else (target_count,)
        )
        if array.dtype != np.bool_ or array.shape != expected_shape:
            raise ModelingError(f"design audit {name} shape or dtype changed")
    if not np.array_equal(audit.evidence_support, support_mask):
        raise ModelingError("design audit evidence support changed")
    if np.asarray(audit.target_scales).shape != (target_count,):
        raise ModelingError("design audit target scales changed")
    if np.asarray(audit.nuisance_target_slopes).shape != (
        nuisance_matrix.shape[1],
        target_count,
    ):
        raise ModelingError("design audit nuisance target slopes changed")
    if np.asarray(audit.projection_diagonal).shape != (target_count,):
        raise ModelingError("design audit projection diagonal changed")

    centered_nuisance = within_channel_center(nuisance_matrix, labels)
    (
        centered_scales,
        residual_scales,
        nuisance_pseudoinverse,
        target_slopes,
    ) = _sparse_target_statistics(
        targets,
        centered_nuisance,
        labels,
    )
    scale_reference = max(
        1.0,
        float(np.max(centered_scales, initial=0.0)),
        float(np.max(residual_scales, initial=0.0)),
    )
    scale_tolerance = np.finfo(np.float64).eps * max(targets.shape) * scale_reference
    zero_variance = centered_scales <= scale_tolerance
    nuisance_aliased = (~zero_variance) & (residual_scales <= scale_tolerance)
    fitted = support_mask & ~zero_variance & ~nuisance_aliased
    target_scales = np.full(target_count, np.nan, dtype=np.float64)
    target_scales[fitted] = residual_scales[fitted]
    for name, observed, expected in (
        ("zero_variance", audit.zero_variance, zero_variance),
        ("nuisance_aliased", audit.nuisance_aliased, nuisance_aliased),
        ("fitted", audit.fitted, fitted),
    ):
        if not np.array_equal(observed, expected):
            raise ModelingError(f"design audit {name} changed")
    if not np.allclose(
        audit.target_scales,
        target_scales,
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    ):
        raise ModelingError("design audit target scales changed")
    if not np.allclose(
        audit.nuisance_target_slopes,
        target_slopes,
        rtol=1e-12,
        atol=1e-12,
    ):
        raise ModelingError("design audit nuisance target slopes changed")

    active_indices = np.flatnonzero(fitted)
    nullspace = _as_float_matrix(
        audit.nullspace_basis,
        name="design audit nullspace",
    )
    if nullspace.shape[1] != target_count:
        raise ModelingError("design audit nullspace target axis changed")
    if active_indices.size != audit.effective_rank + nullspace.shape[0]:
        raise ModelingError("design audit effective rank changed")
    inactive = ~fitted
    if nullspace.size and np.any(nullspace[:, inactive] != 0.0):
        raise ModelingError("design audit nullspace leaks into unfitted targets")
    active_nullspace = nullspace[:, active_indices]
    projection_diagonal = np.zeros(target_count, dtype=np.float64)
    projection_diagonal[active_indices] = 1.0 - np.sum(
        active_nullspace * active_nullspace,
        axis=0,
    )
    if not np.allclose(
        audit.projection_diagonal,
        projection_diagonal,
        rtol=0.0,
        atol=1e-12,
    ):
        raise ModelingError("design audit projection diagonal changed")

    active_operator: _SparseResidualizedDesign | None = None
    singular_values = np.asarray(audit.singular_values, dtype=np.float64)
    if singular_values.ndim != 1 or not np.all(np.isfinite(singular_values)):
        raise ModelingError("design audit singular values changed")
    if active_indices.size:
        if not singular_values.size or singular_values[0] <= 0.0:
            raise ModelingError("design audit largest singular value is invalid")
        expected_retained = (
            singular_values / singular_values[0] > IDENTIFIABLE_RELATIVE_TOLERANCE
        )
        if not np.array_equal(
            audit.retained_singular_directions,
            expected_retained,
        ):
            raise ModelingError("design audit retained singular directions changed")
        if np.count_nonzero(~expected_retained) != active_nullspace.shape[0]:
            raise ModelingError("design audit nullspace dimension changed")
        active_operator = _SparseResidualizedDesign(
            target_design=targets[:, active_indices].tocsr(),
            channels=labels,
            centered_nuisance=centered_nuisance,
            nuisance_pseudoinverse=nuisance_pseudoinverse,
            target_scales=residual_scales[active_indices],
        )
        _validate_sparse_nullspace(
            dataclasses.replace(
                active_operator,
                target_scales=np.ones(active_indices.size, dtype=np.float64),
            ),
            active_nullspace,
            largest=float(singular_values[0]),
            threshold=float(singular_values[0]) * IDENTIFIABLE_RELATIVE_TOLERANCE,
        )
        if not np.isfinite(audit.condition_number) or audit.condition_number < 1.0:
            raise ModelingError("design audit condition number is invalid")
    elif (
        singular_values.size
        or nullspace.size
        or audit.effective_rank != 0
        or not np.isinf(audit.condition_number)
    ):
        raise ModelingError("empty design audit geometry changed")

    return _ValidatedDesignGeometry(
        audit=audit,
        centered_nuisance=centered_nuisance,
        nuisance_pseudoinverse=nuisance_pseudoinverse,
        active_indices=active_indices,
        active_operator=active_operator,
        active_nullspace=active_nullspace,
    )


def audit_production_design(
    binary_target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: Sequence[object],
    guide_evidence: GuideEvidence,
    *,
    folds: Sequence[Sequence[object]] = FROZEN_HELD_CHANNEL_FOLDS,
) -> ProductionDesignAudit:
    """Compute every support, alias, rank, and condition diagnostic without Y."""

    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=targets.shape[0],
    )
    labels = _as_channels(channels, targets.shape[0])
    evidence = _validate_guide_evidence(
        guide_evidence,
        rows=targets.shape[0],
        target_count=targets.shape[1],
    )
    _validate_frozen_production_channels(labels)
    _validate_evidence_matches_design(evidence, targets)
    normalized_folds = tuple(
        _integer_channel_tuple(fold, name="held fold") for fold in folds
    )
    if normalized_folds != FROZEN_HELD_CHANNEL_FOLDS:
        raise ModelingError(f"production folds must equal {FROZEN_HELD_CHANNEL_FOLDS}")

    fold_support: list[TargetSupportEvidence] = []
    fold_geometry: list[DesignGeometryAudit] = []
    for held_channels in normalized_folds:
        training_rows = ~np.isin(labels, held_channels)
        training_channels = tuple(
            channel for channel in FROZEN_CHANNELS if channel not in held_channels
        )
        support = _derive_target_support_validated(
            evidence,
            binary_target_design=targets,
            channels=labels,
            training_channels=training_channels,
            target_count=targets.shape[1],
        )
        computed = _compute_sparse_design_geometry(
            targets[training_rows].tocsr(),
            nuisance_matrix[training_rows],
            labels[training_rows],
            support=support.supported,
            training_channels=training_channels,
        )
        fold_support.append(support)
        fold_geometry.append(computed.audit)

    final_support = _derive_target_support_validated(
        evidence,
        binary_target_design=targets,
        channels=labels,
        training_channels=FROZEN_CHANNELS,
        target_count=targets.shape[1],
    )
    final_geometry = _compute_sparse_design_geometry(
        targets,
        nuisance_matrix,
        labels,
        support=final_support.supported,
        training_channels=FROZEN_CHANNELS,
    ).audit
    return ProductionDesignAudit(
        held_channel_folds=normalized_folds,
        fold_support=tuple(fold_support),
        fold_geometry=tuple(fold_geometry),
        final_support=final_support,
        final_geometry=final_geometry,
    )


def validate_production_design_audit(
    audit: ProductionDesignAudit,
    binary_target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: Sequence[object],
    guide_evidence: GuideEvidence,
) -> _ValidatedProductionDesignAudit:
    """Validate and materialize a response-free audit before any response fit."""

    if not isinstance(audit, ProductionDesignAudit):
        raise ModelingError("production design audit has the wrong type")
    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=targets.shape[0],
    )
    labels = _as_channels(channels, targets.shape[0])
    evidence = _validate_guide_evidence(
        guide_evidence,
        rows=targets.shape[0],
        target_count=targets.shape[1],
    )
    _validate_frozen_production_channels(labels)
    _validate_evidence_matches_design(evidence, targets)
    if audit.held_channel_folds != FROZEN_HELD_CHANNEL_FOLDS:
        raise ModelingError("production design audit fold definition changed")
    if len(audit.fold_support) != len(FROZEN_HELD_CHANNEL_FOLDS) or len(
        audit.fold_geometry
    ) != len(FROZEN_HELD_CHANNEL_FOLDS):
        raise ModelingError("production design audit fold count changed")

    validated_folds: list[_ValidatedDesignGeometry] = []
    for index, held_channels in enumerate(FROZEN_HELD_CHANNEL_FOLDS):
        training_rows = ~np.isin(labels, held_channels)
        training_channels = tuple(
            channel for channel in FROZEN_CHANNELS if channel not in held_channels
        )
        observed_support = _derive_target_support_validated(
            evidence,
            binary_target_design=targets,
            channels=labels,
            training_channels=training_channels,
            target_count=targets.shape[1],
        )
        expected_support = audit.fold_support[index]
        for field in dataclasses.fields(TargetSupportEvidence):
            if not np.array_equal(
                getattr(expected_support, field.name),
                getattr(observed_support, field.name),
            ):
                raise ModelingError(
                    f"production design audit fold {index} support changed at "
                    f"{field.name}"
                )
        validated_folds.append(
            _validate_sparse_design_geometry(
                audit.fold_geometry[index],
                targets[training_rows].tocsr(),
                nuisance_matrix[training_rows],
                labels[training_rows],
                support=expected_support.supported,
                training_channels=training_channels,
            )
        )

    observed_final_support = _derive_target_support_validated(
        evidence,
        binary_target_design=targets,
        channels=labels,
        training_channels=FROZEN_CHANNELS,
        target_count=targets.shape[1],
    )
    for field in dataclasses.fields(TargetSupportEvidence):
        if not np.array_equal(
            getattr(audit.final_support, field.name),
            getattr(observed_final_support, field.name),
        ):
            raise ModelingError(
                f"production design audit final support changed at {field.name}"
            )
    validated_final = _validate_sparse_design_geometry(
        audit.final_geometry,
        targets,
        nuisance_matrix,
        labels,
        support=audit.final_support.supported,
        training_channels=FROZEN_CHANNELS,
    )
    return _ValidatedProductionDesignAudit(
        audit=audit,
        fold_geometry=tuple(validated_folds),
        final_geometry=validated_final,
    )


def _batched_spd_ridge_solve(
    operator: _SparseResidualizedDesign,
    ridge_rhs: np.ndarray,
    ridge_lambda: float,
    *,
    initial: np.ndarray | None = None,
    relative_tolerance: float = 1e-8,
    absolute_tolerance: float = 1e-10,
    max_iterations: int = 1000,
) -> tuple[np.ndarray, RidgeSolveDiagnostics]:
    """Solve all ridge right-hand sides with deterministic batched CG."""

    right = _as_float_matrix(
        ridge_rhs,
        name="ridge_rhs",
        rows=operator.target_design.shape[1],
    )
    if ridge_lambda <= 0.0 or not np.isfinite(ridge_lambda):
        raise ModelingError("batched ridge_lambda must be finite and positive")
    if relative_tolerance > 1e-8 or absolute_tolerance > 1e-10:
        raise ModelingError("batched ridge tolerances cannot be weakened")
    if relative_tolerance <= 0.0 or absolute_tolerance < 0.0:
        raise ModelingError("batched ridge tolerances must be positive")
    if max_iterations < 1:
        raise ModelingError("max_iterations must be positive")
    if initial is None:
        solution = np.zeros_like(right)
        residual = right.copy()
        operator_applications = 0
        warm_started = False
    else:
        solution = _as_float_matrix(
            initial,
            name="initial ridge coefficients",
            rows=right.shape[0],
        ).copy()
        if solution.shape != right.shape:
            raise ModelingError("initial ridge coefficients have the wrong shape")
        residual = right - (
            operator.rmatmat(operator.matmat(solution)) + ridge_lambda * solution
        )
        operator_applications = 1
        warm_started = True

    rhs_norms = np.linalg.norm(right, axis=0)
    thresholds = absolute_tolerance + relative_tolerance * rhs_norms
    residual_squared = np.sum(residual * residual, axis=0)
    residual_norms = np.sqrt(residual_squared)
    converged = residual_norms <= thresholds
    directions = residual.copy()
    directions[:, converged] = 0.0
    iterations = 0

    while not np.all(converged) and iterations < max_iterations:
        image = (
            operator.rmatmat(operator.matmat(directions)) + ridge_lambda * directions
        )
        operator_applications += 1
        active = ~converged
        denominators = np.sum(directions * image, axis=0)
        if np.any(~np.isfinite(denominators[active])) or np.any(
            denominators[active] <= 0.0
        ):
            raise ModelingError("batched ridge operator is not numerically SPD")
        step = np.zeros(right.shape[1], dtype=np.float64)
        step[active] = residual_squared[active] / denominators[active]
        solution += directions * step
        residual -= image * step
        next_squared = np.sum(residual * residual, axis=0)
        next_norms = np.sqrt(next_squared)
        next_converged = next_norms <= thresholds
        continuing = ~next_converged
        conjugacy = np.zeros(right.shape[1], dtype=np.float64)
        conjugacy[continuing] = next_squared[continuing] / residual_squared[continuing]
        directions = residual + directions * conjugacy
        directions[:, next_converged] = 0.0
        residual_squared = next_squared
        converged = next_converged
        iterations += 1

    final_residual = right - (
        operator.rmatmat(operator.matmat(solution)) + ridge_lambda * solution
    )
    operator_applications += 1
    residual_norms = np.linalg.norm(final_residual, axis=0)
    converged = residual_norms <= thresholds
    relative_residuals = residual_norms / np.maximum(
        rhs_norms,
        np.finfo(np.float64).tiny,
    )
    diagnostics = RidgeSolveDiagnostics(
        ridge_lambda=float(ridge_lambda),
        iterations=iterations,
        operator_applications=operator_applications,
        warm_started=warm_started,
        converged=converged,
        rhs_norms=rhs_norms,
        residual_norms=residual_norms,
        relative_residual_norms=relative_residuals,
        absolute_tolerance=float(absolute_tolerance),
        relative_tolerance=float(relative_tolerance),
    )
    if not np.all(converged):
        worst = float(np.max(relative_residuals[~converged]))
        raise ModelingError(
            "batched ridge solver failed strict convergence: "
            f"{np.count_nonzero(~converged)} components, worst relative "
            f"residual {worst:.3e}, iterations {iterations}"
        )
    return solution, diagnostics


def _prepare_sparse_fit_geometry(
    response: np.ndarray,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    *,
    support: np.ndarray,
    rank: int,
    design_geometry: _ValidatedDesignGeometry | None = None,
) -> _PreparedSparseFitGeometry:
    """Attach response-only quantities to prevalidated design geometry."""

    if design_geometry is None:
        design_geometry = _compute_sparse_design_geometry(
            target_design,
            nuisance,
            channels,
            support=support,
            training_channels=_stable_integer_data_channels(channels),
        )
    design = design_geometry.audit
    centered_nuisance = design_geometry.centered_nuisance
    centered_response = within_channel_center(response, channels)
    response_slopes = _nuisance_slopes(centered_nuisance, centered_response)
    residual_response = centered_response - _matmul(
        centered_nuisance,
        response_slopes,
        name="sparse training response nuisance prediction",
    )
    target_count = target_design.shape[1]
    response_count = response.shape[1]
    response_basis = _canonical_response_basis_from_covariance(
        residual_response,
        rank,
    )
    ridge_rhs = np.zeros(
        (design_geometry.active_indices.size, response_basis.shape[0]),
        dtype=np.float64,
    )

    if design_geometry.active_operator is not None and response_basis.shape[0]:
        response_scores = _matmul(
            residual_response,
            response_basis.T,
            name="sparse training response scores",
        )
        ridge_rhs = design_geometry.active_operator.rmatmat(response_scores)
    return _PreparedSparseFitGeometry(
        requested_response_rank=rank,
        nuisance_response_slopes=response_slopes,
        nuisance_target_slopes=design.nuisance_target_slopes,
        response_basis=response_basis,
        evidence_support=design.evidence_support.copy(),
        zero_variance=design.zero_variance,
        nuisance_aliased=design.nuisance_aliased,
        fitted=design.fitted,
        active_indices=design_geometry.active_indices,
        target_scales=design.target_scales,
        active_operator=design_geometry.active_operator,
        ridge_rhs=ridge_rhs,
        active_nullspace=design_geometry.active_nullspace,
        full_nullspace=design.nullspace_basis,
        projection_diagonal=design.projection_diagonal,
        singular_values=design.singular_values,
        retained_singular_directions=design.retained_singular_directions,
        effective_rank=design.effective_rank,
        condition_number=design.condition_number,
        target_count=target_count,
        response_count=response_count,
    )


def _fit_prepared_sparse_geometry(
    geometry: _PreparedSparseFitGeometry,
    *,
    ridge_lambda: float,
    initial_scaled_coefficients: np.ndarray | None = None,
) -> tuple[FittedEstimator, np.ndarray]:
    """Solve one lambda while reusing all prepared fold geometry."""

    coefficients = np.zeros(
        (geometry.target_count, geometry.response_count),
        dtype=np.float64,
    )
    projected = np.zeros_like(coefficients)
    scaled_coefficients = np.zeros_like(geometry.ridge_rhs)
    projection_postconditions = (0.0, 0.0, 0.0)
    if geometry.active_operator is None:
        solve_diagnostics = RidgeSolveDiagnostics(
            ridge_lambda=float(ridge_lambda),
            iterations=0,
            operator_applications=0,
            warm_started=initial_scaled_coefficients is not None,
            converged=np.ones(geometry.ridge_rhs.shape[1], dtype=bool),
            rhs_norms=np.zeros(geometry.ridge_rhs.shape[1], dtype=np.float64),
            residual_norms=np.zeros(geometry.ridge_rhs.shape[1], dtype=np.float64),
            relative_residual_norms=np.zeros(
                geometry.ridge_rhs.shape[1],
                dtype=np.float64,
            ),
            absolute_tolerance=1e-10,
            relative_tolerance=1e-8,
        )
    else:
        scaled_coefficients, solve_diagnostics = _batched_spd_ridge_solve(
            geometry.active_operator,
            geometry.ridge_rhs,
            ridge_lambda,
            initial=initial_scaled_coefficients,
        )
        active_coefficients = (
            _matmul(
                scaled_coefficients,
                geometry.response_basis,
                name="sparse binary-unit coefficient reconstruction",
            )
            / geometry.target_scales[geometry.active_indices, None]
        )
        if geometry.active_nullspace.size:
            active_projected = active_coefficients - _matmul(
                geometry.active_nullspace.T,
                _matmul(
                    geometry.active_nullspace,
                    active_coefficients,
                    name="sparse nullspace coefficient coordinates",
                ),
                name="sparse identifiable projection",
            )
        else:
            active_projected = active_coefficients
        projection_postconditions = validate_production_projection(
            active_coefficients,
            active_projected,
            geometry.active_nullspace,
        )
        coefficients[geometry.active_indices] = active_coefficients
        projected[geometry.active_indices] = active_projected

    transferable = geometry.fitted & (
        geometry.projection_diagonal >= TRANSFERABLE_PROJECTION_DIAGONAL
    )
    exported = np.zeros_like(projected)
    exported[transferable] = projected[transferable]
    diagnostics = SupportDiagnostics(
        evidence_support=geometry.evidence_support.copy(),
        zero_variance=geometry.zero_variance,
        nuisance_aliased=geometry.nuisance_aliased,
        fitted=geometry.fitted,
        projection_diagonal=geometry.projection_diagonal,
        transferable=transferable,
        target_scales=geometry.target_scales,
        singular_values=geometry.singular_values,
        retained_singular_directions=geometry.retained_singular_directions,
        effective_rank=geometry.effective_rank,
        condition_number=geometry.condition_number,
        relative_singular_tolerance=IDENTIFIABLE_RELATIVE_TOLERANCE,
        projection_diagonal_threshold=TRANSFERABLE_PROJECTION_DIAGONAL,
        nullspace_annihilation_error=projection_postconditions[0],
        projection_idempotence_error=projection_postconditions[1],
        gauge_shift_invariance_error=projection_postconditions[2],
    )
    estimator = FittedEstimator(
        ridge_lambda=ridge_lambda,
        requested_response_rank=geometry.requested_response_rank,
        effective_response_rank=geometry.response_basis.shape[0],
        nuisance_response_slopes=geometry.nuisance_response_slopes,
        nuisance_target_slopes=geometry.nuisance_target_slopes,
        response_basis=geometry.response_basis,
        coefficients_binary_units=coefficients,
        identifiable_projection=None,
        nullspace_basis=geometry.full_nullspace,
        projected_coefficients=projected,
        exported_effects=exported,
        diagnostics=diagnostics,
        ridge_solve_diagnostics=solve_diagnostics,
    )
    return estimator, scaled_coefficients


def _fit_sparse_estimator(
    response: np.ndarray,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    *,
    support: np.ndarray,
    ridge_lambda: float,
    rank: int,
) -> FittedEstimator:
    geometry = _prepare_sparse_fit_geometry(
        response,
        target_design,
        nuisance,
        channels,
        support=support,
        rank=rank,
    )
    estimator, _ = _fit_prepared_sparse_geometry(
        geometry,
        ridge_lambda=ridge_lambda,
    )
    return estimator


def fit_estimator(
    response: np.ndarray,
    binary_target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: Sequence[object],
    guide_evidence: GuideEvidence,
    *,
    ridge_lambda: float,
    rank: int = FROZEN_RESPONSE_RANK,
) -> FittedEstimator:
    """Fit the final all-channel estimator under the frozen production contract."""

    penalty, requested_rank = _validate_frozen_fit_controls(ridge_lambda, rank)
    response_matrix = _as_float_matrix(response, name="response")
    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
        rows=response_matrix.shape[0],
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=response_matrix.shape[0],
    )
    labels = _as_channels(channels, response_matrix.shape[0])
    evidence = _validate_guide_evidence(
        guide_evidence,
        rows=response_matrix.shape[0],
        target_count=targets.shape[1],
    )
    _validate_frozen_production_channels(labels)
    _validate_evidence_matches_design(evidence, targets)
    support = _derive_target_support_validated(
        evidence,
        binary_target_design=targets,
        channels=labels,
        training_channels=FROZEN_CHANNELS,
        target_count=targets.shape[1],
    )
    return _fit_sparse_estimator(
        response_matrix,
        targets,
        nuisance_matrix,
        labels,
        ridge_lambda=penalty,
        rank=requested_rank,
        support=support.supported,
    )


def _predict_sparse_held(
    estimator: FittedEstimator,
    target_design: sparse.csr_matrix,
    channels: np.ndarray,
    centered_nuisance: np.ndarray,
) -> np.ndarray:
    raw_prediction = np.asarray(
        target_design @ estimator.exported_effects,
        dtype=np.float64,
    )
    prediction = within_channel_center(raw_prediction, channels)
    if centered_nuisance.shape[1]:
        prediction -= _matmul(
            centered_nuisance,
            _matmul(
                estimator.nuisance_target_slopes,
                estimator.exported_effects,
                name="sparse held target nuisance coefficients",
            ),
            name="sparse held target nuisance prediction",
        )
    return prediction


def _evaluate_sparse_held(
    estimator: FittedEstimator,
    response: np.ndarray,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
) -> HeldChannelResult:
    centered_response = within_channel_center(response, channels)
    centered_nuisance = within_channel_center(nuisance, channels)
    residual_response = centered_response - _matmul(
        centered_nuisance,
        estimator.nuisance_response_slopes,
        name="sparse held response nuisance prediction",
    )
    prediction = _predict_sparse_held(
        estimator,
        target_design,
        channels,
        centered_nuisance,
    )
    unmodeled = np.asarray(
        target_design[:, estimator.diagnostics.unsupported].sum(axis=1)
    ).ravel()
    return HeldChannelResult(
        centered_response=centered_response,
        centered_target_design=None,
        centered_nuisance=centered_nuisance,
        residual_response=residual_response,
        residual_target_design=None,
        predicted_residual_response=prediction,
        unmodeled_target_burden=unmodeled,
    )


@dataclasses.dataclass(frozen=True)
class ChannelFoldResult:
    """A training-only fit plus held transforms from complete held channels."""

    estimator: FittedEstimator
    target_support: TargetSupportEvidence
    held_channels: tuple[int, ...]
    training_rows: np.ndarray
    held_rows: np.ndarray
    held: HeldChannelResult


def fit_channel_fold(
    response: np.ndarray,
    binary_target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: Sequence[object],
    guide_evidence: GuideEvidence,
    *,
    held_channels: Sequence[object],
    ridge_lambda: float,
    rank: int = FROZEN_RESPONSE_RANK,
) -> ChannelFoldResult:
    """Fit one frozen complete-channel fold with training-only support."""

    penalty, requested_rank = _validate_frozen_fit_controls(ridge_lambda, rank)
    response_matrix = _as_float_matrix(response, name="response")
    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
        rows=response_matrix.shape[0],
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=response_matrix.shape[0],
    )
    labels = _as_channels(channels, response_matrix.shape[0])
    evidence = _validate_guide_evidence(
        guide_evidence,
        rows=response_matrix.shape[0],
        target_count=targets.shape[1],
    )
    _validate_frozen_production_channels(labels)
    _validate_evidence_matches_design(evidence, targets)
    held_channel_tuple = _integer_channel_tuple(
        held_channels,
        name="held_channels",
    )
    if held_channel_tuple not in FROZEN_HELD_CHANNEL_FOLDS:
        raise ModelingError(f"held_channels must be one of {FROZEN_HELD_CHANNEL_FOLDS}")
    held_rows = np.isin(labels, held_channel_tuple)
    train = ~held_rows
    training_channels = tuple(
        channel for channel in FROZEN_CHANNELS if channel not in held_channel_tuple
    )
    support = _derive_target_support_validated(
        evidence,
        binary_target_design=targets,
        channels=labels,
        training_channels=training_channels,
        target_count=targets.shape[1],
    )

    estimator = _fit_sparse_estimator(
        response_matrix[train],
        targets[train].tocsr(),
        nuisance_matrix[train],
        labels[train],
        ridge_lambda=penalty,
        rank=requested_rank,
        support=support.supported,
    )
    held = _evaluate_sparse_held(
        estimator,
        response_matrix[held_rows],
        targets[held_rows].tocsr(),
        nuisance_matrix[held_rows],
        labels[held_rows],
    )
    return ChannelFoldResult(
        estimator=estimator,
        target_support=support,
        held_channels=held_channel_tuple,
        training_rows=train.copy(),
        held_rows=held_rows.copy(),
        held=held,
    )


def _held_response_inverse_variance_weights(
    residual_response: np.ndarray,
) -> np.ndarray:
    variances = np.var(residual_response, axis=0, ddof=0)
    reference = max(1.0, float(np.max(variances, initial=0.0)))
    tolerance = np.finfo(np.float64).eps * max(residual_response.shape) * reference
    informative = variances > tolerance
    if not np.any(informative):
        raise ModelingError("held response has no non-zero-variance coordinate")
    weights = np.zeros(variances.shape, dtype=np.float64)
    weights[informative] = 1.0 / variances[informative]
    return weights


def _mean_held_cell_variance_weighted_sse_gain(
    residual_response: np.ndarray,
    prediction: np.ndarray,
    weights: np.ndarray,
) -> float:
    if residual_response.shape != prediction.shape:
        raise ModelingError("held prediction and response shapes differ")
    if weights.shape != (residual_response.shape[1],):
        raise ModelingError("held response weights have the wrong shape")
    denominator = float(np.sum(weights))
    if not np.isfinite(denominator) or denominator <= 0.0:
        raise ModelingError("held response weights have no positive mass")
    residual = residual_response - prediction
    per_cell_gain = (
        (residual_response * residual_response - residual * residual) * weights
    ).sum(axis=1) / denominator
    if not np.all(np.isfinite(per_cell_gain)):
        raise ModelingError("held-cell SSE gain is non-finite")
    return float(np.mean(per_cell_gain))


def _largest_frozen_lambda_within_relative_gain(
    mean_gains: np.ndarray,
) -> float:
    gains = np.asarray(mean_gains, dtype=np.float64)
    if gains.shape != (len(FROZEN_LAMBDA_GRID),):
        raise ModelingError("mean_gains must have one value per frozen lambda")
    if not np.all(np.isfinite(gains)):
        raise ModelingError("mean_gains contains non-finite values")
    best = float(np.max(gains))
    tolerance = LAMBDA_RELATIVE_GAIN_TOLERANCE * max(
        abs(best),
        np.finfo(np.float64).tiny,
    )
    eligible = gains >= best - tolerance
    return max(
        penalty
        for penalty, is_eligible in zip(
            FROZEN_LAMBDA_GRID,
            eligible,
            strict=True,
        )
        if is_eligible
    )


@dataclasses.dataclass(frozen=True)
class FoldGeometryDiagnostics:
    """Lambda-invariant geometry and held-population diagnostics for one fold."""

    zero_variance: np.ndarray
    nuisance_aliased: np.ndarray
    projection_diagonal: np.ndarray
    transferable: np.ndarray
    effective_rank: int
    condition_number: float
    held_unmodeled_target_burden: np.ndarray


@dataclasses.dataclass(frozen=True)
class ProductionSelectionResult:
    """Frozen four-fold lambda selection and all-channel refit."""

    selected_lambda: float
    estimator: FittedEstimator
    lambda_grid: np.ndarray
    held_channel_folds: tuple[tuple[int, ...], ...]
    fold_mean_gains: np.ndarray
    mean_held_cell_gains: np.ndarray
    held_cell_counts: np.ndarray
    fold_support: tuple[TargetSupportEvidence, ...]
    fold_solver_diagnostics: tuple[tuple[RidgeSolveDiagnostics, ...], ...]
    fold_geometry_diagnostics: tuple[FoldGeometryDiagnostics, ...]
    fold_selected_effects: np.ndarray
    final_support: TargetSupportEvidence


def select_production_estimator(
    response: np.ndarray,
    binary_target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: Sequence[object],
    guide_evidence: GuideEvidence,
    *,
    lambda_grid: Sequence[float] = FROZEN_LAMBDA_GRID,
    rank: int = FROZEN_RESPONSE_RANK,
    folds: Sequence[Sequence[object]] = FROZEN_HELD_CHANNEL_FOLDS,
    design_audit: ProductionDesignAudit | _ValidatedProductionDesignAudit | None = None,
) -> ProductionSelectionResult:
    """Run the frozen selector on a response-free audited design geometry."""

    penalties, requested_rank, held_folds = _validate_frozen_selector_controls(
        lambda_grid,
        rank,
        folds,
    )
    response_matrix = _as_float_matrix(response, name="response")
    targets = _as_csr_binary_design(
        binary_target_design,
        name="binary_target_design",
        rows=response_matrix.shape[0],
    )
    nuisance_matrix = _as_float_matrix(
        nuisance,
        name="nuisance",
        rows=response_matrix.shape[0],
    )
    labels = _as_channels(channels, response_matrix.shape[0])
    evidence = _validate_guide_evidence(
        guide_evidence,
        rows=response_matrix.shape[0],
        target_count=targets.shape[1],
    )
    _validate_frozen_production_channels(labels)
    _validate_evidence_matches_design(evidence, targets)
    if design_audit is None:
        design_audit = audit_production_design(
            targets,
            nuisance_matrix,
            labels,
            GuideEvidence(
                guide_presence=evidence.guide_presence,
                guide_target_indices=evidence.guide_target_indices,
            ),
            folds=held_folds,
        )
    validated_audit = (
        design_audit
        if isinstance(design_audit, _ValidatedProductionDesignAudit)
        else validate_production_design_audit(
            design_audit,
            targets,
            nuisance_matrix,
            labels,
            GuideEvidence(
                guide_presence=evidence.guide_presence,
                guide_target_indices=evidence.guide_target_indices,
            ),
        )
    )
    if validated_audit.audit.held_channel_folds != held_folds:
        raise ModelingError("validated production design audit folds changed")

    fold_mean_gains = np.empty(
        (len(held_folds), len(penalties)),
        dtype=np.float64,
    )
    held_cell_counts = np.empty(len(held_folds), dtype=np.int64)
    fold_support: list[TargetSupportEvidence] = []
    fold_solver_diagnostics: list[tuple[RidgeSolveDiagnostics, ...]] = []
    fold_geometry_diagnostics: list[FoldGeometryDiagnostics] = []
    fold_effect_candidates: list[tuple[np.ndarray, ...]] = []
    for fold_index, held_channels in enumerate(held_folds):
        held_rows = np.isin(labels, held_channels)
        training_rows = ~held_rows
        support = validated_audit.audit.fold_support[fold_index]
        fold_support.append(support)
        held_cell_counts[fold_index] = int(np.count_nonzero(held_rows))
        geometry = _prepare_sparse_fit_geometry(
            response_matrix[training_rows],
            targets[training_rows].tocsr(),
            nuisance_matrix[training_rows],
            labels[training_rows],
            support=support.supported,
            rank=requested_rank,
            design_geometry=validated_audit.fold_geometry[fold_index],
        )
        held_centered_response = within_channel_center(
            response_matrix[held_rows],
            labels[held_rows],
        )
        held_centered_nuisance = within_channel_center(
            nuisance_matrix[held_rows],
            labels[held_rows],
        )
        held_residual_response = held_centered_response - _matmul(
            held_centered_nuisance,
            geometry.nuisance_response_slopes,
            name="selector held response nuisance prediction",
        )
        held_weights = _held_response_inverse_variance_weights(held_residual_response)
        warm_start: np.ndarray | None = None
        lambda_diagnostics: list[RidgeSolveDiagnostics] = []
        lambda_effects: list[np.ndarray] = []
        final_fold_estimator: FittedEstimator | None = None
        for lambda_index, penalty in enumerate(penalties):
            estimator, warm_start = _fit_prepared_sparse_geometry(
                geometry,
                ridge_lambda=penalty,
                initial_scaled_coefficients=warm_start,
            )
            if estimator.ridge_solve_diagnostics is None:
                raise ModelingError("production ridge diagnostics are missing")
            final_fold_estimator = estimator
            lambda_diagnostics.append(estimator.ridge_solve_diagnostics)
            lambda_effects.append(estimator.exported_effects.astype(np.float32))
            held_prediction = _predict_sparse_held(
                estimator,
                targets[held_rows].tocsr(),
                labels[held_rows],
                held_centered_nuisance,
            )
            fold_mean_gains[fold_index, lambda_index] = (
                _mean_held_cell_variance_weighted_sse_gain(
                    held_residual_response,
                    held_prediction,
                    held_weights,
                )
            )
        if final_fold_estimator is None:
            raise ModelingError("frozen lambda grid produced no fold estimator")
        fold_solver_diagnostics.append(tuple(lambda_diagnostics))
        fold_effect_candidates.append(tuple(lambda_effects))
        held_unmodeled = np.asarray(
            targets[held_rows][:, final_fold_estimator.diagnostics.unsupported].sum(
                axis=1
            )
        ).ravel()
        fold_geometry_diagnostics.append(
            FoldGeometryDiagnostics(
                zero_variance=final_fold_estimator.diagnostics.zero_variance.copy(),
                nuisance_aliased=(
                    final_fold_estimator.diagnostics.nuisance_aliased.copy()
                ),
                projection_diagonal=(
                    final_fold_estimator.diagnostics.projection_diagonal.copy()
                ),
                transferable=final_fold_estimator.diagnostics.transferable.copy(),
                effective_rank=final_fold_estimator.diagnostics.effective_rank,
                condition_number=final_fold_estimator.diagnostics.condition_number,
                held_unmodeled_target_burden=held_unmodeled,
            )
        )

    if int(np.sum(held_cell_counts)) != response_matrix.shape[0]:
        raise ModelingError("frozen folds must hold every cell exactly once")
    mean_gains = np.average(
        fold_mean_gains,
        axis=0,
        weights=held_cell_counts,
    )
    selected_lambda = _largest_frozen_lambda_within_relative_gain(mean_gains)
    selected_lambda_index = penalties.index(selected_lambda)
    fold_selected_effects = np.stack(
        [candidates[selected_lambda_index] for candidates in fold_effect_candidates],
        axis=0,
    )
    del fold_effect_candidates
    final_support = validated_audit.audit.final_support
    final_geometry = _prepare_sparse_fit_geometry(
        response_matrix,
        targets,
        nuisance_matrix,
        labels,
        support=final_support.supported,
        rank=requested_rank,
        design_geometry=validated_audit.final_geometry,
    )
    estimator, _ = _fit_prepared_sparse_geometry(
        final_geometry,
        ridge_lambda=selected_lambda,
    )
    return ProductionSelectionResult(
        selected_lambda=selected_lambda,
        estimator=estimator,
        lambda_grid=np.asarray(penalties, dtype=np.float64),
        held_channel_folds=held_folds,
        fold_mean_gains=fold_mean_gains,
        mean_held_cell_gains=mean_gains,
        held_cell_counts=held_cell_counts,
        fold_support=tuple(fold_support),
        fold_solver_diagnostics=tuple(fold_solver_diagnostics),
        fold_geometry_diagnostics=tuple(fold_geometry_diagnostics),
        fold_selected_effects=fold_selected_effects,
        final_support=final_support,
    )


def exact_unsupported_fallback(
    incumbent: np.ndarray,
    source: np.ndarray,
    supported: np.ndarray,
) -> np.ndarray:
    """Use source values only on supported entries.

    The output starts as a direct incumbent copy, and unsupported entries are
    never assigned. Their bytes are therefore exactly the incumbent bytes.
    """

    incumbent_array = np.asarray(incumbent)
    source_array = np.asarray(source)
    if incumbent_array.shape != source_array.shape:
        raise ModelingError("incumbent and source must have identical shapes")
    if incumbent_array.ndim != 2:
        raise ModelingError("incumbent and source must be two-dimensional")
    mask = np.asarray(supported, dtype=bool)
    if mask.shape == (incumbent_array.shape[0],):
        mask = mask[:, None]
    try:
        mask = np.broadcast_to(mask, incumbent_array.shape)
    except ValueError as error:
        raise ModelingError("supported mask is not broadcastable to effects") from error
    output = incumbent_array.copy()
    output[mask] = source_array[mask]
    return output
