#!/usr/bin/env python3
"""Strict, versioned artifact contracts for experiment 110.

The writers operate only on caller-supplied arrays. They do not know how to
locate STATE outputs, source truth, held truth, or challenge data.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from arbitration import (
    CONTROL_BANKS,
    DEVELOPMENT_CONTEXTS,
    FEATURE_NAMES,
    FIXED_GRID,
    NULL_REPETITIONS,
    FullSelectionNull,
    LocoSelection,
    apply_row_cosine_distance_cap,
    apply_row_cap_scales,
    cosine_distance_rows,
    fixed_grid_array,
    make_development_context,
    sha256_array,
    sparse_state_replacement,
)


SCHEMA_VERSION = 1
EXPERIMENT = "110-state-disagreement-arbitration"
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEVELOPMENT_CONTEXT_STAGE = "development-context-freeze"
DEVELOPMENT_SELECTION_STAGE = "development-selection-freeze"
HELD_PREDICTION_STAGE = "held-prediction-freeze"
ARTIFACT_NAME = "artifact.npz"
MANIFEST_NAME = "manifest.json"
CONFIG_WIDTH = 4
SHA256_LENGTH = 64


class ArtifactContractError(ValueError):
    """Raised when an artifact fails a schema or provenance contract."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def _is_sha256(value: object) -> bool:
    if not isinstance(value, str) or len(value) != SHA256_LENGTH:
        return False
    return all(character in "0123456789abcdef" for character in value)


def _require_keys(
    mapping: Mapping[str, object],
    expected: set[str],
    *,
    label: str,
) -> None:
    observed = set(mapping)
    if observed != expected:
        missing = sorted(expected - observed)
        extra = sorted(observed - expected)
        raise ArtifactContractError(
            f"{label} keys differ; missing={missing}, extra={extra}"
        )


def _array(
    arrays: Mapping[str, np.ndarray],
    name: str,
    *,
    ndim: int,
    kind: str,
    finite: bool = True,
) -> np.ndarray:
    value = np.asarray(arrays[name])
    if value.ndim != ndim:
        raise ArtifactContractError(f"{name} must have rank {ndim}")
    if value.dtype.kind not in kind:
        raise ArtifactContractError(
            f"{name} dtype kind {value.dtype.kind!r} not in {kind!r}"
        )
    if finite and value.dtype.kind == "f" and not np.isfinite(value).all():
        raise ArtifactContractError(f"{name} contains non-finite values")
    return value


def _grid_index(config: np.ndarray) -> int:
    values = np.asarray(config, dtype=np.float64)
    if values.shape != (CONFIG_WIDTH,):
        raise ArtifactContractError("candidate configuration must have four values")
    matches = np.flatnonzero(
        np.all(np.isclose(fixed_grid_array(), values, rtol=0, atol=0), axis=1)
    )
    if len(matches) != 1:
        raise ArtifactContractError(
            f"configuration is outside the frozen grid: {values}"
        )
    return int(matches[0])


def selection_arrays(
    selection: LocoSelection,
    null: FullSelectionNull,
) -> dict[str, np.ndarray]:
    """Convert numerical selection outputs to the strict NPZ schema."""

    model = selection.final_model
    return {
        "feature_names": np.asarray(FEATURE_NAMES, dtype="U"),
        "development_contexts": np.asarray(DEVELOPMENT_CONTEXTS, dtype="U"),
        "grid": fixed_grid_array(),
        "rotation_scores": np.asarray(
            selection.grid_scores.rotation_scores, dtype=np.float64
        ),
        "rotation_winner_indices": np.asarray(
            selection.grid_scores.rotation_winner_indices, dtype=np.int64
        ),
        "rotation_unique_winner": np.asarray(
            selection.grid_scores.rotation_unique_winner, dtype=bool
        ),
        "selected_config_index": np.asarray([selection.selected_index], dtype=np.int64),
        "selected_config": selection.selected_config.as_array(),
        "scaler_mean": np.asarray(model.standardizer.mean, dtype=np.float64),
        "scaler_scale": np.asarray(model.standardizer.scale, dtype=np.float64),
        "coefficient": np.asarray(model.coefficient, dtype=np.float64),
        "intercept": np.asarray([model.intercept], dtype=np.float64),
        "null_selected_config_indices": np.asarray(
            null.selected_config_indices, dtype=np.int64
        ),
        "null_selected_configs": np.asarray(null.selected_configs, dtype=np.float64),
        "null_selected_macro_scores": np.asarray(
            null.selected_macro_scores, dtype=np.float64
        ),
        "null_rotation_scores": np.asarray(null.rotation_scores, dtype=np.float64),
        "null_rotation_winner_indices": np.asarray(
            null.rotation_winner_indices, dtype=np.int64
        ),
        "null_rotation_unique_winner": np.asarray(
            null.rotation_unique_winner, dtype=bool
        ),
        "null_common_rotation_winner": np.asarray(
            null.common_rotation_winner, dtype=bool
        ),
        "null_label_sha256": np.asarray(null.label_sha256, dtype="U64"),
        "null_scaler_mean": np.asarray(null.scaler_mean, dtype=np.float64),
        "null_scaler_scale": np.asarray(null.scaler_scale, dtype=np.float64),
        "null_coefficient": np.asarray(null.coefficient, dtype=np.float64),
        "null_intercept": np.asarray(null.intercept, dtype=np.float64),
    }


DEVELOPMENT_CONTEXT_KEYS = {
    "context",
    "feature_names",
    "genes",
    "targets",
    "target_blocks_by_row",
    "incumbent",
    "state_banks",
    "state_mean",
    "a14_cross_donor_agreement",
    "control_log_cpm",
    "reference_lfc",
    "reference_significant",
    "own_target_coordinate_mask",
    "state_supported",
    "forbidden_gene_mask",
    "geometry_mask",
    "direction_mask",
    "mutable_mask",
    "pair_rows",
    "pair_cols",
    "pair_target_blocks",
    "features",
    "labels",
}


def validate_development_context_arrays(
    arrays: Mapping[str, np.ndarray],
) -> None:
    """Recompute source-development features and labels from frozen arrays."""

    _require_keys(arrays, DEVELOPMENT_CONTEXT_KEYS, label=DEVELOPMENT_CONTEXT_STAGE)
    context_array = _array(arrays, "context", ndim=1, kind="US")
    feature_names = _array(arrays, "feature_names", ndim=1, kind="US")
    genes = _array(arrays, "genes", ndim=1, kind="US")
    targets = _array(arrays, "targets", ndim=1, kind="US")
    target_blocks = _array(arrays, "target_blocks_by_row", ndim=1, kind="US")
    incumbent = _array(arrays, "incumbent", ndim=2, kind="f")
    state_banks = _array(arrays, "state_banks", ndim=3, kind="f")
    state_mean = _array(arrays, "state_mean", ndim=2, kind="f")
    cross_donor = _array(arrays, "a14_cross_donor_agreement", ndim=2, kind="f")
    control_log_cpm = _array(arrays, "control_log_cpm", ndim=1, kind="f")
    reference = _array(arrays, "reference_lfc", ndim=2, kind="f", finite=False)
    significant = _array(arrays, "reference_significant", ndim=2, kind="b")
    own_target = _array(arrays, "own_target_coordinate_mask", ndim=2, kind="b")
    supported = _array(arrays, "state_supported", ndim=1, kind="b")
    forbidden = _array(arrays, "forbidden_gene_mask", ndim=1, kind="b")
    geometry = _array(arrays, "geometry_mask", ndim=1, kind="b")
    direction = _array(arrays, "direction_mask", ndim=2, kind="b")
    mutable = _array(arrays, "mutable_mask", ndim=2, kind="b")
    rows = _array(arrays, "pair_rows", ndim=1, kind="iu")
    cols = _array(arrays, "pair_cols", ndim=1, kind="iu")
    pair_blocks = _array(arrays, "pair_target_blocks", ndim=1, kind="US")
    features = _array(arrays, "features", ndim=2, kind="f")
    labels = _array(arrays, "labels", ndim=1, kind="b")

    if context_array.shape != (1,) or str(context_array[0]) not in DEVELOPMENT_CONTEXTS:
        raise ArtifactContractError("development context identity is invalid")
    context = str(context_array[0])
    if not np.array_equal(feature_names.astype(str), np.asarray(FEATURE_NAMES)):
        raise ArtifactContractError("development feature names changed")
    if incumbent.shape != (len(targets), len(genes)):
        raise ArtifactContractError("development axes do not match incumbent matrix")
    if target_blocks.shape != targets.shape or not np.array_equal(
        target_blocks.astype(str), targets.astype(str)
    ):
        raise ArtifactContractError(
            "development target blocks must be the canonical target axis"
        )
    if state_banks.shape != (CONTROL_BANKS, *incumbent.shape):
        raise ArtifactContractError("development artifact lacks four STATE banks")
    if state_mean.shape != incumbent.shape or not np.allclose(
        state_mean, state_banks.mean(axis=0), rtol=0, atol=1e-12
    ):
        raise ArtifactContractError("development STATE mean differs from bank mean")
    if cross_donor.shape != incumbent.shape:
        raise ArtifactContractError("development cross-donor matrix is misaligned")
    if control_log_cpm.shape != (len(genes),):
        raise ArtifactContractError("development control vector is misaligned")
    for name, value in (
        ("reference_lfc", reference),
        ("reference_significant", significant),
        ("own_target_coordinate_mask", own_target),
        ("direction_mask", direction),
        ("mutable_mask", mutable),
    ):
        if value.shape != incumbent.shape:
            raise ArtifactContractError(f"{name} is misaligned")
    for name, value in (
        ("state_supported", supported),
        ("forbidden_gene_mask", forbidden),
        ("geometry_mask", geometry),
    ):
        if value.shape != (len(genes),):
            raise ArtifactContractError(f"{name} is misaligned")
    if not np.array_equal(geometry, ~forbidden):
        raise ArtifactContractError("development geometry is not exactly non-forbidden")

    rebuilt = make_development_context(
        name=context,
        incumbent_lfc=incumbent,
        state_bank_lfc=state_banks,
        a14_cross_donor_agreement=cross_donor,
        control_log_cpm=control_log_cpm,
        reference_lfc=reference,
        reference_significant=significant,
        state_supported_gene_mask=supported,
        forbidden_gene_mask=forbidden,
        own_target_coordinate_mask=own_target,
        target_blocks_by_row=target_blocks,
    )
    expected = {
        "state_mean": rebuilt.state,
        "direction_mask": rebuilt.direction_mask,
        "mutable_mask": rebuilt.mutable_mask,
        "pair_rows": rebuilt.pair_rows,
        "pair_cols": rebuilt.pair_cols,
        "pair_target_blocks": rebuilt.target_blocks,
        "features": rebuilt.features,
        "labels": rebuilt.labels,
    }
    observed = {
        "state_mean": state_mean,
        "direction_mask": direction,
        "mutable_mask": mutable,
        "pair_rows": rows,
        "pair_cols": cols,
        "pair_target_blocks": pair_blocks,
        "features": features,
        "labels": labels,
    }
    for name in expected:
        if not np.array_equal(observed[name], expected[name]):
            raise ArtifactContractError(
                f"development {name} differs from pure reconstruction"
            )


SELECTION_KEYS = {
    "feature_names",
    "development_contexts",
    "grid",
    "rotation_scores",
    "rotation_winner_indices",
    "rotation_unique_winner",
    "selected_config_index",
    "selected_config",
    "scaler_mean",
    "scaler_scale",
    "coefficient",
    "intercept",
    "null_selected_config_indices",
    "null_selected_configs",
    "null_selected_macro_scores",
    "null_rotation_scores",
    "null_rotation_winner_indices",
    "null_rotation_unique_winner",
    "null_common_rotation_winner",
    "null_label_sha256",
    "null_scaler_mean",
    "null_scaler_scale",
    "null_coefficient",
    "null_intercept",
}


def validate_selection_arrays(arrays: Mapping[str, np.ndarray]) -> None:
    """Validate the complete development-selection and permutation-null freeze."""

    _require_keys(arrays, SELECTION_KEYS, label=DEVELOPMENT_SELECTION_STAGE)
    feature_names = _array(arrays, "feature_names", ndim=1, kind="US")
    contexts = _array(arrays, "development_contexts", ndim=1, kind="US")
    grid = _array(arrays, "grid", ndim=2, kind="f")
    rotation_scores = _array(arrays, "rotation_scores", ndim=2, kind="f")
    rotation_winners = _array(arrays, "rotation_winner_indices", ndim=1, kind="iu")
    rotation_unique = _array(arrays, "rotation_unique_winner", ndim=1, kind="b")
    selected_index = _array(arrays, "selected_config_index", ndim=1, kind="iu")
    selected_config = _array(arrays, "selected_config", ndim=1, kind="f")
    scaler_mean = _array(arrays, "scaler_mean", ndim=1, kind="f")
    scaler_scale = _array(arrays, "scaler_scale", ndim=1, kind="f")
    coefficient = _array(arrays, "coefficient", ndim=1, kind="f")
    intercept = _array(arrays, "intercept", ndim=1, kind="f")

    expected_grid = fixed_grid_array()
    if not np.array_equal(feature_names.astype(str), np.asarray(FEATURE_NAMES)):
        raise ArtifactContractError("selector feature names changed")
    if not np.array_equal(contexts.astype(str), np.asarray(DEVELOPMENT_CONTEXTS)):
        raise ArtifactContractError("development context order changed")
    if not np.array_equal(grid, expected_grid):
        raise ArtifactContractError("candidate grid changed")
    expected_rotations = len(DEVELOPMENT_CONTEXTS)
    if rotation_scores.shape != (expected_rotations, len(FIXED_GRID)):
        raise ArtifactContractError("rotation score matrix has the wrong shape")
    if rotation_winners.shape != (expected_rotations,) or rotation_unique.shape != (
        expected_rotations,
    ):
        raise ArtifactContractError("rotation winner arrays have the wrong shape")
    if selected_index.shape != (1,) or not 0 <= selected_index[0] < len(FIXED_GRID):
        raise ArtifactContractError("selected grid index is invalid")
    if not np.all(rotation_unique):
        raise ArtifactContractError("a development rotation lacks a unique winner")
    if not np.all(rotation_winners == selected_index[0]):
        raise ArtifactContractError("development rotations do not share one winner")
    for rotation, claimed in zip(rotation_scores, rotation_winners, strict=True):
        maximum = float(rotation.max())
        actual = np.flatnonzero(np.isclose(rotation, maximum, rtol=0, atol=1e-12))
        if len(actual) != 1 or int(actual[0]) != int(claimed):
            raise ArtifactContractError(
                "rotation winner is not the unique grid-score maximum"
            )
    if not np.array_equal(selected_config, grid[selected_index[0]]):
        raise ArtifactContractError("selected configuration differs from grid row")
    _grid_index(selected_config)
    for name, value in (
        ("scaler_mean", scaler_mean),
        ("scaler_scale", scaler_scale),
        ("coefficient", coefficient),
    ):
        if value.shape != (len(FEATURE_NAMES),):
            raise ArtifactContractError(f"{name} has the wrong feature width")
    if np.any(scaler_scale <= 0) or intercept.shape != (1,):
        raise ArtifactContractError("invalid selected-model normalization or intercept")

    null_indices = _array(arrays, "null_selected_config_indices", ndim=1, kind="iu")
    null_configs = _array(arrays, "null_selected_configs", ndim=2, kind="f")
    null_scores = _array(arrays, "null_selected_macro_scores", ndim=1, kind="f")
    null_rotation_scores = _array(arrays, "null_rotation_scores", ndim=3, kind="f")
    null_winners = _array(arrays, "null_rotation_winner_indices", ndim=2, kind="iu")
    null_unique = _array(arrays, "null_rotation_unique_winner", ndim=2, kind="b")
    null_common = _array(arrays, "null_common_rotation_winner", ndim=1, kind="b")
    null_hashes = _array(arrays, "null_label_sha256", ndim=2, kind="US")
    null_mean = _array(arrays, "null_scaler_mean", ndim=2, kind="f")
    null_scale = _array(arrays, "null_scaler_scale", ndim=2, kind="f")
    null_coef = _array(arrays, "null_coefficient", ndim=2, kind="f")
    null_intercept = _array(arrays, "null_intercept", ndim=1, kind="f")
    shapes = {
        "null_selected_config_indices": (NULL_REPETITIONS,),
        "null_selected_configs": (NULL_REPETITIONS, CONFIG_WIDTH),
        "null_selected_macro_scores": (NULL_REPETITIONS,),
        "null_rotation_scores": (
            NULL_REPETITIONS,
            expected_rotations,
            len(FIXED_GRID),
        ),
        "null_rotation_winner_indices": (
            NULL_REPETITIONS,
            expected_rotations,
        ),
        "null_rotation_unique_winner": (
            NULL_REPETITIONS,
            expected_rotations,
        ),
        "null_common_rotation_winner": (NULL_REPETITIONS,),
        "null_label_sha256": (NULL_REPETITIONS, expected_rotations),
        "null_scaler_mean": (NULL_REPETITIONS, len(FEATURE_NAMES)),
        "null_scaler_scale": (NULL_REPETITIONS, len(FEATURE_NAMES)),
        "null_coefficient": (NULL_REPETITIONS, len(FEATURE_NAMES)),
        "null_intercept": (NULL_REPETITIONS,),
    }
    observed = {
        "null_selected_config_indices": null_indices.shape,
        "null_selected_configs": null_configs.shape,
        "null_selected_macro_scores": null_scores.shape,
        "null_rotation_scores": null_rotation_scores.shape,
        "null_rotation_winner_indices": null_winners.shape,
        "null_rotation_unique_winner": null_unique.shape,
        "null_common_rotation_winner": null_common.shape,
        "null_label_sha256": null_hashes.shape,
        "null_scaler_mean": null_mean.shape,
        "null_scaler_scale": null_scale.shape,
        "null_coefficient": null_coef.shape,
        "null_intercept": null_intercept.shape,
    }
    if observed != shapes:
        raise ArtifactContractError(
            f"null artifact shapes differ: observed={observed}, expected={shapes}"
        )
    if np.any(null_indices < 0) or np.any(null_indices >= len(FIXED_GRID)):
        raise ArtifactContractError("null selected grid index is invalid")
    if not np.array_equal(null_configs, grid[null_indices]):
        raise ArtifactContractError("null configurations differ from indexed grid rows")
    for repetition in range(NULL_REPETITIONS):
        calculated_winners = []
        calculated_unique = []
        for rotation in null_rotation_scores[repetition]:
            maximum = float(rotation.max())
            actual = np.flatnonzero(np.isclose(rotation, maximum, rtol=0, atol=1e-12))
            calculated_winners.append(int(actual[0]))
            calculated_unique.append(len(actual) == 1)
        if not np.array_equal(
            null_winners[repetition],
            np.asarray(calculated_winners, dtype=np.int64),
        ):
            raise ArtifactContractError("null rotation winner differs from scores")
        if not np.array_equal(
            null_unique[repetition],
            np.asarray(calculated_unique, dtype=bool),
        ):
            raise ArtifactContractError("null rotation uniqueness differs from scores")
        expected_common = bool(
            all(calculated_unique) and len(set(calculated_winners)) == 1
        )
        if bool(null_common[repetition]) != expected_common:
            raise ArtifactContractError("null common-winner flag differs from scores")
        macro = null_rotation_scores[repetition].mean(axis=0)
        maximum = float(macro.max())
        expected_index = int(
            np.flatnonzero(np.isclose(macro, maximum, rtol=0, atol=1e-12))[0]
        )
        if int(null_indices[repetition]) != expected_index:
            raise ArtifactContractError("null selected index differs from macro scores")
        if float(null_scores[repetition]) != float(macro[expected_index]):
            raise ArtifactContractError("null selected score differs from macro scores")
    if np.any(null_scale <= 0):
        raise ArtifactContractError("null standardization contains nonpositive scales")
    if not all(_is_sha256(value) for value in null_hashes.ravel().astype(str)):
        raise ArtifactContractError("null label hashes are malformed")


HELD_KEYS = {
    "feature_names",
    "genes",
    "targets",
    "incumbent",
    "state_banks",
    "state_mean",
    "state_supported",
    "forbidden_gene_mask",
    "geometry_mask",
    "pair_rows",
    "pair_cols",
    "features",
    "selected_probability",
    "replacement_mask",
    "row_cap_scale",
    "candidate",
    "selected_config",
    "null_selected_configs",
    "null_pair_selected",
    "null_row_cap_scale",
}


def validate_held_prediction_arrays(arrays: Mapping[str, np.ndarray]) -> None:
    """Validate a truth-free observed/null prediction freeze."""

    _require_keys(arrays, HELD_KEYS, label=HELD_PREDICTION_STAGE)
    feature_names = _array(arrays, "feature_names", ndim=1, kind="US")
    genes = _array(arrays, "genes", ndim=1, kind="US")
    targets = _array(arrays, "targets", ndim=1, kind="US")
    incumbent = _array(arrays, "incumbent", ndim=2, kind="f")
    state_banks = _array(arrays, "state_banks", ndim=3, kind="f")
    state_mean = _array(arrays, "state_mean", ndim=2, kind="f")
    supported = _array(arrays, "state_supported", ndim=1, kind="b")
    forbidden = _array(arrays, "forbidden_gene_mask", ndim=1, kind="b")
    geometry = _array(arrays, "geometry_mask", ndim=1, kind="b")
    rows = _array(arrays, "pair_rows", ndim=1, kind="iu")
    cols = _array(arrays, "pair_cols", ndim=1, kind="iu")
    features = _array(arrays, "features", ndim=2, kind="f")
    probabilities = _array(arrays, "selected_probability", ndim=1, kind="f")
    replacement = _array(arrays, "replacement_mask", ndim=2, kind="b")
    row_scales = _array(arrays, "row_cap_scale", ndim=1, kind="f")
    candidate = _array(arrays, "candidate", ndim=2, kind="f")
    selected_config = _array(arrays, "selected_config", ndim=1, kind="f")
    null_configs = _array(arrays, "null_selected_configs", ndim=2, kind="f")
    null_selected = _array(arrays, "null_pair_selected", ndim=2, kind="b")
    null_scales = _array(arrays, "null_row_cap_scale", ndim=2, kind="f")

    if not np.array_equal(feature_names.astype(str), np.asarray(FEATURE_NAMES)):
        raise ArtifactContractError("held feature names changed")
    if len(genes) != incumbent.shape[1] or len(targets) != incumbent.shape[0]:
        raise ArtifactContractError("held axes do not match incumbent matrix")
    matrix_shape = incumbent.shape
    if state_banks.shape != (CONTROL_BANKS, *matrix_shape):
        raise ArtifactContractError("held artifact lacks four aligned STATE banks")
    if state_mean.shape != matrix_shape or candidate.shape != matrix_shape:
        raise ArtifactContractError("held effect matrices are misaligned")
    if not np.allclose(state_mean, state_banks.mean(axis=0), rtol=0, atol=1e-12):
        raise ArtifactContractError("STATE aggregate differs from the four-bank mean")
    gene_shape = (matrix_shape[1],)
    if (
        supported.shape != gene_shape
        or forbidden.shape != gene_shape
        or geometry.shape != gene_shape
    ):
        raise ArtifactContractError("held gene masks are misaligned")
    if np.any(geometry & forbidden):
        raise ArtifactContractError("forbidden genes entered PDS geometry")
    if np.any(incumbent[:, forbidden] != 0):
        raise ArtifactContractError("incumbent is nonzero on forbidden coordinates")
    pair_shape = (len(rows),)
    if (
        cols.shape != pair_shape
        or probabilities.shape != pair_shape
        or features.shape != (len(rows), len(FEATURE_NAMES))
    ):
        raise ArtifactContractError("held disagreement-pair arrays are misaligned")
    if (
        np.any(rows < 0)
        or np.any(rows >= matrix_shape[0])
        or np.any(cols < 0)
        or np.any(cols >= matrix_shape[1])
    ):
        raise ArtifactContractError("held pair coordinates are out of bounds")
    linear = rows * matrix_shape[1] + cols
    if len(np.unique(linear)) != len(linear):
        raise ArtifactContractError("held pair coordinates contain duplicates")
    mutable = np.broadcast_to(supported & ~forbidden, matrix_shape)
    if not np.all(mutable[rows, cols]):
        raise ArtifactContractError("held pair list contains an immutable coordinate")
    if np.any(np.sign(incumbent[rows, cols]) == np.sign(state_mean[rows, cols])):
        raise ArtifactContractError("held pair list contains a sign agreement")
    if np.any((probabilities < 0) | (probabilities > 1)):
        raise ArtifactContractError("held selector probabilities leave [0, 1]")

    config_index = _grid_index(selected_config)
    config = FIXED_GRID[config_index]
    expected_replacement = np.zeros(matrix_shape, dtype=bool)
    selected_pairs = probabilities >= config.probability_threshold
    expected_replacement[rows[selected_pairs], cols[selected_pairs]] = True
    if not np.array_equal(replacement, expected_replacement):
        raise ArtifactContractError("held replacement mask differs from thresholding")
    if replacement.shape != matrix_shape or np.any(replacement & ~mutable):
        raise ArtifactContractError("held replacement mask violates support")
    if row_scales.shape != (matrix_shape[0],) or np.any(
        (row_scales < 0) | (row_scales > 1)
    ):
        raise ArtifactContractError("held row cap scales are invalid")
    uncapped = sparse_state_replacement(
        incumbent, state_mean, replacement, config.state_blend_fraction
    )
    reconstructed = apply_row_cap_scales(incumbent, uncapped, row_scales)
    reconstructed[~mutable] = incumbent[~mutable]
    reconstructed[:, forbidden] = 0.0
    if not np.array_equal(candidate, reconstructed):
        raise ArtifactContractError("held candidate is not exactly reconstructible")
    expected_candidate, expected_scales = apply_row_cosine_distance_cap(
        incumbent,
        uncapped,
        geometry,
        config.maximum_cosine_distance,
    )
    expected_candidate[~mutable] = incumbent[~mutable]
    expected_candidate[:, forbidden] = 0.0
    if not np.array_equal(row_scales, expected_scales) or not np.array_equal(
        candidate, expected_candidate
    ):
        raise ArtifactContractError("held candidate differs from the frozen cap rule")
    distances = cosine_distance_rows(incumbent, candidate, geometry)
    if np.any(distances > config.maximum_cosine_distance + 5e-13):
        raise ArtifactContractError("held candidate exceeds its cosine-distance cap")
    if not np.array_equal(candidate[~mutable], incumbent[~mutable]):
        raise ArtifactContractError("held unsupported fallback is not exact")
    if np.any(candidate[:, forbidden] != 0):
        raise ArtifactContractError("held forbidden coordinate is not exactly zero")

    if null_configs.shape != (NULL_REPETITIONS, CONFIG_WIDTH):
        raise ArtifactContractError("held null configuration matrix has wrong shape")
    if null_selected.shape != (NULL_REPETITIONS, len(rows)):
        raise ArtifactContractError("held null pair-selection matrix has wrong shape")
    if null_scales.shape != (NULL_REPETITIONS, matrix_shape[0]):
        raise ArtifactContractError("held null cap-scale matrix has wrong shape")
    if np.any((null_scales < 0) | (null_scales > 1)):
        raise ArtifactContractError("held null cap scales leave [0, 1]")
    for repetition in range(NULL_REPETITIONS):
        null_index = _grid_index(null_configs[repetition])
        null_config = FIXED_GRID[null_index]
        null_replacement = np.zeros(matrix_shape, dtype=bool)
        chosen = null_selected[repetition]
        null_replacement[rows[chosen], cols[chosen]] = True
        null_uncapped = sparse_state_replacement(
            incumbent,
            state_mean,
            null_replacement,
            null_config.state_blend_fraction,
        )
        null_candidate = apply_row_cap_scales(
            incumbent, null_uncapped, null_scales[repetition]
        )
        null_candidate[~mutable] = incumbent[~mutable]
        null_candidate[:, forbidden] = 0.0
        null_distances = cosine_distance_rows(incumbent, null_candidate, geometry)
        if np.any(null_distances > null_config.maximum_cosine_distance + 5e-13):
            raise ArtifactContractError(
                f"held null {repetition} exceeds its cosine-distance cap"
            )
        expected_null, expected_null_scales = apply_row_cosine_distance_cap(
            incumbent,
            null_uncapped,
            geometry,
            null_config.maximum_cosine_distance,
        )
        expected_null[~mutable] = incumbent[~mutable]
        expected_null[:, forbidden] = 0.0
        if not np.array_equal(
            null_scales[repetition], expected_null_scales
        ) or not np.array_equal(null_candidate, expected_null):
            raise ArtifactContractError(
                f"held null {repetition} differs from the frozen cap rule"
            )


def _probability_from_parameters(
    features: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    coefficient: np.ndarray,
    intercept: float,
) -> np.ndarray:
    standardized = (features - mean[None, :]) / scale[None, :]
    logits = standardized @ coefficient + intercept
    probabilities = np.empty_like(logits, dtype=np.float64)
    positive = logits >= 0
    probabilities[positive] = 1.0 / (1.0 + np.exp(-logits[positive]))
    exp_logits = np.exp(logits[~positive])
    probabilities[~positive] = exp_logits / (1.0 + exp_logits)
    return probabilities


def validate_held_against_selection(
    held_arrays: Mapping[str, np.ndarray],
    selection_arrays_value: Mapping[str, np.ndarray],
) -> None:
    """Bind held observed/null masks to the frozen selector parameters."""

    validate_held_prediction_arrays(held_arrays)
    validate_selection_arrays(selection_arrays_value)
    features = np.asarray(held_arrays["features"], dtype=np.float64)
    selected_config = np.asarray(held_arrays["selected_config"], dtype=np.float64)
    if not np.array_equal(selected_config, selection_arrays_value["selected_config"]):
        raise ArtifactContractError(
            "held selected configuration differs from selection freeze"
        )
    probabilities = _probability_from_parameters(
        features,
        np.asarray(selection_arrays_value["scaler_mean"], dtype=np.float64),
        np.asarray(selection_arrays_value["scaler_scale"], dtype=np.float64),
        np.asarray(selection_arrays_value["coefficient"], dtype=np.float64),
        float(np.asarray(selection_arrays_value["intercept"])[0]),
    )
    if not np.allclose(
        probabilities,
        np.asarray(held_arrays["selected_probability"], dtype=np.float64),
        rtol=0,
        atol=1e-12,
    ):
        raise ArtifactContractError(
            "held probabilities differ from the frozen selected model"
        )

    null_configs = np.asarray(
        selection_arrays_value["null_selected_configs"], dtype=np.float64
    )
    if not np.array_equal(null_configs, held_arrays["null_selected_configs"]):
        raise ArtifactContractError(
            "held null configurations differ from selection freeze"
        )
    null_selected = np.asarray(held_arrays["null_pair_selected"], dtype=bool)
    for repetition in range(NULL_REPETITIONS):
        null_probability = _probability_from_parameters(
            features,
            np.asarray(selection_arrays_value["null_scaler_mean"])[repetition],
            np.asarray(selection_arrays_value["null_scaler_scale"])[repetition],
            np.asarray(selection_arrays_value["null_coefficient"])[repetition],
            float(np.asarray(selection_arrays_value["null_intercept"])[repetition]),
        )
        threshold = float(null_configs[repetition, 1])
        if not np.array_equal(null_selected[repetition], null_probability >= threshold):
            raise ArtifactContractError(
                f"held null {repetition} mask differs from frozen selector"
            )


COMMON_MANIFEST_KEYS = {
    "schema_version",
    "experiment",
    "stage",
    "artifact",
    "artifact_sha256",
    "protocol",
    "protocol_sha256",
    "implementation_sha256",
    "input_artifacts",
    "held_truth_read",
}
DEVELOPMENT_CONTEXT_MANIFEST_KEYS = COMMON_MANIFEST_KEYS | {
    "context",
    "contains_source_reference_truth",
    "source_truth_scope",
    "state_training_exclusion",
    "in_context_trained_state_forbidden",
    "control_banks",
    "gene_axis_sha256",
    "target_axis_sha256",
    "state_checkpoint_manifest",
    "state_checkpoint_manifest_sha256",
}
SELECTION_MANIFEST_KEYS = COMMON_MANIFEST_KEYS | {
    "development_contexts",
    "feature_names",
    "grid_sha256",
    "null_repetitions",
    "checkpoint_policy",
}
HELD_MANIFEST_KEYS = COMMON_MANIFEST_KEYS | {
    "surface",
    "contains_reference_truth",
    "sealed_before_held_truth",
    "control_banks",
    "gene_axis_sha256",
    "target_axis_sha256",
    "selection_manifest",
    "selection_manifest_sha256",
}


def validate_manifest(
    manifest: Mapping[str, object],
    *,
    stage: str,
    artifact_sha256: str,
) -> None:
    """Validate exact manifest keys and truth-sealing declarations."""

    expected_keys = {
        DEVELOPMENT_CONTEXT_STAGE: DEVELOPMENT_CONTEXT_MANIFEST_KEYS,
        DEVELOPMENT_SELECTION_STAGE: SELECTION_MANIFEST_KEYS,
        HELD_PREDICTION_STAGE: HELD_MANIFEST_KEYS,
    }.get(stage)
    if expected_keys is None:
        raise ArtifactContractError(f"unknown artifact stage: {stage}")
    _require_keys(manifest, expected_keys, label=f"{stage} manifest")
    expected = {
        "schema_version": SCHEMA_VERSION,
        "experiment": EXPERIMENT,
        "stage": stage,
        "artifact": ARTIFACT_NAME,
        "artifact_sha256": artifact_sha256,
        "protocol": "preregistered_protocol.json",
        "held_truth_read": False,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ArtifactContractError(
                f"manifest {key}={manifest.get(key)!r}, expected {value!r}"
            )
    if not _is_sha256(manifest.get("protocol_sha256")):
        raise ArtifactContractError("protocol_sha256 is malformed")
    if manifest.get("protocol_sha256") != sha256_file(
        HERE / "preregistered_protocol.json"
    ):
        raise ArtifactContractError("frozen preregistered protocol changed")
    implementation = manifest.get("implementation_sha256")
    expected_implementation = {
        "arbitration.py": sha256_file(HERE / "arbitration.py"),
        "artifact_contract.py": sha256_file(HERE / "artifact_contract.py"),
    }
    if implementation != expected_implementation:
        raise ArtifactContractError("implementation source hashes changed")
    inputs = manifest.get("input_artifacts")
    if not isinstance(inputs, list) or not inputs:
        raise ArtifactContractError("input_artifacts must be a nonempty list")
    for item in inputs:
        if (
            not isinstance(item, dict)
            or set(item) != {"path", "sha256", "role"}
            or not isinstance(item["path"], str)
            or not isinstance(item["role"], str)
            or not _is_sha256(item["sha256"])
        ):
            raise ArtifactContractError("input_artifacts entry is malformed")

    if stage == DEVELOPMENT_CONTEXT_STAGE:
        context = manifest.get("context")
        if context not in DEVELOPMENT_CONTEXTS:
            raise ArtifactContractError("development manifest context is invalid")
        expected_development = {
            "contains_source_reference_truth": True,
            "source_truth_scope": "public-development-only",
            "state_training_exclusion": context,
            "in_context_trained_state_forbidden": True,
            "control_banks": CONTROL_BANKS,
        }
        for key, value in expected_development.items():
            if manifest.get(key) != value:
                raise ArtifactContractError(f"development manifest {key} changed")
        for key in (
            "gene_axis_sha256",
            "target_axis_sha256",
            "state_checkpoint_manifest_sha256",
        ):
            if not _is_sha256(manifest.get(key)):
                raise ArtifactContractError(f"{key} is malformed")
        if not isinstance(manifest.get("state_checkpoint_manifest"), str):
            raise ArtifactContractError("STATE checkpoint manifest path is missing")
    elif stage == DEVELOPMENT_SELECTION_STAGE:
        if manifest.get("development_contexts") != list(DEVELOPMENT_CONTEXTS):
            raise ArtifactContractError("selection manifest context order changed")
        if manifest.get("feature_names") != list(FEATURE_NAMES):
            raise ArtifactContractError("selection manifest feature names changed")
        if manifest.get("grid_sha256") != sha256_array(fixed_grid_array()):
            raise ArtifactContractError("selection manifest grid hash changed")
        if manifest.get("null_repetitions") != NULL_REPETITIONS:
            raise ArtifactContractError("selection manifest null count changed")
        checkpoint = manifest.get("checkpoint_policy")
        expected_checkpoint = {
            "name": "last.ckpt",
            "global_step": 20_000,
            "best_ckpt_forbidden": True,
            "strict_context_holdout_required": True,
        }
        if checkpoint != expected_checkpoint:
            raise ArtifactContractError("selection checkpoint policy changed")
    else:
        if manifest.get("contains_reference_truth") is not False:
            raise ArtifactContractError("held prediction manifest permits truth")
        if manifest.get("sealed_before_held_truth") is not True:
            raise ArtifactContractError("held prediction artifact was not sealed")
        if manifest.get("control_banks") != CONTROL_BANKS:
            raise ArtifactContractError("held prediction control-bank count changed")
        for key in (
            "gene_axis_sha256",
            "target_axis_sha256",
            "selection_manifest_sha256",
        ):
            if not _is_sha256(manifest.get(key)):
                raise ArtifactContractError(f"{key} is malformed")
        if not isinstance(manifest.get("surface"), str):
            raise ArtifactContractError("held surface is missing")
        if not isinstance(manifest.get("selection_manifest"), str):
            raise ArtifactContractError("selection manifest path is missing")


def verify_manifest_inputs(
    manifest: Mapping[str, object],
    provenance_root: Path,
) -> None:
    """Verify every bound input path and hash without following paths outside root."""

    root = provenance_root.resolve()
    inputs = manifest["input_artifacts"]
    assert isinstance(inputs, list)
    resolved: dict[str, str] = {}
    for item in inputs:
        assert isinstance(item, dict)
        relative = Path(str(item["path"]))
        if relative.is_absolute():
            raise ArtifactContractError("input artifact paths must be root-relative")
        path = (root / relative).resolve()
        if root not in path.parents:
            raise ArtifactContractError(f"input artifact escapes root: {relative}")
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256_file(path)
        if observed != item["sha256"]:
            raise ArtifactContractError(f"input artifact hash changed: {relative}")
        resolved[relative.as_posix()] = observed
    if manifest["stage"] == HELD_PREDICTION_STAGE:
        selection_relative = Path(str(manifest["selection_manifest"]))
        selection_key = selection_relative.as_posix()
        if resolved.get(selection_key) != manifest["selection_manifest_sha256"]:
            raise ArtifactContractError(
                "held selection manifest is not identically input-bound"
            )
    elif manifest["stage"] == DEVELOPMENT_CONTEXT_STAGE:
        checkpoint_relative = Path(str(manifest["state_checkpoint_manifest"]))
        checkpoint_key = checkpoint_relative.as_posix()
        if resolved.get(checkpoint_key) != manifest["state_checkpoint_manifest_sha256"]:
            raise ArtifactContractError(
                "development STATE checkpoint manifest is not identically input-bound"
            )
    elif manifest["stage"] == DEVELOPMENT_SELECTION_STAGE:
        expected_roles = {
            f"development-context-{context}" for context in DEVELOPMENT_CONTEXTS
        }
        observed_roles = {
            str(item["role"]) for item in inputs if isinstance(item, dict)
        }
        if observed_roles != expected_roles:
            raise ArtifactContractError(
                "selection inputs are not exactly the three development contexts"
            )


def validate_arrays(stage: str, arrays: Mapping[str, np.ndarray]) -> None:
    if stage == DEVELOPMENT_CONTEXT_STAGE:
        validate_development_context_arrays(arrays)
    elif stage == DEVELOPMENT_SELECTION_STAGE:
        validate_selection_arrays(arrays)
    elif stage == HELD_PREDICTION_STAGE:
        validate_held_prediction_arrays(arrays)
    else:
        raise ArtifactContractError(f"unknown artifact stage: {stage}")


def write_artifact(
    output_dir: Path,
    *,
    stage: str,
    arrays: Mapping[str, np.ndarray],
    manifest: Mapping[str, object],
    provenance_root: Path = ROOT,
) -> Path:
    """Atomically write an exact-schema artifact and refuse overwrite."""

    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    validate_arrays(stage, arrays)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        artifact_path = staging / ARTIFACT_NAME
        with artifact_path.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
        artifact_hash = sha256_file(artifact_path)
        final_manifest = dict(manifest)
        final_manifest["artifact_sha256"] = artifact_hash
        validate_manifest(
            final_manifest,
            stage=stage,
            artifact_sha256=artifact_hash,
        )
        verify_manifest_inputs(final_manifest, provenance_root)
        manifest_path = staging / MANIFEST_NAME
        manifest_path.write_text(
            json.dumps(final_manifest, indent=2, sort_keys=True) + "\n"
        )
        os.replace(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output_dir / MANIFEST_NAME


def load_artifact(
    directory: Path,
    *,
    expected_stage: str,
    provenance_root: Path = ROOT,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Load and fully revalidate a frozen experiment-110 artifact."""

    directory = directory.resolve()
    artifact_path = directory / ARTIFACT_NAME
    manifest_path = directory / MANIFEST_NAME
    if not artifact_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(f"incomplete experiment-110 artifact: {directory}")
    manifest = json.loads(manifest_path.read_text())
    artifact_hash = sha256_file(artifact_path)
    validate_manifest(
        manifest,
        stage=expected_stage,
        artifact_sha256=artifact_hash,
    )
    verify_manifest_inputs(manifest, provenance_root)
    with np.load(artifact_path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    validate_arrays(expected_stage, arrays)
    return arrays, manifest


def load_selection_with_contexts(
    selection_directory: Path,
    context_directories: Mapping[str, Path],
    *,
    provenance_root: Path = ROOT,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, object],
    dict[str, tuple[dict[str, np.ndarray], dict[str, object]]],
]:
    """Load a selection only with all three schema-valid source contexts."""

    if set(context_directories) != set(DEVELOPMENT_CONTEXTS):
        raise ArtifactContractError(
            "selection validation requires all three development contexts"
        )
    selection_arrays_value, selection_manifest = load_artifact(
        selection_directory,
        expected_stage=DEVELOPMENT_SELECTION_STAGE,
        provenance_root=provenance_root,
    )
    inputs = {
        str(item["role"]): item
        for item in selection_manifest["input_artifacts"]
        if isinstance(item, dict)
    }
    loaded_contexts = {}
    for context in DEVELOPMENT_CONTEXTS:
        arrays, manifest = load_artifact(
            context_directories[context],
            expected_stage=DEVELOPMENT_CONTEXT_STAGE,
            provenance_root=provenance_root,
        )
        if manifest["context"] != context:
            raise ArtifactContractError(
                f"development artifact has wrong context: {manifest['context']}"
            )
        role = f"development-context-{context}"
        expected = inputs[role]
        expected_path = (provenance_root.resolve() / str(expected["path"])).resolve()
        actual_path = (context_directories[context].resolve() / MANIFEST_NAME).resolve()
        if expected_path != actual_path:
            raise ArtifactContractError(
                f"selection points to a different {context} manifest"
            )
        if sha256_file(actual_path) != expected["sha256"]:
            raise ArtifactContractError(
                f"selection-bound {context} manifest hash changed"
            )
        loaded_contexts[context] = (arrays, manifest)
    return selection_arrays_value, selection_manifest, loaded_contexts


def load_held_with_selection(
    held_directory: Path,
    selection_directory: Path,
    *,
    provenance_root: Path = ROOT,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, object],
    dict[str, np.ndarray],
    dict[str, object],
]:
    """Load a held freeze only when its exact selection freeze also validates."""

    selection_arrays_value, selection_manifest = load_artifact(
        selection_directory,
        expected_stage=DEVELOPMENT_SELECTION_STAGE,
        provenance_root=provenance_root,
    )
    held_arrays, held_manifest = load_artifact(
        held_directory,
        expected_stage=HELD_PREDICTION_STAGE,
        provenance_root=provenance_root,
    )
    expected_selection_path = (
        provenance_root.resolve() / str(held_manifest["selection_manifest"])
    ).resolve()
    actual_selection_path = (selection_directory.resolve() / MANIFEST_NAME).resolve()
    if expected_selection_path != actual_selection_path:
        raise ArtifactContractError(
            "held artifact points to a different selection manifest"
        )
    if sha256_file(actual_selection_path) != held_manifest["selection_manifest_sha256"]:
        raise ArtifactContractError("held artifact selection manifest hash changed")
    validate_held_against_selection(held_arrays, selection_arrays_value)
    return held_arrays, held_manifest, selection_arrays_value, selection_manifest
