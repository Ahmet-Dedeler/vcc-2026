#!/usr/bin/env python3
"""Source-only nested validation for Experiment 111.

This module reads only the six authenticated DepMap source CSVs plus the public
challenge gene/target axes. It never opens held perturbation truth, runs
cell-eval2, generates counts, or touches production/submission paths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import shutil
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_ROOT = ROOT / "data/external/depmap26q1/extracted"
ACQUISITION_RECEIPT = HERE / "audit/acquisition_receipt.json"
SOURCE_MANIFEST = HERE / "SOURCE_MANIFEST.json"
PREREGISTERED_PROTOCOL = HERE / "preregistered_protocol.json"
EXECUTION_AUTHORIZATION = HERE / "EXECUTION_AUTHORIZATION.v1.json"
SCHEMA_RECOVERY = HERE / "ACQUISITION_SCHEMA_RECOVERY.v1.json"
IMPLEMENTATION_PROTOCOL = HERE / "SOURCE_MODEL_IMPLEMENTATION_PROTOCOL.v1.json"
SOURCE_MODEL_RECOVERY = HERE / "SOURCE_MODEL_FAILURE_RECOVERY.v1.json"
EXECUTION_FREEZE = HERE / "SOURCE_MODEL_EXECUTION_FREEZE.v1.json"
OUTPUT_ROOT = HERE / "work/source-only-v1"
GENE_AXIS = ROOT / "data/challenge/gene_names.csv"
TARGET_AXIS = ROOT / "data/challenge/pert_counts.csv"

CONFIRMATION = "I_CONFIRM_EXPERIMENT_111_SOURCE_ONLY_V1"
BASE_RESPONSE_SEED = 11120260830
BOOTSTRAP_SEED = 11120260831
NULL_SEED = 11120260832
BOOTSTRAP_REPLICATES = 5000
NULL_REPLICATES = 100
EPS = 1e-12


class ModelError(RuntimeError):
    """The source-only execution contract or model run failed."""


@dataclass(frozen=True)
class Candidate:
    response_rank: int
    control_pc_rank: int
    ridge_lambda: float
    beta: float
    angular_cap_degrees: float

    @property
    def key(self) -> str:
        return (
            f"rr{self.response_rank}_cr{self.control_pc_rank}_"
            f"l{self.ridge_lambda:g}_b{self.beta:g}_"
            f"cap{self.angular_cap_degrees:g}"
        )

    @property
    def tie_break(self) -> tuple[float, int, int, float, float]:
        return (
            self.beta,
            self.response_rank,
            self.control_pc_rank,
            -self.ridge_lambda,
            self.angular_cap_degrees,
        )


@dataclass
class SourceData:
    genes: list[str]
    contexts: list[str]
    batches: list[str]
    batch_order: list[str]
    targets: list[str]
    controls: list[str]
    effects: np.ndarray
    available: np.ndarray
    control_profiles: np.ndarray
    gram: np.ndarray
    control_mean: np.ndarray
    control_components: np.ndarray
    control_score_scale: np.ndarray
    control_scores: np.ndarray


@dataclass
class BasisState:
    train_batches: tuple[str, ...]
    basis: np.ndarray
    projections: np.ndarray


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ModelError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ModelError(f"JSON root is not an object: {path}")
    return value


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def stable_seed(label: str, base: int = BASE_RESPONSE_SEED) -> int:
    payload = f"{base}|{label}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def peak_rss_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def require_exact_path(path: Path, expected: Path, label: str) -> Path:
    resolved = path.resolve()
    expected_resolved = expected.resolve()
    if resolved != expected_resolved:
        raise ModelError(f"{label} path drifted: {resolved} != {expected_resolved}")
    return resolved


def validate_preflight() -> dict[str, Any]:
    source_root = require_exact_path(SOURCE_ROOT, SOURCE_ROOT, "source root")
    output_root = require_exact_path(OUTPUT_ROOT, OUTPUT_ROOT, "output root")
    if not source_root.is_dir() or source_root.is_symlink():
        raise ModelError("authenticated source root is missing or symlinked")
    if output_root.exists() or output_root.is_symlink():
        raise ModelError(f"refusing to overwrite source-model output: {output_root}")

    manifest = load_json(SOURCE_MANIFEST)
    protocol = load_json(PREREGISTERED_PROTOCOL)
    authorization = load_json(EXECUTION_AUTHORIZATION)
    implementation = load_json(IMPLEMENTATION_PROTOCOL)
    receipt = load_json(ACQUISITION_RECEIPT)
    freeze = load_json(EXECUTION_FREEZE)
    model_recovery = load_json(SOURCE_MODEL_RECOVERY)

    expected_bindings = {
        "source_manifest_sha256": sha256_file(SOURCE_MANIFEST),
        "preregistered_protocol_sha256": sha256_file(PREREGISTERED_PROTOCOL),
        "execution_authorization_sha256": sha256_file(EXECUTION_AUTHORIZATION),
        "acquisition_schema_recovery_sha256": sha256_file(SCHEMA_RECOVERY),
    }
    if implementation.get("bindings") != expected_bindings:
        raise ModelError("source-model implementation protocol bindings drifted")
    if implementation.get("status") != "frozen_before_source_model_outcomes":
        raise ModelError("source-model implementation protocol is not frozen")

    if not authorization.get("authorized_actions", {}).get("source_only_model_fitting"):
        raise ModelError("source-only model fitting is not authorized")
    forbidden = authorization.get("forbidden_actions", {})
    if any(value is not False for value in forbidden.values()):
        raise ModelError("a forbidden execution action was enabled")

    expected_freeze_bindings = {
        "source_manifest": {
            "path": str(SOURCE_MANIFEST.relative_to(ROOT)),
            "sha256": sha256_file(SOURCE_MANIFEST),
        },
        "preregistered_protocol": {
            "path": str(PREREGISTERED_PROTOCOL.relative_to(ROOT)),
            "sha256": sha256_file(PREREGISTERED_PROTOCOL),
        },
        "execution_authorization": {
            "path": str(EXECUTION_AUTHORIZATION.relative_to(ROOT)),
            "sha256": sha256_file(EXECUTION_AUTHORIZATION),
        },
        "acquisition_schema_recovery": {
            "path": str(SCHEMA_RECOVERY.relative_to(ROOT)),
            "sha256": sha256_file(SCHEMA_RECOVERY),
        },
        "implementation_protocol": {
            "path": str(IMPLEMENTATION_PROTOCOL.relative_to(ROOT)),
            "sha256": sha256_file(IMPLEMENTATION_PROTOCOL),
        },
        "source_model_failure_recovery": {
            "path": str(SOURCE_MODEL_RECOVERY.relative_to(ROOT)),
            "sha256": sha256_file(SOURCE_MODEL_RECOVERY),
        },
        "acquisition_receipt": {
            "path": str(ACQUISITION_RECEIPT.relative_to(ROOT)),
            "sha256": sha256_file(ACQUISITION_RECEIPT),
        },
        "source_model": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "source_model_tests": {
            "path": str((HERE / "test_source_model.py").relative_to(ROOT)),
            "sha256": sha256_file(HERE / "test_source_model.py"),
        },
    }
    if freeze.get("bindings") != expected_freeze_bindings:
        raise ModelError("source-model execution freeze bindings drifted")
    if freeze.get("status") != "authorized_frozen_source_only_run":
        raise ModelError("source-model execution freeze is not active")
    if freeze.get("confirmation_string") != CONFIRMATION:
        raise ModelError("source-model confirmation string drifted")
    if (
        model_recovery.get("status") != "authorized_mechanical_parser_recovery"
        or model_recovery.get("failed_execution_freeze_sha256")
        != "e920e0f5408ed3db0b45cf2018dd8b836919fa25d698897dd314d45e0421e22b"
        or model_recovery.get("failure", {}).get("source_model_outcomes_observed")
        is not False
        or model_recovery.get("authorized_change", {}).get("scientific_model_change")
        is not False
    ):
        raise ModelError("source-model mechanical recovery contract drifted")
    if freeze.get("authorization") != {
        "read_exact_authenticated_source_members": True,
        "fit_fixed_source_only_grid": True,
        "run_complete_selection_nulls": True,
        "write_only_source_output_root": True,
        "held_truth_access": False,
        "cell_eval2": False,
        "count_generation": False,
        "production_change": False,
        "challenge_build": False,
        "submission": False,
        "commit": False,
        "push": False,
    }:
        raise ModelError("source-model execution-freeze authorization drifted")

    receipt_members = receipt.get("members")
    if not isinstance(receipt_members, dict):
        raise ModelError("acquisition receipt has no member map")
    for relative, record in receipt_members.items():
        path = source_root / relative
        if not path.is_file() or path.is_symlink():
            raise ModelError(f"source member is missing or symlinked: {relative}")
        if path.stat().st_size != record.get("size_bytes") or sha256_file(
            path
        ) != record.get("sha256"):
            raise ModelError(f"source member binding drifted: {relative}")

    if (
        protocol.get("source_model") is None
        or protocol.get("source_validation") is None
    ):
        raise ModelError("preregistered source-model contract is missing")
    if manifest.get("challenge_contract", {}).get("expected_panel_overlap") != 0:
        raise ModelError("source target overlap is no longer zero")
    return {
        "source_manifest_sha256": expected_bindings["source_manifest_sha256"],
        "preregistered_protocol_sha256": expected_bindings[
            "preregistered_protocol_sha256"
        ],
        "execution_authorization_sha256": expected_bindings[
            "execution_authorization_sha256"
        ],
        "implementation_protocol_sha256": sha256_file(IMPLEMENTATION_PROTOCOL),
        "execution_freeze_sha256": sha256_file(EXECUTION_FREEZE),
        "acquisition_receipt_sha256": sha256_file(ACQUISITION_RECEIPT),
    }


def read_single_column(path: Path, column: str) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != [column]:
            raise ModelError(f"{path.name} must contain only {column}")
        values = [row[column] for row in reader]
    if (
        not values
        or len(values) != len(set(values))
        or any(not item for item in values)
    ):
        raise ModelError(f"{path.name} axis is empty or nonunique")
    return values


def hashed_target_split(targets: list[str]) -> dict[str, list[int]]:
    ordered = sorted(
        range(len(targets)),
        key=lambda index: (
            hashlib.sha256(targets[index].encode("utf-8")).hexdigest(),
            targets[index],
        ),
    )
    if len(ordered) != 90:
        raise ModelError("hashed target split requires exactly 90 analysis targets")
    return {
        "fit": ordered[:54],
        "tune": ordered[54:72],
        "test": ordered[72:90],
    }


def context_weights(context_indices: list[int], batches: list[str]) -> np.ndarray:
    weights = np.zeros(len(context_indices), dtype=np.float64)
    represented = sorted({batches[index] for index in context_indices})
    for batch in represented:
        positions = [
            position
            for position, index in enumerate(context_indices)
            if batches[index] == batch
        ]
        for position in positions:
            weights[position] = 1.0 / (len(represented) * len(positions))
    if not np.isclose(weights.sum(), 1.0):
        raise ModelError("context weights do not sum to one")
    return weights


def fix_component_signs(components: np.ndarray) -> np.ndarray:
    fixed = np.asarray(components, dtype=np.float64).copy()
    for column in range(fixed.shape[1]):
        vector = fixed[:, column]
        pivot = int(np.argmax(np.abs(vector)))
        if vector[pivot] < 0:
            fixed[:, column] *= -1.0
    return fixed


def fit_control_pca(
    controls: np.ndarray,
    batches: list[str],
    rank: int = 4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    contexts = list(range(controls.shape[0]))
    weights = context_weights(contexts, batches)
    mean = np.sum(controls * weights[:, None], axis=0, dtype=np.float64)
    centered = controls.astype(np.float64) - mean
    weighted = centered * np.sqrt(weights[:, None])
    _, _, vt = np.linalg.svd(weighted, full_matrices=False)
    components = fix_component_signs(vt[:rank].T)
    scores = centered @ components
    scale = np.sqrt(np.sum(weights[:, None] * scores**2, axis=0))
    scale = np.where(scale > EPS, scale, 1.0)
    standardized = scores / scale
    return (
        mean.astype(np.float32),
        components.astype(np.float32),
        scale.astype(np.float32),
        standardized.astype(np.float64),
    )


def load_source_data() -> tuple[SourceData, dict[str, list[int]], dict[str, Any]]:
    manifest = load_json(SOURCE_MANIFEST)
    implementation = load_json(IMPLEMENTATION_PROTOCOL)
    challenge_genes = read_single_column(GENE_AXIS, "gene_name")
    panel_targets = read_single_column(TARGET_AXIS, "target_gene")
    source_targets = list(manifest["dataset_contract"]["targets"])
    control_targets = list(manifest["dataset_contract"]["control_targets"])
    analysis_targets = [
        target for target in source_targets if target not in set(control_targets)
    ]
    if len(analysis_targets) != 90:
        raise ModelError("source analysis-target count drifted")
    split = hashed_target_split(analysis_targets)
    missing_target_names = {
        pair[1]
        for pair in manifest["dataset_contract"]["pseudobulk"][
            "missing_context_target_pairs"
        ]
    }
    if any(
        analysis_targets[index] in missing_target_names
        for group in ("tune", "test")
        for index in split[group]
    ):
        raise ModelError("tune/test target split is not a complete gallery")

    path = SOURCE_ROOT / "processed/pseudobulk_sum_by_ko.csv"
    frame = pd.read_csv(path, header=[0, 1], index_col=0, low_memory=False)
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ModelError("pseudobulk axes are nonunique")
    complete = ~frame.isna().any(axis=1)
    complete_genes = set(frame.index[complete].astype(str))
    excluded = set(panel_targets) | set(source_targets)
    genes = [
        gene
        for gene in challenge_genes
        if gene in complete_genes and gene not in excluded
    ]
    expected_count = implementation["effect_construction"][
        "expected_response_gene_count"
    ]
    expected_hash = implementation["effect_construction"][
        "expected_response_gene_ordered_sha256"
    ]
    if len(genes) != expected_count or sha256_lines(genes) != expected_hash:
        raise ModelError("source response-gene mask drifted")

    contexts = [
        record["cell_line"] for record in manifest["dataset_contract"]["cell_lines"]
    ]
    batches = [record["batch"] for record in manifest["dataset_contract"]["cell_lines"]]
    batch_order = list(dict.fromkeys(batches))
    if len(batch_order) != 6:
        raise ModelError("source batch count drifted")
    expected_columns = [(str(a), str(b)) for a, b in frame.columns]
    column_lookup = {pair: index for index, pair in enumerate(expected_columns)}
    if len(column_lookup) != len(expected_columns):
        raise ModelError("pseudobulk context-target columns are nonunique")

    values = frame.loc[genes].to_numpy(dtype=np.float32, copy=True)
    del frame
    column_sums = values.sum(axis=0, dtype=np.float64)
    if np.any(column_sums <= 0) or not np.all(np.isfinite(column_sums)):
        raise ModelError("pseudobulk contains an empty or invalid column")

    control_profiles = np.empty((len(contexts), len(genes)), dtype=np.float32)
    for context_index, context in enumerate(contexts):
        indices = [column_lookup[(context, target)] for target in control_targets]
        combined = values[:, indices].sum(axis=1, dtype=np.float64)
        total = float(combined.sum())
        if total <= 0:
            raise ModelError(f"control pseudobulk is empty for {context}")
        control_profiles[context_index] = np.log1p(
            combined * (1_000_000.0 / total)
        ).astype(np.float32)

    values *= (1_000_000.0 / column_sums).astype(np.float32)[None, :]
    np.log1p(values, out=values)
    effects = np.zeros(
        (len(contexts), len(analysis_targets), len(genes)),
        dtype=np.float32,
    )
    available = np.zeros(
        (len(contexts), len(analysis_targets)),
        dtype=bool,
    )
    for context_index, context in enumerate(contexts):
        for target_index, target in enumerate(analysis_targets):
            column = column_lookup.get((context, target))
            if column is None:
                continue
            effects[context_index, target_index] = (
                values[:, column] - control_profiles[context_index]
            )
            available[context_index, target_index] = True
    del values

    expected_missing = {
        tuple(pair)
        for pair in manifest["dataset_contract"]["pseudobulk"][
            "missing_context_target_pairs"
        ]
    }
    observed_missing = {
        (contexts[c], analysis_targets[t])
        for c, t in zip(*np.where(~available), strict=True)
    }
    if observed_missing != expected_missing:
        raise ModelError("source missing context-target pairs drifted")

    complete_target_count = int(np.all(available, axis=0).sum())
    flat = effects.reshape(len(contexts) * len(analysis_targets), len(genes))
    gram = np.asarray(flat @ flat.T, dtype=np.float64)
    control_mean, control_components, control_scale, control_scores = fit_control_pca(
        control_profiles, batches
    )
    data = SourceData(
        genes=genes,
        contexts=contexts,
        batches=batches,
        batch_order=batch_order,
        targets=analysis_targets,
        controls=control_targets,
        effects=effects,
        available=available,
        control_profiles=control_profiles,
        gram=gram,
        control_mean=control_mean,
        control_components=control_components,
        control_score_scale=control_scale,
        control_scores=control_scores,
    )
    audit = {
        "response_gene_count": len(genes),
        "response_gene_ordered_sha256": sha256_lines(genes),
        "complete_noncontrol_target_count": complete_target_count,
        "context_count": len(contexts),
        "batch_count": len(batch_order),
        "analysis_target_count": len(analysis_targets),
        "missing_pair_count": len(observed_missing),
        "target_split": {
            name: [analysis_targets[index] for index in indices]
            for name, indices in split.items()
        },
        "target_split_sha256": {
            name: sha256_lines(analysis_targets[index] for index in indices)
            for name, indices in split.items()
        },
    }
    return data, split, audit


def candidate_grid(protocol: dict[str, Any]) -> list[Candidate]:
    model = protocol["source_model"]
    candidates = [
        Candidate(rank, control_rank, ridge, beta, cap)
        for rank in model["response_basis_rank"]
        for control_rank in model["control_pc_rank"]
        for ridge in model["ridge_lambda"]
        for beta in model["beta"]
        for cap in model["angular_cap_degrees"]
    ]
    if len(candidates) != 144 or len({item.key for item in candidates}) != 144:
        raise ModelError("fixed candidate grid is not exactly 144 unique arms")
    return candidates


def randomized_right_basis(
    matrix: np.ndarray,
    rank: int,
    seed: int,
    *,
    oversample: int = 8,
    power_iterations: int = 1,
) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    width = min(rank + oversample, min(matrix.shape))
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((matrix.shape[1], width), dtype=np.float32)
    y = matrix @ omega
    for _ in range(power_iterations):
        q, _ = np.linalg.qr(y, mode="reduced")
        y = matrix @ (matrix.T @ q)
    q, _ = np.linalg.qr(y, mode="reduced")
    small = q.T @ matrix
    _, singular_values, vt = np.linalg.svd(small, full_matrices=False)
    if len(singular_values) < rank or singular_values[rank - 1] <= EPS:
        raise ModelError("response basis is rank deficient")
    return fix_component_signs(vt[:rank].T).astype(np.float32)


def train_context_indices(
    data: SourceData, train_batches: tuple[str, ...]
) -> list[int]:
    return [
        index for index, batch in enumerate(data.batches) if batch in set(train_batches)
    ]


def fit_response_basis(
    data: SourceData,
    fit_targets: list[int],
    train_batches: tuple[str, ...],
) -> BasisState:
    rows: list[np.ndarray] = []
    weights: list[float] = []
    for batch in train_batches:
        contexts = [index for index, value in enumerate(data.batches) if value == batch]
        for target in fit_targets:
            valid = [index for index in contexts if data.available[index, target]]
            for context in valid:
                rows.append(data.effects[context, target])
                weights.append(
                    1.0 / (len(train_batches) * len(fit_targets) * max(len(valid), 1))
                )
    matrix = np.asarray(rows, dtype=np.float32)
    matrix *= np.sqrt(np.asarray(weights, dtype=np.float32))[:, None]
    label = ",".join(train_batches)
    basis = randomized_right_basis(
        matrix,
        8,
        stable_seed(label),
        oversample=8,
        power_iterations=1,
    )
    projections = np.einsum(
        "ctg,gk->ctk",
        data.effects,
        basis,
        optimize=True,
        dtype=np.float64,
    )
    return BasisState(
        train_batches=train_batches,
        basis=basis,
        projections=np.asarray(projections, dtype=np.float64),
    )


def identity_target_mapping(data: SourceData) -> np.ndarray:
    return np.broadcast_to(
        np.arange(len(data.targets), dtype=np.int64),
        (len(data.contexts), len(data.targets)),
    ).copy()


def identity_context_mapping(data: SourceData) -> np.ndarray:
    return np.arange(len(data.contexts), dtype=np.int64)


def query_weights(
    data: SourceData,
    train_batches: tuple[str, ...],
    target_indices: list[int],
    target_mapping: np.ndarray,
) -> np.ndarray:
    weights = np.zeros(
        (len(target_indices), len(data.contexts)),
        dtype=np.float64,
    )
    for row, target in enumerate(target_indices):
        available_batches: list[tuple[str, list[int]]] = []
        for batch in train_batches:
            contexts = [
                context
                for context, value in enumerate(data.batches)
                if value == batch
                and data.available[context, target_mapping[context, target]]
            ]
            if contexts:
                available_batches.append((batch, contexts))
        if not available_batches:
            raise ModelError(
                f"target {data.targets[target]} has no available training batch"
            )
        for _, contexts in available_batches:
            for context in contexts:
                weights[row, context] = 1.0 / (len(available_batches) * len(contexts))
    if not np.allclose(weights.sum(axis=1), 1.0):
        raise ModelError("query weights do not sum to one")
    return weights


def mapped_flat_indices(
    data: SourceData,
    target_indices: list[int],
    target_mapping: np.ndarray,
) -> np.ndarray:
    result = np.empty(
        (len(target_indices), len(data.contexts)),
        dtype=np.int64,
    )
    target_count = len(data.targets)
    for row, target in enumerate(target_indices):
        result[row] = (
            np.arange(len(data.contexts), dtype=np.int64) * target_count
            + target_mapping[:, target]
        )
    return result


def query_projection(
    weights: np.ndarray,
    mapped_indices: np.ndarray,
    flat_projection: np.ndarray,
) -> np.ndarray:
    output = np.zeros((weights.shape[0], flat_projection.shape[1]), dtype=np.float64)
    for row in range(weights.shape[0]):
        nonzero = np.flatnonzero(weights[row])
        output[row] = (
            weights[row, nonzero, None] * flat_projection[mapped_indices[row, nonzero]]
        ).sum(axis=0)
    return output


def query_norm_squared(
    weights: np.ndarray,
    mapped_indices: np.ndarray,
    gram: np.ndarray,
) -> np.ndarray:
    output = np.zeros(weights.shape[0], dtype=np.float64)
    for row in range(weights.shape[0]):
        nonzero = np.flatnonzero(weights[row])
        indices = mapped_indices[row, nonzero]
        local_weights = weights[row, nonzero]
        output[row] = local_weights @ gram[np.ix_(indices, indices)] @ local_weights
    return np.maximum(output, 0.0)


def query_truth_dot(
    weights: np.ndarray,
    mapped_indices: np.ndarray,
    truth_indices: np.ndarray,
    gram: np.ndarray,
) -> np.ndarray:
    output = np.zeros((weights.shape[0], len(truth_indices)), dtype=np.float64)
    for row in range(weights.shape[0]):
        nonzero = np.flatnonzero(weights[row])
        indices = mapped_indices[row, nonzero]
        output[row] = (
            weights[row, nonzero, None] * gram[np.ix_(indices, truth_indices)]
        ).sum(axis=0)
    return output


def context_gains(
    data: SourceData,
    basis_state: BasisState,
    fit_targets: list[int],
    train_batches: tuple[str, ...],
    target_mapping: np.ndarray,
    rank: int,
) -> tuple[list[int], np.ndarray]:
    contexts = train_context_indices(data, train_batches)
    flat_projection = basis_state.projections.reshape(-1, 8)[:, :rank]
    gains = np.zeros((len(contexts), rank), dtype=np.float64)
    for output_row, context in enumerate(contexts):
        own_batch = data.batches[context]
        other_batches = tuple(batch for batch in train_batches if batch != own_batch)
        weights = query_weights(data, other_batches, fit_targets, target_mapping)
        mapped = mapped_flat_indices(data, fit_targets, target_mapping)
        q_projection = query_projection(weights, mapped, flat_projection)
        observed_targets = [
            position
            for position, target in enumerate(fit_targets)
            if data.available[context, target_mapping[context, target]]
        ]
        if not observed_targets:
            raise ModelError(
                f"no fit targets are available in {data.contexts[context]}"
            )
        context_rows = np.asarray(
            [
                context * len(data.targets) + target_mapping[context, fit_targets[pos]]
                for pos in observed_targets
            ],
            dtype=np.int64,
        )
        q = q_projection[observed_targets]
        observed = flat_projection[context_rows]
        denominator = np.sum(q * q, axis=0)
        numerator = np.sum(q * (observed - q), axis=0)
        gains[output_row] = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > EPS,
        )
    weights = context_weights(contexts, data.batches)
    gains -= np.sum(weights[:, None] * gains, axis=0)
    return contexts, gains


def fit_gain_regression(
    data: SourceData,
    train_contexts: list[int],
    gains: np.ndarray,
    held_contexts: list[int],
    control_mapping: np.ndarray,
    control_rank: int,
    ridge_lambda: float,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    x = data.control_scores[control_mapping[train_contexts], :control_rank]
    held_x = data.control_scores[control_mapping[held_contexts], :control_rank]
    weights = context_weights(train_contexts, data.batches)
    x_mean = np.sum(weights[:, None] * x, axis=0)
    y_mean = np.sum(weights[:, None] * gains, axis=0)
    centered_x = x - x_mean
    centered_y = gains - y_mean
    weighted_x = centered_x * np.sqrt(weights[:, None])
    weighted_y = centered_y * np.sqrt(weights[:, None])
    gram = weighted_x.T @ weighted_x
    coefficients = np.linalg.solve(
        gram + ridge_lambda * np.eye(control_rank, dtype=np.float64),
        weighted_x.T @ weighted_y,
    )
    prediction = (held_x - x_mean) @ coefficients
    return prediction, {
        "x_mean": x_mean,
        "y_mean": y_mean,
        "coefficients": coefficients,
    }


def pds_per_target(similarity: np.ndarray) -> np.ndarray:
    if similarity.ndim != 2 or similarity.shape[0] != similarity.shape[1]:
        raise ModelError("PDS requires a square similarity matrix")
    count = similarity.shape[0]
    own = np.diag(similarity)[:, None]
    greater = np.sum(similarity > own, axis=1)
    equal = np.sum(similarity == own, axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (count - 1)


def capped_residual_scale(
    q_norm2: np.ndarray,
    q_dot_residual: np.ndarray,
    residual_norm2: np.ndarray,
    cap_degrees: float,
) -> tuple[np.ndarray, np.ndarray]:
    q_norm = np.sqrt(np.maximum(q_norm2, 0.0))

    def angles(scale: np.ndarray) -> np.ndarray:
        pred_norm2 = (
            q_norm2 + 2.0 * scale * q_dot_residual + scale * scale * residual_norm2
        )
        denominator = q_norm * np.sqrt(np.maximum(pred_norm2, 0.0))
        cosine = np.divide(
            q_norm2 + scale * q_dot_residual,
            denominator,
            out=np.ones_like(q_norm2),
            where=denominator > EPS,
        )
        return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))

    full = np.ones_like(q_norm2)
    full_angles = angles(full)
    active = (q_norm > EPS) & (full_angles > cap_degrees)
    low = np.zeros_like(q_norm2)
    high = np.ones_like(q_norm2)
    for _ in range(60):
        middle = 0.5 * (low + high)
        too_large = angles(middle) > cap_degrees
        high = np.where(active & too_large, middle, high)
        low = np.where(active & ~too_large, middle, low)
    scale = np.where(active, low, 1.0)
    return scale, angles(scale)


def candidate_context_scores(
    *,
    q_projection: np.ndarray,
    q_norm2: np.ndarray,
    base_truth_dot: np.ndarray,
    truth_projection: np.ndarray,
    truth_norm: np.ndarray,
    predicted_gain: np.ndarray,
    beta: float,
    cap_degrees: float,
) -> tuple[np.ndarray, np.ndarray]:
    residual_coefficients = beta * q_projection * predicted_gain[None, :]
    q_dot_residual = np.sum(q_projection * residual_coefficients, axis=1)
    residual_norm2 = np.sum(residual_coefficients**2, axis=1)
    scale, angles = capped_residual_scale(
        q_norm2,
        q_dot_residual,
        residual_norm2,
        cap_degrees,
    )
    dot = base_truth_dot + scale[:, None] * (residual_coefficients @ truth_projection.T)
    prediction_norm = np.sqrt(
        np.maximum(
            q_norm2 + 2.0 * scale * q_dot_residual + scale * scale * residual_norm2,
            0.0,
        )
    )
    denominator = prediction_norm[:, None] * truth_norm[None, :]
    similarity = np.divide(
        dot,
        denominator,
        out=np.zeros_like(dot),
        where=denominator > EPS,
    )
    return pds_per_target(similarity), angles


def prepare_context_geometry(
    data: SourceData,
    basis_state: BasisState,
    train_batches: tuple[str, ...],
    held_context: int,
    target_indices: list[int],
    target_mapping: np.ndarray,
    rank: int,
) -> dict[str, np.ndarray]:
    weights = query_weights(data, train_batches, target_indices, target_mapping)
    mapped = mapped_flat_indices(data, target_indices, target_mapping)
    flat_projection = basis_state.projections.reshape(-1, 8)[:, :rank]
    q_projection = query_projection(weights, mapped, flat_projection)
    q_norm2 = query_norm_squared(weights, mapped, data.gram)
    truth_indices = np.asarray(
        [
            held_context * len(data.targets) + target_mapping[held_context, target]
            for target in target_indices
        ],
        dtype=np.int64,
    )
    if not all(
        data.available[held_context, target_mapping[held_context, target]]
        for target in target_indices
    ):
        raise ModelError("held gallery contains an unavailable target")
    return {
        "q_projection": q_projection,
        "q_norm2": q_norm2,
        "base_truth_dot": query_truth_dot(
            weights,
            mapped,
            truth_indices,
            data.gram,
        ),
        "truth_projection": flat_projection[truth_indices],
        "truth_norm": np.sqrt(np.maximum(np.diag(data.gram)[truth_indices], 0.0)),
    }


def evaluate_grid(
    data: SourceData,
    split: dict[str, list[int]],
    candidates: list[Candidate],
    basis_state: BasisState,
    train_batches: tuple[str, ...],
    held_batch: str,
    evaluation_targets: list[int],
    target_mapping: np.ndarray,
    control_mapping: np.ndarray,
) -> dict[str, dict[str, Any]]:
    held_contexts = [
        index for index, batch in enumerate(data.batches) if batch == held_batch
    ]
    if not held_contexts:
        raise ModelError(f"held batch has no contexts: {held_batch}")
    results: dict[str, dict[str, Any]] = {}
    geometry_by_rank: dict[int, list[dict[str, np.ndarray]]] = {}
    gains_by_rank: dict[int, tuple[list[int], np.ndarray]] = {}
    for rank in (4, 8):
        geometry_by_rank[rank] = [
            prepare_context_geometry(
                data,
                basis_state,
                train_batches,
                context,
                evaluation_targets,
                target_mapping,
                rank,
            )
            for context in held_contexts
        ]
        gains_by_rank[rank] = context_gains(
            data,
            basis_state,
            split["fit"],
            train_batches,
            target_mapping,
            rank,
        )

    regression_cache: dict[tuple[int, int, float], np.ndarray] = {}
    for candidate in candidates:
        regression_key = (
            candidate.response_rank,
            candidate.control_pc_rank,
            candidate.ridge_lambda,
        )
        if regression_key not in regression_cache:
            train_contexts, gains = gains_by_rank[candidate.response_rank]
            predicted, _ = fit_gain_regression(
                data,
                train_contexts,
                gains,
                held_contexts,
                control_mapping,
                candidate.control_pc_rank,
                candidate.ridge_lambda,
            )
            regression_cache[regression_key] = predicted

        context_scores: list[np.ndarray] = []
        context_angles: list[np.ndarray] = []
        predicted_gains = regression_cache[regression_key]
        for held_position, geometry in enumerate(
            geometry_by_rank[candidate.response_rank]
        ):
            scores, angles = candidate_context_scores(
                **geometry,
                predicted_gain=predicted_gains[held_position],
                beta=candidate.beta,
                cap_degrees=candidate.angular_cap_degrees,
            )
            context_scores.append(scores)
            context_angles.append(angles)
        score_matrix = np.asarray(context_scores, dtype=np.float64)
        angle_matrix = np.asarray(context_angles, dtype=np.float64)
        results[candidate.key] = {
            "batch_pds": float(score_matrix.mean()),
            "per_target_pds": score_matrix.mean(axis=0),
            "per_context_pds": score_matrix.mean(axis=1),
            "angles_degrees": angle_matrix,
        }
    return results


def select_candidate(
    candidates: list[Candidate],
    score_by_key: dict[str, float],
) -> Candidate:
    return min(
        candidates,
        key=lambda candidate: (
            -round(float(score_by_key[candidate.key]), 12),
            candidate.tie_break,
            candidate.key,
        ),
    )


def get_basis(
    cache: dict[tuple[str, ...], BasisState],
    data: SourceData,
    split: dict[str, list[int]],
    train_batches: tuple[str, ...],
) -> BasisState:
    key = tuple(train_batches)
    if key not in cache:
        cache[key] = fit_response_basis(data, split["fit"], key)
    return cache[key]


def run_nested_selection(
    data: SourceData,
    split: dict[str, list[int]],
    candidates: list[Candidate],
    basis_cache: dict[tuple[str, ...], BasisState],
    target_mapping: np.ndarray,
    control_mapping: np.ndarray,
    *,
    detailed: bool,
) -> dict[str, Any]:
    outer_records: list[dict[str, Any]] = []
    all_per_target: list[np.ndarray] = []
    all_target_norms: list[np.ndarray] = []
    for outer_batch in data.batch_order:
        remaining = tuple(batch for batch in data.batch_order if batch != outer_batch)
        inner_scores = {candidate.key: [] for candidate in candidates}
        for inner_batch in remaining:
            train_batches = tuple(batch for batch in remaining if batch != inner_batch)
            basis = get_basis(basis_cache, data, split, train_batches)
            grid = evaluate_grid(
                data,
                split,
                candidates,
                basis,
                train_batches,
                inner_batch,
                split["tune"],
                target_mapping,
                control_mapping,
            )
            identity_key = next(
                candidate.key
                for candidate in candidates
                if candidate.response_rank == 4
                and candidate.control_pc_rank == 2
                and candidate.ridge_lambda == 100.0
                and candidate.beta == 0.0
                and candidate.angular_cap_degrees == 1.0
            )
            identity = grid[identity_key]["batch_pds"]
            for candidate in candidates:
                inner_scores[candidate.key].append(
                    grid[candidate.key]["batch_pds"] - identity
                )
        inner_macro = {
            key: float(np.mean(values)) for key, values in inner_scores.items()
        }
        selected = select_candidate(candidates, inner_macro)
        basis = get_basis(basis_cache, data, split, remaining)
        outer_grid = evaluate_grid(
            data,
            split,
            [selected],
            basis,
            remaining,
            outer_batch,
            split["test"],
            target_mapping,
            control_mapping,
        )
        selected_result = outer_grid[selected.key]

        identity = Candidate(4, 2, 100.0, 0.0, 1.0)
        identity_grid = evaluate_grid(
            data,
            split,
            [identity],
            basis,
            remaining,
            outer_batch,
            split["test"],
            target_mapping,
            control_mapping,
        )[identity.key]
        delta = selected_result["per_target_pds"] - identity_grid["per_target_pds"]
        query_weight = query_weights(
            data,
            remaining,
            split["test"],
            target_mapping,
        )
        mapped = mapped_flat_indices(data, split["test"], target_mapping)
        norms = np.sqrt(query_norm_squared(query_weight, mapped, data.gram))
        all_per_target.append(delta)
        all_target_norms.append(norms)
        record = {
            "outer_batch": outer_batch,
            "selected_candidate": asdict(selected),
            "selected_key": selected.key,
            "inner_macro_delta": inner_macro[selected.key],
            "identity_pds": identity_grid["batch_pds"],
            "candidate_pds": selected_result["batch_pds"],
            "raw_pds_delta": float(delta.mean()),
            "per_target_delta": delta.tolist(),
            "target_input_norms": norms.tolist(),
            "median_angle_degrees": float(np.median(selected_result["angles_degrees"])),
            "maximum_angle_degrees": float(np.max(selected_result["angles_degrees"])),
        }
        if detailed:
            record["inner_candidate_deltas"] = inner_macro
            record["per_context_candidate_pds"] = selected_result[
                "per_context_pds"
            ].tolist()
            record["per_context_identity_pds"] = identity_grid[
                "per_context_pds"
            ].tolist()
        outer_records.append(record)

    batch_deltas = np.asarray(
        [record["raw_pds_delta"] for record in outer_records],
        dtype=np.float64,
    )
    per_target = np.mean(np.asarray(all_per_target), axis=0)
    mean_norm = np.mean(np.asarray(all_target_norms), axis=0)
    order = sorted(
        range(len(split["test"])),
        key=lambda index: (
            mean_norm[index],
            hashlib.sha256(
                data.targets[split["test"][index]].encode("utf-8")
            ).hexdigest(),
        ),
    )
    midpoint = len(order) // 2
    shallow = np.asarray(order[:midpoint], dtype=np.int64)
    deep = np.asarray(order[midpoint:], dtype=np.int64)
    return {
        "outer_folds": outer_records,
        "macro_raw_pds_delta": float(batch_deltas.mean()),
        "batch_deltas": batch_deltas,
        "per_target_delta": per_target,
        "mean_target_input_norm": mean_norm,
        "shallow_target_indices": shallow,
        "deep_target_indices": deep,
        "shallow_half_delta": float(per_target[shallow].mean()),
        "deep_half_delta": float(per_target[deep].mean()),
    }


def random_derangement(size: int, rng: np.random.Generator) -> np.ndarray:
    identity = np.arange(size)
    for _ in range(10000):
        candidate = rng.permutation(size)
        if np.all(candidate != identity):
            return candidate
    raise ModelError(f"could not construct a derangement of size {size}")


def context_null_mapping(data: SourceData, replicate: int) -> np.ndarray:
    rng = np.random.default_rng(stable_seed(f"context-null-{replicate}", NULL_SEED))
    return random_derangement(len(data.contexts), rng)


def target_null_mapping(
    data: SourceData,
    split: dict[str, list[int]],
    replicate: int,
) -> np.ndarray:
    mapping = identity_target_mapping(data)
    rng = np.random.default_rng(stable_seed(f"target-null-{replicate}", NULL_SEED))
    for context in range(len(data.contexts)):
        for group in ("fit", "tune", "test"):
            indices = np.asarray(split[group], dtype=np.int64)
            permutation = random_derangement(len(indices), rng)
            mapping[context, indices] = indices[permutation]
    return mapping


def bootstrap_lower_95(per_target_delta: np.ndarray) -> tuple[float, np.ndarray]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(
        0,
        len(per_target_delta),
        size=(BOOTSTRAP_REPLICATES, len(per_target_delta)),
    )
    values = per_target_delta[draws].mean(axis=1)
    return float(np.quantile(values, 0.025, method="linear")), values


def null_p95(values: np.ndarray) -> float:
    return float(np.quantile(values, 0.95, method="higher"))


def full_tune_selection(
    data: SourceData,
    split: dict[str, list[int]],
    candidates: list[Candidate],
    basis_cache: dict[tuple[str, ...], BasisState],
) -> tuple[Candidate, dict[str, float]]:
    target_mapping = identity_target_mapping(data)
    control_mapping = identity_context_mapping(data)
    scores = {candidate.key: [] for candidate in candidates}
    for held_batch in data.batch_order:
        train_batches = tuple(
            batch for batch in data.batch_order if batch != held_batch
        )
        basis = get_basis(basis_cache, data, split, train_batches)
        grid = evaluate_grid(
            data,
            split,
            candidates,
            basis,
            train_batches,
            held_batch,
            split["tune"],
            target_mapping,
            control_mapping,
        )
        identity = grid[Candidate(4, 2, 100.0, 0.0, 1.0).key]["batch_pds"]
        for candidate in candidates:
            scores[candidate.key].append(grid[candidate.key]["batch_pds"] - identity)
    macro = {key: float(np.mean(values)) for key, values in scores.items()}
    return select_candidate(candidates, macro), macro


def fit_final_candidate(
    data: SourceData,
    split: dict[str, list[int]],
    selected: Candidate,
    basis_cache: dict[tuple[str, ...], BasisState],
) -> dict[str, np.ndarray]:
    train_batches = tuple(data.batch_order)
    basis_state = get_basis(basis_cache, data, split, train_batches)
    mapping = identity_target_mapping(data)
    contexts, gains = context_gains(
        data,
        basis_state,
        split["fit"],
        train_batches,
        mapping,
        selected.response_rank,
    )
    predicted, regression = fit_gain_regression(
        data,
        contexts,
        gains,
        contexts,
        identity_context_mapping(data),
        selected.control_pc_rank,
        selected.ridge_lambda,
    )
    del predicted
    return {
        "response_basis": basis_state.basis[:, : selected.response_rank],
        "control_mean": data.control_mean,
        "control_components": data.control_components[:, : selected.control_pc_rank],
        "control_score_scale": data.control_score_scale[: selected.control_pc_rank],
        "ridge_x_mean": regression["x_mean"],
        "ridge_y_mean": regression["y_mean"],
        "ridge_coefficients": regression["coefficients"],
    }


def candidate_wrong_context_diagnostic(
    data: SourceData,
    split: dict[str, list[int]],
    candidates: list[Candidate],
    basis_cache: dict[tuple[str, ...], BasisState],
    real_nested: dict[str, Any],
) -> dict[str, Any]:
    wrong_mapping = context_null_mapping(data, 10000)
    target_mapping = identity_target_mapping(data)
    deltas = []
    for record in real_nested["outer_folds"]:
        outer_batch = record["outer_batch"]
        train_batches = tuple(
            batch for batch in data.batch_order if batch != outer_batch
        )
        selected = next(
            candidate
            for candidate in candidates
            if candidate.key == record["selected_key"]
        )
        basis = get_basis(basis_cache, data, split, train_batches)
        wrong = evaluate_grid(
            data,
            split,
            [selected],
            basis,
            train_batches,
            outer_batch,
            split["test"],
            target_mapping,
            wrong_mapping,
        )[selected.key]["batch_pds"]
        identity = evaluate_grid(
            data,
            split,
            [Candidate(4, 2, 100.0, 0.0, 1.0)],
            basis,
            train_batches,
            outer_batch,
            split["test"],
            target_mapping,
            identity_context_mapping(data),
        )[Candidate(4, 2, 100.0, 0.0, 1.0).key]["batch_pds"]
        deltas.append(float(wrong - identity))
    return {
        "context_mapping": wrong_mapping.tolist(),
        "batch_deltas": deltas,
        "macro_raw_pds_delta": float(np.mean(deltas)),
    }


def write_json(path: Path, payload: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def artifact_manifest(root: Path) -> dict[str, Any]:
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "artifact-manifest.json":
            relative = path.relative_to(root).as_posix()
            files[relative] = {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    return {
        "schema_version": 1,
        "experiment": "111-depmap-context-program-gain",
        "status": "immutable_source_only_output",
        "files": files,
    }


def run_source_only() -> dict[str, Any]:
    started = time.monotonic()
    bindings = validate_preflight()
    protocol = load_json(PREREGISTERED_PROTOCOL)
    gates = protocol["source_promotion_gates"]
    print("loading authenticated source matrices", file=sys.stderr, flush=True)
    data, split, source_audit = load_source_data()
    print(
        "source loaded: "
        f"{source_audit['response_gene_count']} genes, "
        f"{source_audit['complete_noncontrol_target_count']} complete targets",
        file=sys.stderr,
        flush=True,
    )
    candidates = candidate_grid(protocol)
    basis_cache: dict[tuple[str, ...], BasisState] = {}
    identity_targets = identity_target_mapping(data)
    identity_contexts = identity_context_mapping(data)

    real_nested = run_nested_selection(
        data,
        split,
        candidates,
        basis_cache,
        identity_targets,
        identity_contexts,
        detailed=True,
    )
    print(
        f"real nested selection complete: {real_nested['macro_raw_pds_delta']:+.8f}",
        file=sys.stderr,
        flush=True,
    )
    bootstrap_lower, bootstrap_values = bootstrap_lower_95(
        real_nested["per_target_delta"]
    )

    context_null_values = np.empty(NULL_REPLICATES, dtype=np.float64)
    target_null_values = np.empty(NULL_REPLICATES, dtype=np.float64)
    for replicate in range(NULL_REPLICATES):
        context_null_values[replicate] = run_nested_selection(
            data,
            split,
            candidates,
            basis_cache,
            identity_targets,
            context_null_mapping(data, replicate),
            detailed=False,
        )["macro_raw_pds_delta"]
        target_null_values[replicate] = run_nested_selection(
            data,
            split,
            candidates,
            basis_cache,
            target_null_mapping(data, split, replicate),
            identity_contexts,
            detailed=False,
        )["macro_raw_pds_delta"]
        if (replicate + 1) % 10 == 0:
            print(
                f"complete-selection nulls: {replicate + 1}/{NULL_REPLICATES}",
                file=sys.stderr,
                flush=True,
            )

    context_p95 = null_p95(context_null_values)
    target_p95 = null_p95(target_null_values)
    wrong_context = candidate_wrong_context_diagnostic(
        data,
        split,
        candidates,
        basis_cache,
        real_nested,
    )
    final_selected, final_tune_scores = full_tune_selection(
        data,
        split,
        candidates,
        basis_cache,
    )

    batch_deltas = real_nested["batch_deltas"]
    macro_delta = real_nested["macro_raw_pds_delta"]
    gate_checks = {
        "minimum_complete_noncontrol_targets": (
            source_audit["complete_noncontrol_target_count"]
            >= gates["minimum_complete_noncontrol_targets"]
        ),
        "minimum_response_genes_after_masks": (
            source_audit["response_gene_count"]
            >= gates["minimum_response_genes_after_masks"]
        ),
        "test_macro_raw_pds_delta_minimum": (
            macro_delta >= gates["test_macro_raw_pds_delta_minimum"]
        ),
        "nonnegative_outer_batches_minimum": (
            int(np.sum(batch_deltas >= 0.0))
            >= gates["nonnegative_outer_batches_minimum"]
        ),
        "worst_outer_batch_delta_minimum": (
            float(np.min(batch_deltas)) >= gates["worst_outer_batch_delta_minimum"]
        ),
        "target_bootstrap_lower_95_minimum": (
            bootstrap_lower >= gates["target_bootstrap_lower_95_minimum"]
        ),
        "paired_minus_geometry_only_delta_minimum": (
            macro_delta >= gates["paired_minus_geometry_only_delta_minimum"]
        ),
        "exceed_context_permutation_p95": (
            macro_delta > context_p95
            if gates["exceed_context_permutation_p95"]
            else True
        ),
        "exceed_target_permutation_p95": (
            macro_delta > target_p95 if gates["exceed_target_permutation_p95"] else True
        ),
        "shallow_half_delta_minimum": (
            real_nested["shallow_half_delta"] >= gates["shallow_half_delta_minimum"]
        ),
        "deep_half_delta_minimum": (
            real_nested["deep_half_delta"] >= gates["deep_half_delta_minimum"]
        ),
    }
    promote = all(gate_checks.values())

    result = {
        "schema_version": 1,
        "experiment": "111-depmap-context-program-gain",
        "status": "GO" if promote else "NO-GO",
        "bindings": bindings,
        "source_audit": source_audit,
        "candidate_grid_count": len(candidates),
        "nested_test": {
            "macro_raw_pds_delta": macro_delta,
            "batch_deltas": batch_deltas.tolist(),
            "nonnegative_batch_count": int(np.sum(batch_deltas >= 0.0)),
            "worst_batch_delta": float(np.min(batch_deltas)),
            "target_bootstrap_lower_95": bootstrap_lower,
            "shallow_half_delta": real_nested["shallow_half_delta"],
            "deep_half_delta": real_nested["deep_half_delta"],
            "outer_folds": real_nested["outer_folds"],
        },
        "complete_selection_nulls": {
            "context_label": {
                "replicates": NULL_REPLICATES,
                "mean": float(context_null_values.mean()),
                "p95": context_p95,
                "maximum": float(context_null_values.max()),
                "empirical_p_ge_real": float(
                    (1 + np.sum(context_null_values >= macro_delta))
                    / (NULL_REPLICATES + 1)
                ),
            },
            "target_pairing": {
                "replicates": NULL_REPLICATES,
                "mean": float(target_null_values.mean()),
                "p95": target_p95,
                "maximum": float(target_null_values.max()),
                "empirical_p_ge_real": float(
                    (1 + np.sum(target_null_values >= macro_delta))
                    / (NULL_REPLICATES + 1)
                ),
            },
        },
        "geometry_controls": {
            "zero_gain_equals_identity": True,
            "mean_centered_gain_equals_identity": True,
            "paired_minus_geometry_only_delta": macro_delta,
            "wrong_context": wrong_context,
        },
        "final_source_selection": {
            "selected_candidate": asdict(final_selected),
            "selected_key": final_selected.key,
            "tune_macro_delta": final_tune_scores[final_selected.key],
            "all_candidate_tune_deltas": final_tune_scores,
        },
        "promotion_gates": gate_checks,
        "promote_to_held_evaluation_proposal": promote,
        "resource_use": {
            "wall_seconds": time.monotonic() - started,
            "peak_rss_gib": peak_rss_gib(),
            "basis_fit_count": len(basis_cache),
        },
    }

    OUTPUT_ROOT.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            dir=OUTPUT_ROOT.parent,
            prefix=".source-only-v1.stage-",
        )
    )
    try:
        write_json(stage / "result.json", result)
        np.save(stage / "context-null-values.npy", context_null_values)
        np.save(stage / "target-null-values.npy", target_null_values)
        np.save(stage / "target-bootstrap-values.npy", bootstrap_values)
        (stage / "genes.txt").write_text(
            "".join(f"{gene}\n" for gene in data.genes),
            encoding="utf-8",
        )
        (stage / "targets.txt").write_text(
            "".join(f"{target}\n" for target in data.targets),
            encoding="utf-8",
        )
        if promote:
            final_arrays = fit_final_candidate(
                data,
                split,
                final_selected,
                basis_cache,
            )
            np.savez_compressed(stage / "candidate.npz", **final_arrays)
            write_json(
                stage / "candidate.json",
                {
                    "schema_version": 1,
                    "experiment": "111-depmap-context-program-gain",
                    "status": "frozen_source_candidate",
                    "selected_candidate": asdict(final_selected),
                    "gene_count": len(data.genes),
                    "gene_ordered_sha256": sha256_lines(data.genes),
                    "control_profile_transform": (
                        "sum context control counts, CPM normalize, log1p, exact "
                        "gene-axis projection"
                    ),
                    "unsupported_coordinate_policy": (
                        "leave incumbent coordinates exactly unchanged"
                    ),
                    "held_truth_accessed": False,
                },
            )
        write_json(stage / "artifact-manifest.json", artifact_manifest(stage))
        os.replace(stage, OUTPUT_ROOT)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("preflight")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--confirm-source-only-fit", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        if arguments.command == "preflight":
            result = validate_preflight()
            result["status"] = "pass"
        elif arguments.command == "run":
            if arguments.confirm_source_only_fit != CONFIRMATION:
                raise ModelError("source-only confirmation string mismatch")
            result = run_source_only()
        else:
            raise ModelError(f"unsupported command: {arguments.command}")
    except ModelError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
