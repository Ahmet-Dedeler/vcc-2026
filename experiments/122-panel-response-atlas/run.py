#!/usr/bin/env python3
"""Run Experiment 122 Phase 1 using only frozen public source effects."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

import freeze_plans


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
INVENTORY_PATH = HERE / "SOURCE_INVENTORY.json"
AUTHORIZATION_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
IMPLEMENTATION_FREEZE_PATH = HERE / "IMPLEMENTATION_FREEZE.json"
PLAN_DIR = HERE / "work" / "preregistered-freeze-v1"
DEFAULT_OUTPUT = HERE / "work" / "source-only-phase1-v1"

FAMILY_NAMES = ("replogle_k562", "orion", "cd4")
NULL_DRAWS = 200


class ExperimentError(RuntimeError):
    """A source, preregistration, selection, or output invariant failed."""


@dataclass(frozen=True)
class Candidate:
    name: str
    location: str
    confidence: str
    minimum_family_support: int


@dataclass(frozen=True)
class Matrix:
    name: str
    effect: np.ndarray
    valid: np.ndarray


CANDIDATES = (
    Candidate("identity", "zero", "zero", 0),
    Candidate("equal_mean_magnitude", "mean", "magnitude", 1),
    Candidate("equal_mean_agreement", "mean", "agreement", 1),
    Candidate("equal_median_magnitude", "median", "magnitude", 1),
    Candidate("equal_median_agreement", "median", "agreement", 1),
    Candidate("sign_consensus_median_min2", "sign_consensus_median", "agreement", 2),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ExperimentError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(value.dtype).encode())
    digest.update(b"\0")
    digest.update(json.dumps(list(value.shape), separators=(",", ":")).encode())
    digest.update(b"\0")
    digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(str(row[column]) for row in reader)
    validate_axis(values, str(path))
    return values


def read_text_axis(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    validate_axis(values, str(path), require_unique=require_unique)
    return values


def validate_axis(
    values: Sequence[str],
    label: str,
    *,
    require_unique: bool = True,
) -> None:
    require(bool(values), f"{label}: empty axis")
    require(
        all(value and value == value.strip() for value in values),
        f"{label}: empty or padded labels",
    )
    if require_unique:
        require(len(values) == len(set(values)), f"{label}: duplicate labels")


def parse_bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ExperimentError(f"unexpected boolean value: {value!r}")


def validate_implementation_freeze() -> dict[str, Any]:
    require(
        IMPLEMENTATION_FREEZE_PATH.is_file(),
        "IMPLEMENTATION_FREEZE.json must be created before effect loading",
    )
    freeze = read_json(IMPLEMENTATION_FREEZE_PATH)
    if freeze.get("created_before_effect_matrix_loading") is not True:
        require(
            freeze.get("scientific_spec_frozen_before_effect_matrix_loading") is True,
            "implementation freeze timing declaration changed",
        )
        amendment = freeze.get("implementation_amendment")
        require(isinstance(amendment, dict), "implementation amendment missing")
        amendment_path = ROOT / amendment["path"]
        require(
            amendment_path.is_file()
            and sha256_file(amendment_path) == amendment["sha256"],
            "implementation amendment hash mismatch",
        )
    for relative, expected in freeze["files"].items():
        path = ROOT / relative
        require(path.is_file() and not path.is_symlink(), f"missing code file: {relative}")
        require(sha256_file(path) == expected, f"implementation changed: {relative}")
    return freeze


def validate_authorization(*, hash_large_inputs: bool) -> dict[str, Any]:
    freeze_plans.verify(PLAN_DIR)
    validate_implementation_freeze()
    authorization = read_json(AUTHORIZATION_PATH)
    protocol = read_json(PROTOCOL_PATH)
    inventory = read_json(INVENTORY_PATH)
    require(authorization["source_only_authorized"] is True, "source-only not authorized")
    require(
        authorization["protocol_sha256"] == sha256_file(PROTOCOL_PATH),
        "protocol authorization hash mismatch",
    )
    require(
        authorization["source_inventory_sha256"] == sha256_file(INVENTORY_PATH),
        "inventory authorization hash mismatch",
    )
    require(
        protocol["status"] == "written_before_any_authorized_effect_matrix_was_loaded",
        "protocol timing status changed",
    )
    require(
        inventory["created_before_authorized_effect_matrix_loading"] is True,
        "inventory timing status changed",
    )
    denied = {
        "read_a549",
        "read_pan_h1",
        "read_dld1_experiment_116_matrix_values",
        "read_vcc_held_harness_or_outcomes",
        "read_challenge_context_expression_values",
        "run_cell_eval2",
        "generate_counts",
        "build_challenge_artifact",
        "submit",
        "commit",
        "push",
    }
    for key in denied:
        require(
            authorization["authorization"].get(key) is False,
            f"forbidden authorization changed: {key}",
        )
    for relative, expected in sorted(authorization["allowed_input_files"].items()):
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT), f"input escapes repository: {relative}")
        require(path.is_file() and not path.is_symlink(), f"missing input: {relative}")
        if hash_large_inputs or path.suffix != ".npy":
            require(sha256_file(path) == expected, f"input hash changed: {relative}")
    return authorization


def panel_axes() -> tuple[tuple[str, ...], tuple[str, ...]]:
    genes = read_csv_axis(ROOT / "data/challenge/gene_names.csv", "gene_name")
    panel = read_csv_axis(ROOT / "data/challenge/pert_counts.csv", "target_gene")
    require(len(genes) == 18533, "canonical gene count changed")
    require(len(panel) == 300, "panel target count changed")
    require(set(panel) <= set(genes), "panel is not a subset of canonical genes")
    return genes, panel


def panel_indices(axis: Sequence[str], panel: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    positions: dict[str, list[int]] = {}
    for index, label in enumerate(axis):
        positions.setdefault(label, []).append(index)
    duplicated_panel = [label for label in panel if len(positions.get(label, [])) > 1]
    require(not duplicated_panel, f"duplicated panel coordinates: {duplicated_panel[:5]}")
    panel_positions: list[int] = []
    source_positions: list[int] = []
    for panel_index, label in enumerate(panel):
        if label in positions:
            panel_positions.append(panel_index)
            source_positions.append(positions[label][0])
    return (
        np.asarray(panel_positions, dtype=np.int64),
        np.asarray(source_positions, dtype=np.int64),
    )


def load_compact_member(name: str, root: Path, panel: Sequence[str]) -> Matrix:
    perts = read_text_axis(root / "perts.txt")
    genes = read_text_axis(root / "genes.txt", require_unique=False)
    target_dest, target_source = panel_indices(perts, panel)
    gene_dest, gene_source = panel_indices(genes, panel)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    require(delta.shape == (len(perts), len(genes)), f"{name}: delta shape changed")
    effect = np.zeros((len(panel), len(panel)), dtype=np.float64)
    valid = np.zeros_like(effect, dtype=bool)
    if target_dest.size and gene_dest.size:
        values = np.asarray(
            delta[np.ix_(target_source, gene_source)],
            dtype=np.float64,
        )
        destination = np.ix_(target_dest, gene_dest)
        effect[destination] = np.where(np.isfinite(values), values, 0.0)
        valid[destination] = np.isfinite(values)
    diagonal = np.arange(len(panel))
    effect[diagonal, diagonal] = 0.0
    valid[diagonal, diagonal] = False
    return Matrix(name=name, effect=effect, valid=valid)


def cd4_target_qc(path: Path) -> dict[str, bool]:
    output: dict[str, bool] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            target = str(row["target_gene"])
            require(target not in output, f"{path}: duplicate target {target}")
            output[target] = (
                parse_bool(row["ontarget_significant"])
                and not parse_bool(row["neighboring_gene_KD"])
                and not parse_bool(row["distal_offtarget_flag"])
                and not parse_bool(row["low_target_gex"])
                and float(row["n_guides"]) >= 2.0
            )
    return output


def load_cd4_member(name: str, root: Path, panel: Sequence[str]) -> Matrix:
    perts = read_text_axis(root / "perts.txt")
    genes = read_text_axis(root / "genes.txt", require_unique=False)
    target_dest, target_source = panel_indices(perts, panel)
    gene_dest, gene_source = panel_indices(genes, panel)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    finite_flag = np.load(root / "valid_log_fc.npy", mmap_mode="r", allow_pickle=False)
    expected = (len(perts), len(genes))
    require(delta.shape == expected, f"{name}: delta shape changed")
    require(finite_flag.shape == expected, f"{name}: validity shape changed")
    effect = np.zeros((len(panel), len(panel)), dtype=np.float64)
    valid = np.zeros_like(effect, dtype=bool)
    values = np.asarray(delta[np.ix_(target_source, gene_source)], dtype=np.float64)
    flags = np.asarray(finite_flag[np.ix_(target_source, gene_source)], dtype=bool)
    qc_by_target = cd4_target_qc(root / "obs.csv")
    target_qc = np.asarray(
        [qc_by_target.get(panel[index], False) for index in target_dest],
        dtype=bool,
    )
    flags &= target_qc[:, None] & np.isfinite(values)
    destination = np.ix_(target_dest, gene_dest)
    effect[destination] = np.where(flags, values, 0.0)
    valid[destination] = flags
    diagonal = np.arange(len(panel))
    effect[diagonal, diagonal] = 0.0
    valid[diagonal, diagonal] = False
    return Matrix(name=name, effect=effect, valid=valid)


def nanmedian(values: np.ndarray, axis: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def aggregate_family(
    name: str,
    members: Sequence[Matrix],
    *,
    minimum_member_support: int,
) -> Matrix:
    require(bool(members), f"{name}: no members")
    effects = np.stack([member.effect for member in members])
    valid = np.stack([member.valid for member in members])
    support = valid.sum(axis=0)
    masked = np.where(valid, effects, np.nan)
    effect = nanmedian(masked, axis=0)
    family_valid = support >= minimum_member_support
    effect = np.where(family_valid & np.isfinite(effect), effect, 0.0)
    diagonal = np.arange(effect.shape[0])
    effect[diagonal, diagonal] = 0.0
    family_valid[diagonal, diagonal] = False
    require(np.all(np.isfinite(effect)), f"{name}: nonfinite family effect")
    require(np.all(effect[~family_valid] == 0.0), f"{name}: invalid effect nonzero")
    return Matrix(name=name, effect=effect, valid=family_valid)


def load_families(panel: Sequence[str]) -> tuple[Matrix, ...]:
    k562 = load_compact_member(
        "k562_gwps",
        ROOT / "data/donors/k562_gwps",
        panel,
    )
    hct116 = load_compact_member(
        "orion_hct116",
        ROOT / "data/donors/orion_hct116",
        panel,
    )
    hek293t = load_compact_member(
        "orion_hek293t",
        ROOT / "data/donors/orion_hek293t",
        panel,
    )
    cd4_root = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
    cd4_members = (
        load_cd4_member("cd4_rest", cd4_root / "cd4_rest", panel),
        load_cd4_member("cd4_stim8hr", cd4_root / "cd4_stim8hr", panel),
        load_cd4_member("cd4_stim48hr", cd4_root / "cd4_stim48hr", panel),
    )
    families = (
        aggregate_family("replogle_k562", (k562,), minimum_member_support=1),
        aggregate_family("orion", (hct116, hek293t), minimum_member_support=1),
        aggregate_family("cd4", cd4_members, minimum_member_support=2),
    )
    require(tuple(family.name for family in families) == FAMILY_NAMES, "family order drift")
    return families


def transform_family(
    family: Matrix,
    *,
    target_permutation: np.ndarray | None,
    response_permutation: np.ndarray | None,
) -> Matrix:
    n = family.effect.shape[0]
    rows = np.arange(n) if target_permutation is None else target_permutation
    columns = np.arange(n) if response_permutation is None else response_permutation
    effect = np.asarray(family.effect[np.ix_(rows, columns)], dtype=np.float64).copy()
    valid = np.asarray(family.valid[np.ix_(rows, columns)], dtype=bool).copy()
    diagonal = np.arange(n)
    effect[diagonal, diagonal] = 0.0
    valid[diagonal, diagonal] = False
    return Matrix(family.name, effect, valid)


def row_scales(family: Matrix) -> np.ndarray:
    scales = np.zeros(family.effect.shape[0], dtype=np.float64)
    for row in range(family.effect.shape[0]):
        values = np.abs(family.effect[row, family.valid[row]])
        values = values[np.isfinite(values) & (values > 0)]
        if values.size:
            scales[row] = float(np.median(values))
    return scales


def build_candidate(
    families: Sequence[Matrix],
    candidate: Candidate,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = families[0].effect.shape[0]
    if candidate.name == "identity":
        zeros = np.zeros((n, n), dtype=np.float64)
        support = np.stack([family.valid for family in families]).sum(axis=0)
        diagonal = np.arange(n)
        support[diagonal, diagonal] = 0
        return zeros.copy(), zeros.copy(), support.astype(np.int16)

    scales = np.stack([row_scales(family) for family in families])
    valid = np.stack([family.valid for family in families])
    valid &= scales[:, :, None] > 0
    normalized = np.zeros_like(valid, dtype=np.float64)
    for index, family in enumerate(families):
        np.divide(
            family.effect,
            scales[index, :, None],
            out=normalized[index],
            where=valid[index],
        )
    masked = np.where(valid, normalized, np.nan)
    support = valid.sum(axis=0).astype(np.int16)
    if candidate.location == "mean":
        numerator = np.where(valid, normalized, 0.0).sum(axis=0)
        location = np.divide(
            numerator,
            support,
            out=np.zeros((n, n), dtype=np.float64),
            where=support > 0,
        )
    elif candidate.location in {"median", "sign_consensus_median"}:
        location = nanmedian(masked, axis=0)
        location = np.where(np.isfinite(location), location, 0.0)
    else:
        raise ExperimentError(f"unknown location: {candidate.location}")

    nonzero = valid & (normalized != 0)
    nonzero_count = nonzero.sum(axis=0)
    sign_sum = np.where(nonzero, np.sign(normalized), 0.0).sum(axis=0)
    sign_agreement = np.divide(
        np.abs(sign_sum),
        nonzero_count,
        out=np.zeros((n, n), dtype=np.float64),
        where=nonzero_count > 0,
    )
    allowed = support >= candidate.minimum_family_support
    if candidate.location == "sign_consensus_median":
        allowed &= np.abs(sign_sum) > (0.5 * nonzero_count)
        allowed &= nonzero_count >= candidate.minimum_family_support
    location = np.where(allowed, location, 0.0)

    scale_values = np.broadcast_to(scales[:, :, None], valid.shape)
    output_scale = nanmedian(np.where(valid, scale_values, np.nan), axis=0)
    output_scale = np.where(np.isfinite(output_scale), output_scale, 0.0)
    residual = location * output_scale
    if candidate.confidence == "magnitude":
        confidence = np.abs(location)
    elif candidate.confidence == "agreement":
        confidence = (
            np.abs(location)
            * (support.astype(np.float64) / len(families))
            * sign_agreement
        )
    else:
        raise ExperimentError(f"unknown confidence: {candidate.confidence}")
    residual = np.where(allowed, residual, 0.0)
    confidence = np.where(allowed, confidence, 0.0)
    diagonal = np.arange(n)
    residual[diagonal, diagonal] = 0.0
    confidence[diagonal, diagonal] = 0.0
    support[diagonal, diagonal] = 0
    require(np.all(np.isfinite(residual)), f"{candidate.name}: nonfinite residual")
    require(np.all(np.isfinite(confidence)), f"{candidate.name}: nonfinite confidence")
    require(np.all(confidence >= 0), f"{candidate.name}: negative confidence")
    require(np.all(np.diag(residual) == 0.0), f"{candidate.name}: self effect")
    return residual, confidence, support


def per_family_metrics(
    prediction: np.ndarray,
    confidence: np.ndarray,
    support: np.ndarray,
    held: Matrix,
    panel: Sequence[str],
    *,
    include_targets: bool,
) -> dict[str, Any]:
    credits: list[float] = []
    precisions: list[float] = []
    target_records: list[dict[str, Any]] = []
    n = len(panel)
    canonical_columns = np.arange(n)
    for target in range(n):
        held_valid = held.valid[target].copy()
        held_valid[target] = False
        indices = np.flatnonzero(held_valid)
        if indices.size < 20 or not np.any(support[target] > 0):
            continue
        truth = held.effect[target, indices]
        denominator = float(np.sum(np.abs(truth)))
        if not np.isfinite(denominator) or denominator <= 0:
            continue
        pred = prediction[target, indices]
        credit = float(np.sum(np.sign(truth) * pred) / denominator)
        top_n = max(1, int(math.ceil(0.20 * indices.size)))
        order = np.lexsort(
            (canonical_columns[indices], -confidence[target, indices])
        )[:top_n]
        chosen = indices[order]
        weights = np.abs(held.effect[target, chosen])
        weight_sum = float(np.sum(weights))
        if not np.isfinite(weight_sum) or weight_sum <= 0:
            continue
        predicted_sign = np.sign(prediction[target, chosen])
        truth_sign = np.sign(held.effect[target, chosen])
        correct = (
            (predicted_sign != 0)
            & (truth_sign != 0)
            & (predicted_sign == truth_sign)
        )
        precision = float(np.sum(weights * correct) / weight_sum)
        credits.append(credit)
        precisions.append(precision)
        if include_targets:
            target_records.append(
                {
                    "target": panel[target],
                    "held_valid_responses": int(indices.size),
                    "training_supported_responses": int(
                        np.count_nonzero(support[target, indices])
                    ),
                    "signed_credit": credit,
                    "top20_weighted_sign_precision": precision,
                }
            )
    require(credits, f"{held.name}: no evaluable targets")
    return {
        "family": held.name,
        "target_count": len(credits),
        "signed_credit": float(np.mean(credits)),
        "top20_weighted_sign_precision": float(np.mean(precisions)),
        "targets": target_records if include_targets else None,
    }


def evaluate_grid(
    families: Sequence[Matrix],
    panel: Sequence[str],
    *,
    target_permutation: np.ndarray | None = None,
    response_permutation: np.ndarray | None = None,
    include_targets: bool = False,
) -> tuple[Candidate, list[dict[str, Any]]]:
    metrics: dict[str, list[dict[str, Any]]] = {
        candidate.name: [] for candidate in CANDIDATES
    }
    for held_index, held in enumerate(families):
        training = [
            transform_family(
                family,
                target_permutation=target_permutation,
                response_permutation=response_permutation,
            )
            for index, family in enumerate(families)
            if index != held_index
        ]
        for candidate in CANDIDATES:
            prediction, confidence, support = build_candidate(training, candidate)
            metrics[candidate.name].append(
                per_family_metrics(
                    prediction,
                    confidence,
                    support,
                    held,
                    panel,
                    include_targets=include_targets,
                )
            )

    summaries: list[dict[str, Any]] = []
    for candidate in CANDIDATES:
        folds = metrics[candidate.name]
        credits = [float(fold["signed_credit"]) for fold in folds]
        precisions = [
            float(fold["top20_weighted_sign_precision"]) for fold in folds
        ]
        counts = [int(fold["target_count"]) for fold in folds]
        eligible = candidate.name != "identity" and (
            all(count >= 20 for count in counts)
            and all(credit > 0.0 for credit in credits)
            and all(precision >= 0.5 for precision in precisions)
        )
        summaries.append(
            {
                "name": candidate.name,
                "eligible": eligible,
                "worst_family_signed_credit": min(credits),
                "worst_family_top20_weighted_sign_precision": min(precisions),
                "macro_signed_credit": float(np.mean(credits)),
                "macro_top20_weighted_sign_precision": float(np.mean(precisions)),
                "folds": folds,
            }
        )

    eligible = [summary for summary in summaries if summary["eligible"]]
    if not eligible:
        selected = CANDIDATES[0]
    else:
        order = {candidate.name: index for index, candidate in enumerate(CANDIDATES)}
        best = max(
            eligible,
            key=lambda summary: (
                float(summary["worst_family_signed_credit"]),
                float(summary["worst_family_top20_weighted_sign_precision"]),
                float(summary["macro_signed_credit"]),
                float(summary["macro_top20_weighted_sign_precision"]),
                -order[str(summary["name"])],
            ),
        )
        selected = next(
            candidate for candidate in CANDIDATES if candidate.name == best["name"]
        )
    return selected, summaries


def run_null_family(
    families: Sequence[Matrix],
    panel: Sequence[str],
    bank: np.ndarray,
    kind: str,
) -> dict[str, Any]:
    selected_names: list[str] = []
    worst_credit = np.empty(len(bank), dtype=np.float64)
    worst_precision = np.empty(len(bank), dtype=np.float64)
    for draw, permutation in enumerate(bank):
        kwargs: dict[str, np.ndarray] = {}
        if kind == "target_label":
            kwargs["target_permutation"] = permutation
        elif kind == "response_label":
            kwargs["response_permutation"] = permutation
        else:
            raise ExperimentError(f"unknown null family: {kind}")
        selected, summaries = evaluate_grid(
            families,
            panel,
            include_targets=False,
            **kwargs,
        )
        selected_names.append(selected.name)
        summary = next(item for item in summaries if item["name"] == selected.name)
        worst_credit[draw] = float(summary["worst_family_signed_credit"])
        worst_precision[draw] = float(
            summary["worst_family_top20_weighted_sign_precision"]
        )
        if (draw + 1) % 25 == 0:
            print(f"  {kind}: {draw + 1}/{len(bank)}", flush=True)
    return {
        "family": kind,
        "selected_candidate_names": selected_names,
        "selected_worst_family_signed_credit": worst_credit,
        "selected_worst_family_top20_weighted_sign_precision": worst_precision,
        "p95_worst_family_signed_credit": float(
            np.quantile(worst_credit, 0.95, method="linear")
        ),
        "p95_worst_family_top20_weighted_sign_precision": float(
            np.quantile(worst_precision, 0.95, method="linear")
        ),
    }


def expand_panel_matrix(
    panel_values: np.ndarray,
    panel: Sequence[str],
    genes: Sequence[str],
) -> np.ndarray:
    require(panel_values.shape == (len(panel), len(panel)), "panel matrix shape changed")
    gene_position = {gene: index for index, gene in enumerate(genes)}
    columns = np.asarray([gene_position[gene] for gene in panel], dtype=np.int64)
    output = np.zeros((len(panel), len(genes)), dtype=np.float32)
    output[:, columns] = panel_values.astype(np.float32)
    return output


def invariant_report(
    panel_residual: np.ndarray,
    panel_confidence: np.ndarray,
    full_residual: np.ndarray,
    full_confidence: np.ndarray,
    panel: Sequence[str],
    genes: Sequence[str],
) -> dict[str, Any]:
    gene_position = {gene: index for index, gene in enumerate(genes)}
    panel_columns = np.asarray([gene_position[gene] for gene in panel], dtype=np.int64)
    nonpanel = np.ones(len(genes), dtype=bool)
    nonpanel[panel_columns] = False
    diagonal = np.arange(len(panel))
    checks = {
        "panel_residual_shape": panel_residual.shape == (len(panel), len(panel)),
        "panel_confidence_shape": panel_confidence.shape == (len(panel), len(panel)),
        "target_self_residual_exact_zero": bool(
            np.all(panel_residual[diagonal, diagonal] == 0.0)
        ),
        "target_self_confidence_exact_zero": bool(
            np.all(panel_confidence[diagonal, diagonal] == 0.0)
        ),
        "expanded_nonpanel_residual_exact_zero": bool(
            np.all(full_residual[:, nonpanel] == 0.0)
        ),
        "expanded_nonpanel_confidence_exact_zero": bool(
            np.all(full_confidence[:, nonpanel] == 0.0)
        ),
        "expanded_panel_residual_exact": bool(
            np.array_equal(full_residual[:, panel_columns], panel_residual.astype(np.float32))
        ),
        "expanded_panel_confidence_exact": bool(
            np.array_equal(
                full_confidence[:, panel_columns],
                panel_confidence.astype(np.float32),
            )
        ),
        "finite": bool(
            np.all(np.isfinite(panel_residual))
            and np.all(np.isfinite(panel_confidence))
        ),
    }
    checks["all_pass"] = all(checks.values())
    return checks


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(resolved != HERE and resolved.is_relative_to(HERE), "unsafe output path")
    return resolved


def write_npy(path: Path, value: np.ndarray) -> None:
    with path.open("xb") as handle:
        np.save(handle, np.ascontiguousarray(value), allow_pickle=False)


def write_selected_per_target(
    path: Path,
    selected_summary: dict[str, Any],
) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "held_family",
                "target",
                "held_valid_responses",
                "training_supported_responses",
                "signed_credit",
                "top20_weighted_sign_precision",
            )
        )
        for fold in selected_summary["folds"]:
            for row in fold["targets"]:
                writer.writerow(
                    (
                        fold["family"],
                        row["target"],
                        row["held_valid_responses"],
                        row["training_supported_responses"],
                        f"{row['signed_credit']:.17g}",
                        f"{row['top20_weighted_sign_precision']:.17g}",
                    )
                )


def runtime_versions() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }


def run(output: Path = DEFAULT_OUTPUT) -> Path:
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite {output}")
    print("Authenticating frozen inputs...", flush=True)
    authorization = validate_authorization(hash_large_inputs=True)
    genes, panel = panel_axes()
    print("Loading panel slices from authorized source effects...", flush=True)
    families = load_families(panel)
    print("Evaluating preregistered LOFO candidate family...", flush=True)
    selected, summaries = evaluate_grid(
        families,
        panel,
        include_targets=True,
    )
    selected_summary = next(item for item in summaries if item["name"] == selected.name)

    freeze_manifest = freeze_plans.verify(PLAN_DIR)
    target_bank = np.load(
        PLAN_DIR / "target-label-derangements.npy",
        allow_pickle=False,
    )
    response_bank = np.load(
        PLAN_DIR / "response-label-derangements.npy",
        allow_pickle=False,
    )
    print("Running complete-selection target-label null...", flush=True)
    target_null = run_null_family(families, panel, target_bank, "target_label")
    print("Running complete-selection response-label null...", flush=True)
    response_null = run_null_family(families, panel, response_bank, "response_label")

    all_residuals = []
    all_confidences = []
    all_support = []
    for candidate in CANDIDATES:
        residual, confidence, support = build_candidate(families, candidate)
        all_residuals.append(residual)
        all_confidences.append(confidence)
        all_support.append(support)
    residuals = np.stack(all_residuals)
    confidences = np.stack(all_confidences)
    supports = np.stack(all_support)
    selected_index = [candidate.name for candidate in CANDIDATES].index(selected.name)
    panel_residual = residuals[selected_index]
    panel_confidence = confidences[selected_index]
    full_residual = expand_panel_matrix(panel_residual, panel, genes)
    full_confidence = expand_panel_matrix(panel_confidence, panel, genes)
    invariants = invariant_report(
        panel_residual,
        panel_confidence,
        full_residual,
        full_confidence,
        panel,
        genes,
    )
    require(invariants["all_pass"], "PDS-null invariant failed")

    observed_credit = float(selected_summary["worst_family_signed_credit"])
    target_p95 = float(target_null["p95_worst_family_signed_credit"])
    response_p95 = float(response_null["p95_worst_family_signed_credit"])
    selected_nonidentity = selected.name != "identity"
    all_gates = bool(
        selected_nonidentity
        and selected_summary["eligible"]
        and observed_credit > target_p95
        and observed_credit > response_p95
        and invariants["all_pass"]
    )
    results = {
        "schema_version": 1,
        "experiment": "122-panel-response-atlas",
        "phase": "phase_1_source_only",
        "selected_candidate": selected.name,
        "candidate_order": [candidate.name for candidate in CANDIDATES],
        "candidate_summaries": summaries,
        "family_order": list(FAMILY_NAMES),
        "family_coverage": {
            family.name: {
                "valid_target_rows": int(np.count_nonzero(np.any(family.valid, axis=1))),
                "valid_panel_pairs": int(np.count_nonzero(family.valid)),
            }
            for family in families
        },
        "selected_output": {
            "nonzero_panel_pairs": int(np.count_nonzero(panel_residual)),
            "supported_target_rows": int(np.count_nonzero(np.any(panel_residual != 0, axis=1))),
            "panel_residual_array_sha256": sha256_array(panel_residual.astype(np.float32)),
            "panel_confidence_array_sha256": sha256_array(
                panel_confidence.astype(np.float32)
            ),
        },
        "complete_selection_nulls": {
            "target_label": {
                "draws": len(target_bank),
                "p95_worst_family_signed_credit": target_p95,
                "p95_worst_family_top20_weighted_sign_precision": target_null[
                    "p95_worst_family_top20_weighted_sign_precision"
                ],
                "p_value_worst_family_signed_credit": float(
                    (
                        1
                        + np.count_nonzero(
                            target_null["selected_worst_family_signed_credit"]
                            >= observed_credit
                        )
                    )
                    / (len(target_bank) + 1)
                ),
            },
            "response_label": {
                "draws": len(response_bank),
                "p95_worst_family_signed_credit": response_p95,
                "p95_worst_family_top20_weighted_sign_precision": response_null[
                    "p95_worst_family_top20_weighted_sign_precision"
                ],
                "p_value_worst_family_signed_credit": float(
                    (
                        1
                        + np.count_nonzero(
                            response_null["selected_worst_family_signed_credit"]
                            >= observed_credit
                        )
                    )
                    / (len(response_bank) + 1)
                ),
            },
        },
        "pds_null_invariants": invariants,
        "future_dld1_confirmation": {
            "earned": all_gates,
            "dld1_values_accessed": False,
            "requirements": {
                "selected_candidate_nonidentity": selected_nonidentity,
                "all_three_family_eligibility_gates_pass": bool(
                    selected_summary["eligible"]
                ),
                "observed_credit_above_target_null_p95": observed_credit > target_p95,
                "observed_credit_above_response_null_p95": observed_credit > response_p95,
                "exact_pds_null_invariants": bool(invariants["all_pass"]),
            },
        },
        "source_reconstruction": read_json(INVENTORY_PATH)["coverage_summary"],
        "safety": {
            "a549_accessed": False,
            "pan_h1_accessed": False,
            "dld1_values_accessed": False,
            "held_harness_or_outcomes_accessed": False,
            "challenge_truth_accessed": False,
            "cell_eval2_run": False,
            "counts_generated": False,
            "challenge_build_created": False,
            "submission_created_or_sent": False,
            "commit_or_push": False,
        },
        "runtime_versions": runtime_versions(),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.parent / f".{output.name}.staging"
    require(not stage.exists(), f"staging directory exists: {stage}")
    stage.mkdir()
    try:
        write_npy(stage / "family-effects.npy", np.stack([f.effect for f in families]).astype(np.float32))
        write_npy(stage / "family-valid.npy", np.stack([f.valid for f in families]))
        write_npy(stage / "all-candidate-panel-residuals.npy", residuals.astype(np.float32))
        write_npy(
            stage / "all-candidate-panel-confidence.npy",
            confidences.astype(np.float32),
        )
        write_npy(stage / "all-candidate-panel-support.npy", supports)
        write_npy(stage / "selected-panel-residual.npy", panel_residual.astype(np.float32))
        write_npy(
            stage / "selected-panel-confidence.npy",
            panel_confidence.astype(np.float32),
        )
        write_npy(stage / "selected-full-residual.npy", full_residual)
        write_npy(stage / "selected-full-confidence.npy", full_confidence)
        write_npy(
            stage / "target-null-worst-credit.npy",
            target_null["selected_worst_family_signed_credit"],
        )
        write_npy(
            stage / "target-null-worst-precision.npy",
            target_null["selected_worst_family_top20_weighted_sign_precision"],
        )
        write_npy(
            stage / "response-null-worst-credit.npy",
            response_null["selected_worst_family_signed_credit"],
        )
        write_npy(
            stage / "response-null-worst-precision.npy",
            response_null["selected_worst_family_top20_weighted_sign_precision"],
        )
        (stage / "candidate-names.txt").write_text(
            "".join(f"{candidate.name}\n" for candidate in CANDIDATES),
            encoding="utf-8",
        )
        (stage / "family-names.txt").write_text(
            "".join(f"{name}\n" for name in FAMILY_NAMES),
            encoding="utf-8",
        )
        (stage / "panel-targets.txt").write_text(
            "".join(f"{target}\n" for target in panel),
            encoding="utf-8",
        )
        write_selected_per_target(stage / "selected-lofo-per-target.tsv", selected_summary)
        (stage / "results.json").write_bytes(json_bytes(results))
        files = {}
        for path in sorted(stage.iterdir()):
            files[path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        manifest = {
            "schema_version": 1,
            "experiment": "122-panel-response-atlas",
            "phase": "phase_1_source_only",
            "selected_candidate": selected.name,
            "source_authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "implementation_freeze_sha256": sha256_file(IMPLEMENTATION_FREEZE_PATH),
            "plan_freeze_manifest_sha256": sha256_file(
                PLAN_DIR / "freeze-manifest.json"
            ),
            "allowed_input_count": len(authorization["allowed_input_files"]),
            "files": files,
        }
        (stage / "artifact-manifest.json").write_bytes(json_bytes(manifest))
        receipt = {
            "schema_version": 1,
            "experiment": "122-panel-response-atlas",
            "phase": "phase_1_source_only",
            "status": "immutable_source_only_result",
            "artifact_manifest_sha256": sha256_file(stage / "artifact-manifest.json"),
            "selected_candidate": selected.name,
            "future_dld1_confirmation_earned": all_gates,
            "dld1_values_accessed": False,
            "challenge_truth_accessed": False,
            "counts_or_submission_built": False,
        }
        (stage / "result-receipt.json").write_bytes(json_bytes(receipt))
        (stage / "result-receipt.sha256").write_text(
            f"{sha256_file(stage / 'result-receipt.json')}  result-receipt.json\n",
            encoding="utf-8",
        )
        os.sync()
        stage.rename(output)
        for path in output.iterdir():
            path.chmod(0o444)
        output.chmod(0o555)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def audit() -> dict[str, Any]:
    authorization = validate_authorization(hash_large_inputs=False)
    freeze_manifest = freeze_plans.verify(PLAN_DIR)
    return {
        "authorized": True,
        "allowed_input_count": len(authorization["allowed_input_files"]),
        "candidate_names": [candidate.name for candidate in CANDIDATES],
        "family_names": list(FAMILY_NAMES),
        "null_draws_each": freeze_manifest["draws_per_null_family"],
        "output_exists": DEFAULT_OUTPUT.exists(),
        "runtime_versions": runtime_versions(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.audit:
        print(json.dumps(audit(), indent=2, sort_keys=True))
    else:
        print(json.dumps({"output": str(run(args.output))}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
