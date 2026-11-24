#!/usr/bin/env python3
"""Audit source-unit-correct panel transport with an exact integer budget."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
from typing import Any, Sequence

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
FREEZE_PATH = HERE / "IMPLEMENTATION_FREEZE.json"


class AuditError(RuntimeError):
    """Raised when a frozen contract or mechanical invariant fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_int_rows(rows: Sequence[np.ndarray]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(np.ascontiguousarray(row, dtype="<i8").tobytes())
    return digest.hexdigest()


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(str(row[column]) for row in reader)
    require(len(values) == len(set(values)), f"{path}: duplicate labels")
    return values


def read_text_axis(path: Path) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    require(len(values) == len(set(values)), f"{path}: duplicate labels")
    return values


def decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        for value in values
    )


def read_h5ad_index(node: h5py.Dataset | h5py.Group) -> tuple[str, ...]:
    if isinstance(node, h5py.Dataset):
        return decode_strings(np.asarray(node))
    require(
        str(node.attrs.get("encoding-type", "")) == "nullable-string-array",
        "unsupported H5AD index encoding",
    )
    require(set(node.keys()) == {"mask", "values"}, "nullable string fields changed")
    mask = np.asarray(node["mask"], dtype=bool)
    require(not bool(np.any(mask)), "H5AD gene index contains missing values")
    return decode_strings(np.asarray(node["values"]))


def protocol_input_items(protocol: dict[str, Any]) -> list[dict[str, Any]]:
    inputs = protocol["inputs"]
    return [
        *protocol["prior_failed_evidence"].values(),
        inputs["panel_residual"],
        inputs["experiment_122_manifest"],
        inputs["experiment_122_panel_axis"],
        inputs["experiment_122_receipt"],
        inputs["staged_a15_effects"],
        inputs["staged_a15_metadata"],
        inputs["a15_submission_metadata"],
        inputs["production_profile_script"],
        inputs["production_emitter"],
        inputs["gene_names"],
        inputs["panel_targets"],
        *inputs["released_controls"].values(),
    ]


def verify_freeze() -> None:
    require(FREEZE_PATH.is_file(), "implementation freeze is absent")
    freeze = json.loads(FREEZE_PATH.read_text(encoding="utf-8"))
    require(
        freeze["status"] == "frozen_before_experiment_128_released_control_value_access",
        "implementation freeze status changed",
    )
    for relative, expected in freeze["files"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing frozen file: {path}")
        require(sha256_file(path) == expected, f"implementation hash changed: {path}")


def verify_inputs(protocol: dict[str, Any]) -> None:
    for item in protocol_input_items(protocol):
        path = ROOT / item["path"]
        require(path.is_file(), f"missing input: {path}")
        require(sha256_file(path) == item["sha256"], f"input hash changed: {path}")


def verify_provenance(
    protocol: dict[str, Any],
    genes: tuple[str, ...],
    panel: tuple[str, ...],
) -> dict[str, Any]:
    inputs = protocol["inputs"]
    source_manifest = json.loads(
        (ROOT / inputs["experiment_122_manifest"]["path"]).read_text()
    )
    source_receipt = json.loads(
        (ROOT / inputs["experiment_122_receipt"]["path"]).read_text()
    )
    source_axis = read_text_axis(ROOT / inputs["experiment_122_panel_axis"]["path"])
    residual_hash = inputs["panel_residual"]["sha256"]
    require(source_axis == panel, "Experiment 122 panel axis differs from canonical panel")
    require(
        source_manifest["files"]["selected-panel-residual.npy"]["sha256"]
        == residual_hash,
        "Experiment 122 manifest does not bind the selected residual",
    )
    require(
        source_receipt["artifact_manifest_sha256"]
        == inputs["experiment_122_manifest"]["sha256"],
        "Experiment 122 receipt does not bind its artifact manifest",
    )

    staged_meta = json.loads(
        (ROOT / inputs["staged_a15_metadata"]["path"]).read_text()
    )
    a15_meta = json.loads(
        (ROOT / inputs["a15_submission_metadata"]["path"]).read_text()
    )
    require(
        staged_meta["file_sha256"] == inputs["staged_a15_effects"]["sha256"],
        "staged a15 metadata does not bind the effect file",
    )
    require(staged_meta["shape"] == [len(panel), len(genes)], "a15 shape metadata changed")
    config = a15_meta["config"]
    production = protocol["production_semantics"]
    require(float(config["min_cpm"]) == float(production["min_cpm"]), "a15 CPM gate changed")
    require(
        float(config["effect_scale"]) == float(production["effect_scale"]),
        "a15 effect scale changed",
    )
    require(
        int(config["cells_per_target_context"])
        == int(production["cells_per_perturbation"]),
        "a15 cell count changed",
    )
    implementation = a15_meta["hashes"]["implementation"]
    require(
        implementation["scripts/make_submission_pm_cd4.py"]
        == inputs["production_profile_script"]["sha256"],
        "a15 metadata does not bind the production profile script",
    )
    require(
        implementation["src/vcc2026/donor_transfer.py"]
        == inputs["production_emitter"]["sha256"],
        "a15 metadata does not bind the production emitter",
    )
    return {
        "experiment_122_axis_matches_canonical": True,
        "experiment_122_manifest_binds_residual": True,
        "experiment_122_receipt_binds_manifest": True,
        "staged_a15_metadata_binds_effects": True,
        "a15_metadata_binds_production_implementation": True,
    }


def control_profile_and_libraries(
    path: Path,
    expected_genes: tuple[str, ...],
    *,
    ts: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Read a CSR h5ad in bounded row chunks."""
    with h5py.File(path, "r") as handle:
        require("X" in handle, f"{path}: X missing")
        matrix = handle["X"]
        require(
            str(matrix.attrs.get("encoding-type", "")) == "csr_matrix",
            f"{path}: X is not CSR",
        )
        shape = tuple(map(int, matrix.attrs["shape"]))
        require(shape[1] == len(expected_genes), f"{path}: gene dimension changed")
        require(
            read_h5ad_index(handle["var"]["_index"]) == expected_genes,
            f"{path}: gene order changed",
        )
        data = matrix["data"]
        indices = matrix["indices"]
        indptr = np.asarray(matrix["indptr"], dtype=np.int64)
        require(indptr.shape == (shape[0] + 1,), f"{path}: CSR indptr changed")
        gene_totals = np.zeros(shape[1], dtype=np.float64)
        libraries = np.zeros(shape[0], dtype=np.float64)
        row_chunk = 1024
        for row_start in range(0, shape[0], row_chunk):
            row_stop = min(row_start + row_chunk, shape[0])
            value_start = int(indptr[row_start])
            value_stop = int(indptr[row_stop])
            values = np.asarray(data[value_start:value_stop], dtype=np.float64)
            columns = np.asarray(indices[value_start:value_stop], dtype=np.int64)
            require(bool(np.all(np.isfinite(values))), f"{path}: nonfinite count")
            require(bool(np.all(values >= 0.0)), f"{path}: negative count")
            gene_totals += np.bincount(
                columns,
                weights=values,
                minlength=shape[1],
            )
            local = indptr[row_start : row_stop + 1] - value_start
            nonempty = local[1:] > local[:-1]
            if np.any(nonempty):
                starts = local[:-1][nonempty]
                libraries[row_start:row_stop][nonempty] = np.add.reduceat(
                    values,
                    starts,
                )
    require(float(gene_totals.sum()) > 0.0, f"{path}: zero total")
    require(bool(np.all(libraries > 0.0)), f"{path}: empty control cell")
    require(
        bool(np.all(libraries == np.floor(libraries))),
        f"{path}: fractional control counts",
    )
    b_ctrl = np.log1p(ts * gene_totals / gene_totals.sum())
    return b_ctrl, libraries.astype(np.int64)


def profile_from_effect(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    *,
    ts: float,
    min_cpm: float,
    effect_scale: float,
) -> np.ndarray:
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > min_cpm) & (b_ctrl > 0.0)
    gated = np.where(gate, effect, 0.0)
    if effect_scale != 1.0:
        gated = gated * effect_scale
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / ts
    total = float(profile.sum())
    require(total > 0.0, "degenerate baseline profile")
    return profile / total


def mutable_panel_mask(
    baseline_profile: np.ndarray,
    b_ctrl: np.ndarray,
    residual_row: np.ndarray,
    panel_columns: np.ndarray,
    target_index: int,
    *,
    min_cpm: float,
) -> np.ndarray:
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    mutable = (
        (cpm[panel_columns] > min_cpm)
        & (b_ctrl[panel_columns] > 0.0)
        & (baseline_profile[panel_columns] > 0.0)
        & np.isfinite(residual_row)
    )
    mutable[target_index] = False
    return mutable


def source_coordinate_transform(
    baseline_panel: np.ndarray,
    residual_row: np.ndarray,
    mutable: np.ndarray,
    eta: float,
    *,
    ts: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return clipped source-unit weights and the preclip-negative mask."""
    weights = np.asarray(baseline_panel, dtype=np.float64).copy()
    coordinate = np.log1p(ts * weights[mutable]) + eta * residual_row[mutable]
    with np.errstate(over="raise", invalid="raise"):
        raw = np.expm1(coordinate) / ts
    require(bool(np.all(np.isfinite(raw))), "source transform produced nonfinite values")
    preclip_negative = raw < 0.0
    weights[mutable] = np.clip(raw, 0.0, None)
    return weights, preclip_negative


def largest_remainder(
    budget: int,
    weights: np.ndarray,
    canonical_order: np.ndarray,
    *,
    fallback: np.ndarray,
) -> tuple[np.ndarray, bool]:
    """Allocate an integer budget with canonical-order tie breaking."""
    values = np.asarray(weights, dtype=np.float64)
    order_key = np.asarray(canonical_order, dtype=np.int64)
    fallback_values = np.asarray(fallback, dtype=np.int64)
    require(values.ndim == 1, "largest-remainder weights must be one-dimensional")
    require(order_key.shape == values.shape, "canonical-order shape changed")
    require(fallback_values.shape == values.shape, "fallback shape changed")
    require(budget >= 0, "negative mutable budget")
    require(bool(np.all(np.isfinite(values))), "nonfinite allocation weight")
    require(bool(np.all(values >= 0.0)), "negative allocation weight")
    require(int(fallback_values.sum(dtype=np.int64)) == budget, "fallback budget changed")
    if budget == 0:
        return np.zeros(values.size, dtype=np.int64), False
    weight_sum = math.fsum(map(float, values))
    if weight_sum <= 0.0:
        return fallback_values.copy(), True
    ideal = budget * values / weight_sum
    allocated = np.floor(ideal).astype(np.int64)
    remaining = budget - int(allocated.sum(dtype=np.int64))
    require(0 <= remaining <= values.size, "largest-remainder residual is invalid")
    if remaining:
        fractional = ideal - allocated
        ranking = np.lexsort((order_key, -fractional))
        allocated[ranking[:remaining]] += 1
    require(int(allocated.sum(dtype=np.int64)) == budget, "integer budget not conserved")
    require(bool(np.all(allocated >= 0)), "integer allocation became negative")
    return allocated, False


def integer_budget_transport(
    baseline_totals: np.ndarray,
    baseline_profile: np.ndarray,
    b_ctrl: np.ndarray,
    residual_row: np.ndarray,
    panel_columns: np.ndarray,
    panel_mask: np.ndarray,
    target_index: int,
    eta: float,
    *,
    ts: float,
    min_cpm: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    baseline = np.asarray(baseline_totals, dtype=np.int64)
    candidate = baseline.copy()
    mutable = mutable_panel_mask(
        baseline_profile,
        b_ctrl,
        residual_row,
        panel_columns,
        target_index,
        min_cpm=min_cpm,
    )
    panel_profile = baseline_profile[panel_columns]
    weights, preclip_negative = source_coordinate_transform(
        panel_profile,
        residual_row,
        mutable,
        eta,
        ts=ts,
    )
    mutable_columns = panel_columns[mutable]
    fallback_used = False
    if mutable_columns.size:
        mutable_baseline = baseline[mutable_columns]
        budget = int(mutable_baseline.sum(dtype=np.int64))
        allocation, fallback_used = largest_remainder(
            budget,
            weights[mutable],
            np.flatnonzero(mutable),
            fallback=mutable_baseline,
        )
        candidate[mutable_columns] = allocation
    else:
        budget = 0

    nonpanel = ~panel_mask
    immutable_panel = ~mutable
    baseline_total = int(baseline.sum(dtype=np.int64))
    candidate_total = int(candidate.sum(dtype=np.int64))
    baseline_panel_total = int(baseline[panel_mask].sum(dtype=np.int64))
    candidate_panel_total = int(candidate[panel_mask].sum(dtype=np.int64))
    baseline_linear = ts * baseline.astype(np.float64) / baseline_total
    candidate_linear = ts * candidate.astype(np.float64) / candidate_total
    baseline_log = np.log1p(baseline_linear)
    candidate_log = np.log1p(candidate_linear)
    baseline_effect = baseline_log - b_ctrl
    candidate_effect = candidate_log - b_ctrl
    target_column = int(panel_columns[target_index])

    invariants = {
        "candidate_int64": candidate.dtype == np.dtype(np.int64),
        "candidate_nonnegative": bool(np.all(candidate >= 0)),
        "exact_nonpanel_integer_totals": bool(
            np.array_equal(candidate[nonpanel], baseline[nonpanel])
        ),
        "exact_target_self_integer_total": bool(
            candidate[target_column] == baseline[target_column]
        ),
        "exact_ineligible_panel_integer_totals": bool(
            np.array_equal(
                candidate[panel_columns[immutable_panel]],
                baseline[panel_columns[immutable_panel]],
            )
        ),
        "exact_mutable_integer_budget": bool(
            int(candidate[mutable_columns].sum(dtype=np.int64)) == budget
        ),
        "exact_panel_integer_total": candidate_panel_total == baseline_panel_total,
        "exact_grand_integer_total": candidate_total == baseline_total,
        "exact_nonpanel_scorer_linear": bool(
            np.array_equal(candidate_linear[nonpanel], baseline_linear[nonpanel])
        ),
        "exact_nonpanel_scorer_log1p": bool(
            np.array_equal(candidate_log[nonpanel], baseline_log[nonpanel])
        ),
        "exact_nonpanel_scorer_effect": bool(
            np.array_equal(candidate_effect[nonpanel], baseline_effect[nonpanel])
        ),
        "exact_target_self_scorer_linear": bool(
            candidate_linear[target_column] == baseline_linear[target_column]
        ),
        "exact_target_self_scorer_log1p": bool(
            candidate_log[target_column] == baseline_log[target_column]
        ),
    }

    nonzero_residual = mutable & (residual_row != 0.0)
    sign_columns = panel_columns[nonzero_residual]
    if sign_columns.size:
        observed_sign = np.sign(
            candidate_log[sign_columns] - baseline_log[sign_columns]
        )
        expected_sign = np.sign(residual_row[nonzero_residual])
        sign_retention = float(np.mean(observed_sign == expected_sign))
    else:
        sign_retention = 1.0
    continuous_sign = np.sign(
        weights[nonzero_residual] - panel_profile[nonzero_residual]
    )
    expected_continuous_sign = np.sign(residual_row[nonzero_residual])
    continuous_sign_retention = (
        float(np.mean(continuous_sign == expected_continuous_sign))
        if expected_continuous_sign.size
        else 1.0
    )
    return candidate, {
        "all_invariants_pass": bool(all(invariants.values())),
        "invariants": invariants,
        "mutable_coordinates": int(mutable_columns.size),
        "mutable_integer_budget": budget,
        "preclip_negative_coordinates": int(np.sum(preclip_negative)),
        "preclip_negative_row": bool(np.any(preclip_negative)),
        "continuous_source_sign_retention": continuous_sign_retention,
        "integer_transformed_sign_retention": sign_retention,
        "fallback_used": fallback_used,
        "changed_mutable_coordinates": int(
            np.count_nonzero(candidate[mutable_columns] != baseline[mutable_columns])
        ),
        "mutable_l1_count_change": int(
            np.abs(
                candidate[mutable_columns].astype(np.int64)
                - baseline[mutable_columns].astype(np.int64)
            ).sum(dtype=np.int64)
        ),
    }


def summarize(values: Sequence[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
    }


def run(output: Path) -> dict[str, Any]:
    verify_freeze()
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    verify_inputs(protocol)

    inputs = protocol["inputs"]
    production = protocol["production_semantics"]
    ts = float(production["ts"])
    min_cpm = float(production["min_cpm"])
    effect_scale = float(production["effect_scale"])
    cells = int(production["cells_per_perturbation"])
    genes = read_csv_axis(ROOT / inputs["gene_names"]["path"], "gene_name")
    panel = read_csv_axis(ROOT / inputs["panel_targets"]["path"], "target_gene")
    require(len(genes) == 18_533, "official gene count changed")
    require(len(panel) == 300, "official panel count changed")
    require(set(panel) <= set(genes), "panel is not a subset of the gene axis")
    provenance = verify_provenance(protocol, genes, panel)

    gene_position = {gene: index for index, gene in enumerate(genes)}
    panel_columns = np.asarray([gene_position[gene] for gene in panel], dtype=np.int64)
    panel_mask = np.zeros(len(genes), dtype=bool)
    panel_mask[panel_columns] = True
    effects = np.load(
        ROOT / inputs["staged_a15_effects"]["path"],
        allow_pickle=False,
    ).astype(np.float64)
    residual = np.load(
        ROOT / inputs["panel_residual"]["path"],
        allow_pickle=False,
    ).astype(np.float64)
    require(effects.shape == (len(panel), len(genes)), "a15 effect shape changed")
    require(residual.shape == (len(panel), len(panel)), "panel residual shape changed")
    require(bool(np.all(np.isfinite(residual))), "panel residual became nonfinite")
    require(bool(np.all(np.diag(residual) == 0.0)), "panel residual diagonal changed")

    contexts = tuple(production["context_order"])
    doses = tuple(map(float, protocol["fixed_doses"]))
    fixed_libraries = tuple(map(int, production["representative_libraries"]))
    scheme_names = tuple(f"fixed_L_{value}" for value in fixed_libraries) + (
        "actual_production_pick_L",
    )
    actual_libraries = np.zeros((len(contexts), len(panel)), dtype=np.int64)
    sign_retention = np.zeros(
        (len(contexts), len(doses), len(scheme_names), len(panel)),
        dtype=np.float64,
    )
    preclip_negative = np.zeros(
        (len(contexts), len(doses), len(panel)),
        dtype=np.int16,
    )
    invariant_pass = np.zeros(
        (len(contexts), len(doses), len(scheme_names), len(panel)),
        dtype=bool,
    )
    mutable_counts = np.zeros((len(contexts), len(panel)), dtype=np.int16)
    context_results: dict[str, Any] = {}
    pick_rng = np.random.default_rng(int(production["library_pick_seed"]))
    all_pass = True

    for context_index, context in enumerate(contexts):
        item = inputs["released_controls"][context]
        b_ctrl, control_libraries = control_profile_and_libraries(
            ROOT / item["path"],
            genes,
            ts=ts,
        )
        baseline_profiles = np.stack(
            [
                profile_from_effect(
                    row,
                    b_ctrl,
                    ts=ts,
                    min_cpm=min_cpm,
                    effect_scale=effect_scale,
                )
                for row in effects
            ]
        )
        for target_index in range(len(panel)):
            picks = pick_rng.choice(
                len(control_libraries),
                cells,
                replace=len(control_libraries) < cells,
            )
            actual_libraries[context_index, target_index] = int(
                min(
                    np.median(control_libraries[picks]),
                    int(production["max_counts_per_cell"]),
                )
            )
            mutable_counts[context_index, target_index] = int(
                np.sum(
                    mutable_panel_mask(
                        baseline_profiles[target_index],
                        b_ctrl,
                        residual[target_index],
                        panel_columns,
                        target_index,
                        min_cpm=min_cpm,
                    )
                )
            )

        context_doses: dict[str, Any] = {}
        for dose_index, eta in enumerate(doses):
            scheme_results: dict[str, Any] = {}
            preclip_counts: list[float] = []
            preclip_rows = 0
            continuous_signs: list[float] = []
            for scheme_index, scheme_name in enumerate(scheme_names):
                rows_pass: list[bool] = []
                signs: list[float] = []
                changed_coordinates: list[float] = []
                l1_changes: list[float] = []
                mutable_budgets: list[float] = []
                fallbacks = 0
                baseline_rows: list[np.ndarray] = []
                candidate_rows: list[np.ndarray] = []
                for target_index in range(len(panel)):
                    library = (
                        fixed_libraries[scheme_index]
                        if scheme_index < len(fixed_libraries)
                        else int(actual_libraries[context_index, target_index])
                    )
                    baseline_totals = np.floor(
                        cells * library * baseline_profiles[target_index]
                    ).astype(np.int64)
                    candidate_totals, audit = integer_budget_transport(
                        baseline_totals,
                        baseline_profiles[target_index],
                        b_ctrl,
                        residual[target_index],
                        panel_columns,
                        panel_mask,
                        target_index,
                        eta,
                        ts=ts,
                        min_cpm=min_cpm,
                    )
                    rows_pass.append(bool(audit["all_invariants_pass"]))
                    signs.append(float(audit["integer_transformed_sign_retention"]))
                    changed_coordinates.append(float(audit["changed_mutable_coordinates"]))
                    l1_changes.append(float(audit["mutable_l1_count_change"]))
                    mutable_budgets.append(float(audit["mutable_integer_budget"]))
                    fallbacks += int(audit["fallback_used"])
                    baseline_rows.append(baseline_totals)
                    candidate_rows.append(candidate_totals)
                    sign_retention[
                        context_index,
                        dose_index,
                        scheme_index,
                        target_index,
                    ] = audit["integer_transformed_sign_retention"]
                    invariant_pass[
                        context_index,
                        dose_index,
                        scheme_index,
                        target_index,
                    ] = audit["all_invariants_pass"]
                    if scheme_index == 0:
                        preclip_negative[
                            context_index,
                            dose_index,
                            target_index,
                        ] = audit["preclip_negative_coordinates"]
                        preclip_counts.append(
                            float(audit["preclip_negative_coordinates"])
                        )
                        preclip_rows += int(audit["preclip_negative_row"])
                        continuous_signs.append(
                            float(audit["continuous_source_sign_retention"])
                        )
                scheme_pass = bool(all(rows_pass))
                all_pass &= scheme_pass
                scheme_results[scheme_name] = {
                    "all_rows_pass": scheme_pass,
                    "failed_rows": int(len(rows_pass) - sum(rows_pass)),
                    "integer_transformed_sign_retention": summarize(signs),
                    "changed_mutable_coordinates": summarize(changed_coordinates),
                    "mutable_l1_count_change": summarize(l1_changes),
                    "mutable_integer_budget": summarize(mutable_budgets),
                    "degenerate_fallback_rows": fallbacks,
                    "baseline_integer_totals_sha256": sha256_int_rows(baseline_rows),
                    "candidate_integer_totals_sha256": sha256_int_rows(candidate_rows),
                }
            context_doses[f"{eta:.2f}"] = {
                "diagnostic_only": True,
                "preclip_negative_coordinates": {
                    "total": int(sum(preclip_counts)),
                    "rows": preclip_rows,
                    "per_row": summarize(preclip_counts),
                },
                "continuous_source_sign_retention": summarize(continuous_signs),
                "library_schemes": scheme_results,
            }

        context_results[context] = {
            "control_cells": int(control_libraries.size),
            "control_library_size": summarize(control_libraries.tolist()),
            "replayed_actual_production_L": {
                **summarize(actual_libraries[context_index].tolist()),
                "unique_values": int(np.unique(actual_libraries[context_index]).size),
            },
            "mutable_panel_coordinates": summarize(
                mutable_counts[context_index].tolist()
            ),
            "doses": context_doses,
        }

    result = {
        "schema_version": 1,
        "experiment": "128-panel-integer-budget-transport",
        "status": "mechanical_audit_complete",
        "source_units": "delta in log1p(50000 * relative_abundance)",
        "doses_are_diagnostics_only": True,
        "provenance": provenance,
        "axes": {
            "contexts": list(contexts),
            "doses": list(doses),
            "library_schemes": list(scheme_names),
            "targets": len(panel),
            "genes": len(genes),
        },
        "contexts": context_results,
        "decision": {
            "all_integer_and_scorer_space_invariants_pass": bool(all_pass),
            "production_promotion_authorized": False,
            "biological_confirmation_required": True,
        },
        "forbidden_actions_respected": {
            "held_outcomes_accessed": False,
            "challenge_perturbation_truth_accessed": False,
            "cell_eval2_run": False,
            "count_h5ad_or_vcc_built": False,
            "submission_sent": False,
            "commit": False,
            "push": False,
        },
    }
    require(all_pass, "integer-budget or scorer-space invariant failed")

    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    stage = output.parent / f".{output.name}.tmp-{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    files = {
        "actual-library-medians.npy": actual_libraries,
        "integer-sign-retention.npy": sign_retention,
        "preclip-negative-counts.npy": preclip_negative,
        "invariant-pass.npy": invariant_pass,
        "mutable-coordinate-counts.npy": mutable_counts,
    }
    for name, value in files.items():
        with (stage / name).open("xb") as handle:
            np.save(handle, np.ascontiguousarray(value), allow_pickle=False)
    result_path = stage / "results.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    manifest_files: dict[str, Any] = {}
    for path in sorted(stage.iterdir()):
        manifest_files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    manifest = {
        "schema_version": 1,
        "experiment": "128-panel-integer-budget-transport",
        "files": manifest_files,
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "implementation_freeze_sha256": sha256_file(FREEZE_PATH),
    }
    (stage / "artifact-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    os.replace(stage, output)
    for path in output.iterdir():
        path.chmod(0o444)
    output.chmod(0o555)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "work/integration-audit-v1",
    )
    args = parser.parse_args()
    require(args.audit, "run requires --audit")
    result = run(args.output.resolve())
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
