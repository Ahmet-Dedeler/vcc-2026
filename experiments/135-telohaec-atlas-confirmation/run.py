#!/usr/bin/env python3
"""One-shot TeloHAEC confirmation of the frozen Experiment 122 panel atlas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "preregistered_protocol.json"
FREEZE = HERE / "IMPLEMENTATION_FREEZE.json"
OUTPUT = HERE / "work/telohaec-confirmation-v1"
H5AD = ROOT / "data/external/igvf-telohaec/IGVFFI8645YMGX.h5ad"
GUIDE_QC = ROOT / "data/external/igvf-telohaec/guide_qc.csv"
GENES = ROOT / "data/challenge/gene_names.csv"
PANEL = ROOT / "data/challenge/pert_counts.csv"
EXP128 = ROOT / "experiments/128-panel-integer-budget-transport/run.py"
EXP128_FREEZE = (
    ROOT / "experiments/128-panel-integer-budget-transport/IMPLEMENTATION_FREEZE.json"
)
EXP122_ROOT = (
    ROOT / "experiments/122-panel-response-atlas/work/source-only-phase1-v1"
)
RESIDUAL = EXP122_ROOT / "selected-panel-residual.npy"
CONFIDENCE = EXP122_ROOT / "selected-panel-confidence.npy"
RESIDUAL_BANK = EXP122_ROOT / "all-candidate-panel-residuals.npy"
CONFIDENCE_BANK = EXP122_ROOT / "all-candidate-panel-confidence.npy"
SUPPORT_BANK = EXP122_ROOT / "all-candidate-panel-support.npy"
CANDIDATE_NAMES = EXP122_ROOT / "candidate-names.txt"
CANDIDATE_AXIS = EXP122_ROOT / "panel-targets.txt"
CANDIDATE_RESULTS = EXP122_ROOT / "results.json"
CANDIDATE_MANIFEST = EXP122_ROOT / "artifact-manifest.json"
CANDIDATE_RECEIPT = EXP122_ROOT / "result-receipt.json"
INCUMBENT_EFFECTS = (
    ROOT
    / "experiments/103-residual-only-production/work/staged-a15-linux/a15-effects.npy"
)
INCUMBENT_METADATA = INCUMBENT_EFFECTS.with_suffix(".npy.json")
TS = 50_000.0
CELLS = 400
MIN_CPM = 10.0
MAX_LIBRARY = 1_000_000


class ExperimentError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ExperimentError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(np.asarray(contiguous.shape, dtype="<i8").tobytes())
    digest.update(contiguous.tobytes())
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(row[column].strip().upper() for row in reader)
    require(len(values) == len(set(values)), f"{path}: duplicate labels")
    return values


def decoded(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def summarize(values: np.ndarray | list[float] | list[int]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    require(array.size > 0 and np.all(np.isfinite(array)), "invalid summary input")
    return {
        "minimum": float(np.min(array)),
        "median": float(np.median(array)),
        "maximum": float(np.max(array)),
    }


def load_exp128() -> Any:
    spec = importlib.util.spec_from_file_location("exp128_transport", EXP128)
    require(
        spec is not None and spec.loader is not None, "cannot import Experiment 128"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_protocol_constants(protocol: dict[str, Any]) -> None:
    axes = protocol["axes"]
    transport = protocol["transport"]
    require(int(axes["expected_genes"]) == 18_533, "protocol gene count drift")
    require(
        int(axes["expected_panel_targets"]) == 300,
        "protocol panel count drift",
    )
    require(float(transport["source_ts"]) == TS, "protocol source TS drift")
    require(int(transport["cells"]) == CELLS, "protocol cell count drift")
    require(
        float(transport["minimum_control_cpm"]) == MIN_CPM,
        "protocol CPM gate drift",
    )
    require(
        int(transport["maximum_library"]) == MAX_LIBRARY,
        "protocol library cap drift",
    )
    require(float(protocol["candidate"]["dose"]) == 1.0, "candidate dose drift")


def verify_freeze(*, hash_large_source: bool) -> dict[str, Any]:
    require(FREEZE.is_file(), "implementation freeze is absent")
    freeze = read_json(FREEZE)
    require(
        freeze.get("status")
        == "frozen_before_experiment_135_telohaec_panel_response_aggregation_or_analysis",
        "wrong freeze status",
    )
    for relative, record in freeze["files"].items():
        path = ROOT / relative
        require(
            path.is_file() and not path.is_symlink(), f"missing frozen file: {relative}"
        )
        require(path.stat().st_size == record["size_bytes"], f"size drift: {relative}")
        require(sha256_file(path) == record["sha256"], f"hash drift: {relative}")
    protocol = read_json(PROTOCOL)
    verify_protocol_constants(protocol)
    require(
        H5AD.stat().st_size == protocol["source"]["expected_matrix_size_bytes"],
        "source size drift",
    )
    if hash_large_source:
        require(
            sha256_file(H5AD) == protocol["source"]["expected_matrix_sha256"],
            "source SHA-256 drift",
        )
    require(
        sha256_file(EXP128) == protocol["transport"]["implementation_sha256"],
        "Experiment 128 implementation drift",
    )
    require(
        sha256_file(EXP128_FREEZE) == protocol["transport"]["freeze_sha256"],
        "Experiment 128 freeze drift",
    )
    return protocol


def verify_candidate_provenance(protocol: dict[str, Any]) -> dict[str, bool]:
    candidate = protocol["candidate"]
    expected = candidate["artifact_sha256"]
    paths = {
        "selected-panel-residual.npy": RESIDUAL,
        "selected-panel-confidence.npy": CONFIDENCE,
        "all-candidate-panel-residuals.npy": RESIDUAL_BANK,
        "all-candidate-panel-confidence.npy": CONFIDENCE_BANK,
        "all-candidate-panel-support.npy": SUPPORT_BANK,
        "candidate-names.txt": CANDIDATE_NAMES,
        "panel-targets.txt": CANDIDATE_AXIS,
        "results.json": CANDIDATE_RESULTS,
    }
    require(
        sha256_file(CANDIDATE_MANIFEST) == candidate["artifact_manifest_sha256"],
        "Experiment 122 artifact manifest hash drift",
    )
    require(
        sha256_file(CANDIDATE_RECEIPT) == candidate["result_receipt_sha256"],
        "Experiment 122 result receipt hash drift",
    )
    manifest = read_json(CANDIDATE_MANIFEST)
    receipt = read_json(CANDIDATE_RECEIPT)
    results = read_json(CANDIDATE_RESULTS)
    require(
        manifest.get("selected_candidate") == candidate["name"],
        "Experiment 122 manifest selected candidate drift",
    )
    require(
        receipt.get("selected_candidate") == candidate["name"],
        "Experiment 122 receipt selected candidate drift",
    )
    require(
        results.get("selected_candidate") == candidate["name"],
        "Experiment 122 results selected candidate drift",
    )
    require(
        receipt.get("artifact_manifest_sha256")
        == candidate["artifact_manifest_sha256"],
        "Experiment 122 receipt does not bind its manifest",
    )
    for name, path in paths.items():
        require(sha256_file(path) == expected[name], f"candidate hash drift: {name}")
        require(
            manifest["files"][name]["sha256"] == expected[name],
            f"Experiment 122 manifest does not bind {name}",
        )
    return {
        "manifest_hash_matches": True,
        "receipt_hash_matches": True,
        "receipt_binds_manifest": True,
        "manifest_binds_all_candidate_artifacts": True,
        "selected_candidate_matches_across_manifest_receipt_and_results": True,
    }


def load_guide_contract(
    panel: tuple[str, ...], protocol: dict[str, Any]
) -> dict[str, Any]:
    panel_set = set(panel)
    controls: list[str] = []
    target_spacers: dict[str, list[str]] = {target: [] for target in panel}
    expected_cells: dict[str, int] = {}
    with GUIDE_QC.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "spacer",
            "clean_target",
            "eligible",
            "eligibility_note",
            "observed_cells",
        }
        require(
            reader.fieldnames is not None and required <= set(reader.fieldnames),
            "guide_qc schema drift",
        )
        for row in reader:
            if row["eligible"] != "1":
                continue
            spacer = row["spacer"].strip().upper()
            target = row["clean_target"].strip().upper()
            require(spacer not in expected_cells, f"duplicate eligible spacer: {spacer}")
            expected_cells[spacer] = int(row["observed_cells"])
            if row["eligibility_note"] == "eligible_control":
                controls.append(spacer)
            elif target in panel_set:
                target_spacers[target].append(spacer)
    retained = {target: values for target, values in target_spacers.items() if values}
    minimum_cells_per_guide = int(protocol["truth"]["minimum_cells_per_guide"])
    minimum_cells_per_target = int(
        protocol["truth"]["minimum_eligible_cells_per_target"]
    )
    minimum_guides_per_target = int(
        protocol["truth"]["minimum_eligible_guides_per_target"]
    )
    require(
        len(controls) == protocol["source"]["expected_controls"],
        f"expected {protocol['source']['expected_controls']} controls, found {len(controls)}",
    )
    require(
        len(retained) == protocol["source"]["expected_targets"],
        f"expected {protocol['source']['expected_targets']} targets, found {len(retained)}",
    )
    require(
        all(len(values) >= minimum_guides_per_target for values in retained.values()),
        "target has too few eligible guides",
    )
    require(
        all(expected_cells[spacer] >= minimum_cells_per_guide for spacer in controls),
        "control guide falls below the frozen cell floor",
    )
    for target, spacers in retained.items():
        require(
            all(
                expected_cells[spacer] >= minimum_cells_per_guide
                for spacer in spacers
            ),
            f"{target}: guide falls below the frozen cell floor",
        )
        require(
            sum(expected_cells[spacer] for spacer in spacers)
            >= minimum_cells_per_target,
            f"{target}: target falls below the frozen cell floor",
        )
    all_spacers = controls + [
        spacer for target in panel for spacer in retained.get(target, ())
    ]
    require(len(all_spacers) == len(set(all_spacers)), "eligible spacer overlap")
    return {
        "controls": tuple(controls),
        "targets": retained,
        "expected_cells": expected_cells,
    }


def feature_maps(
    raw_names: list[str], panel: tuple[str, ...]
) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    unique: list[str] = []
    lookup: dict[str, int] = {}
    raw_to_unique = np.empty(len(raw_names), dtype=np.int32)
    panel_lookup = {gene: index for index, gene in enumerate(panel)}
    raw_to_panel = np.full(len(raw_names), -1, dtype=np.int16)
    for index, raw in enumerate(raw_names):
        require(":" in raw, f"unexpected feature label: {raw}")
        symbol, ensembl = raw.split(":", 1)
        symbol = symbol.strip().upper()
        require(
            symbol and ensembl.startswith("ENSG"), f"unexpected feature label: {raw}"
        )
        if symbol not in lookup:
            lookup[symbol] = len(unique)
            unique.append(symbol)
        raw_to_unique[index] = lookup[symbol]
        if symbol in panel_lookup:
            raw_to_panel[index] = panel_lookup[symbol]
    require(set(panel) <= set(unique), "TeloHAEC lacks a panel response symbol")
    return tuple(unique), raw_to_unique, raw_to_panel


def barcode_spacer(value: object) -> str:
    parts = decoded(value).split(":", 2)
    require(len(parts) == 3, f"unexpected barcode label: {decoded(value)!r}")
    return parts[1].upper()


def map_cells(
    barcodes: h5py.Dataset, spacer_to_row: dict[str, int], chunk: int = 4096
) -> np.ndarray:
    output = np.full(barcodes.shape[0], -1, dtype=np.int32)
    for start in range(0, barcodes.shape[0], chunk):
        stop = min(start + chunk, barcodes.shape[0])
        output[start:stop] = [
            spacer_to_row.get(barcode_spacer(value), -1)
            for value in barcodes[start:stop]
        ]
    return output


def scan_controls(
    handle: h5py.File,
    control_spacers: tuple[str, ...],
    expected_cells: dict[str, int],
    raw_to_unique: np.ndarray,
    unique_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    spacer_to_row = {spacer: index for index, spacer in enumerate(control_spacers)}
    cell_rows = map_cells(handle["obs/barcodes"], spacer_to_row)
    sums = np.zeros((len(control_spacers), unique_count), dtype=np.float64)
    totals = np.zeros(len(control_spacers), dtype=np.float64)
    x = handle["X"]
    indptr = np.asarray(x["indptr"], dtype=np.int64)
    control_cells = np.flatnonzero(cell_rows >= 0)
    observed_counts = np.bincount(
        cell_rows[control_cells], minlength=len(control_spacers)
    )
    expected_counts = np.asarray(
        [expected_cells[spacer] for spacer in control_spacers], dtype=np.int64
    )
    require(
        np.array_equal(observed_counts, expected_counts),
        "control guide cell counts drifted",
    )
    control_cell_libraries = np.empty(len(control_cells), dtype=np.float64)
    for output_row, cell in enumerate(control_cells):
        lo, hi = int(indptr[cell]), int(indptr[cell + 1])
        values = np.asarray(x["data"][lo:hi], dtype=np.float64)
        columns = np.asarray(x["indices"][lo:hi], dtype=np.int64)
        require(
            np.all(np.isfinite(values))
            and np.all(values >= 0)
            and np.all(values == np.floor(values)),
            "invalid control count value",
        )
        guide = int(cell_rows[cell])
        np.add.at(sums[guide], raw_to_unique[columns], values)
        library = float(values.sum(dtype=np.float64))
        totals[guide] += library
        control_cell_libraries[output_row] = library
    require(np.all(totals > 0), "empty eligible control guide")
    require(np.all(control_cell_libraries > 0), "empty eligible control cell")
    return sums, totals, control_cell_libraries


def scan_targets(
    handle: h5py.File,
    target_order: tuple[str, ...],
    target_spacers: dict[str, list[str]],
    expected_cells: dict[str, int],
    raw_to_panel: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...], dict[str, list[int]]]:
    spacers: list[str] = []
    target_rows: dict[str, list[int]] = {}
    for target in target_order:
        rows: list[int] = []
        for spacer in target_spacers[target]:
            rows.append(len(spacers))
            spacers.append(spacer)
        target_rows[target] = rows
    spacer_to_row = {spacer: index for index, spacer in enumerate(spacers)}
    cell_rows = map_cells(handle["obs/barcodes"], spacer_to_row)
    panel_sums = np.zeros((len(spacers), 300), dtype=np.float64)
    totals = np.zeros(len(spacers), dtype=np.float64)
    x = handle["X"]
    indptr = np.asarray(x["indptr"], dtype=np.int64)
    target_cells = np.flatnonzero(cell_rows >= 0)
    observed_counts = np.bincount(cell_rows[target_cells], minlength=len(spacers))
    expected_counts = np.asarray(
        [expected_cells[spacer] for spacer in spacers], dtype=np.int64
    )
    require(
        np.array_equal(observed_counts, expected_counts),
        "target guide cell counts drifted",
    )
    for cell in target_cells:
        lo, hi = int(indptr[cell]), int(indptr[cell + 1])
        values = np.asarray(x["data"][lo:hi], dtype=np.float64)
        columns = np.asarray(x["indices"][lo:hi], dtype=np.int64)
        require(
            np.all(np.isfinite(values))
            and np.all(values >= 0)
            and np.all(values == np.floor(values)),
            "invalid targeting count value",
        )
        panel_columns = raw_to_panel[columns]
        keep = panel_columns >= 0
        guide = int(cell_rows[cell])
        if np.any(keep):
            np.add.at(panel_sums[guide], panel_columns[keep], values[keep])
        totals[guide] += float(values.sum(dtype=np.float64))
    require(np.all(totals > 0), "empty eligible targeting guide")
    return panel_sums, totals, tuple(spacers), target_rows


def response_bins(control_panel_abundance: np.ndarray) -> np.ndarray:
    require(control_panel_abundance.shape == (300,), "response abundance shape drift")
    order = np.lexsort((np.arange(300), control_panel_abundance))
    bins = np.empty(300, dtype=np.int16)
    bins[order] = np.arange(300, dtype=np.int16) // 30
    require(set(map(int, bins)) == set(range(10)), "response bins changed")
    require(
        all(np.count_nonzero(bins == index) == 30 for index in range(10)),
        "response bins became unbalanced",
    )
    return bins


def freeze_nulls(
    target_panel_indices: np.ndarray,
    bins: np.ndarray,
    draws: int,
    *,
    target_seed: int,
    response_seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Freeze structurally exchangeable target and response null mappings."""
    target_panel_indices = np.asarray(target_panel_indices, dtype=np.int16)
    target_rng = np.random.default_rng(target_seed)
    response_rng = np.random.default_rng(response_seed)
    target_base = np.arange(target_panel_indices.size, dtype=np.int16)
    target_rows = np.stack(
        [target_rng.permutation(target_base) for _ in range(draws)]
    )
    target_responses = np.tile(np.arange(300, dtype=np.int16), (draws, 1))
    for draw, row_permutation in enumerate(target_rows):
        target_responses[draw, target_panel_indices] = target_panel_indices[
            row_permutation
        ]

    response = np.empty(
        (draws, target_panel_indices.size, 300), dtype=np.int16
    )
    identity = np.arange(300, dtype=np.int16)
    for draw in range(draws):
        for row, target_index in enumerate(target_panel_indices):
            permutation = identity.copy()
            for bin_id in range(10):
                members = np.flatnonzero(bins == bin_id).astype(np.int16)
                movable = members[members != target_index]
                permutation[movable] = response_rng.permutation(movable)
            require(permutation[target_index] == target_index, "response null moved self")
            response[draw, row] = permutation
    return target_rows, target_responses, response


def load_candidate(
    panel: tuple[str, ...], protocol: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, bool]]:
    provenance = verify_candidate_provenance(protocol)
    names = tuple(CANDIDATE_NAMES.read_text(encoding="utf-8").splitlines())
    candidate_index = int(protocol["candidate"]["candidate_index"])
    require(
        names[candidate_index] == protocol["candidate"]["name"],
        "candidate index drift",
    )
    require(
        tuple(
            value.strip().upper()
            for value in CANDIDATE_AXIS.read_text(encoding="utf-8").splitlines()
        )
        == panel,
        "candidate panel axis drift",
    )
    residual = np.load(RESIDUAL, allow_pickle=False).astype(np.float64)
    confidence = np.load(CONFIDENCE, allow_pickle=False).astype(np.float64)
    residual_bank = np.load(RESIDUAL_BANK, allow_pickle=False)
    confidence_bank = np.load(CONFIDENCE_BANK, allow_pickle=False)
    support_bank = np.load(SUPPORT_BANK, allow_pickle=False)
    support = support_bank[candidate_index] > 0
    require(
        residual.shape == confidence.shape == support.shape == (300, 300),
        "candidate shape drift",
    )
    require(
        residual_bank.shape == confidence_bank.shape == support_bank.shape,
        "candidate bank shape drift",
    )
    require(
        np.array_equal(residual, residual_bank[candidate_index]),
        "selected residual differs from candidate-bank row",
    )
    require(
        np.array_equal(confidence, confidence_bank[candidate_index]),
        "selected confidence differs from candidate-bank row",
    )
    require(np.all(np.isfinite(residual)), "nonfinite residual")
    require(
        np.all(np.isfinite(confidence)) and np.all(confidence >= 0),
        "invalid confidence",
    )
    require(
        np.all(np.diag(residual) == 0) and np.all(np.diag(confidence) == 0),
        "candidate self coordinate drift",
    )
    return residual, confidence, support, provenance


def load_incumbent_effects(
    canonical_genes: tuple[str, ...],
    unique_genes: tuple[str, ...],
    panel: tuple[str, ...],
    protocol: dict[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    incumbent = protocol["incumbent"]
    require(
        sha256_file(INCUMBENT_EFFECTS) == incumbent["effect_sha256"],
        "staged-a15 effect hash drift",
    )
    require(
        sha256_file(INCUMBENT_METADATA) == incumbent["metadata_sha256"],
        "staged-a15 metadata hash drift",
    )
    metadata = read_json(INCUMBENT_METADATA)
    require(
        metadata["file_sha256"] == incumbent["effect_sha256"],
        "staged-a15 metadata does not bind the effect array",
    )
    require(
        metadata["shape"] == [len(panel), len(canonical_genes)],
        "staged-a15 metadata shape drift",
    )
    effects = np.load(INCUMBENT_EFFECTS, allow_pickle=False).astype(np.float64)
    require(
        effects.shape == (len(panel), len(canonical_genes)),
        "staged-a15 effect shape drift",
    )
    require(np.all(np.isfinite(effects)), "staged-a15 contains nonfinite effects")
    canonical_index = {gene: index for index, gene in enumerate(canonical_genes)}
    mapped = np.zeros((len(panel), len(unique_genes)), dtype=np.float64)
    covered = np.zeros(len(unique_genes), dtype=bool)
    for unique_index, gene in enumerate(unique_genes):
        source_index = canonical_index.get(gene)
        if source_index is not None:
            mapped[:, unique_index] = effects[:, source_index]
            covered[unique_index] = True
    return mapped, {
        "canonical_genes": len(canonical_genes),
        "telohaec_unique_genes": len(unique_genes),
        "mapped_unique_genes": int(covered.sum()),
        "unmapped_unique_genes_zero_filled": int((~covered).sum()),
        "all_panel_response_genes_mapped": bool(set(panel) <= set(unique_genes)),
    }


def replay_libraries(
    control_cell_libraries: np.ndarray,
    target_count: int,
    *,
    cells: int,
    seed: int,
    maximum: int,
) -> np.ndarray:
    libraries = np.asarray(control_cell_libraries, dtype=np.float64)
    require(libraries.size >= cells, "not enough control cells for no-replacement replay")
    require(np.all(libraries > 0), "invalid control-cell library")
    rng = np.random.default_rng(seed)
    output = np.empty(target_count, dtype=np.int64)
    for target_index in range(target_count):
        picks = rng.choice(libraries.size, cells, replace=False)
        output[target_index] = int(
            min(math.floor(float(np.median(libraries[picks]))), maximum)
        )
    require(np.all(output > 0), "replayed library is empty")
    return output


def build_incumbent_rows(
    exp128: Any,
    mapped_effects: np.ndarray,
    target_panel_indices: np.ndarray,
    b_ctrl: np.ndarray,
    panel_columns: np.ndarray,
    replayed_libraries: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    profiles = np.empty(
        (target_panel_indices.size, mapped_effects.shape[1]), dtype=np.float64
    )
    totals = np.empty_like(profiles, dtype=np.int64)
    predictions = np.empty((target_panel_indices.size, 300), dtype=np.float64)
    control_panel = b_ctrl[panel_columns]
    for row, target_index in enumerate(target_panel_indices):
        profile = exp128.profile_from_effect(
            mapped_effects[target_index],
            b_ctrl,
            ts=TS,
            min_cpm=MIN_CPM,
            effect_scale=1.0,
        )
        integer_totals = np.floor(
            CELLS * int(replayed_libraries[target_index]) * profile
        ).astype(np.int64)
        total = int(integer_totals.sum(dtype=np.int64))
        require(total > 0, "incumbent integerization produced an empty row")
        profiles[row] = profile
        totals[row] = integer_totals
        predictions[row] = (
            np.log1p(TS * integer_totals[panel_columns].astype(np.float64) / total)
            - control_panel
        )
    return profiles, totals, predictions


def transport_prediction(
    exp128: Any,
    incumbent_profile: np.ndarray,
    incumbent_totals: np.ndarray,
    b_ctrl: np.ndarray,
    panel_columns: np.ndarray,
    panel_mask: np.ndarray,
    residual: np.ndarray,
    target_index: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    candidate_totals, audit = exp128.integer_budget_transport(
        incumbent_totals,
        incumbent_profile,
        b_ctrl,
        residual,
        panel_columns,
        panel_mask,
        target_index,
        1.0,
        ts=TS,
        min_cpm=MIN_CPM,
    )
    total = int(candidate_totals.sum(dtype=np.int64))
    require(total > 0, "candidate integer transport produced an empty row")
    prediction = (
        np.log1p(TS * candidate_totals[panel_columns].astype(np.float64) / total)
        - b_ctrl[panel_columns]
    )
    audit["sign_retention"] = audit["integer_transformed_sign_retention"]
    return prediction, audit


def target_metrics(
    prediction: np.ndarray, truth: np.ndarray, target_index: int
) -> tuple[float, float]:
    valid = np.isfinite(truth)
    valid[target_index] = False
    indices = np.flatnonzero(valid)
    require(indices.size == 299, "truth response count changed")
    denominator = float(np.abs(truth[indices]).sum())
    require(denominator > 0, "zero truth magnitude")
    credit = float(np.sum(np.sign(truth[indices]) * prediction[indices]) / denominator)
    top_n = int(math.ceil(0.20 * indices.size))
    chosen = indices[
        np.lexsort((indices, -np.abs(prediction[indices])))[:top_n]
    ]
    weights = np.abs(truth[chosen])
    require(float(weights.sum()) > 0, "zero top-20 truth magnitude")
    correct = (
        (np.sign(prediction[chosen]) != 0)
        & (np.sign(truth[chosen]) != 0)
        & (np.sign(prediction[chosen]) == np.sign(truth[chosen]))
    )
    precision = float(np.sum(weights * correct) / np.sum(weights))
    return credit, precision


def baseline_metrics(
    target_panel_indices: np.ndarray,
    truth: np.ndarray,
    incumbent_predictions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    credits = np.empty(target_panel_indices.size, dtype=np.float64)
    precisions = np.empty(target_panel_indices.size, dtype=np.float64)
    for row, target_index in enumerate(target_panel_indices):
        credits[row], precisions[row] = target_metrics(
            incumbent_predictions[row], truth[row], int(target_index)
        )
    return credits, precisions


def evaluate_rows(
    exp128: Any,
    target_order: tuple[str, ...],
    panel: tuple[str, ...],
    target_panel_indices: np.ndarray,
    truth: np.ndarray,
    residual: np.ndarray,
    confidence: np.ndarray,
    incumbent_profiles: np.ndarray,
    incumbent_totals: np.ndarray,
    b_ctrl: np.ndarray,
    panel_columns: np.ndarray,
    panel_mask: np.ndarray,
    *,
    target_mapping: np.ndarray | None = None,
    response_mapping: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    credits = np.empty(len(target_order), dtype=np.float64)
    precisions = np.empty(len(target_order), dtype=np.float64)
    audits: list[dict[str, Any]] = []
    for row, target in enumerate(target_order):
        source_row = row if target_mapping is None else int(target_mapping[row])
        candidate_target = int(target_panel_indices[source_row])
        r = residual[candidate_target].copy()
        c = confidence[candidate_target].copy()
        if response_mapping is not None:
            mapping = (
                response_mapping
                if response_mapping.ndim == 1
                else response_mapping[row]
            )
            r = r[mapping]
            c = c[mapping]
        destination = int(target_panel_indices[row])
        require(r[destination] == 0.0, "candidate/null moved the structural self-zero")
        require(c[destination] == 0.0, "confidence/null moved the structural self-zero")
        prediction, audit = transport_prediction(
            exp128,
            incumbent_profiles[row],
            incumbent_totals[row],
            b_ctrl,
            panel_columns,
            panel_mask,
            r,
            destination,
        )
        credits[row], precisions[row] = target_metrics(
            prediction, truth[row], destination
        )
        audit["target"] = target
        audit["candidate_source_target"] = panel[candidate_target]
        audits.append(audit)
    return credits, precisions, audits


def guide_truth_and_diagnostics(
    target_order: tuple[str, ...],
    target_panel_indices: np.ndarray,
    target_rows: dict[str, list[int]],
    guide_panel_profiles: np.ndarray,
    control_panel: np.ndarray,
    expected_cells: dict[str, int],
    spacers: tuple[str, ...],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    truth = np.empty((len(target_order), 300), dtype=np.float64)
    guide_counts = np.empty(len(target_order), dtype=np.int16)
    target_cell_counts = np.empty(len(target_order), dtype=np.int32)
    guide_sign_agreement = np.empty(len(target_order), dtype=np.float64)
    held_out_sign_min = np.empty(len(target_order), dtype=np.float64)
    held_out_sign_median = np.empty(len(target_order), dtype=np.float64)
    held_out_nmae_max = np.empty(len(target_order), dtype=np.float64)
    for row, target in enumerate(target_order):
        indices = np.asarray(target_rows[target], dtype=np.int64)
        guide_effects = guide_panel_profiles[indices] - control_panel[None, :]
        target_truth = np.median(guide_effects, axis=0)
        truth[row] = target_truth
        guide_counts[row] = len(indices)
        target_cell_counts[row] = sum(
            expected_cells[spacers[index]] for index in indices
        )
        valid = np.ones(300, dtype=bool)
        valid[int(target_panel_indices[row])] = False
        held_out_sign: list[float] = []
        held_out_nmae: list[float] = []
        for held in range(len(indices)):
            keep = np.arange(len(indices)) != held
            remaining_median = np.median(guide_effects[keep], axis=0)
            informative = valid & (np.sign(remaining_median) != 0)
            require(
                np.any(informative),
                f"{target}: no informative leave-one-guide coordinates",
            )
            held_out_sign.append(
                float(
                    np.mean(
                        np.sign(guide_effects[held, informative])
                        == np.sign(remaining_median[informative])
                    )
                )
            )
            denominator = float(np.abs(remaining_median[valid]).sum())
            require(
                denominator > 0,
                f"{target}: zero leave-one-guide reference magnitude",
            )
            held_out_nmae.append(
                float(
                    np.abs(
                        guide_effects[held, valid] - remaining_median[valid]
                    ).sum()
                    / denominator
                )
            )
        guide_sign_agreement[row] = float(np.mean(held_out_sign))
        held_out_sign_min[row] = min(held_out_sign)
        held_out_sign_median[row] = float(np.median(held_out_sign))
        held_out_nmae_max[row] = max(held_out_nmae)
    return truth, {
        "guide_counts": guide_counts,
        "target_cell_counts": target_cell_counts,
        "guide_sign_agreement": guide_sign_agreement,
        "held_out_guide_sign_agreement_min": held_out_sign_min,
        "held_out_guide_sign_agreement_median": held_out_sign_median,
        "held_out_guide_normalized_mae_max": held_out_nmae_max,
    }


def add_one_p(null: np.ndarray, observed: float) -> float:
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def metric_null_summary(
    credit_null: np.ndarray,
    precision_null: np.ndarray,
    observed_credit_delta: float,
    observed_precision_delta: float,
) -> dict[str, float]:
    return {
        "signed_credit_delta_p95": float(
            np.quantile(credit_null, 0.95, method="linear")
        ),
        "precision_delta_p95": float(
            np.quantile(precision_null, 0.95, method="linear")
        ),
        "signed_credit_delta_p": add_one_p(
            credit_null, observed_credit_delta
        ),
        "precision_delta_p": add_one_p(
            precision_null, observed_precision_delta
        ),
    }


def run_evaluation() -> dict[str, Any]:
    protocol = verify_freeze(hash_large_source=True)
    panel = read_axis(PANEL, "target_gene")
    canonical_genes = read_axis(GENES, "gene_name")
    require(len(panel) == 300 and len(canonical_genes) == 18533, "canonical axis drift")
    guide_contract = load_guide_contract(panel, protocol)
    target_order = tuple(
        target for target in panel if target in guide_contract["targets"]
    )
    require(
        len(target_order) == protocol["source"]["expected_targets"],
        "target order drift",
    )
    panel_index = {target: index for index, target in enumerate(panel)}
    target_panel_indices = np.asarray(
        [panel_index[target] for target in target_order], dtype=np.int16
    )
    residual, confidence, support, provenance = load_candidate(panel, protocol)
    exp128 = load_exp128()

    staging = OUTPUT.with_name(f".{OUTPUT.name}.staging")
    require(
        not OUTPUT.exists() and not staging.exists(), "refusing to overwrite output"
    )
    staging.mkdir(parents=True)
    try:
        with h5py.File(H5AD, "r") as handle:
            raw_names = [decoded(value) for value in handle["var/_index"][:]]
            unique_genes, raw_to_unique, raw_to_panel = feature_maps(raw_names, panel)
            control_sums, control_totals, control_cell_libraries = scan_controls(
                handle,
                guide_contract["controls"],
                guide_contract["expected_cells"],
                raw_to_unique,
                len(unique_genes),
            )
            pooled_control = control_sums.sum(axis=0)
            require(
                float(pooled_control.sum()) == float(control_totals.sum()),
                "pooled control total drift",
            )
            b_ctrl = np.log1p(TS * pooled_control / pooled_control.sum())
            control_abundance = np.expm1(b_ctrl) / TS
            control_abundance /= control_abundance.sum()
            unique_index = {gene: index for index, gene in enumerate(unique_genes)}
            panel_columns = np.asarray(
                [unique_index[gene] for gene in panel], dtype=np.int64
            )
            panel_mask = np.zeros(len(unique_genes), dtype=bool)
            panel_mask[panel_columns] = True
            mapped_effects, incumbent_mapping = load_incumbent_effects(
                canonical_genes, unique_genes, panel, protocol
            )
            replayed_libraries = replay_libraries(
                control_cell_libraries,
                len(panel),
                cells=int(protocol["transport"]["cells"]),
                seed=int(protocol["transport"]["library_seed"]),
                maximum=int(protocol["transport"]["maximum_library"]),
            )
            incumbent_profiles, incumbent_totals, incumbent_predictions = (
                build_incumbent_rows(
                    exp128,
                    mapped_effects,
                    target_panel_indices,
                    b_ctrl,
                    panel_columns,
                    replayed_libraries,
                )
            )
            bins = response_bins(control_abundance[panel_columns])
            draws = int(protocol["nulls"]["draws_each"])
            target_null, target_response_null, response_null = freeze_nulls(
                target_panel_indices,
                bins,
                draws,
                target_seed=int(protocol["nulls"]["target_label_seed"]),
                response_seed=int(protocol["nulls"]["response_label_seed"]),
            )
            np.save(staging / "response-bins.npy", bins, allow_pickle=False)
            np.save(
                staging / "target-null-row-permutations.npy",
                target_null,
                allow_pickle=False,
            )
            np.save(
                staging / "target-null-response-permutations.npy",
                target_response_null,
                allow_pickle=False,
            )
            np.save(
                staging / "response-null-permutations.npy",
                response_null,
                allow_pickle=False,
            )
            np.save(
                staging / "replayed-libraries-all-300.npy",
                replayed_libraries,
                allow_pickle=False,
            )
            np.save(
                staging / "incumbent-panel-predictions.npy",
                incumbent_predictions,
                allow_pickle=False,
            )
            frozen_names = (
                "response-bins.npy",
                "target-null-row-permutations.npy",
                "target-null-response-permutations.npy",
                "response-null-permutations.npy",
                "replayed-libraries-all-300.npy",
                "incumbent-panel-predictions.npy",
            )
            null_freeze = {
                "files": {
                    name: sha256_file(staging / name) for name in frozen_names
                },
                "experiment_135_targeting_panel_aggregation_started": False,
                "prior_raw_chunk_streaming_disclosed": True,
                "telohaec_panel_responses_previously_aggregated_or_analyzed": False,
            }
            (staging / "pretarget-null-freeze.json").write_text(
                json.dumps(null_freeze, indent=2, sort_keys=True) + "\n"
            )

        access_receipt = {
            "experiment": "135-telohaec-atlas-confirmation",
            "telohaec_panel_response_aggregation_started": True,
            "pretarget_null_freeze_sha256": sha256_file(
                staging / "pretarget-null-freeze.json"
            ),
            "source_sha256": protocol["source"]["expected_matrix_sha256"],
            "prior_raw_chunk_streaming_disclosed": True,
        }
        (staging / "target-value-access-receipt.json").write_text(
            json.dumps(access_receipt, indent=2, sort_keys=True) + "\n"
        )
        with h5py.File(H5AD, "r") as handle:
            panel_sums, target_totals, target_spacers, target_rows = scan_targets(
                handle,
                target_order,
                guide_contract["targets"],
                guide_contract["expected_cells"],
                raw_to_panel,
            )

        control_panel = b_ctrl[panel_columns]
        guide_panel_profiles = np.log1p(TS * panel_sums / target_totals[:, None])
        truth, guide_diagnostics = guide_truth_and_diagnostics(
            target_order,
            target_panel_indices,
            target_rows,
            guide_panel_profiles,
            control_panel,
            guide_contract["expected_cells"],
            target_spacers,
        )
        require(
            int(guide_diagnostics["target_cell_counts"].min())
            >= int(protocol["truth"]["minimum_eligible_cells_per_target"]),
            "observed target-cell total fell below the frozen floor",
        )

        baseline_credit, baseline_precision = baseline_metrics(
            target_panel_indices, truth, incumbent_predictions
        )
        candidate_credit, candidate_precision, observed_audits = evaluate_rows(
            exp128,
            target_order,
            panel,
            target_panel_indices,
            truth,
            residual,
            confidence,
            incumbent_profiles,
            incumbent_totals,
            b_ctrl,
            panel_columns,
            panel_mask,
        )
        credit_delta = candidate_credit - baseline_credit
        precision_delta = candidate_precision - baseline_precision
        observed_credit_delta = float(credit_delta.mean())
        observed_precision_delta = float(precision_delta.mean())

        target_credit_null = np.empty(draws, dtype=np.float64)
        target_precision_null = np.empty(draws, dtype=np.float64)
        response_credit_null = np.empty(draws, dtype=np.float64)
        response_precision_null = np.empty(draws, dtype=np.float64)
        all_null_invariants = True
        for draw in range(draws):
            credit, precision, audits = evaluate_rows(
                exp128,
                target_order,
                panel,
                target_panel_indices,
                truth,
                residual,
                confidence,
                incumbent_profiles,
                incumbent_totals,
                b_ctrl,
                panel_columns,
                panel_mask,
                target_mapping=target_null[draw],
                response_mapping=target_response_null[draw],
            )
            target_credit_null[draw] = float(np.mean(credit - baseline_credit))
            target_precision_null[draw] = float(
                np.mean(precision - baseline_precision)
            )
            all_null_invariants &= all(
                record["all_invariants_pass"] for record in audits
            )
            credit, precision, audits = evaluate_rows(
                exp128,
                target_order,
                panel,
                target_panel_indices,
                truth,
                residual,
                confidence,
                incumbent_profiles,
                incumbent_totals,
                b_ctrl,
                panel_columns,
                panel_mask,
                response_mapping=response_null[draw],
            )
            response_credit_null[draw] = float(np.mean(credit - baseline_credit))
            response_precision_null[draw] = float(
                np.mean(precision - baseline_precision)
            )
            all_null_invariants &= all(
                record["all_invariants_pass"] for record in audits
            )

        bootstrap_rng = np.random.default_rng(int(protocol["bootstrap"]["seed"]))
        bootstrap_indices = bootstrap_rng.integers(
            0,
            len(credit_delta),
            size=(
                int(protocol["bootstrap"]["replicates"]),
                len(credit_delta),
            ),
        )
        bootstrap_credit_delta = credit_delta[bootstrap_indices].mean(axis=1)
        bootstrap_precision_delta = precision_delta[bootstrap_indices].mean(axis=1)
        credit_lower = float(
            np.quantile(bootstrap_credit_delta, 0.025, method="linear")
        )
        precision_lower = float(
            np.quantile(bootstrap_precision_delta, 0.025, method="linear")
        )
        nulls = {
            "target_label": metric_null_summary(
                target_credit_null,
                target_precision_null,
                observed_credit_delta,
                observed_precision_delta,
            ),
            "response_label": metric_null_summary(
                response_credit_null,
                response_precision_null,
                observed_credit_delta,
                observed_precision_delta,
            ),
        }
        sign_retention = np.asarray(
            [record["sign_retention"] for record in observed_audits],
            dtype=np.float64,
        )
        thresholds = protocol["promotion_gates"]
        provenance_pass = all(provenance.values())
        all_observed_invariants = all(
            record["all_invariants_pass"] for record in observed_audits
        )
        gates = {
            "exact_candidate_and_dependency_hashes": provenance_pass,
            "evaluable_targets_exactly_65": len(target_order) == 65,
            "truth_responses_per_target_exactly_299": truth.shape == (65, 300),
            "candidate_absolute_signed_credit_strictly_positive": float(
                candidate_credit.mean()
            )
            > 0.0,
            "candidate_absolute_precision_at_least_0p5": float(
                candidate_precision.mean()
            )
            >= float(thresholds["candidate_absolute_precision_at_least"]),
            "macro_signed_credit_delta_at_least_0p005": observed_credit_delta
            >= float(thresholds["macro_signed_credit_delta_at_least"]),
            "macro_precision_delta_at_least_0p01": observed_precision_delta
            >= float(thresholds["macro_precision_delta_at_least"]),
            "signed_credit_delta_bootstrap_lower_above_zero": credit_lower
            > float(
                thresholds[
                    "signed_credit_delta_bootstrap_lower_95_strictly_above"
                ]
            ),
            "precision_delta_bootstrap_lower_above_zero": precision_lower
            > float(
                thresholds["precision_delta_bootstrap_lower_95_strictly_above"]
            ),
            "signed_credit_delta_above_both_null_p95": all(
                observed_credit_delta > record["signed_credit_delta_p95"]
                for record in nulls.values()
            ),
            "precision_delta_above_both_null_p95": all(
                observed_precision_delta > record["precision_delta_p95"]
                for record in nulls.values()
            ),
            "all_null_p_values_at_most_0p05": all(
                record[key]
                <= float(thresholds["both_add_one_null_p_values_at_most"])
                for record in nulls.values()
                for key in ("signed_credit_delta_p", "precision_delta_p")
            ),
            "median_transport_sign_retention_at_least_0p9": float(
                np.median(sign_retention)
            )
            >= float(thresholds["median_transport_sign_retention_at_least"]),
            "all_integer_transport_invariants": bool(
                all_observed_invariants and all_null_invariants
            ),
        }
        result = {
            "schema_version": 2,
            "experiment": "135-telohaec-atlas-confirmation",
            "candidate": protocol["candidate"]["name"],
            "comparator": protocol["incumbent"]["name"],
            "dose": protocol["candidate"]["dose"],
            "target_count": len(target_order),
            "panel_response_count": 300,
            "control_guide_count": len(guide_contract["controls"]),
            "candidate_support": {
                "supported_coordinates": int(support.sum()),
                "unsupported_finite_coordinates": int(
                    np.count_nonzero((~support) & np.isfinite(residual))
                ),
                "used_as_transport_gate": False,
            },
            "provenance": provenance,
            "incumbent_mapping": incumbent_mapping,
            "replayed_library": summarize(replayed_libraries),
            "baseline": {
                "macro_signed_credit": float(baseline_credit.mean()),
                "macro_top20_weighted_sign_precision": float(
                    baseline_precision.mean()
                ),
            },
            "candidate_metrics": {
                "macro_signed_credit": float(candidate_credit.mean()),
                "macro_top20_weighted_sign_precision": float(
                    candidate_precision.mean()
                ),
            },
            "candidate_minus_baseline": {
                "macro_signed_credit": observed_credit_delta,
                "macro_top20_weighted_sign_precision": observed_precision_delta,
            },
            "guide_diagnostics": {
                key: summarize(value) for key, value in guide_diagnostics.items()
            },
            "bootstrap_scope": (
                "Stability within the availability-selected 65 TeloHAEC targets only; "
                "not an IID panel-wide generalization interval."
            ),
            "bootstrap_delta": {
                "signed_credit_lower_95": credit_lower,
                "signed_credit_upper_95": float(
                    np.quantile(bootstrap_credit_delta, 0.975, method="linear")
                ),
                "precision_lower_95": precision_lower,
                "precision_upper_95": float(
                    np.quantile(bootstrap_precision_delta, 0.975, method="linear")
                ),
            },
            "nulls": nulls,
            "transport": {
                "median_sign_retention": float(np.median(sign_retention)),
                "minimum_sign_retention": float(np.min(sign_retention)),
                "maximum_sign_retention": float(np.max(sign_retention)),
                "preclip_negative_coordinates": int(
                    sum(
                        record["preclip_negative_coordinates"]
                        for record in observed_audits
                    )
                ),
                "fallback_rows": int(
                    sum(record["fallback_used"] for record in observed_audits)
                ),
                "observed_invariants_pass": bool(all_observed_invariants),
                "all_null_invariants_pass": bool(all_null_invariants),
            },
            "gates": gates,
            "promotion_passes": bool(all(gates.values())),
            "authorization_if_passed": protocol["authorization_if_passed"],
            "no_challenge_build_or_submission": True,
            "post_outcome_value_release": "aggregate_only",
            "sealed_value_hashes": {
                "truth_panel": sha256_array(truth.astype("<f4")),
                "incumbent_signed_credit": sha256_array(baseline_credit),
                "candidate_signed_credit": sha256_array(candidate_credit),
                "signed_credit_delta": sha256_array(credit_delta),
                "incumbent_precision": sha256_array(baseline_precision),
                "candidate_precision": sha256_array(candidate_precision),
                "precision_delta": sha256_array(precision_delta),
                "target_null_credit_delta": sha256_array(target_credit_null),
                "target_null_precision_delta": sha256_array(
                    target_precision_null
                ),
                "response_null_credit_delta": sha256_array(
                    response_credit_null
                ),
                "response_null_precision_delta": sha256_array(
                    response_precision_null
                ),
                "bootstrap_credit_delta": sha256_array(
                    bootstrap_credit_delta
                ),
                "bootstrap_precision_delta": sha256_array(
                    bootstrap_precision_delta
                ),
                **{
                    key: sha256_array(value)
                    for key, value in guide_diagnostics.items()
                },
            },
        }
        (staging / "target-order.txt").write_text(
            "\n".join(target_order) + "\n", encoding="utf-8"
        )
        (staging / "results.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n"
        )
        manifest = {}
        for path in sorted(staging.iterdir()):
            if path.name == "artifact-manifest.json":
                continue
            manifest[path.name] = {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        (staging / "artifact-manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "experiment": "135-telohaec-atlas-confirmation",
                    "files": manifest,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        os.replace(staging, OUTPUT)
        for path in OUTPUT.iterdir():
            path.chmod(0o444)
        OUTPUT.chmod(0o555)
        return result
    except Exception:
        if staging.exists():
            failure = staging / "FAILED.txt"
            if not failure.exists():
                failure.write_text(
                    "Experiment 135 stopped after staging began. Preserve this "
                    "directory; do not rerun under the same experiment number.\n"
                )
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--preflight", action="store_true")
    actions.add_argument("--evaluate", action="store_true")
    parser.add_argument("--confirm-value-access", default="")
    args = parser.parse_args()
    if args.preflight:
        protocol = verify_freeze(hash_large_source=False)
        panel = read_axis(PANEL, "target_gene")
        genes = read_axis(GENES, "gene_name")
        contract = load_guide_contract(panel, protocol)
        provenance = verify_candidate_provenance(protocol)
        require(
            sha256_file(INCUMBENT_EFFECTS)
            == protocol["incumbent"]["effect_sha256"],
            "incumbent effect hash drift",
        )
        require(
            sha256_file(INCUMBENT_METADATA)
            == protocol["incumbent"]["metadata_sha256"],
            "incumbent metadata hash drift",
        )
        print(
            json.dumps(
                {
                    "status": "READY",
                    "candidate": protocol["candidate"]["name"],
                    "comparator": protocol["incumbent"]["name"],
                    "dose": protocol["candidate"]["dose"],
                    "targets": len(contract["targets"]),
                    "controls": len(contract["controls"]),
                    "canonical_genes": len(genes),
                    "candidate_provenance": provenance,
                    "telohaec_panel_response_aggregation_or_analysis": False,
                    "prior_raw_chunk_streaming_disclosed": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    require(
        args.confirm_value_access
        == "I_CONFIRM_EXPERIMENT_135_TELOHAEC_PANEL_V1",
        "exact value-access confirmation is required",
    )
    print(json.dumps(run_evaluation(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
