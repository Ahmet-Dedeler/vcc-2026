#!/usr/bin/env python3
"""Metadata preflight and fail-closed value boundary for Experiment 130."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "preregistered_protocol.json"
RECEIPT = HERE / "SOURCE_METADATA_RECEIPT.v1.json"
FREEZE = HERE / "IMPLEMENTATION_FREEZE.json"
BINDING = HERE / "CANDIDATE_BINDING_AMENDMENT.json"
BOUNDARY = HERE / "VALUE_ACCESS_BOUNDARY.v1.json"
SCHEMA = HERE / "CANDIDATE_BINDING_AMENDMENT.schema.json"
TELO_H5AD = ROOT / "data/external/igvf-telohaec/IGVFFI8645YMGX.h5ad"
GUIDE_QC = ROOT / "data/external/igvf-telohaec/guide_qc.csv"
GENES = ROOT / "data/challenge/gene_names.csv"
PANEL = ROOT / "data/challenge/pert_counts.csv"
EXP105_PROVENANCE = ROOT / "experiments/105-telohaec-residual/PROVENANCE.json"
EXP128_FREEZE = ROOT / "experiments/128-panel-integer-budget-transport/IMPLEMENTATION_FREEZE.json"
EXP128_RUN = ROOT / "experiments/128-panel-integer-budget-transport/run.py"


class ExperimentError(RuntimeError):
    pass


class ValueAccessLocked(ExperimentError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ExperimentError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(row[column].strip() for row in reader)
    require(len(values) == len(set(values)), f"{path}: duplicate axis labels")
    return values


def h5ad_structure(path: Path) -> dict[str, Any]:
    """Inspect group/dataset metadata only; never slice X, layers, or obs values."""
    with h5py.File(path, "r") as handle:
        require(set(handle.keys()) == {"X", "obs", "var"}, "H5AD root structure changed")
        x = handle["X"]
        require(isinstance(x, h5py.Group), "X is no longer a sparse group")
        return {
            "root_keys": sorted(handle.keys()),
            "x": {
                "encoding_type": str(x.attrs["encoding-type"]),
                "shape": [int(v) for v in x.attrs["shape"]],
                "datasets": {
                    key: {"shape": list(x[key].shape), "dtype": str(x[key].dtype)}
                    for key in sorted(x.keys())
                },
            },
            "obs": {
                "keys": sorted(handle["obs"].keys()),
                "column_order": [str(v) for v in handle["obs"].attrs["column-order"]],
            },
            "var": {
                "keys": sorted(handle["var"].keys()),
                "column_order": [str(v) for v in handle["var"].attrs["column-order"]],
            },
            "layers_present": "layers" in handle,
        }


def guide_mapping(panel: Sequence[str]) -> dict[str, Any]:
    panel_set = set(panel)
    targets: dict[str, list[dict[str, Any]]] = {target: [] for target in panel}
    controls: list[dict[str, Any]] = []
    with GUIDE_QC.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "guide_id", "spacer", "guide_type", "clean_target",
            "observed_cells", "target_eligible_cells", "eligible", "eligibility_note"
        }
        require(reader.fieldnames is not None and required.issubset(reader.fieldnames), "guide_qc schema changed")
        for row in reader:
            if row["eligible"] != "1":
                continue
            item = {
                "guide_id": row["guide_id"],
                "spacer": row["spacer"],
                "guide_type": row["guide_type"],
                "observed_cells": int(row["observed_cells"]),
            }
            target = row["clean_target"]
            if target in panel_set:
                targets[target].append(item)
            elif row["eligibility_note"] == "eligible_control":
                controls.append(item)
    retained = {target: guides for target, guides in targets.items() if guides}
    require(len(retained) == 65, f"expected 65 retained targets, found {len(retained)}")
    for target, guides in retained.items():
        require(len(guides) >= 2, f"{target}: fewer than two eligible guides")
        require(sum(g["observed_cells"] for g in guides) >= 20, f"{target}: fewer than 20 eligible cells")
    require(len(controls) == 549, f"expected 549 controls, found {len(controls)}")
    return {"targets": retained, "controls": controls}


def build_metadata_receipt() -> dict[str, Any]:
    genes = csv_axis(GENES, "gene_name")
    panel = csv_axis(PANEL, "target_gene")
    require(len(genes) == 18533, "canonical gene count changed")
    require(len(panel) == 300, "canonical panel count changed")
    provenance = read_json(EXP105_PROVENANCE)
    source = provenance["local_files"]["IGVFFI8645YMGX.h5ad"]
    stat = TELO_H5AD.stat()
    require(stat.st_size == source["size_bytes_observed"], "TeloHAEC H5AD size changed")
    mapping = guide_mapping(panel)
    retained = [target for target in panel if target in mapping["targets"]]
    guide_counts = [len(mapping["targets"][target]) for target in retained]
    cell_counts = [sum(g["observed_cells"] for g in mapping["targets"][target]) for target in retained]
    receipt = {
        "schema_version": 1,
        "experiment": "130-telohaec-panel-confirmation",
        "scope": "metadata_only_no_X_or_layer_values_read",
        "source_file": {
            "path": str(TELO_H5AD.relative_to(ROOT)),
            "size_bytes": stat.st_size,
            "md5_from_verified_experiment_105_receipt": source["md5_observed"],
            "sha256_from_verified_experiment_105_receipt": source["sha256_observed"],
            "content_hash_not_recomputed_in_experiment_130": True
        },
        "structure": h5ad_structure(TELO_H5AD),
        "bindings": {
            "experiment_105_provenance": {"path": str(EXP105_PROVENANCE.relative_to(ROOT)), "sha256": sha256_file(EXP105_PROVENANCE)},
            "guide_qc": {"path": str(GUIDE_QC.relative_to(ROOT)), "sha256": sha256_file(GUIDE_QC)},
            "canonical_genes": {"path": str(GENES.relative_to(ROOT)), "sha256": sha256_file(GENES), "count": len(genes)},
            "canonical_panel": {"path": str(PANEL.relative_to(ROOT)), "sha256": sha256_file(PANEL), "count": len(panel)},
            "experiment_128_freeze": {"path": str(EXP128_FREEZE.relative_to(ROOT)), "sha256": sha256_file(EXP128_FREEZE)},
            "experiment_128_run": {"path": str(EXP128_RUN.relative_to(ROOT)), "sha256": sha256_file(EXP128_RUN)}
        },
        "panel_target_order": retained,
        "guide_mapping": mapping,
        "summary": {
            "targets": len(retained),
            "controls": len(mapping["controls"]),
            "eligible_guides_min": min(guide_counts),
            "eligible_guides_median": float(np.median(guide_counts)),
            "eligible_guides_max": max(guide_counts),
            "eligible_cells_min": min(cell_counts),
            "eligible_cells_median": float(np.median(cell_counts)),
            "eligible_cells_max": max(cell_counts)
        },
        "biological_value_access": False
    }
    return receipt


def write_receipt() -> None:
    require(not RECEIPT.exists(), "refusing to overwrite metadata receipt")
    RECEIPT.write_text(json.dumps(build_metadata_receipt(), indent=2, sort_keys=True) + "\n")
    RECEIPT.chmod(0o444)


def validate_freeze() -> None:
    require(FREEZE.is_file(), "Experiment 130 implementation freeze missing")
    freeze = read_json(FREEZE)
    for relative, digest in freeze["files"].items():
        path = ROOT / relative
        require(path.is_file() and not path.is_symlink(), f"missing frozen file: {relative}")
        require(sha256_file(path) == digest, f"frozen file changed: {relative}")


def validate_candidate_binding() -> dict[str, Any]:
    if not BINDING.is_file():
        raise ValueAccessLocked("candidate-binding amendment is absent; biological value access remains locked")
    validate_freeze()
    binding = read_json(BINDING)
    require(binding.get("experiment") == "130-telohaec-panel-confirmation", "wrong binding experiment")
    require(binding.get("created_before_telohaec_value_access") is True, "binding timing not attested")
    e125 = binding["experiment_125"]
    result_path = ROOT / e125["results_path"]
    manifest_path = ROOT / e125["artifact_manifest_path"]
    receipt_path = ROOT / e125["result_receipt_path"]
    for path, digest in (
        (result_path, e125["results_sha256"]),
        (manifest_path, e125["artifact_manifest_sha256"]),
        (receipt_path, e125["result_receipt_sha256"]),
    ):
        require(path.is_file() and not path.is_symlink(), f"missing immutable Experiment 125 artifact: {path}")
        require(sha256_file(path) == digest, f"Experiment 125 hash mismatch: {path}")
    result = read_json(result_path)
    require(result.get("decision") == "GO" and result.get("promotion_passes") is True, "Experiment 125 did not pass")
    require(result.get("gates") and all(result["gates"].values()), "Experiment 125 gate failure")
    counts = result["selected_candidate_counts_across_15_outer_target_folds"]
    order = result["candidate_order"]
    selected = max(order, key=lambda name: (int(counts[name]), -order.index(name)))
    selection = binding["selection"]
    require(selection["candidate_name"] == selected, "binding violates frozen modal selection rule")
    require(int(selection["selection_count"]) == int(counts[selected]) >= 8, "candidate lacks 8/15 selection support")
    definition = next(item for item in result["candidate_definitions"] if item["name"] == selected)
    require(definition["transport"] == selection["transport"] == "production_log1p", "candidate transport is not production_log1p")
    attest = binding["attestations"]
    require(all(attest.values()), "binding attestations incomplete")
    for key in ("residual", "confidence", "support", "panel_axis"):
        path = ROOT / binding["candidate_artifacts"][f"{key}_path"]
        digest = binding["candidate_artifacts"][f"{key}_sha256"]
        require(path.is_file() and not path.is_symlink(), f"missing candidate {key}")
        require(sha256_file(path) == digest, f"candidate {key} hash mismatch")
    return binding


def open_telohaec_values() -> h5py.File:
    """The only authorized biological-value entry point."""
    validate_candidate_binding()
    return h5py.File(TELO_H5AD, "r")


def target_metrics(prediction: np.ndarray, confidence: np.ndarray, truth: np.ndarray) -> tuple[float, float]:
    require(prediction.shape == confidence.shape == truth.shape, "metric shapes differ")
    valid = np.isfinite(truth)
    require(int(valid.sum()) >= 20, "too few truth responses")
    indices = np.flatnonzero(valid)
    denominator = float(np.abs(truth[indices]).sum())
    require(denominator > 0, "zero truth magnitude")
    credit = float(np.sum(np.sign(truth[indices]) * prediction[indices]) / denominator)
    top_n = int(math.ceil(0.20 * len(indices)))
    chosen = indices[np.lexsort((indices, -confidence[indices]))[:top_n]]
    weights = np.abs(truth[chosen])
    require(float(weights.sum()) > 0, "zero top-20 truth magnitude")
    correct = (np.sign(prediction[chosen]) != 0) & (np.sign(prediction[chosen]) == np.sign(truth[chosen]))
    precision = float(np.sum(weights * correct) / np.sum(weights))
    return credit, precision


def derangement(rng: np.random.Generator, n: int) -> np.ndarray:
    require(n > 1, "derangement requires n > 1")
    base = np.arange(n)
    while True:
        draw = rng.permutation(n)
        if np.all(draw != base):
            return draw


def main() -> None:
    parser = argparse.ArgumentParser()
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--metadata-preflight", action="store_true")
    actions.add_argument("--check-lock", action="store_true")
    actions.add_argument("--evaluate", action="store_true")
    args = parser.parse_args()
    if args.metadata_preflight:
        write_receipt()
        print(json.dumps({"receipt": str(RECEIPT), "biological_value_access": False}, indent=2))
    elif args.check_lock:
        try:
            validate_candidate_binding()
        except ValueAccessLocked as exc:
            print(json.dumps({"status": "LOCKED", "reason": str(exc)}, indent=2))
            return
        raise ExperimentError("unexpected unlocked state")
    else:
        handle = open_telohaec_values()
        handle.close()
        raise ExperimentError("value gate opened, but biological evaluation is intentionally not executed by this metadata-only task")


if __name__ == "__main__":
    main()
