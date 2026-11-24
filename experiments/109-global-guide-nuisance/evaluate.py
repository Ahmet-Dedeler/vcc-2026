#!/usr/bin/env python3
"""Held-context evaluation scaffold for experiment 109.

Only contract checks and pure adapter/gate helpers are implemented. There is
no command that opens held truth or runs scoring in this source-only scaffold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping
from pathlib import Path

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import fit_freeze
import model
import numpy as np

HERE = Path(__file__).resolve().parent
EXPERIMENT = "109-global-guide-nuisance"


def canonicalize_panel_rows(
    values: np.ndarray,
    target_ids: tuple[str, ...],
    canonical_targets: tuple[str, ...],
    *,
    allow_remap: bool = False,
) -> np.ndarray:
    """Require canonical row order unless the caller explicitly requests remapping."""
    rows = np.asarray(values)
    targets = tuple(map(str, target_ids))
    canonical = tuple(map(str, canonical_targets))
    if len(targets) != len(set(targets)) or set(targets) != set(canonical):
        raise ValueError("target rows do not exactly cover the canonical VCC panel")
    if rows.ndim < 1 or rows.shape[0] != len(targets):
        raise ValueError("values do not align with the supplied target axis")
    if targets == canonical:
        return rows
    if not allow_remap:
        raise ValueError(
            "target rows are not in canonical VCC panel order; "
            "explicit remapping is required"
        )
    source_row = {target: index for index, target in enumerate(targets)}
    return rows[
        np.asarray([source_row[target] for target in canonical], dtype=np.int64)
    ]


def verify_frozen_authorization(frozen_dir: Path) -> dict[str, object]:
    """Refuse held evaluation until its separately reviewed adapter exists."""
    lifecycle = fit_freeze.read_lifecycle(frozen_dir)
    if lifecycle is not None and lifecycle.get("state") == "held_truth_opened":
        raise ValueError("held truth was already opened; evaluation is not resumable")
    raise RuntimeError(
        "held evaluation is deliberately disabled in Experiment 109 pre-source "
        "hardening; no held truth, scorer, or production adapter may be opened"
    )


def _verify_frozen_authorization_disabled(frozen_dir: Path) -> dict[str, object]:
    """Retained unreachable verifier for the later held-evaluation review."""
    manifest_path = frozen_dir / "freeze-manifest.json"
    authorization_path = frozen_dir / "held-evaluation-authorized.json"
    manifest = json.loads(manifest_path.read_text())
    authorization = json.loads(authorization_path.read_text())
    if authorization.get("authorized") is not True:
        raise ValueError("held evaluation is not authorized")
    if authorization.get("experiment") != EXPERIMENT:
        raise ValueError("authorization names the wrong experiment")
    if authorization.get("purpose") != "held_context_evaluation_only":
        raise ValueError("authorization purpose is not held-context evaluation")
    if authorization.get("no_held_truth_opened_during_fit") is not True:
        raise ValueError("fit phase did not attest a closed held set")
    if manifest.get("experiment") != EXPERIMENT:
        raise ValueError("wrong frozen experiment")
    if manifest.get("status") != "frozen_before_held_truth":
        raise ValueError("freeze status is not pre-held")
    if manifest.get("post_held_refit_possible") is not False:
        raise ValueError("post-held refit guard is absent")
    if manifest.get("candidate_grid_frozen_before_expression") is not True:
        raise ValueError("candidate grid was not frozen before expression")
    if manifest.get("held_context_truth_opened") is not False:
        raise ValueError("manifest does not attest a closed held set")
    if authorization.get("freeze_manifest_sha256") != fit_freeze.sha256_file(
        manifest_path
    ):
        raise ValueError("freeze manifest changed after authorization")
    if manifest.get("protocol_sha256") != fit_freeze.sha256_file(
        fit_freeze.PROTOCOL_PATH
    ):
        raise ValueError("preregistered protocol changed after freeze")
    fit_freeze.validate_source_test(manifest.get("source_test", {}))
    source_test_path = frozen_dir / "source-test.json"
    if json.loads(source_test_path.read_text()) != manifest["source_test"]:
        raise ValueError("source-test artifact differs from freeze manifest")

    source_hashes = manifest.get("source_code_sha256")
    if not isinstance(source_hashes, Mapping) or set(source_hashes) != set(
        fit_freeze.SOURCE_FILES
    ):
        raise ValueError("source-code hash set is incomplete")
    for name, path in fit_freeze.SOURCE_FILES.items():
        if source_hashes[name] != fit_freeze.sha256_file(path):
            raise ValueError(f"source code changed after freeze: {name}")

    protocol = fit_freeze.validate_preregistration()
    splits = manifest.get("splits")
    if not isinstance(splits, Mapping) or set(splits) != set(
        fit_freeze.REQUIRED_SPLITS
    ):
        raise ValueError("frozen split set is incomplete")
    expected_split_counts = protocol["frozen_source_audit_expectations"]["split_counts"]
    for name in fit_freeze.REQUIRED_SPLITS:
        split = splits[name]
        if int(split["count"]) != int(expected_split_counts[name]):
            raise ValueError(f"{name} split count changed")
        if not fit_freeze.is_sha256(split["targets_sha256"]):
            raise ValueError(f"{name} split hash is invalid")
    split_values = {
        name: tuple(map(str, splits[name]["targets"]))
        for name in fit_freeze.REQUIRED_SPLITS
    }

    artifacts = manifest.get("frozen_artifacts")
    if not isinstance(artifacts, Mapping) or set(artifacts) != set(
        fit_freeze.REQUIRED_FROZEN_ARTIFACTS
    ):
        raise ValueError("frozen artifact set is incomplete")
    for filename, expected in artifacts.items():
        if not fit_freeze.is_sha256(expected):
            raise ValueError(f"invalid frozen artifact hash: {filename}")
        if fit_freeze.sha256_file(frozen_dir / filename) != expected:
            raise ValueError(f"frozen artifact changed: {filename}")
    fit_freeze.validate_source_evidence(
        frozen_dir,
        manifest["source_test"],
        split_values,
    )

    axes = manifest.get("axes")
    if not isinstance(axes, Mapping):
        raise TypeError("frozen axis contracts are missing")
    expected_axis_names = {
        "centering_targets",
        "gene_axis",
        "panel_targets",
        "development_targets",
        "response_axis",
        "perturbmatch_targets",
        "perturbmatch_genes",
    }
    if set(axes) != expected_axis_names:
        raise ValueError("frozen axis contract set is incomplete")
    contracts = {
        name: model.AxisContract(
            count=int(axes[name]["count"]),
            sha256=str(axes[name]["sha256"]),
        )
        for name in expected_axis_names
    }
    if not all(fit_freeze.is_sha256(value.sha256) for value in contracts.values()):
        raise ValueError("invalid frozen axis hash")
    expectations = protocol["frozen_source_audit_expectations"]
    if contracts["centering_targets"].count != int(
        expectations["source_sequence_supported_singletons"]
    ):
        raise ValueError("singleton centering population count changed")
    if contracts["panel_targets"].count != 300:
        raise ValueError("VCC panel target count changed")
    if contracts["development_targets"].count != int(
        expectations["identifiable_nonpanel_multiguide_targets"]
    ):
        raise ValueError("development target count changed")

    axis_payload = json.loads((frozen_dir / "axes.json").read_text())
    if set(axis_payload) != {
        "centering_targets",
        "gene_axis",
        "panel_targets",
        "development_targets",
        "response_genes",
        "perturbmatch_targets",
        "perturbmatch_genes",
    }:
        raise ValueError("axes.json field set is incomplete")
    model.validate_axis(
        axis_payload["centering_targets"],
        contracts["centering_targets"],
        label="singleton centering axis",
        require_sorted=True,
    )
    response_axis = model.authenticate_response_axis(
        axis_payload["gene_axis"],
        axis_payload["panel_targets"],
        axis_payload["development_targets"],
        gene_axis_contract=contracts["gene_axis"],
        panel_targets_contract=contracts["panel_targets"],
        development_targets_contract=contracts["development_targets"],
        response_axis_contract=contracts["response_axis"],
    )
    response_genes = [
        gene
        for gene, keep in zip(
            response_axis.gene_axis,
            response_axis.response_mask,
            strict=True,
        )
        if keep
    ]
    if response_genes != axis_payload["response_genes"]:
        raise ValueError("axes.json response axis is not derived from exclusions")
    canonical = fit_freeze.load_canonical_source_axes(protocol)
    canonical_to_frozen = {
        "singleton_centering_targets": "centering_targets",
        "gene_axis": "gene_axis",
        "panel_targets": "panel_targets",
        "development_targets": "development_targets",
        "response_axis": "response_genes",
        "perturbmatch_targets": "perturbmatch_targets",
        "perturbmatch_genes": "perturbmatch_genes",
    }
    for canonical_name, frozen_name in canonical_to_frozen.items():
        if canonical["axes"][canonical_name]["values"] != axis_payload[frozen_name]:
            raise ValueError(f"{frozen_name} differs from canonical source axes")

    expected_splits = fit_freeze.development_split(
        response_axis.development_targets,
        response_axis.panel_targets,
    )
    for name in fit_freeze.REQUIRED_SPLITS:
        if list(expected_splits[name]) != splits[name].get("targets"):
            raise ValueError(f"{name} split differs from canonical deterministic split")
        if (
            fit_freeze.sha256_lines(expected_splits[name])
            != splits[name]["targets_sha256"]
        ):
            raise ValueError(f"{name} split hash differs from target list")

    frozen_model = fit_freeze.load_frozen_model(
        frozen_dir / "selected-model.npz",
        expected_role="deployment_model",
    )
    fit_freeze.validate_frozen_model_structure(
        frozen_model,
        response_axis=response_axis,
        expected_fit_targets_sha256=contracts["development_targets"].sha256,
    )
    return manifest


def map_response_correction_to_source_slot(
    fitted: model.StableSubspaceModel,
    response_correction: np.ndarray,
    target_ids: tuple[str, ...],
    gene_axis: tuple[str, ...],
    *,
    remap_targets_to_canonical: bool = False,
) -> np.ndarray:
    """Return canonical-panel rows, rejecting or explicitly remapping input order."""
    source_axes = fit_freeze.load_canonical_source_axes()
    _, gene_contract = fit_freeze.canonical_axis("gene_axis", source_axes)
    panel_targets, _ = fit_freeze.canonical_axis("panel_targets", source_axes)
    pm_targets, _ = fit_freeze.canonical_axis("perturbmatch_targets", source_axes)
    pm_genes, _ = fit_freeze.canonical_axis("perturbmatch_genes", source_axes)
    observed_genes = model.validate_axis(
        gene_axis,
        gene_contract,
        label="K562 source gene axis",
    )
    targets = tuple(map(str, target_ids))
    correction = np.asarray(response_correction)
    if correction.shape != (len(targets), fitted.response_width):
        raise ValueError("response correction shape differs from supplied target axis")
    correction = canonicalize_panel_rows(
        correction,
        targets,
        panel_targets,
        allow_remap=remap_targets_to_canonical,
    )
    targets = panel_targets
    if fitted.gene_axis_sha256 != gene_contract.sha256:
        raise ValueError("frozen model and K562 source gene axes differ")
    full = model.expand_response_corrections(fitted, correction)
    if full.shape != (len(targets), len(observed_genes)):
        raise ValueError("expanded correction shape differs from canonical axes")
    pm_target_set = set(pm_targets)
    pm_gene_set = set(pm_genes)
    return model.apply_k562_slot_retention(
        full,
        np.asarray([target in pm_target_set for target in targets], dtype=bool),
        np.asarray([gene in pm_gene_set for gene in observed_genes], dtype=bool),
    )


def promotion_decision(
    metrics: Mapping[str, float | bool],
    protocol: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Apply the frozen scalar held-context gates without tuning."""
    frozen = (
        fit_freeze.validate_preregistration() if protocol is None else dict(protocol)
    )
    gates = frozen["promotion_gates"]
    checks = {
        "source_test": bool(metrics["source_test_passed"]),
        "source_bootstrap": bool(metrics["source_bootstrap_passed"]),
        "source_nulls": bool(metrics["source_nulls_passed"]),
        "correction_mean": float(metrics["correction_mean_abs_max"])
        <= float(gates["correction_mean_over_source_singleton_support_abs_max"]),
        "unsupported_zero": float(
            metrics["correction_outside_supported_singletons_abs_max"]
        )
        <= float(gates["correction_outside_supported_singletons_abs_max"]),
        "feature_nullspace_invariance": float(
            metrics["feature_nullspace_prediction_change_abs_max"]
        )
        <= float(gates["feature_nullspace_prediction_change_abs_max"]),
        "excluded_coordinates": float(metrics["excluded_coordinate_abs_max"])
        <= float(gates["excluded_coordinate_abs_max"]),
        "eta_zero_identity": bool(metrics["eta_zero_exact_identity"]),
        "h1_pds": float(metrics["h1_raw_pds_delta"])
        >= float(gates["h1_raw_pds_delta_minimum"]),
        "hct116_pds": float(metrics["hct116_raw_pds_delta"])
        >= float(gates["hct116_raw_pds_delta_minimum"]),
        "hek293t_pds": float(metrics["hek293t_raw_pds_delta"])
        >= float(gates["hek293t_raw_pds_delta_minimum"]),
        "h1_directional": float(metrics["h1_directional_credit_delta"])
        >= float(gates["h1_directional_credit_delta_minimum"]),
        "hek293t_directional": float(metrics["hek293t_directional_credit_delta"])
        >= float(gates["hek293t_directional_credit_delta_minimum"]),
        "macro_pds": float(metrics["macro_raw_pds_delta"])
        >= float(gates["equal_donor_macro_raw_pds_delta_minimum"]),
        "singleton_macro_pds": float(metrics["singleton_macro_raw_pds_delta"])
        >= float(gates["singleton_only_macro_raw_pds_delta_minimum"]),
        "above_label_null": bool(metrics["above_frozen_label_null_p95"]),
        "baseline_exact": bool(metrics["exact_baseline_reconstruction"]),
        "coverage_unchanged": bool(metrics["unchanged_target_coverage"]),
    }
    if "kolf_raw_pds_delta" in metrics:
        checks["kolf_diagnostic"] = float(metrics["kolf_raw_pds_delta"]) >= float(
            gates["kolf_raw_pds_delta_minimum_if_reported"]
        )
    return {
        "decision": "PROMOTE" if all(checks.values()) else "REJECT",
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-contract",
        action="store_true",
        help="validate preregistration and evaluator source only",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.check_contract:
        raise SystemExit(
            "source-only scaffold: held evaluation is unavailable; use --check-contract"
        )
    protocol = fit_freeze.validate_preregistration()
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(
        json.dumps(
            {
                "experiment": EXPERIMENT,
                "status": protocol["status"],
                "evaluator_sha256": source_hash,
                "held_truth_opened": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
