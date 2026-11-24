#!/usr/bin/env python3
"""Synthetic tests for experiment 112's frozen correction contracts."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import correction as C  # noqa: E402


def test_control_recognition_is_exact_and_negr1_is_not_control() -> None:
    assert C.corrected_is_control("Non-Targeting")
    for label in (
        "non-targeting",
        "Non-Targeting ",
        "NTC",
        "control",
        "negative-control",
        "neg_1",
        "NEGR1",
        None,
    ):
        assert not C.corrected_is_control(label)
    assert C.legacy_is_control("Non-Targeting")
    assert C.legacy_is_control("NEGR1")


def test_expression_filters_are_exact_and_limited_to_two_source_groups() -> None:
    assert C.expression_filter("Non-Targeting") == (
        ("pass_guide_filter", "=", 1),
        ("gene_target", "=", "Non-Targeting"),
    )
    assert C.expression_filter("NEGR1") == (
        ("pass_guide_filter", "=", 1),
        ("gene_target", "=", "NEGR1"),
    )
    try:
        C.expression_filter("ADNP")
    except C.ContractError:
        pass
    else:
        raise AssertionError("unrelated perturbation expression was authorized")


def test_metadata_counts_identify_only_the_legacy_negr1_excess() -> None:
    counts = C.count_metadata_rows(
        [
            "Non-Targeting",
            "NEGR1",
            "NEGATIVE-1",
            "ADNP",
            "Non-Targeting",
            "NEGR1",
        ],
        [1, 1, 1, 1, 0, 0],
    )
    assert counts == C.ControlCounts(
        exact_control=1,
        negr1=1,
        legacy_control=3,
    )
    for line, expected in C.EXPECTED_LINES.items():
        C.validate_control_counts(
            line,
            C.ControlCounts(
                exact_control=expected["exact_control_count"],
                negr1=expected["negr1_count"],
                legacy_control=expected["legacy_control_count"],
            ),
        )


def test_metadata_count_contract_rejects_extra_legacy_control_labels() -> None:
    expected = C.EXPECTED_LINES["HCT116"]
    try:
        C.validate_control_counts(
            "HCT116",
            C.ControlCounts(
                exact_control=expected["exact_control_count"],
                negr1=expected["negr1_count"],
                legacy_control=expected["legacy_control_count"] + 1,
            ),
        )
    except C.ContractError:
        pass
    else:
        raise AssertionError("an unexplained legacy control row was accepted")


def test_translation_algebra_preserves_target_log_profiles_exactly() -> None:
    old_delta = np.asarray(
        [
            [1.0, -2.0, 0.5],
            [-0.25, 0.75, 4.0],
        ],
        dtype=np.float64,
    )
    old_baseline = np.asarray([8.0, 4.0, 2.0], dtype=np.float64)
    corrected_baseline = np.asarray([4.0, 2.0, 1.0], dtype=np.float64)
    corrected = C.translate_deltas(
        old_delta,
        old_baseline,
        corrected_baseline,
        output_dtype=np.float64,
    )
    incumbent_profiles = C.reconstruct_target_log_profiles(
        old_delta,
        old_baseline,
    )
    corrected_profiles = C.reconstruct_target_log_profiles(
        corrected,
        corrected_baseline,
    )
    np.testing.assert_array_equal(corrected_profiles, incumbent_profiles)


def test_exact_baseline_reconstruction_matches_legacy_builder() -> None:
    exact_control_sum = np.asarray([10, 30, 60, 0], dtype=np.float64)
    negr1_sum = np.asarray([2, 8, 10, 0], dtype=np.float64)
    old_baseline = C.log1p_cpm_from_sum(
        exact_control_sum + negr1_sum,
        dtype=np.float32,
    )
    reconstructed_old, corrected = C.reconstruct_control_baselines(
        exact_control_sum,
        negr1_sum,
    )
    C.require_exact_legacy_baseline(old_baseline, reconstructed_old)
    np.testing.assert_array_equal(
        corrected,
        C.log1p_cpm_from_sum(exact_control_sum, dtype=np.float32),
    )

    changed = reconstructed_old.copy()
    changed[0] = np.nextafter(changed[0], np.float32(np.inf))
    try:
        C.require_exact_legacy_baseline(old_baseline, changed)
    except C.ContractError:
        pass
    else:
        raise AssertionError("a non-bit-identical legacy baseline was accepted")


def test_target_coverage_and_axes_must_remain_unchanged() -> None:
    genes = ["G0", "G1", "G2"]
    targets = ["ADNP", "MYRF"]
    C.validate_axis_identity(genes, list(genes), targets, list(targets))

    for changed_targets in (
        ["MYRF", "ADNP"],
        ["ADNP"],
        ["ADNP", "MYRF", "NEGR1"],
    ):
        try:
            C.validate_axis_identity(genes, genes, targets, changed_targets)
        except C.ContractError:
            pass
        else:
            raise AssertionError(f"target-axis drift was accepted: {changed_targets}")


def fake_source_audit() -> dict[str, object]:
    lines = {}
    for index, (line, expected) in enumerate(C.EXPECTED_LINES.items()):
        lines[line] = {
            "raw_parquet_files": [
                {
                    "path": f"data/external/orion/data/{line}_Batch1.parquet",
                    "size_bytes": 1000 + index,
                    "row_groups": expected["row_groups"],
                    "sha256": str(index + 1) * 64,
                }
            ],
            "row_groups": expected["row_groups"],
            "accepted_label_counts": {
                C.EXACT_CONTROL_LABEL: expected["exact_control_count"]
            },
            "negr1_false_control_count": expected["negr1_count"],
            "legacy_control_count": expected["legacy_control_count"],
            "donor_perturbations_contain_negr1": False,
            "legacy_artifacts": {
                name: chr(ord("a") + index) * 64 for name in C.LEGACY_ARTIFACT_NAMES
            },
        }
    return {
        "schema_version": 1,
        "experiment": C.EXPERIMENT,
        "stage": "metadata-only-source-audit",
        "expression_column_read": False,
        "columns_read": list(C.METADATA_ONLY_COLUMNS),
        "challenge_panel": {
            "sha256": "f" * 64,
            "contains_negr1": False,
        },
        "lines": lines,
    }


def test_source_audit_binds_provenance_without_expression() -> None:
    audit = fake_source_audit()
    C.validate_source_metadata_audit(audit)

    changed = json.loads(json.dumps(audit))
    changed["expression_column_read"] = True
    try:
        C.validate_source_metadata_audit(changed)
    except C.ContractError:
        pass
    else:
        raise AssertionError("a source audit that opened expression was accepted")

    changed = json.loads(json.dumps(audit))
    changed["lines"]["HEK293T"]["negr1_false_control_count"] += 1
    try:
        C.validate_source_metadata_audit(changed)
    except C.ContractError:
        pass
    else:
        raise AssertionError("source count drift was accepted")


def test_preregistered_protocol_keeps_execution_and_submission_closed() -> None:
    protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
    assert protocol["comparison"]["binary_only"] is True
    assert protocol["comparison"]["candidate_grid"] == []
    assert protocol["future_expression_projection"]["authorized_now"] is False
    assert protocol["held_evaluation"]["hct116_direction_metrics"].startswith(
        "diagnostic_only"
    )
    assert protocol["current_authorization"]["read_gene_expression"] is False
    assert protocol["current_authorization"]["open_held_truth"] is False
    assert protocol["current_authorization"]["fit_whitening_operator"] is False
    assert protocol["current_authorization"]["edit_shared_scripts"] is False
    assert protocol["current_authorization"]["challenge_build"] is False
    assert protocol["current_authorization"]["submission"] is False
    assert protocol["current_authorization"]["commit"] is False
    assert protocol["current_authorization"]["push"] is False


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
