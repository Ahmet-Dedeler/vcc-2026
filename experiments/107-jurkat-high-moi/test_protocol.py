from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    assert protocol["status"] == (
        "acquisition_complete_and_protocol_frozen_before_"
        "transcriptome_coordinate_values_are_parsed"
    )

    for record in protocol["immutable_inputs"].values():
        path = record.get("path")
        expected = record.get("sha256")
        if path is None:
            continue
        assert sha256_file(ROOT / path) == expected

    assignment = protocol["cell_and_assignment_freeze"]
    cell_call = assignment["transcriptome_cell_call"]
    assert cell_call["method"] == "streaming_python_port_of_dropletutils_emptydrops"
    assert cell_call["semantic_reference"].endswith("version_1.31.2")
    assert cell_call["parameters"] == {
        "lower": 100,
        "niters": 10000,
        "test.ambient": False,
        "ignore": None,
        "alpha": "Inf",
        "retain": (
            "barcodeRanks(matrix, lower=100, exclude.from=50, window=1, "
            "gradient.threshold=-1)$knee"
        ),
        "round": False,
        "BPPARAM": "SerialParam",
    }
    assert cell_call["accepted_cell"] == "FDR_le_0.01_including_forced_retain_barcodes"
    assert cell_call["rng_contract"]["numpy_version"] == "2.5.2"
    assert cell_call["rng_contract"]["bit_generator"] == "PCG64"
    assert cell_call["rng_contract"]["chunking_must_not_change_results"] is True
    assert cell_call["rng_contract"]["r_dqrng_bit_identity_claimed"] is False
    assert (
        cell_call["port_validation_before_real_coordinates"][
            "real_data_execution_blocked_until_tests_pass"
        ]
        is True
    )

    guide_call = assignment["guide_assignment"]
    assert guide_call["adjusted_p_value_maximum"] == 0.05
    assert guide_call["test"].startswith("scipy.stats.binom.sf(observed_umi - 1")
    assert "exactly_all_83401_guides" in guide_call["multiple_testing"]
    assert guide_call["scipy_version"] == "1.18.0"

    hash_call = assignment["hash_assignment"]
    assert hash_call["source"].startswith("deMULTIplex_1.0.2")
    assert hash_call["q_values"] == {
        "start": 0.01,
        "stop": 0.99,
        "step": 0.02,
        "inclusive": True,
    }
    assert "for_rounds_1_through_3" in hash_call["round_state_machine"]
    assert hash_call["early_stop"].startswith("if_a_round_produces_no_new")
    assert assignment["hash_assignment"]["negative_rescue"] is False
    assert (
        "expression_qc_retention_by_guide_target_channel_and_condition"
        in assignment["required_reports_before_expression_fit"]
    )

    response_axis = protocol["response_axis"]
    assert response_axis["required_unique_source_symbols"] == 374
    assert response_axis["estimated_vcc_visible_coordinates"] == 349

    estimator = protocol["estimator"]
    assert estimator["primary_condition"] == "untreated"
    assert estimator["activated_condition_use"].startswith("diagnostic_only")
    assert estimator["design_unit"] == "target_not_guide"
    assert "channel_one_hot" not in estimator["nuisance_covariates"]
    assert estimator["channel_effect_elimination"].startswith(
        "before_any_fold_fit_subtract_each_channels_own_mean"
    )
    assert estimator["nuisance_residualization"].startswith(
        "fit_training_only_float64_least_squares_to_both_response"
    )
    assert estimator["target_feature_scaling"].startswith(
        "training_only_unit_standard_deviation"
    )
    assert (
        "training_channels_only"
        in estimator["fold_local_estimability"]["support_thresholds"]
    )
    assert "P_ident" in estimator["fold_local_estimability"]["identifiable_projection"]
    assert estimator["fold_local_estimability"]["target_transferability"].endswith(
        "at_least_0.99"
    )
    assert estimator["effect_gauge"].startswith(
        "export_only_the_identifiable_row_space_component"
    )
    assert (
        "null_of_the_training_residualized_design"
        in estimator["nullspace_invariance_test"]
    )
    assert (
        "hard_error"
        in assignment["guide_assignment"]["zero_total_or_zero_probability_policy"]
    )
    assert "guide" not in estimator["selection_metric"]

    held = protocol["independent_source_validation"]
    assert "frozen_to_a14" in held["operator_rule"]
    assert held["source_fit_rule"].startswith(
        "fit_the_GSE249595_source_model_once_without_Nadig"
    )

    integration = protocol["integration"]
    assert integration["stage"] == "post_whitening_masked_residual"
    assert integration["ordinary_fourth_donor_insertion"] == "forbidden"
    assert integration["candidate_formula_unsupported"] == "incumbent_bit_exact"
    assert (
        "supported_by_supported_principal_submatrix" in integration["source_whitening"]
    )
    assert integration["context_a_only_route"].startswith("forbidden")
    assert 0.0 in integration["blend_grid"]

    gates = protocol["promotion_gates"]
    assert gates["nested_blend_macro_pds_delta_minimum"] >= 0.003
    assert gates["nested_blend_pds_delta_each_gallery_minimum"] >= 0.0
    assert gates["target_label_permutations"] >= 99
    assert gates["exact_fallback_required"] is True

    forbidden = set(protocol["forbidden_before_all_source_gates_pass"])
    assert {
        "challenge_build",
        "challenge_submission",
        "count_level_scoring",
    } <= forbidden
    print("experiment 107 preregistration contract passes")


if __name__ == "__main__":
    main()
