from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_run():
    spec = importlib.util.spec_from_file_location("experiment_114_run", HERE / "run.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


RUN = load_run()


class Frame:
    def __init__(self, columns: list[str], values: np.ndarray):
        self.columns = columns
        self._values = values

    def to_numpy(self, dtype):
        return np.asarray(self._values, dtype=dtype)


def test_checked_in_protocol_and_manifest_pass() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    source = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    RUN.validate_protocol(protocol)
    RUN.validate_source_manifest(source)
    report = RUN.run_manifest()
    assert report["eligible_count"] == 41
    assert report["incumbent_uncovered_count"] == 28
    assert report["low_confidence_count"] == 13
    assert report["eligible_ordered_sha256"] == (
        "508daf6356aeb62fd30e9048482bdaca9b1b001ec80e8aef9dabe5b32bc0dcb0"
    )
    assert report["depmap_large_download_authorized"]
    assert report["source_feature_build_authorized"]
    assert not report["source_candidate_freeze_authorized"]
    assert not report["held_truth_access_authorized"]
    assert report["approved_role_count"] == 0
    assert not report["freeze_boundary_authorized"]
    assert report["depmap_feature_shape"] == [17916, 1178]
    assert report["string_feature_shape"] == [21688, 64]


def test_schema_amendment_binds_original_acquisition_and_new_shape() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    amendment = RUN.load_json(RUN.SOURCE_SCHEMA_AMENDMENT_PATH)
    RUN.validate_source_schema_amendment(amendment)
    assert amendment["acquisition_protocol_sha256"] == (RUN.ACQUISITION_PROTOCOL_SHA256)
    assert (
        protocol["representations"]["depmap_24q4"]["rows_cell_lines"],
        protocol["representations"]["depmap_24q4"]["columns_genes"],
    ) == RUN.EXPECTED_DEPMAP_SHAPE
    assert protocol["source_schema_amendment"]["sha256"] == (
        RUN.EXPECTED_SOURCE_SCHEMA_AMENDMENT_SHA256
    )


def test_protocol_is_frozen_to_five_named_arrays_and_four_experimental_arms() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    assert tuple(protocol["arms"]) == RUN.ARMS
    assert protocol["primary_arm"] == "string_depmap_blend_w1"
    assert protocol["secondary_arms"] == ["string_neighbor", "depmap_neighbor"]
    assert protocol["controls"] == ["identity", "current_fallback"]


def test_depmap_transform_matches_center_nan_zero_and_l2_contract() -> None:
    frame = Frame(
        ["A (1)", "B (2)", "C (3)", "D (4)"],
        np.asarray(
            [
                [1.0, 5.0, np.nan, 2.0],
                [2.0, 5.0, np.nan, np.nan],
                [3.0, 5.0, 7.0, np.nan],
            ],
            dtype=np.float32,
        ),
    )
    features = RUN.depmap_features_from_frame(frame)
    assert features.genes.tolist() == ["A", "B"]
    expected = np.asarray([-1.0, 0.0, 1.0]) / np.sqrt(2.0)
    np.testing.assert_allclose(features.matrix[0], expected, atol=1e-7)
    np.testing.assert_array_equal(features.matrix[1], np.zeros(3))


def test_depmap_transform_is_bit_exact_to_pinned_upstream_float32_path() -> None:
    values = np.asarray(
        [
            [1.0, 2.0, np.nan, np.inf],
            [3.0, 4.0, 8.0, -np.inf],
            [5.0, 6.0, 9.0, 7.0],
        ],
        dtype=np.float32,
    )
    frame = Frame(["A (1)", "B (2)", "C (3)", "D (4)"], values)
    observed = RUN.depmap_features_from_frame(frame)

    matrix = values.T
    seen = np.isfinite(matrix).sum(axis=1)
    keep = seen >= 2
    expected = np.where(np.isfinite(matrix[keep]), matrix[keep], np.nan)
    expected = np.nan_to_num(
        expected - np.nanmean(expected, axis=1, keepdims=True),
        nan=0.0,
    )
    norms = np.linalg.norm(expected, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    expected = expected / norms

    assert observed.genes.tolist() == ["A", "B", "C"]
    np.testing.assert_array_equal(observed.matrix, expected)


def test_json_duplicate_keys_fail_closed() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "duplicate.json"
        path.write_text('{"a": 1, "a": 2}\\n')
        try:
            RUN.load_json(path)
        except RUN.ContractError as error:
            assert "duplicate JSON key" in str(error)
        else:
            raise AssertionError("duplicate JSON key was accepted")


def test_depmap_duplicate_symbol_prefix_fails_closed() -> None:
    frame = Frame(
        ["A (1)", "A (2)"],
        np.asarray([[1.0, 2.0], [2.0, 3.0]], dtype=np.float32),
    )
    try:
        RUN.depmap_features_from_frame(frame)
    except RUN.ContractError as error:
        assert "not unique" in str(error)
    else:
        raise AssertionError("duplicate DepMap symbols were accepted")


def test_weight_one_blend_is_equal_cosine_share_when_both_present() -> None:
    left = RUN.GeneFeatures(
        np.asarray(["A", "B"]),
        RUN.row_normalize(np.asarray([[1.0, 0.0], [0.0, 1.0]])),
    )
    right = RUN.GeneFeatures(
        np.asarray(["A", "B"]),
        RUN.row_normalize(np.asarray([[1.0, 1.0], [1.0, -1.0]])),
    )
    blend = RUN.blend_features(left, right, 1.0)
    similarity = float(blend.matrix[0] @ blend.matrix[1])
    expected = 0.5 * float(left.matrix[0] @ left.matrix[1]) + 0.5 * float(
        right.matrix[0] @ right.matrix[1]
    )
    np.testing.assert_allclose(similarity, expected, atol=1e-7)


def test_neighbour_weights_select_exactly_25_under_cutoff_ties() -> None:
    similarity = np.full(30, 0.5)
    weights = RUN.neighbour_weights(similarity)
    assert np.count_nonzero(weights) == 25
    np.testing.assert_array_equal(np.flatnonzero(weights), np.arange(25))


def synthetic_bank() -> RUN.SourceBank:
    return RUN.SourceBank(
        targets=np.asarray(["Q", "A", "A", "B", "C"]),
        genes=np.asarray(["PANEL", "g1", "g2", "g3"]),
        effects=np.asarray(
            [
                [0.0, 50.0, 50.0, 50.0],
                [0.0, 2.0, 4.0, 0.0],
                [0.0, 4.0, 0.0, 8.0],
                [0.0, -2.0, 2.0, 1.0],
                [0.0, 1.0, -3.0, 2.0],
            ]
        ),
        covered=np.asarray(
            [
                [True, True, True, True],
                [True, True, True, False],
                [True, True, False, True],
                [True, True, True, True],
                [True, True, True, True],
            ]
        ),
        source_ids=np.asarray(["s1", "s1", "s2", "s1", "s1"]),
    )


def test_consensus_is_equal_donor_and_coordinate_coverage_aware() -> None:
    consensus = RUN.build_consensus(synthetic_bank())
    row = int(np.flatnonzero(consensus.targets == "A")[0])
    np.testing.assert_allclose(consensus.effects[row], [0.0, 3.0, 4.0, 8.0])
    assert consensus.carrier_count[row] == 2
    assert consensus.covered[row].all()


def test_build_arm_is_leave_target_out_masked_and_exact_outside_scope() -> None:
    consensus = RUN.build_consensus(synthetic_bank())
    targets = np.asarray(["Q", "INELIGIBLE", "MISSING"])
    genes = consensus.genes
    incumbent = np.asarray(
        [
            [0.0, 3.0, 4.0, 0.0],
            [0.0, 7.0, 0.0, 0.0],
            [0.0, 0.0, 5.0, 0.0],
        ],
        dtype=np.float64,
    )
    features = RUN.GeneFeatures(
        np.asarray(["Q", "A", "B", "C"]),
        RUN.row_normalize(
            np.asarray(
                [
                    [1.0, 0.0],
                    [1.0, 0.0],
                    [0.5, 0.5],
                    [0.0, 1.0],
                ]
            )
        ),
    )
    effect, _, audit = RUN.build_arm(
        incumbent_targets=targets,
        incumbent_genes=genes,
        incumbent_effects=incumbent,
        consensus=consensus,
        features=features,
        eligible_targets={"Q", "MISSING"},
        panel_targets={"PANEL"},
    )
    assert not np.array_equal(effect[0], incumbent[0])
    assert np.array_equal(effect[1], incumbent[1])
    assert np.array_equal(effect[2], incumbent[2])
    assert effect[:, 0].tolist() == [0.0, 0.0, 0.0]
    assert audit["neighbor_audit"]["Q"]["same_target_source_rows"] == 0
    assert "MISSING" in audit["unsupported_targets"]
    np.testing.assert_allclose(
        np.linalg.norm(effect[0, 1:]),
        np.linalg.norm(incumbent[0, 1:]),
    )


def test_candidate_confidence_is_exact_abs_of_frozen_effect_dtype() -> None:
    consensus = RUN.build_consensus(synthetic_bank())
    features = RUN.GeneFeatures(
        np.asarray(["Q", "A", "B", "C"]),
        RUN.row_normalize(
            np.asarray([[1.0, 0.0], [0.91, 0.09], [0.5, 0.5], [0.0, 1.0]])
        ),
    )
    incumbent = np.asarray([[0.0, 0.31415927, 0.27182818, 0.0]], dtype=np.float32)
    effect, confidence, _ = RUN.build_arm(
        incumbent_targets=["Q"],
        incumbent_genes=consensus.genes,
        incumbent_effects=incumbent,
        consensus=consensus,
        features=features,
        eligible_targets={"Q"},
        panel_targets={"PANEL"},
    )
    np.testing.assert_array_equal(
        confidence,
        np.abs(effect.astype(np.float64)),
    )


def test_same_target_signature_cannot_influence_candidate() -> None:
    bank = synthetic_bank()
    consensus = RUN.build_consensus(bank)
    features = RUN.GeneFeatures(
        np.asarray(["Q", "A", "B", "C"]),
        RUN.row_normalize(np.asarray([[1, 0], [0.9, 0.1], [0.5, 0.5], [0, 1]])),
    )
    incumbent = np.asarray([[0.0, 1.0, 0.0, 0.0]])
    first, _, _ = RUN.build_arm(
        incumbent_targets=["Q"],
        incumbent_genes=consensus.genes,
        incumbent_effects=incumbent,
        consensus=consensus,
        features=features,
        eligible_targets={"Q"},
        panel_targets={"PANEL"},
    )
    changed_bank = RUN.SourceBank(
        bank.targets,
        bank.genes,
        bank.effects.copy(),
        bank.covered,
        bank.source_ids,
    )
    changed_bank.effects[0, 1:] = 1e9
    second, _, _ = RUN.build_arm(
        incumbent_targets=["Q"],
        incumbent_genes=consensus.genes,
        incumbent_effects=incumbent,
        consensus=RUN.build_consensus(changed_bank),
        features=features,
        eligible_targets={"Q"},
        panel_targets={"PANEL"},
    )
    np.testing.assert_array_equal(first, second)


def test_permutation_null_is_deterministic_and_stays_within_strata() -> None:
    fallback = ["A", "B", "C", "D"]
    low = ["E", "F", "G"]
    first = RUN.permuted_query_labels(fallback, low, 7)
    second = RUN.permuted_query_labels(fallback, low, 7)
    assert first == second
    assert {first[target] for target in fallback} == set(fallback)
    assert {first[target] for target in low} == set(low)
    assert RUN.permuted_query_labels(fallback, low, 8) != first


def test_pds_uses_raw_uncentered_complete_gallery() -> None:
    truth = np.asarray(
        [
            [-2.25014117, 0.38636960, -0.58164084],
            [0.10927970, -0.07570153, 0.20211440],
            [0.69417194, -0.75836975, 1.42098202],
            [0.72609379, 0.84373266, 1.16486398],
        ]
    )
    prediction = truth.copy()
    score, per_target = RUN.raw_uncentered_pds(prediction, truth)
    assert score == 1.0
    np.testing.assert_array_equal(per_target, np.ones(4))
    centered_prediction = prediction - prediction.mean(axis=0, keepdims=True)
    centered_score, _ = RUN.raw_uncentered_pds(centered_prediction, truth)
    assert centered_score < score


def test_directional_precision_selection_is_prediction_only() -> None:
    prediction = np.asarray([[100.0, 2.0, -1.0, 0.5, 0.1]])
    confidence = np.asarray([[100.0, 4.0, 3.0, 2.0, 1.0]])
    truth = np.asarray([[-1, 1, -1, -1, 1]], dtype=np.int8)
    score, per_target = RUN.source_blind_directional_precision(
        prediction,
        confidence,
        truth,
        fraction=0.5,
        allowed_mask=np.asarray([False, True, True, True, True]),
    )
    assert score == 1.0
    assert per_target.tolist() == [1.0]


def result(pds: float, precision: float) -> dict[str, object]:
    return {"raw_pds": pds, "directional_precision": precision}


def test_promotion_requires_all_surfaces_macro_precision_and_null() -> None:
    surfaces = {
        "clean_a": {
            "current_fallback": result(0.60, 0.55),
            RUN.PRIMARY_ARM: result(0.604, 0.56),
        },
        "clean_b": {
            "current_fallback": result(0.70, 0.60),
            RUN.PRIMARY_ARM: result(0.704, 0.61),
        },
    }
    decision = RUN.promotion_decision(surfaces, [0.001] * RUN.NULL_DRAWS)
    assert decision["promote"] is True
    negative = json.loads(json.dumps(surfaces))
    negative["clean_b"][RUN.PRIMARY_ARM]["raw_pds"] = 0.699
    assert (
        RUN.promotion_decision(negative, [0.001] * RUN.NULL_DRAWS)["promote"] is False
    )
    assert (
        RUN.promotion_decision(surfaces, [0.005] * RUN.NULL_DRAWS)["promote"] is False
    )


def test_mutating_outputs_are_locked() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        fake = Path(directory) / "fake.csv"
        fake.write_text("x\n")
        commands = [
            (
                [
                    "freeze",
                    "--input-manifest",
                    str(fake),
                    "--out",
                    str(Path(directory) / "frozen"),
                ],
                "source_candidate_freeze=false",
            ),
            (
                [
                    "score",
                    "--frozen",
                    str(fake),
                    "--truth-manifest",
                    str(fake),
                    "--out",
                    str(Path(directory) / "score.json"),
                ],
                "held_truth_access=false",
            ),
        ]
        for arguments, message in commands:
            completed = subprocess.run(
                [sys.executable, str(HERE / "run.py"), *arguments],
                text=True,
                capture_output=True,
            )
            assert completed.returncode != 0
            assert message in completed.stderr


def test_reviewed_boundary_independently_blocks_freeze() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    protocol["current_authorization"]["source_candidate_freeze"] = True
    try:
        RUN.freeze_from_manifest(
            HERE / "not-approved.json",
            HERE / "work" / "not-approved",
            protocol=protocol,
        )
    except RUN.ContractError as error:
        assert "reviewed boundary" in str(error)
    else:
        raise AssertionError("freeze bypassed the independent boundary lock")


def test_download_requires_explicit_confirmation() -> None:
    completed = subprocess.run(
        [sys.executable, str(HERE / "run.py"), "download-depmap"],
        text=True,
        capture_output=True,
    )
    assert completed.returncode != 0
    assert "large download requires --confirm-large-download" in completed.stderr


def test_source_feature_build_rejects_noncanonical_paths() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        fake = Path(directory) / "fake.csv"
        fake.write_text("x\n")
        protocol = RUN.load_json(RUN.PROTOCOL_PATH)
        source = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
        try:
            RUN.build_depmap_features(
                fake,
                RUN.DEPMAP_FEATURES_PATH,
                protocol=protocol,
                source_manifest=source,
            )
        except RUN.ContractError as error:
            assert "frozen experiment source path" in str(error)
        else:
            raise AssertionError("noncanonical source path was accepted")


def test_exact_eligibility_rejects_caller_supplied_row_selection() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    source = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    targets = np.asarray(RUN.read_csv_column(RUN.PANEL_PATH, "target_gene"))
    donor_sets = {
        name: set(RUN.read_lines(RUN.ROOT / f"data/donors/{name}/perts.txt"))
        for name in protocol["incumbent_contract"]["real_donors"]
    }
    perturbmatch = set(
        RUN.read_lines(RUN.ROOT / "data/external/perturbmatch/perts.txt")
    )
    carriers = np.asarray(
        [sum(target in values for values in donor_sets.values()) for target in targets]
    )
    fallback = np.asarray(
        [
            target not in donor_sets["k562_gwps"] and target not in perturbmatch
            for target in targets
        ]
    )
    incumbent = {
        "targets": targets,
        "incumbent_fallback": fallback.copy(),
        "source_carrier_count": carriers,
        "expected_source_count": np.full(len(targets), 3, dtype=np.int16),
    }
    assert (
        RUN.validate_exact_eligibility(incumbent, protocol, source)["eligible_targets"]
        == protocol["eligibility"]["eligible_targets"]
    )
    incumbent["incumbent_fallback"][0] = not incumbent["incumbent_fallback"][0]
    try:
        RUN.validate_exact_eligibility(incumbent, protocol, source)
    except RUN.ContractError as error:
        assert "fallback flags" in str(error)
    else:
        raise AssertionError("caller-supplied eligibility was accepted")


def test_source_policy_rejects_held_surface_source() -> None:
    bank = RUN.SourceBank(
        targets=np.asarray(["A", "B"]),
        genes=np.asarray(["g1"]),
        effects=np.asarray([[1.0], [2.0]]),
        covered=np.ones((2, 1), dtype=bool),
        source_ids=np.asarray(["source", "held"]),
    )
    approved = {
        "allowed_source_ids": ["source", "held"],
        "forbidden_source_ids": ["held"],
        "held_surface_source_id": "held",
    }
    receipt = {
        "source_ids": ["source", "held"],
        "forbidden_source_ids": ["held"],
        "held_surface_source_id": "held",
    }
    try:
        RUN.validate_source_id_policy(bank, receipt, approved)
    except RUN.ContractError as error:
        assert "forbidden held-surface" in str(error)
    else:
        raise AssertionError("held-surface source was accepted")


def test_role_manifest_must_be_hash_pinned_by_boundary() -> None:
    try:
        RUN.approved_role_for_manifest(HERE / "test_run.py", {"approved_roles": {}})
    except RUN.ContractError as error:
        assert "not hash-pinned" in str(error)
    else:
        raise AssertionError("unapproved role manifest was accepted")


def test_synthetic_freeze_writes_only_sealed_source_artifacts(monkeypatch) -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        monkeypatch.setattr(RUN, "HERE", root)
        role = "synthetic"
        role_dir = root / "work" / "roles" / role
        role_dir.mkdir(parents=True)
        targets = np.asarray(["Q", "L", "SAFE"])
        genes = np.asarray(["Q", "L", "SAFE", "g1", "g2"])
        incumbent_path = role_dir / "incumbent.npz"
        np.savez_compressed(
            incumbent_path,
            targets=targets,
            genes=genes,
            effects=np.asarray(
                [
                    [0.0, 0.0, 0.0, 2.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0, 1.0],
                    [0.0, 0.0, 0.0, 0.0, 2.0],
                ],
                dtype=np.float32,
            ),
            incumbent_fallback=np.asarray([True, False, False]),
            source_carrier_count=np.asarray([1, 2, 3]),
            expected_source_count=np.asarray(3),
        )
        source_path = role_dir / "source-bank.npz"
        np.savez_compressed(
            source_path,
            targets=np.asarray(["Q", "A", "B"]),
            genes=genes,
            effects=np.asarray(
                [
                    [0.0, 0.0, 0.0, 100.0, 100.0],
                    [0.0, 0.0, 0.0, 1.0, 2.0],
                    [0.0, 0.0, 0.0, -1.0, 2.0],
                ]
            ),
            covered=np.asarray(
                [
                    [False, False, False, True, True],
                    [False, False, False, True, True],
                    [False, False, False, True, True],
                ]
            ),
            source_ids=np.asarray(["s1", "s1", "s2"]),
        )
        depmap_features = RUN.GeneFeatures(
            np.asarray(["Q", "L", "A", "B"]),
            RUN.row_normalize(np.asarray([[1, 0], [0, 1], [1, 0], [0, 1]])),
        )
        string_features = RUN.GeneFeatures(
            np.asarray(["Q", "L", "A", "B"]),
            RUN.row_normalize(np.asarray([[1, 0], [0, 1], [0.8, 0.2], [0.2, 0.8]])),
        )

        def record(path: Path) -> dict[str, str]:
            return {
                "path": path.resolve().relative_to(RUN.ROOT.resolve()).as_posix(),
                "sha256": RUN.sha256_file(path),
            }

        protocol = RUN.load_json(RUN.PROTOCOL_PATH)
        protocol["current_authorization"]["source_candidate_freeze"] = True
        fallback = ["Q"]
        low = ["L"]
        eligible = ["Q", "L"]
        protocol["eligibility"] = {
            "incumbent_uncovered_targets": fallback,
            "incumbent_uncovered_ordered_sha256": RUN.sha256_lines(fallback),
            "low_confidence_targets": low,
            "low_confidence_ordered_sha256": RUN.sha256_lines(low),
            "eligible_targets": eligible,
            "eligible_ordered_sha256": RUN.sha256_lines(eligible),
        }
        protocol_path = root / "protocol.json"
        protocol_path.write_text(json.dumps(protocol))
        monkeypatch.setattr(RUN, "PROTOCOL_PATH", protocol_path)

        builder = record(HERE / "test_run.py")
        role_source_manifest = record(HERE / "ROLE_SOURCE_MANIFEST.json")
        source_inputs = [
            {
                **builder,
                "source_only": True,
                "contains_held_truth": False,
            }
        ]
        incumbent = RUN.load_incumbent(incumbent_path)
        source_bank = RUN.load_source_bank(source_path)
        incumbent_receipt_path = role_dir / "incumbent-receipt.json"
        incumbent_receipt_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "experiment": "114-depmap-coessentiality",
                    "stage": "role_incumbent_source_only",
                    "role": role,
                    "surface": "synthetic",
                    "source_only": True,
                    "contains_held_truth": False,
                    "sealed_before_held_truth": True,
                    "artifact_path": record(incumbent_path)["path"],
                    "artifact_sha256": record(incumbent_path)["sha256"],
                    "builder": builder,
                    "role_source_manifest": role_source_manifest,
                    "source_inputs": source_inputs,
                    "endpoint": "a14-shaped",
                    "array_sha256": {
                        name: RUN.sha256_array(incumbent[name])
                        for name in (
                            "targets",
                            "genes",
                            "effects",
                            "incumbent_fallback",
                            "source_carrier_count",
                            "expected_source_count",
                        )
                    },
                }
            )
        )
        source_receipt_path = role_dir / "source-bank-receipt.json"
        source_receipt_path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "experiment": "114-depmap-coessentiality",
                    "stage": "role_source_bank_source_only",
                    "role": role,
                    "surface": "synthetic",
                    "source_only": True,
                    "contains_held_truth": False,
                    "sealed_before_held_truth": True,
                    "artifact_path": record(source_path)["path"],
                    "artifact_sha256": record(source_path)["sha256"],
                    "builder": builder,
                    "role_source_manifest": role_source_manifest,
                    "source_inputs": source_inputs,
                    "effect_space": "synthetic-final",
                    "source_ids": ["s1", "s2"],
                    "forbidden_source_ids": ["held"],
                    "held_surface_source_id": "held",
                    "array_sha256": {
                        "targets": RUN.sha256_array(source_bank.targets),
                        "genes": RUN.sha256_array(source_bank.genes),
                        "effects": RUN.sha256_array(source_bank.effects),
                        "covered": RUN.sha256_array(source_bank.covered),
                        "source_ids": RUN.sha256_array(source_bank.source_ids),
                    },
                }
            )
        )
        input_manifest = role_dir / "manifest.json"
        input_manifest.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "experiment": "114-depmap-coessentiality",
                    "stage": "source_only_role_inputs",
                    "role": role,
                    "surface": "synthetic",
                    "role_directory": role_dir.resolve()
                    .relative_to(RUN.ROOT.resolve())
                    .as_posix(),
                    "contains_held_truth": False,
                    "source_only": True,
                    "sealed_before_held_truth": True,
                    "role_source_manifest": role_source_manifest,
                    "incumbent": {
                        "artifact": record(incumbent_path),
                        "receipt": record(incumbent_receipt_path),
                    },
                    "source_bank": {
                        "artifact": record(source_path),
                        "receipt": record(source_receipt_path),
                    },
                }
            )
        )
        boundary_path = root / "boundary.json"
        boundary_path.write_text("{}")
        monkeypatch.setattr(RUN, "SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH", boundary_path)
        output = role_dir / "frozen"
        approved = {
            "manifest_path": record(input_manifest)["path"],
            "manifest_sha256": record(input_manifest)["sha256"],
            "surface": "synthetic",
            "role_directory": role_dir.resolve()
            .relative_to(RUN.ROOT.resolve())
            .as_posix(),
            "output_directory": output.resolve()
            .relative_to(RUN.ROOT.resolve())
            .as_posix(),
            "incumbent_endpoint": "a14-shaped",
            "source_effect_space": "synthetic-final",
            "allowed_source_ids": ["s1", "s2"],
            "forbidden_source_ids": ["held"],
            "held_surface_source_id": "held",
        }
        boundary = {
            "source_candidate_freeze_authorized": True,
            "approved_roles": {role: approved},
            "features": {
                "depmap": {
                    "artifact_sha256": "synthetic-depmap",
                    "receipt_sha256": "synthetic-depmap-receipt",
                },
                "string": {"artifact_sha256": "synthetic-string"},
            },
        }
        monkeypatch.setattr(RUN, "load_freeze_boundary", lambda: boundary)
        monkeypatch.setattr(
            RUN,
            "validate_canonical_axes",
            lambda boundary, source: (targets, genes),
        )
        monkeypatch.setattr(
            RUN,
            "validate_depmap_feature_artifact",
            lambda boundary, protocol, source: depmap_features,
        )
        monkeypatch.setattr(
            RUN,
            "validate_string_feature_artifact",
            lambda boundary, source: string_features,
        )
        synthetic_eligibility = {
            "incumbent_uncovered_targets": fallback,
            "low_confidence_targets": low,
            "eligible_targets": eligible,
            "incumbent_uncovered_ordered_sha256": RUN.sha256_lines(fallback),
            "low_confidence_ordered_sha256": RUN.sha256_lines(low),
            "eligible_ordered_sha256": RUN.sha256_lines(eligible),
        }
        monkeypatch.setattr(
            RUN,
            "validate_exact_eligibility",
            lambda incumbent, protocol, source: synthetic_eligibility,
        )
        manifest = RUN.freeze_from_manifest(
            input_manifest,
            output,
            protocol=protocol,
        )
        assert manifest["contains_held_truth"] is False
        assert manifest["sealed_before_held_truth"] is True
        assert manifest["identity_current_fallback_bit_identical"] is True
        assert manifest["permutation_draws"] == 100
        assert manifest["eligibility"]["eligible_targets"] == ["Q", "L"]
        assert set(path.name for path in output.iterdir()) == {
            "arms.npz",
            "manifest.json",
            "permutation_nulls.npz",
        }
        loaded_manifest, arrays = RUN.load_frozen_directory(output)
        assert loaded_manifest["surface"] == "synthetic"
        np.testing.assert_array_equal(
            arrays["effect__identity"],
            arrays["effect__current_fallback"],
        )
        np.testing.assert_array_equal(
            arrays["confidence__string_depmap_blend_w1"],
            np.abs(arrays["effect__string_depmap_blend_w1"].astype(np.float64)),
        )


def test_output_confinement_fails_closed() -> None:
    try:
        RUN.require_under_experiment(Path("/tmp/experiment-114-escape"))
    except RUN.ContractError as error:
        assert "escapes experiment directory" in str(error)
    else:
        raise AssertionError("outside output path was accepted")
