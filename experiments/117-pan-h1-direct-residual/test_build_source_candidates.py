"""Synthetic tests for Experiment 117's source-only candidate freeze."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "experiment_117_build_source_candidates",
    HERE / "build_source_candidates.py",
)
assert SPEC is not None and SPEC.loader is not None
BUILD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BUILD
SPEC.loader.exec_module(BUILD)


def synthetic_role() -> object:
    genes = np.asarray(["PANEL", "A", "HELD", "B", "NO_OPERATOR", "C"])
    targets = np.asarray(["T_NEW", "T_OLD", "T_MISSING", "T_ZERO"])
    incumbent = np.asarray(
        [
            [0.0, 1.0, 0.0, 2.0, 3.0, 4.0],
            [0.0, 2.0, 0.0, 1.0, 4.0, 3.0],
            [0.0, 3.0, 0.0, 4.0, 5.0, 2.0],
            [0.0, 4.0, 0.0, 3.0, 2.0, 1.0],
        ],
        dtype=np.float64,
    )
    columns = np.asarray([1, 3, 5], dtype=np.int64)
    return BUILD.RoleBundle(
        name="synthetic",
        targets=targets,
        genes=genes,
        score_indices=np.asarray([1, 3, 4, 5], dtype=np.int64),
        incumbent=incumbent,
        eigenvalues=np.ones(3, dtype=np.float64),
        eigenvectors=np.eye(3, dtype=np.float64),
        operator_columns=columns,
        input_bindings=(),
    )


def synthetic_compact(directory: Path) -> object:
    genes = ("PANEL", "A", "HELD", "B", "NO_OPERATOR", "C")
    targets = ("T_OLD", "T_NEW", "T_ZERO")
    delta = np.asarray(
        [
            [99.0, 2.0, 88.0, 1.0, 77.0, -2.0],
            [99.0, 1.0, 88.0, -1.0, 77.0, 2.0],
            [99.0, 0.0, 88.0, 0.0, 77.0, 0.0],
        ],
        dtype=np.float32,
    )
    measured = np.ones_like(delta, dtype=bool)
    measured[1, 3] = False
    return BUILD.CompactDonor(
        directory=directory,
        genes=genes,
        targets=targets,
        delta=delta,
        measured=measured,
        receipt_sha256="a" * 64,
        artifact_manifest_sha256="b" * 64,
    )


def test_protocol_requires_two_non_h1_cross_context_gates() -> None:
    protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
    roles = protocol["validation_design"]["roles"]
    assert roles["h1"]["promotion_authorizing"] is False
    assert roles["orion_hek293t"]["promotion_authorizing"] is True
    assert roles["nadig_jurkat"]["promotion_authorizing"] is True
    assert roles["nadig_jurkat"]["cannot_authorize_without_hek293t"] is True
    gate = protocol["promotion_gate"]
    assert gate["h1_same_context_result_cannot_replace_any_cross_context_gate"] is True
    assert gate["nadig_alone_cannot_authorize_promotion"] is True


def test_protocol_freezes_repaired_evaluation_mechanics() -> None:
    protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
    design = protocol["validation_design"]
    directional = design["directional_credit"]
    assert "scaled_pan_minus_destination_incumbent" in directional["ranking"]
    assert directional["eligible_rows"].startswith("target_effective_novel")
    aggregation = design["cross_context_aggregation"]
    assert aggregation["surface_order"] == [
        "orion_hek293t",
        *BUILD.NADIG_ROLES,
    ]
    assert "divided_by_2" in aggregation["cross_context_macro"]
    null = design["target_label_permutation_null"]
    assert null["bit_generator"] == "numpy_random_PCG64"
    assert "never_permute_or_reuse_the_observed_normalized_residual" in null[
        "null_reconstruction"
    ]
    bootstrap = design["target_cluster_bootstrap"]
    assert bootstrap["bit_generator"] == "numpy_random_PCG64"
    assert "method_linear" in bootstrap["confidence_interval"]


def test_scaffold_authorizes_only_source_candidate_building() -> None:
    protocol, source, boundary = BUILD.validate_scaffold()
    assert source["pan_source"]["real_source_authorized"] is True
    assert boundary["source_candidate_build_authorized"] is True
    assert set(boundary["approved_roles"]) == set(source["roles"])
    assert all(boundary["approved_roles"].values())
    for key in (
        "held_truth_access_authorized",
        "evaluation_authorized",
        "cell_eval2_authorized",
        "count_generation_authorized",
        "challenge_build_authorized",
        "submission_authorized",
    ):
        assert boundary[key] is False
    for key in (
        "held_truth_access",
        "evaluation",
        "cell_eval2",
        "count_generation",
        "challenge_build",
        "submission",
    ):
        assert protocol["current_authorization"][key] is False
        assert source["authorization"].get(key, False) is False
    BUILD.assert_real_build_authorized(source, boundary)


def test_reviewed_compact_donor_authenticates_without_role_outcomes() -> None:
    source = json.loads((HERE / "SOURCE_INPUT_MANIFEST.json").read_text())
    compact = BUILD.authenticate_compact_donor(source)
    assert compact.receipt_sha256 == source["pan_source"]["receipt_sha256"]
    assert (
        compact.artifact_manifest_sha256
        == source["pan_source"]["artifact_manifest_sha256"]
    )


@pytest.mark.parametrize(
    ("drift", "message"),
    [
        ("receipt", "receipt hash drifted"),
        ("artifact_manifest", "artifact manifest hash drifted"),
        ("delta", "delta.npy pinned hash drifted"),
        ("measured_mask", "measured_mask.npy pinned hash drifted"),
    ],
)
def test_compact_pin_validator_rejects_every_reviewed_hash_drift(
    drift: str,
    message: str,
) -> None:
    source = json.loads((HERE / "SOURCE_INPUT_MANIFEST.json").read_text())
    pan = copy.deepcopy(source["pan_source"])
    artifact_hashes = {
        "delta.npy": pan["delta_sha256"],
        "measured_mask.npy": pan["measured_mask_sha256"],
    }
    receipt = pan["receipt_sha256"]
    manifest = pan["artifact_manifest_sha256"]
    if drift == "receipt":
        receipt = "0" * 64
    elif drift == "artifact_manifest":
        manifest = "0" * 64
    elif drift == "delta":
        artifact_hashes["delta.npy"] = "0" * 64
    else:
        artifact_hashes["measured_mask.npy"] = "0" * 64
    with pytest.raises(BUILD.ContractError, match=message):
        BUILD.validate_compact_hash_pins(
            pan,
            receipt_sha256=receipt,
            artifact_manifest_sha256=manifest,
            artifact_hashes=artifact_hashes,
        )


@pytest.mark.parametrize(
    ("key", "message"),
    [
        ("receipt_sha256", "receipt hash drifted"),
        ("artifact_manifest_sha256", "artifact manifest hash drifted"),
        ("delta_sha256", "delta.npy pinned hash drifted"),
        ("measured_mask_sha256", "measured_mask.npy pinned hash drifted"),
    ],
)
def test_source_build_aborts_before_output_on_any_compact_pin_drift(
    key: str,
    message: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protocol, source, boundary = BUILD.validate_scaffold()
    drifted = copy.deepcopy(source)
    drifted["pan_source"][key] = "0" * 64
    monkeypatch.setattr(
        BUILD,
        "validate_scaffold",
        lambda: (protocol, drifted, boundary),
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        output = Path(directory) / "must-not-exist"
        with pytest.raises(BUILD.ContractError, match=message):
            BUILD.build(output)
        assert not output.exists()


def test_build_refuses_to_overwrite_existing_output_before_array_load() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        with pytest.raises(BUILD.ContractError, match="refusing to overwrite"):
            BUILD.build(Path(directory))


def test_exact_panel_coverage_contract_is_64_and_57() -> None:
    source = json.loads((HERE / "SOURCE_INPUT_MANIFEST.json").read_text())
    assert source["coverage"]["vcc_panel_supported_count"] == 64
    assert source["coverage"]["vcc_panel_new_vs_vcc2025_h1_count"] == 57
    assert (
        source["coverage"]["vcc_panel_supported_ordered_sha256"]
        == "f7fee89539444565949d9f7d3ac39e3b5ffa727c165fbc39165d86fa2ee48bc2"
    )
    assert (
        source["coverage"]["vcc_panel_new_vs_vcc2025_h1_ordered_sha256"]
        == "9761fd655ac15dcb16dee00a55b6fe56c977b76acc86246f79b65cb49be87f1d"
    )


def test_row_supported_whitener_matches_explicit_principal_blocks() -> None:
    eigenvalues = np.asarray([1.0, 4.0, 9.0])
    angle = 0.37
    eigenvectors = np.asarray(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    raw = np.asarray([[2.0, 5.0, 7.0], [3.0, 4.0, 8.0]])
    support = np.asarray([[True, False, True], [False, True, True]])
    observed = BUILD.apply_row_supported_whitener(
        raw, support, eigenvalues, eigenvectors
    )
    weights = eigenvalues ** (-BUILD.GAMMA / 2.0)
    full = (eigenvectors * weights) @ eigenvectors.T
    expected = np.zeros_like(raw)
    for row in range(len(raw)):
        positions = np.flatnonzero(support[row])
        expected[row, positions] = (
            raw[row, positions] @ full[np.ix_(positions, positions)]
        )
    np.testing.assert_allclose(observed, expected, rtol=1e-13, atol=1e-13)


def test_candidate_builder_preserves_every_unsupported_cell_bit_exactly() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        role = synthetic_role()
        compact = synthetic_compact(Path(directory))
        arrays, audit = BUILD.build_candidate_arrays(
            role,
            compact,
            panel_targets=["PANEL"],
            old_h1_targets=["T_OLD"],
            enforce_production_coverage=False,
        )
    base = arrays["incumbent_effects"]
    all_candidate = arrays["pan_all_f100"]
    novel_candidate = arrays["pan_novel_f100"]
    support = arrays["coordinate_supported"]
    target_supported = arrays["target_supported"]
    novel = arrays["target_novel_vs_vcc2025_h1"]
    effective_all = arrays["target_effective_all"]
    effective_novel = arrays["target_effective_novel"]

    np.testing.assert_array_equal(target_supported, [True, True, False, True])
    np.testing.assert_array_equal(novel, [True, False, False, True])
    np.testing.assert_array_equal(effective_all, [True, True, False, False])
    np.testing.assert_array_equal(effective_novel, [True, False, False, False])
    assert not np.any(support[:, [0, 2, 4]])
    np.testing.assert_array_equal(all_candidate[:, [0, 2, 4]], base[:, [0, 2, 4]])
    np.testing.assert_array_equal(all_candidate[2], base[2])
    np.testing.assert_array_equal(all_candidate[3], base[3])
    np.testing.assert_array_equal(novel_candidate[1], base[1])
    np.testing.assert_array_equal(novel_candidate[2], base[2])
    np.testing.assert_array_equal(novel_candidate[3], base[3])
    assert not np.array_equal(novel_candidate[0], base[0])
    assert audit["pan_supported_targets"] == 3
    assert audit["pan_novel_targets"] == 2


def test_target_alignment_is_by_symbol_not_source_row_order() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        role = synthetic_role()
        compact = synthetic_compact(Path(directory))
        arrays, _ = BUILD.build_candidate_arrays(
            role,
            compact,
            panel_targets=["PANEL"],
            old_h1_targets=["T_OLD"],
            enforce_production_coverage=False,
        )
    source = arrays["source_effects"]
    assert source[0, 1] == pytest.approx(1.0)
    assert source[1, 1] == pytest.approx(2.0)
    assert source[2].sum() == 0.0


def test_identity_and_zero_fraction_are_bit_exact() -> None:
    role = synthetic_role()
    source = np.ones_like(role.incumbent)
    support = np.ones_like(role.incumbent, dtype=bool)
    rows = np.ones(len(role.targets), dtype=bool)
    candidate, effective = BUILD.blend_row_supported(
        role.incumbent, source, support, rows, 0.0
    )
    np.testing.assert_array_equal(
        candidate.view(np.uint8), role.incumbent.view(np.uint8)
    )
    assert not np.any(effective)


def test_effective_mask_excludes_zero_nonfinite_and_zero_residual_rows() -> None:
    incumbent = np.asarray(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
        ]
    )
    source = np.asarray(
        [
            [0.0, 1.0],
            [0.0, 0.0],
            [np.nan, 1.0],
            [2.0, 0.0],
        ]
    )
    support = np.ones_like(source, dtype=bool)
    residual, effective, effective_support = (
        BUILD.normalized_residual_from_assignment(
            incumbent,
            source,
            support,
            np.ones(4, dtype=bool),
        )
    )
    np.testing.assert_array_equal(effective, [True, False, False, False])
    assert np.all(effective_support[0])
    assert not np.any(effective_support[1:])
    assert np.all(residual[1:] == 0.0)


def test_null_permutations_use_only_effective_rows_and_preserve_novelty() -> None:
    effective = np.asarray([True, True, True, False, False, True])
    novel = np.asarray([True, True, False, False, False, True])
    first = BUILD.deterministic_target_label_permutations(
        effective, novel, role="synthetic", replicates=10
    )
    second = BUILD.deterministic_target_label_permutations(
        effective, novel, role="synthetic", replicates=10
    )
    np.testing.assert_array_equal(first, second)
    assert np.all(first[:, ~effective] == -1)
    for draw in first:
        for stratum in (False, True):
            rows = np.flatnonzero(effective & (novel == stratum))
            np.testing.assert_array_equal(np.sort(draw[rows]), rows)
            if len(rows) > 1:
                assert np.all(draw[rows] != rows)


def test_directional_ranking_uses_actual_normalized_residual_and_effective_rows() -> None:
    source = np.asarray([[0.0, 100.0, 2.0, 0.0], [1.0, 4.0, 0.0, -2.0]])
    residual = np.asarray([[0.0, 1.0, -5.0, 0.0], [9.0, 8.0, 0.0, -7.0]])
    support = np.asarray([[False, True, True, False], [True, False, False, True]])
    effective = np.asarray([True, False])
    ranking = BUILD.directional_rank_indices(residual, support, effective)
    np.testing.assert_array_equal(ranking, [[2, 1], [-1, -1]])
    assert np.argmax(np.abs(source[0])) == 1


def test_null_reconstruction_moves_source_support_and_renormalizes_at_destination() -> None:
    incumbent = np.asarray([[3.0, 4.0, 9.0], [6.0, 8.0, 7.0]])
    source = np.asarray([[0.0, 2.0, 2.0], [5.0, 0.0, 5.0]])
    support = np.asarray([[False, True, True], [True, False, True]])
    novel = np.asarray([True, True])
    mapping = np.asarray([1, 0], dtype=np.int32)
    candidate, residual, effective, effective_support = (
        BUILD.reconstruct_target_label_null(
            incumbent,
            source,
            support,
            novel,
            mapping,
            fraction=0.1,
        )
    )
    np.testing.assert_array_equal(effective, [True, True])
    np.testing.assert_array_equal(
        effective_support,
        [[True, False, True], [False, True, True]],
    )
    assert residual[0, 1] == 0.0
    assert residual[1, 0] == 0.0
    assert candidate[0, 1] == incumbent[0, 1]
    assert candidate[1, 0] == incumbent[1, 0]


def test_null_reconstruction_uses_destination_norm_not_observed_residual() -> None:
    incumbent = np.asarray([[3.0, 4.0], [6.0, 8.0]])
    source = np.asarray([[1.0, 2.0], [5.0, 1.0]])
    support = np.ones_like(source, dtype=bool)
    novel = np.asarray([True, True])
    mapping = np.asarray([1, 0], dtype=np.int32)
    candidate, residual, effective, _ = BUILD.reconstruct_target_label_null(
        incumbent,
        source,
        support,
        novel,
        mapping,
        fraction=0.1,
    )
    np.testing.assert_array_equal(effective, [True, True])
    expected_scaled = np.asarray(
        [
            source[1] * (np.linalg.norm(incumbent[0]) / np.linalg.norm(source[1])),
            source[0] * (np.linalg.norm(incumbent[1]) / np.linalg.norm(source[0])),
        ]
    )
    expected_residual = expected_scaled - incumbent
    np.testing.assert_allclose(residual, expected_residual)
    np.testing.assert_allclose(candidate, incumbent + 0.1 * expected_residual)


def test_role_overlap_counts_are_enforced() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        role = synthetic_role()
        compact = synthetic_compact(Path(directory))
        with pytest.raises(BUILD.ContractError, match="Pan target overlap changed"):
            BUILD.build_candidate_arrays(
                role,
                compact,
                panel_targets=["PANEL"],
                old_h1_targets=["T_OLD"],
                expected_pan_overlap=4,
                expected_novel_overlap=2,
            )
        with pytest.raises(BUILD.ContractError, match="novel Pan target overlap changed"):
            BUILD.build_candidate_arrays(
                role,
                compact,
                panel_targets=["PANEL"],
                old_h1_targets=["T_OLD"],
                expected_pan_overlap=3,
                expected_novel_overlap=1,
            )


def test_cross_context_aggregation_is_equal_context_not_equal_surface() -> None:
    values = {"orion_hek293t": 1.0}
    values.update({name: 0.0 for name in BUILD.NADIG_ROLES})
    result = BUILD.deterministic_cross_context_aggregate(values)
    assert result == {
        "orion_hek293t": 1.0,
        "nadig_jurkat": 0.0,
        "cross_context": 0.5,
    }
    with pytest.raises(BUILD.ContractError, match="incomplete"):
        BUILD.deterministic_cross_context_aggregate({"orion_hek293t": 1.0})


def test_bootstrap_rng_and_type7_quantile_are_frozen() -> None:
    labels = np.asarray([0, 0, 1, 2, 2], dtype=np.int32)
    first = BUILD.target_cluster_bootstrap_weights(
        labels,
        role="orion_hek293t",
        replicate_zero_based=17,
    )
    second = BUILD.target_cluster_bootstrap_weights(
        labels,
        role="orion_hek293t",
        replicate_zero_based=17,
    )
    np.testing.assert_array_equal(first, second)
    assert first[0] == first[1]
    assert first[3] == first[4]
    assert first.sum() > 0
    assert BUILD.type7_quantile(np.asarray([0.0, 10.0]), 0.25) == pytest.approx(2.5)


def test_output_path_must_remain_inside_experiment() -> None:
    with pytest.raises(BUILD.ContractError, match="must be a child"):
        BUILD._safe_output(Path("/tmp/experiment-117-escape"))


def test_npz_writer_is_byte_deterministic() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        arrays = {
            "z": np.asarray([[1.0, 2.0]], dtype=np.float32),
            "a": np.asarray(["x", "y"]),
        }
        first = root / "first.npz"
        second = root / "second.npz"
        BUILD.write_deterministic_npz(first, arrays)
        BUILD.write_deterministic_npz(second, arrays)
        assert first.read_bytes() == second.read_bytes()
        with np.load(first, allow_pickle=False) as archive:
            assert archive.files == ["a", "z"]
            np.testing.assert_array_equal(archive["z"], arrays["z"])


@pytest.mark.parametrize(
    "array",
    [
        np.empty((0,), dtype=np.float64),
        np.empty((0, 3), dtype=np.int16),
        np.empty((2, 0), dtype=np.float32),
        np.empty((0, 0), dtype=np.bool_),
        np.empty((0,), dtype="U7"),
    ],
)
def test_zero_byte_array_hash_uses_header_and_no_payload(array: np.ndarray) -> None:
    expected = hashlib.sha256()
    expected.update(b"exp107-array-v1\0")
    expected.update(array.dtype.str.encode("ascii"))
    expected.update(b"\0")
    expected.update(
        json.dumps(list(array.shape), separators=(",", ":")).encode("ascii")
    )
    expected.update(b"\0")
    assert array.nbytes == 0
    assert BUILD.sha256_array(array) == expected.hexdigest()


def test_empty_array_bindings_are_stable_and_shape_sensitive() -> None:
    first = np.empty((0, 3), dtype=np.float32)
    second = np.empty((0, 3), dtype=np.float32)
    different_shape = np.empty((3, 0), dtype=np.float32)
    expected = {
        "dtype": "<f4",
        "shape": [0, 3],
        "c_order": True,
        "sha256_array": BUILD.sha256_array(first),
    }
    assert BUILD.array_binding(first) == expected
    assert BUILD.array_binding(second) == expected
    assert (
        BUILD.array_binding(different_shape)["sha256_array"]
        != expected["sha256_array"]
    )


def test_manifest_never_names_held_outcome_arrays() -> None:
    text = (HERE / "SOURCE_INPUT_MANIFEST.json").read_text()
    assert "data/donors/nadig_jurkat/delta.npy" not in text
    assert "h1-test.h5ad" not in text
    assert "orion-hek293t-thinned.h5ad" not in text
