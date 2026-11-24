from __future__ import annotations

import copy
import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_builder():
    spec = importlib.util.spec_from_file_location(
        "experiment_114_role_builder",
        HERE / "build_role_inputs.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


BUILDER = load_builder()
RUN = BUILDER.RUN


def test_checked_in_role_source_manifest_and_locks() -> None:
    manifest = BUILDER.load_role_source_manifest()
    BUILDER.validate_checked_in_locks()
    assert tuple(manifest["roles"]) == ("h1", "orion_hek293t")
    assert manifest["canonical_axes"]["targets"]["count"] == 300
    assert manifest["canonical_axes"]["genes"]["count"] == 18_533
    assert manifest["eligibility"]["eligible_count"] == 41
    assert manifest["authorization"]["source_candidate_freeze"] is False
    assert manifest["authorization"]["approved_roles"] == {}


def test_canonical_axis_records_verify_with_axis_metadata() -> None:
    manifest = BUILDER.load_role_source_manifest()
    for name, record in manifest["canonical_axes"].items():
        assert BUILDER.verify_axis_record(record, f"canonical {name} axis").is_file()


def test_role_source_policies_are_exact_and_leave_surface_out() -> None:
    roles = BUILDER.load_role_source_manifest()["roles"]
    assert roles["h1"]["allowed_source_ids"] == [
        "k562_gwps",
        "orion_hct116",
        "orion_hek293t",
    ]
    assert roles["h1"]["held_surface_source_id"] == "h1_2025"
    assert roles["h1"]["forbidden_source_ids"] == ["h1_2025"]
    assert roles["orion_hek293t"]["allowed_source_ids"] == [
        "k562_gwps",
        "orion_hct116",
    ]
    assert roles["orion_hek293t"]["held_surface_source_id"] == "orion_hek293t"
    assert roles["orion_hek293t"]["forbidden_source_ids"] == ["orion_hek293t"]
    for role in roles.values():
        assert not (set(role["allowed_source_ids"]) & set(role["forbidden_source_ids"]))


def test_role_operator_geometry_matches_experiment_107_contract() -> None:
    manifest = BUILDER.load_role_source_manifest()
    upstream = RUN.load_json(
        RUN.ROOT / "experiments/107-jurkat-high-moi/strict_validation_contract.json"
    )
    for role, record in manifest["roles"].items():
        operator = record["operator"]
        surface = upstream["held_context_gate"]["surfaces"][role]
        assert operator["target_exclusion_count"] == surface["target_count"]
        assert (
            operator["target_exclusion_ordered_sha256"] == surface["target_axis_sha256"]
        )
        assert operator["column_count"] == surface["whitener_column_count"]
        assert (
            operator["columns_int64_sha256"] == surface["whitener_columns_int64_sha256"]
        )
        assert (
            operator["basis_rows_by_donor"] == surface["whitener_basis_rows_by_donor"]
        )


def test_role_source_donor_hashes_match_experiment_107_receipts() -> None:
    manifest = BUILDER.load_role_source_manifest()
    donor_records = manifest["source_inputs"]["donors"]
    for donor in donor_records:
        upstream_path = (
            RUN.ROOT
            / "experiments/107-jurkat-high-moi/work/strict-validation"
            / donor
            / "donor.json"
        )
        upstream = RUN.load_json(upstream_path)
        bindings = {
            Path(record["path"]).name: record["sha256"]
            for record in upstream["transitive_inputs"]
            if record["path"].startswith(f"data/donors/{donor}/")
        }
        assert {
            Path(record["path"]).name: record["sha256"]
            for record in donor_records[donor].values()
        } == bindings


def test_exact_a14_setting_drift_fails_closed() -> None:
    manifest = BUILDER.load_role_source_manifest()
    drifted = copy.deepcopy(manifest)
    drifted["a14"]["alpha"] = 0.983
    try:
        BUILDER.validate_role_source_manifest(drifted)
    except RUN.ContractError as error:
        assert "a14 settings" in str(error)
    else:
        raise AssertionError("drifted a14 settings were accepted")


def test_duplicate_gene_policy_keeps_last_occurrence() -> None:
    target_columns, source_columns = BUILDER.last_occurrence_map(
        ["g2", "g1", "g2", "outside"],
        ["g1", "g2", "g3"],
    )
    np.testing.assert_array_equal(target_columns, [0, 1])
    np.testing.assert_array_equal(source_columns, [1, 2])


def test_neighbor_closure_is_leave_target_out_and_exactly_25_under_ties() -> None:
    sources = ["QUERY", *(f"S{index:02d}" for index in range(30))]
    genes = np.asarray(sources)
    matrix = np.ones((len(genes), 1), dtype=np.float32)
    features = RUN.GeneFeatures(genes, matrix)
    selected, audit = BUILDER.source_target_closure(
        sources,
        ["QUERY"],
        {arm: features for arm in BUILDER.SOURCE_BANK_ARMS},
    )
    assert selected == tuple(f"S{index:02d}" for index in range(25))
    for arm in BUILDER.SOURCE_BANK_ARMS:
        assert audit["arms"][arm]["query_neighbor_pair_count"] == 25
        assert audit["arms"][arm]["unsupported_queries"] == []


def synthetic_donor() -> BUILDER.RoleDonor:
    delta = np.asarray(
        [
            [1.0, 2.0],
            [3.0, 6.0],
        ],
        dtype=np.float32,
    )
    mean = delta.mean(axis=0)
    centered = delta - mean
    scale = float(np.median(np.linalg.norm(centered, axis=1)))
    return BUILDER.RoleDonor(
        name="source",
        delta_path=HERE / "test_build_role_inputs.py",
        genes=np.asarray(["g1", "g2"]),
        perts=np.asarray(["A", "B"]),
        mean=mean,
        tcols_all=np.asarray([2, 3], dtype=np.int64),
        dcols_all=np.asarray([0, 1], dtype=np.int64),
        tcols_unique=np.asarray([2, 3], dtype=np.int64),
        dcols_unique=np.asarray([0, 1], dtype=np.int64),
        pert_index={"A": 0, "B": 1},
        scale=scale,
        delta=delta,
    )


def test_source_bank_uses_canonical_axis_and_masks_panel_and_held_coordinates() -> None:
    donor = synthetic_donor()
    operator = BUILDER.RoleOperator(
        evals=np.ones(2, dtype=np.float64),
        evecs=np.eye(2, dtype=np.float64),
        columns=np.asarray([2, 3], dtype=np.int64),
        held_targets=("HELD",),
        artifact_path=HERE / "test_build_role_inputs.py",
        manifest_path=HERE / "test_build_role_inputs.py",
        artifact_sha256=RUN.sha256_file(HERE / "test_build_role_inputs.py"),
        manifest_sha256=RUN.sha256_file(HERE / "test_build_role_inputs.py"),
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        arrays, audit = BUILDER.build_source_bank_arrays(
            donors={"source": donor},
            allowed_source_ids=["source"],
            canonical_genes=np.asarray(["PANEL", "unused", "g1", "g2", "HELD"]),
            canonical_targets=np.asarray(["PANEL"]),
            held_targets=["HELD"],
            operator=operator,
            selected_targets=["B", "A"],
            staging=Path(directory),
            chunk_rows=1,
        )
        assert arrays["targets"].tolist() == ["A", "B"]
        assert arrays["source_ids"].tolist() == ["source", "source"]
        assert arrays["effects"].shape == (2, 5)
        assert arrays["effects"][:, [0, 1, 4]].tolist() == [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
        assert not arrays["covered"][:, [0, 1, 4]].any()
        assert arrays["covered"][:, [2, 3]].all()
        assert audit["panel_coordinates_zero_and_uncovered"] is True
        assert audit["held_surface_coordinates_zero_and_uncovered"] is True
        del arrays


def test_deterministic_npz_has_stable_bytes_and_refuses_overwrite() -> None:
    arrays = {
        "targets": np.asarray(["A", "B"]),
        "effects": np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64),
    }
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        first = root / "first.npz"
        second = root / "second.npz"
        BUILDER.deterministic_npz(first, arrays)
        BUILDER.deterministic_npz(second, arrays)
        assert RUN.sha256_file(first) == RUN.sha256_file(second)
        with np.load(first, allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive["targets"], arrays["targets"])
            np.testing.assert_array_equal(archive["effects"], arrays["effects"])
        try:
            BUILDER.deterministic_npz(first, arrays)
        except RUN.ContractError as error:
            assert "refusing to overwrite" in str(error)
        else:
            raise AssertionError("deterministic NPZ writer overwrote an artifact")


def test_role_bundle_manifest_is_compatible_with_run_validators(monkeypatch) -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        fake_here = Path(directory)
        role = "synthetic"
        output = fake_here / "work" / "roles" / role
        relative_output = output.relative_to(RUN.ROOT.resolve()).as_posix()
        monkeypatch.setattr(RUN, "HERE", fake_here)
        role_spec = {
            "surface": "synthetic",
            "held_surface_source_id": "held",
            "allowed_source_ids": ["source"],
            "forbidden_source_ids": ["held"],
            "incumbent_endpoint": "a14-role-heldout-canonical-300-v1",
            "source_effect_space": (
                "a14_role_donor_deviation_role_whitened_no_global_renorm_v1"
            ),
            "output_directory": relative_output,
        }
        incumbent = {
            "targets": np.asarray(["PANEL"]),
            "genes": np.asarray(["PANEL", "g1"]),
            "effects": np.asarray([[0.0, 1.0]], dtype=np.float64),
            "incumbent_fallback": np.asarray([True]),
            "source_carrier_count": np.asarray([1], dtype=np.int16),
            "expected_source_count": np.asarray([3], dtype=np.int16),
        }
        source = {
            "targets": np.asarray(["A"]),
            "genes": np.asarray(["PANEL", "g1"]),
            "effects": np.asarray([[0.0, 2.0]], dtype=np.float64),
            "covered": np.asarray([[False, True]]),
            "source_ids": np.asarray(["source"]),
        }
        fixture_record = BUILDER.receipt_record(HERE / "test_build_role_inputs.py")
        operator = BUILDER.RoleOperator(
            evals=np.ones(1),
            evecs=np.eye(1),
            columns=np.asarray([1], dtype=np.int64),
            held_targets=("HELD",),
            artifact_path=HERE / "test_build_role_inputs.py",
            manifest_path=HERE / "test_build_role_inputs.py",
            artifact_sha256=fixture_record["sha256"],
            manifest_sha256=fixture_record["sha256"],
        )
        role_manifest = BUILDER.seal_role_bundle(
            role=role,
            role_spec=role_spec,
            incumbent_arrays=incumbent,
            source_arrays=source,
            incumbent_audit={"synthetic": True},
            source_audit={"synthetic": True},
            neighbor_audit={"synthetic": True},
            operator=operator,
            manifest=BUILDER.load_role_source_manifest(),
            output_directory=output,
            incumbent_source_inputs=[fixture_record],
            source_bank_source_inputs=[fixture_record],
        )
        approved = {
            "surface": "synthetic",
            "role_directory": relative_output,
        }
        RUN.validate_input_manifest(role_manifest, role=role, approved=approved)
        role_directory = RUN.validate_role_manifest(
            role_manifest,
            role=role,
            approved=approved,
        )
        incumbent_receipt = RUN.load_json(role_directory / "incumbent-receipt.json")
        source_receipt = RUN.load_json(role_directory / "source-bank-receipt.json")
        RUN.validate_role_receipt(
            incumbent_receipt,
            stage="role_incumbent_source_only",
            role=role,
            surface="synthetic",
            artifact_path=role_directory / "incumbent.npz",
            artifact_sha256=role_manifest["incumbent"]["artifact"]["sha256"],
        )
        RUN.validate_role_receipt(
            source_receipt,
            stage="role_source_bank_source_only",
            role=role,
            surface="synthetic",
            artifact_path=role_directory / "source-bank.npz",
            artifact_sha256=role_manifest["source_bank"]["artifact"]["sha256"],
        )


def test_real_build_refuses_to_overlap_experiment_107(monkeypatch) -> None:
    monkeypatch.setattr(BUILDER, "active_exp107_build_pids", lambda: (12345,))
    try:
        BUILDER.build_role(
            "h1",
            confirmed=True,
            chunk_rows=1,
            max_rss_gib=1.0,
        )
    except RUN.ContractError as error:
        assert "Experiment 107 is still building" in str(error)
    else:
        raise AssertionError("real role build overlapped Experiment 107")


def test_metadata_command_is_read_only_and_reports_unapproved_roles() -> None:
    report = BUILDER.metadata_report()
    assert report["source_candidate_freeze"] is False
    assert report["approved_roles"] == {}
    assert tuple(report["roles"]) == ("h1", "orion_hek293t")
    assert all(
        "build-role --role" in item["future_command"]
        for item in report["roles"].values()
    )
