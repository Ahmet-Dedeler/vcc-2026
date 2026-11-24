from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent


def load_freezer():
    spec = importlib.util.spec_from_file_location(
        "experiment_114_freeze_source_candidates_v1",
        HERE / "freeze_source_candidates_v1.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


FREEZER = load_freezer()


def test_checked_in_approval_and_predecessor_chain_pass() -> None:
    approval = FREEZER.load_approval()
    FREEZER.verify_predecessors(approval)
    assert FREEZER.RUN.sha256_file(FREEZER.APPROVAL_PATH) == (
        FREEZER.EXPECTED_APPROVAL_SHA256
    )
    assert approval["authorization"] == FREEZER.EXPECTED_AUTHORIZATION
    assert tuple(approval["approved_roles"]) == ("h1", "orion_hek293t")


def test_effective_authorization_opens_only_source_candidate_freeze() -> None:
    protocol, boundary = FREEZER.effective_authorization(FREEZER.load_approval())
    enabled = {
        key for key, value in protocol["current_authorization"].items() if value
    }
    assert enabled == {
        "small_public_code_and_metadata_inspection",
        "local_synthetic_tests",
        "source_candidate_freeze",
    }
    assert boundary["source_candidate_freeze_authorized"] is True
    assert tuple(boundary["approved_roles"]) == ("h1", "orion_hek293t")


def test_approval_rejects_any_extra_authorization() -> None:
    approval = copy.deepcopy(FREEZER.load_approval())
    approval["authorization"]["held_truth_access"] = True
    try:
        FREEZER.validate_approval(approval)
    except FREEZER.ContractError as error:
        assert "authorizes more than source_candidate_freeze" in str(error)
    else:
        raise AssertionError("held-truth authorization was accepted")


def test_role_manifest_hash_and_output_paths_are_exact() -> None:
    approval = FREEZER.load_approval()
    for role, expected_hash in FREEZER.EXPECTED_ROLE_MANIFEST_SHA256.items():
        record = approval["approved_roles"][role]
        manifest = FREEZER.RUN.resolve_repository_path(
            record["manifest_path"], f"{role} manifest"
        )
        output = FREEZER.RUN.resolve_repository_path(
            record["output_directory"], f"{role} output"
        )
        assert FREEZER.RUN.sha256_file(manifest) == expected_hash
        assert record["output_directory"] == FREEZER.EXPECTED_OUTPUTS[role]
        assert not output.exists()
        assert not output.is_symlink()


def test_operator_mask_and_zeroing_accepts_exact_synthetic_geometry() -> None:
    genes = np.asarray(["PANEL", "HELD", "g1", "g2"])
    targets = np.asarray(["PANEL"])
    incumbent = {"effects": np.asarray([[0.0, 0.0, 2.0, 3.0]])}
    bank = SimpleNamespace(
        targets=np.asarray(["A", "B"]),
        effects=np.asarray(
            [
                [0.0, 0.0, 1.0, 2.0],
                [0.0, 0.0, -1.0, 3.0],
            ]
        ),
        covered=np.asarray(
            [
                [False, False, True, True],
                [False, False, True, True],
            ]
        ),
    )
    operator = SimpleNamespace(
        columns=np.asarray([2, 3], dtype=np.int64),
        held_targets=("HELD",),
    )
    receipt = {
        "parameters": {
            "artifact": {
                "covered_coordinates_per_row": 2,
                "panel_coordinates_zero_and_uncovered": True,
                "held_surface_coordinates_zero_and_uncovered": True,
            }
        }
    }
    report = FREEZER.validate_operator_mask_and_zeroing(
        role="synthetic",
        incumbent=incumbent,
        source_bank=bank,
        canonical_targets=targets,
        canonical_genes=genes,
        operator=operator,
        source_receipt=receipt,
        chunk_rows=1,
    )
    assert report == {
        "operator_column_count": 2,
        "panel_coordinate_count": 1,
        "held_coordinate_count_on_gene_axis": 1,
    }


def test_operator_mask_and_zeroing_rejects_extra_covered_coordinate() -> None:
    genes = np.asarray(["PANEL", "HELD", "g1"])
    bank = SimpleNamespace(
        targets=np.asarray(["A"]),
        effects=np.asarray([[0.0, 0.0, 1.0]]),
        covered=np.asarray([[True, False, True]]),
    )
    try:
        FREEZER.validate_operator_mask_and_zeroing(
            role="synthetic",
            incumbent={"effects": np.asarray([[0.0, 0.0, 1.0]])},
            source_bank=bank,
            canonical_targets=np.asarray(["PANEL"]),
            canonical_genes=genes,
            operator=SimpleNamespace(
                columns=np.asarray([2], dtype=np.int64),
                held_targets=("HELD",),
            ),
            source_receipt={
                "parameters": {
                    "artifact": {
                        "covered_coordinates_per_row": 1,
                        "panel_coordinates_zero_and_uncovered": True,
                        "held_surface_coordinates_zero_and_uncovered": True,
                    }
                }
            },
        )
    except FREEZER.ContractError as error:
        assert "coverage differs from the operator mask" in str(error)
    else:
        raise AssertionError("extra covered coordinate was accepted")


def test_freeze_requires_explicit_confirmation() -> None:
    try:
        FREEZER.freeze_role("h1", confirmed=False)
    except FREEZER.ContractError as error:
        assert "--confirm-source-only-freeze" in str(error)
    else:
        raise AssertionError("real freeze ran without explicit confirmation")


def test_real_role_locks_and_source_only_preflight_pass() -> None:
    report = FREEZER.audit_all()
    assert report["ready_for_independent_review"] is True
    assert report["real_candidates_frozen"] is False
    assert report["held_truth_access_authorized"] is False
    assert report["evaluation_authorized"] is False
    assert report["count_generation_authorized"] is False
    assert report["submission_authorized"] is False
    assert tuple(report["roles"]) == ("h1", "orion_hek293t")
    for role, role_report in report["roles"].items():
        assert role_report["manifest_sha256"] == (
            FREEZER.EXPECTED_ROLE_MANIFEST_SHA256[role]
        )
        assert role_report["output_absent"] is True
        assert role_report["held_truth_opened"] is False
        assert role_report["eligibility"]["eligible_targets"]
