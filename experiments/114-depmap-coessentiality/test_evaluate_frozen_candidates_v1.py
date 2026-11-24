from __future__ import annotations

import copy
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent


def load_evaluator():
    spec = importlib.util.spec_from_file_location(
        "experiment_114_evaluate_frozen_candidates_v1",
        HERE / "evaluate_frozen_candidates_v1.py",
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


EVALUATOR = load_evaluator()


def test_checked_in_approval_is_exact_and_non_promoting() -> None:
    approval = EVALUATOR.load_approval()
    assert EVALUATOR.RUN.sha256_file(EVALUATOR.APPROVAL_PATH) == (
        EVALUATOR.EXPECTED_APPROVAL_SHA256
    )
    assert approval["authorization"] == EVALUATOR.EXPECTED_AUTHORIZATION
    assert approval["surface"]["id"] == "orion_hek293t_full_outer_holdout"
    assert approval["surface"]["diagnostic_only"] is True
    assert approval["surface"]["direction_valid"] is False
    assert approval["surface"]["metric"]["name"] == "raw_uncentered_pds"
    assert approval["surface"]["metric"]["promotion_authorized"] is False
    assert approval["invariants"]["promotion_authorized"] is False


def test_approval_rejects_promotion_or_extra_capability() -> None:
    approval = EVALUATOR.load_approval()
    promoted = copy.deepcopy(approval)
    promoted["surface"]["metric"]["promotion_authorized"] = True
    with pytest.raises(EVALUATOR.ContractError, match="metric contract drifted"):
        EVALUATOR.validate_approval(promoted)

    expanded = copy.deepcopy(approval)
    expanded["authorization"]["model_fitting"] = True
    with pytest.raises(EVALUATOR.ContractError, match="authorization drifted"):
        EVALUATOR.validate_approval(expanded)


def test_cli_accepts_no_caller_supplied_paths() -> None:
    parsed = EVALUATOR.parser().parse_args(["audit"])
    assert vars(parsed) == {"command": "audit"}
    parsed = EVALUATOR.parser().parse_args(
        ["evaluate", "--confirm-held-evaluation"]
    )
    assert vars(parsed) == {
        "command": "evaluate",
        "confirm_held_evaluation": True,
    }
    with pytest.raises(SystemExit):
        EVALUATOR.parser().parse_args(["evaluate", "--out", "elsewhere"])
    with pytest.raises(SystemExit):
        EVALUATOR.parser().parse_args(
            ["evaluate", "--truth-manifest", "elsewhere"]
        )


def test_pretruth_audit_validates_bundles_without_opening_truth(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    approval = EVALUATOR.load_approval()
    held_paths = {
        EVALUATOR.RUN.resolve_repository_path(
            approval["surface"]["truth"][name]["path"], name
        )
        for name in ("delta", "b_ctrl")
    }
    original_hash = EVALUATOR.RUN.sha256_file
    original_np_load = EVALUATOR.np.load

    def guarded_hash(path: Path, *args, **kwargs):
        resolved = Path(path).resolve()
        if resolved in held_paths:
            raise AssertionError(f"pretruth audit hashed held data: {resolved}")
        return original_hash(path, *args, **kwargs)

    def guarded_np_load(file, *args, **kwargs):
        if isinstance(file, (str, Path)) and Path(file).resolve() in held_paths:
            raise AssertionError(f"pretruth audit opened held data: {file}")
        return original_np_load(file, *args, **kwargs)

    monkeypatch.setattr(EVALUATOR.RUN, "sha256_file", guarded_hash)
    monkeypatch.setattr(EVALUATOR.np, "load", guarded_np_load)
    report = EVALUATOR.audit_pretruth()
    public = EVALUATOR._public_audit(report)
    assert public["all_frozen_bundles_validated"] is True
    assert public["null_mappings_regenerated"] is True
    assert public["changed_target_count"] == 35
    assert public["changed_target_ordered_sha256"] == (
        "e81c113e76a627b4d8d44a6998d33f3fc533ee11e539b413016552a0c35b3587"
    )
    assert public["truth_matrix_hashed"] is False
    assert public["truth_matrix_opened"] is False
    assert public["output_absent"] is True


def test_access_receipt_is_exclusive_and_records_preopen_state() -> None:
    approval = EVALUATOR.load_approval()
    work = HERE / "work"
    with tempfile.TemporaryDirectory(dir=work) as temporary:
        directory = Path(temporary) / "out"
        receipt = directory / "OUTCOME_ACCESS_RECEIPT.v1.json"
        result = directory / "result.json"
        EVALUATOR.create_access_receipt(
            approval,
            directory,
            receipt,
            result,
        )
        payload = json.loads(receipt.read_text())
        assert payload["truth_hashed_before_receipt"] is False
        assert payload["truth_opened_before_receipt"] is False
        assert payload["diagnostic_only"] is True
        assert payload["promotion_authorized"] is False
        with pytest.raises(
            EVALUATOR.ContractError, match="output already exists"
        ):
            EVALUATOR.create_access_receipt(
                approval,
                directory,
                receipt,
                result,
            )


def test_real_evaluation_requires_explicit_confirmation() -> None:
    with pytest.raises(
        EVALUATOR.ContractError, match="--confirm-held-evaluation"
    ):
        EVALUATOR.evaluate(confirmed=False)


def test_receipt_creation_precedes_any_truth_verification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    marker = {"receipt": False}
    approval = EVALUATOR.load_approval()
    arrays = {
        "genes": np.asarray(["g1", "g2"]),
    }
    monkeypatch.setattr(
        EVALUATOR,
        "audit_pretruth",
        lambda **kwargs: {
            "_approval": approval,
            "_loaded": {
                EVALUATOR.EXPECTED_ROLE: (
                    {"role": EVALUATOR.EXPECTED_ROLE},
                    arrays,
                )
            },
            "_mapping": {},
        },
    )
    with tempfile.TemporaryDirectory(dir=HERE / "work") as temporary:
        directory = Path(temporary) / "out"
        receipt = directory / "OUTCOME_ACCESS_RECEIPT.v1.json"
        result = directory / "result.json"
        monkeypatch.setattr(
            EVALUATOR,
            "resolve_output_paths",
            lambda _approval: (directory, receipt, result),
        )

        def fake_receipt(*args, **kwargs):
            marker["receipt"] = True

        class StopAfterReceipt(RuntimeError):
            pass

        def fake_truth(*args, **kwargs):
            assert marker["receipt"] is True
            raise StopAfterReceipt

        monkeypatch.setattr(EVALUATOR, "create_access_receipt", fake_receipt)
        monkeypatch.setattr(
            EVALUATOR, "_verify_truth_bundle_after_receipt", fake_truth
        )
        with pytest.raises(StopAfterReceipt):
            EVALUATOR.evaluate(confirmed=True)


def test_result_contract_cannot_claim_promotion_or_directionality() -> None:
    approval = EVALUATOR.load_approval()
    assert approval["authorization"]["evaluation"] is True
    assert approval["authorization"]["held_truth_access"] is True
    for forbidden in (
        "source_candidate_freeze",
        "model_fitting",
        "count_generation",
        "challenge_build",
        "submission",
        "commit",
        "push",
    ):
        assert approval["authorization"][forbidden] is False
    assert approval["surface"]["direction_valid"] is False
    assert approval["surface"]["diagnostic_only"] is True
