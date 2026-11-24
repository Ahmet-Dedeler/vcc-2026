"""Truth-free synthetic tests for Experiment 117's held evaluator."""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import stat
import sys
from pathlib import Path

import numpy as np
import pytest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "experiment_117_evaluate_held_v1",
    HERE / "evaluate_held_v1.py",
)
assert SPEC is not None and SPEC.loader is not None
EVAL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EVAL
SPEC.loader.exec_module(EVAL)


def bundle(role: str, arrays: dict[str, np.ndarray]) -> object:
    return EVAL.SourceBundle(
        role=role,
        manifest_path=Path(f"/synthetic/{role}/source-candidates.json"),
        npz_path=Path(f"/synthetic/{role}/source-candidates.npz"),
        manifest={},
        arrays=arrays,
    )


def complete_surfaces(value: float = 0.01) -> dict[str, float]:
    return {role: value for role in EVAL.CROSS_CONTEXT_ROLES}


class GuardedArrays(dict[str, np.ndarray]):
    def __init__(
        self,
        values: dict[str, np.ndarray],
        forbidden: set[str],
    ) -> None:
        super().__init__(values)
        self.forbidden = forbidden
        self.accessed: list[str] = []

    def __getitem__(self, key: str) -> np.ndarray:
        if key in self.forbidden:
            raise AssertionError(f"forbidden held arm accessed: {key}")
        self.accessed.append(key)
        return super().__getitem__(key)


def test_pds_uses_midranks_and_zero_prediction_rows() -> None:
    truth = np.asarray([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    prediction = np.asarray([[1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    np.testing.assert_allclose(
        EVAL.pds_per_target(prediction, truth),
        np.asarray([0.75, 0.75, 0.5]),
    )


def test_directional_precision_counts_missing_zero_truth_as_incorrect() -> None:
    base = np.zeros((2, 10), dtype=np.float64)
    base[0, :2] = [-1.0, 1.0]
    base[1, :2] = [-1.0, -1.0]
    candidate = base.copy()
    candidate[0, :2] = [1.0, 1.0]
    candidate[1, :2] = [0.0, 1.0]
    arrays = {
        "identity_f000": base,
        "pan_novel_f100": candidate,
        "target_effective_novel": np.asarray([True, True]),
        "coordinate_supported": np.ones((2, 10), dtype=bool),
        "directional_rank_indices": np.tile(
            np.arange(10, dtype=np.int32), (2, 1)
        ),
    }

    def reference_sign(row: int, selected: np.ndarray) -> np.ndarray:
        assert selected.tolist() == [0, 1]
        return (
            np.asarray([1, 0], dtype=np.int8)
            if row == 0
            else np.asarray([-1, 1], dtype=np.int8)
        )

    result = EVAL.directional_precision_delta(
        bundle("synthetic", arrays),
        reference_sign,
        arm="pan_novel_f100",
        effective_mask="target_effective_novel",
        ranking_name="directional_rank_indices",
    )
    np.testing.assert_allclose(result["candidate"], [0.5, 0.5])
    np.testing.assert_allclose(result["incumbent"], [0.0, 0.5])
    np.testing.assert_allclose(result["delta"], [0.5, 0.0])
    assert result["surface_delta"] == pytest.approx(0.25)


def test_h1_cannot_enter_promotion_and_all_cross_context_roles_are_mandatory() -> None:
    draws = np.full(EVAL.BOOTSTRAP_REPLICATES, 0.01)
    bootstrap = {
        "repetitions": EVAL.BOOTSTRAP_REPLICATES,
        "lower_95": 0.01,
        "upper_95": 0.01,
        "draws": draws,
    }
    report = EVAL.promotion_report(
        complete_surfaces(),
        complete_surfaces(0.1),
        np.zeros(EVAL.NULL_REPLICATES),
        bootstrap,
        invariants_pass=True,
    )
    assert report["passes"] is True
    assert report["h1_used_for_promotion"] is False

    with pytest.raises(EVAL.EvaluationError, match="mandatory"):
        EVAL.promotion_report(
            {**complete_surfaces(), "h1": 1000.0},
            complete_surfaces(0.1),
            np.zeros(EVAL.NULL_REPLICATES),
            bootstrap,
            invariants_pass=True,
        )
    missing = complete_surfaces()
    missing.pop("nadig.gallery_6")
    with pytest.raises(EVAL.EvaluationError, match="mandatory"):
        EVAL.aggregate_cross_context(missing)


def test_cross_context_aggregation_is_equal_gallery_then_equal_context() -> None:
    surfaces = {"orion_hek293t": 1.0}
    surfaces.update({role: 0.0 for role in EVAL.NADIG_ROLES})
    assert EVAL.aggregate_cross_context(surfaces) == {
        "orion_hek293t": 1.0,
        "nadig_jurkat": 0.0,
        "cross_context": 0.5,
    }


@pytest.mark.parametrize(
    ("role", "allowed"),
    [
        ("h1", EVAL.H1_HELD_ARMS),
        ("orion_hek293t", EVAL.AUTHORIZING_HELD_ARMS),
        ("nadig.gallery_1", EVAL.AUTHORIZING_HELD_ARMS),
    ],
)
def test_forbidden_held_arms_are_never_accessed(
    role: str,
    allowed: tuple[str, str],
) -> None:
    all_arms = {
        "identity_f000",
        "pan_novel_f050",
        "pan_novel_f100",
        "pan_novel_f200",
        "pan_all_f050",
        "pan_all_f100",
        "pan_all_f200",
    }
    values = {
        "score_indices": np.asarray([0, 1], dtype=np.int64),
        **{
            arm: np.asarray([[1.0, 0.0], [0.0, 1.0]])
            for arm in all_arms
        },
    }
    guarded = GuardedArrays(values, all_arms - set(allowed))
    scored = EVAL._score_held_arms(
        bundle(role, guarded),
        np.asarray([[1.0, 0.0], [0.0, 1.0]]),
    )
    assert tuple(scored) == allowed
    assert set(guarded.accessed) == {"score_indices", *allowed}


def test_forbidden_held_arms_are_never_persisted() -> None:
    directional = {
        "rows": np.asarray([0]),
        "incumbent": np.asarray([0.0]),
        "candidate": np.asarray([1.0]),
        "delta": np.asarray([1.0]),
    }
    output: dict[str, np.ndarray] = {}
    EVAL.add_role_result_arrays(
        output,
        "orion_hek293t",
        {
            "pds": {
                "identity_f000": np.asarray([0.5]),
                "pan_novel_f100": np.asarray([0.6]),
            },
            "directional": directional,
        },
    )
    assert "orion_hek293t_identity_f000_pds" in output
    assert "orion_hek293t_pan_novel_f100_pds" in output
    assert not any("f050" in name or "f200" in name or "pan_all" in name for name in output)

    forbidden_output: dict[str, np.ndarray] = {}
    with pytest.raises(EVAL.EvaluationError, match="forbidden held arm"):
        EVAL.add_role_result_arrays(
            forbidden_output,
            "orion_hek293t",
            {
                "pds": {
                    "identity_f000": np.asarray([0.5]),
                    "pan_novel_f100": np.asarray([0.6]),
                    "pan_novel_f200": np.asarray([0.7]),
                },
                "directional": directional,
            },
        )
    assert forbidden_output == {}


def test_misordered_nadig_score_columns_fail_closed() -> None:
    role = bundle(
        "nadig.gallery_1",
        {
            "genes": np.asarray(["A", "B", "C"]),
            "score_indices": np.asarray([0, 2], dtype=np.int64),
        },
    )
    passing = {"score_indices_in_nadig": [1, 0]}
    truth_genes = ("C", "A", "B")
    truth_indices, official_indices = EVAL._nadig_score_mapping(
        role,
        passing,
        truth_genes,
    )
    np.testing.assert_array_equal(truth_indices, [1, 0])
    np.testing.assert_array_equal(official_indices, [0, 2])

    with pytest.raises(EVAL.EvaluationError, match="score gene order changed"):
        EVAL._nadig_score_mapping(
            role,
            {"score_indices_in_nadig": [0, 1]},
            truth_genes,
        )


def test_frozen_null_reconstruction_moves_source_and_support_without_mutation() -> None:
    incumbent = np.asarray(
        [[2.0, 4.0, 8.0], [3.0, 6.0, 9.0], [5.0, 10.0, 20.0]],
        dtype=np.float32,
    )
    source = np.asarray(
        [[1.0, 3.0, 50.0], [7.0, 8.0, 9.0], [2.0, 1.0, 25.0]],
        dtype=np.float32,
    )
    support = np.asarray(
        [[True, True, False], [False, False, False], [True, False, True]]
    )
    novel = np.asarray([True, False, True])
    mapping = np.asarray([2, -1, 0], dtype=np.int32)
    originals = [value.tobytes() for value in (incumbent, source, support, novel, mapping)]

    residual, effective, moved_support = EVAL.normalized_residual_from_mapping(
        incumbent, source, support, novel, mapping
    )
    expected0 = source[2, [0, 2]].astype(np.float64)
    expected0 *= np.linalg.norm(
        incumbent[0, [0, 2]].astype(np.float64)
    ) / np.linalg.norm(expected0)
    expected2 = source[0, [0, 1]].astype(np.float64)
    expected2 *= np.linalg.norm(
        incumbent[2, [0, 1]].astype(np.float64)
    ) / np.linalg.norm(expected2)
    np.testing.assert_allclose(residual[0, [0, 2]], expected0 - incumbent[0, [0, 2]])
    np.testing.assert_allclose(residual[2, [0, 1]], expected2 - incumbent[2, [0, 1]])
    np.testing.assert_array_equal(effective, [True, False, True])
    np.testing.assert_array_equal(moved_support[0], support[2])
    np.testing.assert_array_equal(moved_support[2], support[0])
    assert [value.tobytes() for value in (incumbent, source, support, novel, mapping)] == originals

    with pytest.raises(EVAL.EvaluationError, match="out of bounds"):
        EVAL.normalized_residual_from_mapping(
            incumbent, source, support, novel, np.asarray([3, -1, 0])
        )


def test_cluster_bootstrap_is_deterministic_and_uses_type7_quantiles() -> None:
    deltas = {
        role: np.asarray([0.0, 1.0, 2.0, 4.0], dtype=np.float64)
        for role in EVAL.CROSS_CONTEXT_ROLES
    }
    labels = {
        role: np.asarray([0, 0, 2, 3], dtype=np.int32)
        for role in EVAL.CROSS_CONTEXT_ROLES
    }
    first = EVAL.cluster_bootstrap(deltas, labels, repetitions=31)
    second = EVAL.cluster_bootstrap(deltas, labels, repetitions=31)
    np.testing.assert_array_equal(first["draws"], second["draws"])
    assert first["lower_95"] == float(
        np.quantile(first["draws"], 0.025, method="linear")
    )
    assert first["upper_95"] == float(
        np.quantile(first["draws"], 0.975, method="linear")
    )


def test_authentication_never_invokes_truth_verification(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path
    here = root / "experiments/117-pan-h1-direct-residual"
    source_root = here / "work/source-candidates"
    source_root.mkdir(parents=True)
    protocol_path = here / "HELD_EVALUATION_PROTOCOL.v1.json"
    authorization_path = here / "HELD_EVALUATION_AUTHORIZATION.v1.json"
    entrypoint_path = here / "evaluate_held_v1.py"
    top_path = source_root / "source-candidate-inputs.json"
    plan_path = root / "plan.json"
    receipt_path = here / "work/held-evaluation-v1-access-receipt.json"
    result_directory = here / "work/held-evaluation-v1"

    protocol = {
        "schema_version": 1,
        "freeze_revision": 2,
        "experiment": EVAL.EXPERIMENT,
        "status": (
            "frozen_before_experiment_117_held_truth_access_"
            "after_independent_audit_repair"
        ),
        "held_outcome_exposure": EVAL._protocol_held_outcome_exposure(),
    }
    protocol_path.write_text(json.dumps(protocol))
    entrypoint_path.write_text("# synthetic\n")
    plan_path.write_text(
        json.dumps(
            {
                "status": "pretruth_complete",
                "nadig_outcome_opened": False,
                "galleries": [],
            }
        )
    )
    role_bindings = {
        role: {
            "manifest": {
                "path": f"synthetic/{role}/source-candidates.json",
                "bytes": 1,
                "sha256": role.encode().hex().ljust(64, "0")[:64],
            },
            "npz": {
                "path": f"synthetic/{role}/source-candidates.npz",
                "bytes": 1,
                "sha256": role.encode().hex().ljust(64, "1")[:64],
            },
        }
        for role in EVAL.ALL_ROLES
    }
    top = {
        "status": "source_candidates_frozen_before_held_outcome_access",
        "held_truth_opened": False,
        "roles": {
            role: {
                "manifest": role_bindings[role]["manifest"]["path"],
                "sha256": role_bindings[role]["manifest"]["sha256"],
            }
            for role in EVAL.ALL_ROLES
        },
    }
    top_path.write_text(json.dumps(top))

    def binding(path: Path) -> dict[str, object]:
        return {
            "path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": EVAL.sha256_file(path),
        }

    declaration = {"path": "truth/never-open", "bytes": 1, "sha256": "a" * 64}
    authorization = {
        "schema_version": 1,
        "freeze_revision": 2,
        "experiment": EVAL.EXPERIMENT,
        "status": (
            "held_evaluation_v1_repaired_and_refrozen_"
            "pending_independent_review"
        ),
        "authorization": {
            "source_candidate_read": True,
            "held_truth_access_after_receipt": True,
            "held_evaluation_after_receipt": True,
            "model_refit": False,
            "arm_or_fraction_selection": False,
            "cell_eval2": False,
            "count_generation": False,
            "challenge_build": False,
            "submission": False,
            "commit": False,
            "push": False,
        },
        "execution_contract": {
            "access_receipt": receipt_path.relative_to(root).as_posix(),
            "result_directory": result_directory.relative_to(root).as_posix(),
            "receipt_must_be_exclusive_fsynced_and_read_only_before_truth_hash_or_load": True,
            "receipt_survives_failure_and_prevents_replay": True,
            "caller_supplied_paths": False,
            "result_overwrite": False,
            "result_manifest_immutable": True,
            "required_null_replicates": EVAL.NULL_REPLICATES,
            "required_cluster_bootstrap_replicates": EVAL.BOOTSTRAP_REPLICATES,
            "required_cross_context_roles": list(EVAL.CROSS_CONTEXT_ROLES),
            "h1_diagnostic_cannot_promote": True,
        },
        "held_output_contract": EVAL._held_output_contract(),
        "protocol": binding(protocol_path),
        "evaluator": binding(entrypoint_path),
        "predecessors": {"experiment_107_pretruth_plan": binding(plan_path)},
        "source_candidates": {
            "top_manifest": binding(top_path),
            "roles": list(EVAL.ALL_ROLES),
        },
        "role_bindings": role_bindings,
        "role_mappings": {role: {} for role in EVAL.ALL_ROLES},
        "truth_bindings": {
            "h1": {
                "harness": declaration,
                "reference_de": declaration,
                "scored_arms": list(EVAL.H1_HELD_ARMS),
            },
            "orion_hek293t": {
                "harness": declaration,
                "reference_de": declaration,
                "scored_arms": list(EVAL.AUTHORIZING_HELD_ARMS),
            },
            "nadig_jurkat": {
                "files": {
                    "delta": declaration,
                    "genes": declaration,
                    "targets": declaration,
                },
                "scored_arms": list(EVAL.AUTHORIZING_HELD_ARMS),
            },
        },
    }
    authorization_path.write_text(json.dumps(authorization))

    monkeypatch.setattr(EVAL, "ROOT", root)
    monkeypatch.setattr(EVAL, "HERE", here)
    monkeypatch.setattr(EVAL, "PROTOCOL_PATH", protocol_path)
    monkeypatch.setattr(EVAL, "AUTHORIZATION_PATH", authorization_path)
    monkeypatch.setattr(EVAL, "ENTRYPOINT_PATH", entrypoint_path)
    monkeypatch.setattr(EVAL, "EXPECTED_PROTOCOL_SHA256", EVAL.sha256_file(protocol_path))
    monkeypatch.setattr(EVAL, "SOURCE_ROOT", source_root)
    monkeypatch.setattr(EVAL, "SOURCE_TOP_PATH", top_path)
    monkeypatch.setattr(EVAL, "ACCESS_RECEIPT_PATH", receipt_path)
    monkeypatch.setattr(EVAL, "RESULT_DIRECTORY", result_directory)
    monkeypatch.setattr(
        EVAL,
        "load_source_bundle",
        lambda role, role_binding: bundle(role, {}),
    )
    monkeypatch.setattr(EVAL, "_assert_source_invariants", lambda *args: None)
    monkeypatch.setattr(EVAL, "_validate_nadig_mapping", lambda *args: None)
    monkeypatch.setattr(EVAL, "_validate_h5ad_mapping", lambda *args: None)
    monkeypatch.setattr(
        EVAL,
        "_verify_truth_binding_after_receipt",
        lambda *args: pytest.fail("authentication touched truth verification"),
    )

    authenticated = EVAL.authenticate_pretruth()
    assert set(authenticated["bundles"]) == set(EVAL.ALL_ROLES)


def test_receipt_is_fsynced_read_only_and_precedes_truth_callback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path
    here = root / "experiment"
    work = here / "work"
    work.mkdir(parents=True)
    receipt = work / "receipt.json"
    result = work / "result"
    monkeypatch.setattr(EVAL, "ROOT", root)
    monkeypatch.setattr(EVAL, "HERE", here)
    monkeypatch.setattr(EVAL, "AUTHORIZATION_PATH", here / "authorization.json")
    monkeypatch.setattr(EVAL, "PROTOCOL_PATH", here / "protocol.json")
    monkeypatch.setattr(EVAL, "ENTRYPOINT_PATH", here / "evaluator.py")
    monkeypatch.setattr(EVAL, "ACCESS_RECEIPT_PATH", receipt)
    monkeypatch.setattr(EVAL, "RESULT_DIRECTORY", result)
    EVAL.AUTHORIZATION_PATH.write_text("{}")

    fsync_calls: list[int] = []
    original_fsync = os.fsync

    def tracked_fsync(descriptor: int) -> None:
        fsync_calls.append(descriptor)
        original_fsync(descriptor)

    monkeypatch.setattr(EVAL.os, "fsync", tracked_fsync)
    authenticated = {
        "authorization_sha256": "b" * 64,
        "authorization": {
            "protocol": {},
            "evaluator": {},
            "source_candidates": {"top_manifest": {}},
            "truth_bindings": {},
        },
    }

    with EVAL.exclusive_access_receipt(authenticated) as lease:
        assert receipt.is_file()
        assert stat.S_IMODE(receipt.stat().st_mode) == 0o444
        assert len(fsync_calls) >= 2

        def truth_callback(binding: object, label: str) -> Path:
            assert receipt.is_file()
            assert len(fsync_calls) >= 2
            return root / "truth"

        monkeypatch.setattr(EVAL, "_verify_source_binding", truth_callback)
        assert (
            EVAL._verify_truth_binding_after_receipt({}, "synthetic truth", lease)
            == root / "truth"
        )

    with pytest.raises(EVAL.EvaluationError, match="already claimed"):
        with EVAL.exclusive_access_receipt(authenticated):
            pass


def test_cli_exposes_no_caller_supplied_paths() -> None:
    parsed = EVAL.parser().parse_args(["run", "--confirm-open-held-truth"])
    assert vars(parsed) == {"command": "run", "confirm_open_held_truth": True}
    with pytest.raises(SystemExit):
        EVAL.parser().parse_args(["run", "--truth", "/tmp/not-authorized"])


def test_result_manifest_is_hash_bound_read_only_and_never_overwritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path
    here = root / "experiment"
    work = here / "work"
    work.mkdir(parents=True)
    destination = work / "held-evaluation-v1"
    monkeypatch.setattr(EVAL, "ROOT", root)
    monkeypatch.setattr(EVAL, "HERE", here)
    monkeypatch.setattr(EVAL, "RESULT_DIRECTORY", destination)

    manifest = EVAL.write_immutable_results(
        {"status": "synthetic"},
        {"values": np.asarray([1.0, 2.0])},
        [],
    )
    manifest_path = destination / EVAL.RESULT_MANIFEST_NAME
    stored = json.loads(manifest_path.read_text())
    unsigned = copy.deepcopy(stored)
    binding_hash = unsigned.pop("binding_sha256")
    assert binding_hash == EVAL.sha256_json_payload(unsigned)
    assert manifest == stored
    assert stat.S_IMODE(destination.stat().st_mode) == 0o555
    assert stat.S_IMODE(manifest_path.stat().st_mode) == 0o444
    for output in stored["outputs"].values():
        path = root / output["path"]
        assert path.stat().st_size == output["bytes"]
        assert EVAL.sha256_file(path) == output["sha256"]

    with pytest.raises(EVAL.EvaluationError, match="already exists"):
        EVAL.write_immutable_results({}, {}, [])


def test_symlinked_source_component_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "root"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (root / "linked").symlink_to(outside, target_is_directory=True)
    monkeypatch.setattr(EVAL, "ROOT", root)
    with pytest.raises(EVAL.EvaluationError, match="symlink"):
        EVAL._reject_symlink_components(root / "linked/file", "synthetic")


def test_protocol_and_authorization_freeze_required_design() -> None:
    protocol = json.loads((HERE / "HELD_EVALUATION_PROTOCOL.v1.json").read_text())
    authorization = json.loads(
        (HERE / "HELD_EVALUATION_AUTHORIZATION.v1.json").read_text()
    )
    assert protocol["source_candidate_contract"]["primary_arm"] == "pan_novel_f100"
    assert protocol["role_policy"]["h1"]["promotion_authorizing"] is False
    assert protocol["held_outcome_exposure"][
        "h1_scored_and_persisted_arms"
    ] == list(EVAL.H1_HELD_ARMS)
    assert protocol["held_outcome_exposure"][
        "authorizing_surface_scored_and_persisted_arms"
    ] == list(EVAL.AUTHORIZING_HELD_ARMS)
    assert protocol["promotion_gate"][
        "hek293t_and_all_six_nadig_galleries_are_mandatory"
    ]
    assert protocol["null_and_uncertainty"]["target_label_permutation"][
        "replicates"
    ] == 99
    assert protocol["null_and_uncertainty"]["target_cluster_bootstrap"][
        "replicates"
    ] == 20_000
    assert authorization["execution_contract"]["caller_supplied_paths"] is False
    assert authorization["truth_bindings"]["h1"]["scored_arms"] == list(
        EVAL.H1_HELD_ARMS
    )
    assert authorization["truth_bindings"]["orion_hek293t"][
        "scored_arms"
    ] == list(EVAL.AUTHORIZING_HELD_ARMS)
    assert authorization["truth_bindings"]["nadig_jurkat"]["scored_arms"] == list(
        EVAL.AUTHORIZING_HELD_ARMS
    )
