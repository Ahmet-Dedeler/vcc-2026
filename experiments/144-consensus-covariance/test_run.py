#!/usr/bin/env python3
"""Synthetic and static contract tests for Experiment 144 source-only v2."""

from __future__ import annotations

import contextlib
import importlib.util
import inspect
import io
import json
import os
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_PATH = HERE / "run.py"
MEMORY_GUARD_PATH = HERE / "memory_guard.py"


def load_run():
    spec = importlib.util.spec_from_file_location("experiment144_v2_tested", RUN_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


M = load_run()


def load_memory_guard():
    spec = importlib.util.spec_from_file_location(
        "experiment144_external_memory_guard_tested",
        MEMORY_GUARD_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


G = load_memory_guard()


class RecordingGuard:
    def __init__(self) -> None:
        self.labels: list[str] = []

    def check(self, label: str) -> None:
        self.labels.append(label)


def write_synthetic_donor(
    root: Path,
    targets: tuple[str, ...],
    genes: tuple[str, ...],
    values: np.ndarray,
) -> None:
    root.mkdir(parents=True)
    (root / "perts.txt").write_text("\n".join(targets) + "\n")
    (root / "genes.txt").write_text("\n".join(genes) + "\n")
    np.save(root / "delta.npy", np.asarray(values, dtype=np.float32))


class Experiment144V2Tests(unittest.TestCase):
    def test_candidate_family_is_exact_and_frozen(self) -> None:
        self.assertEqual(
            M.candidate_configs(),
            (
                {"name": "identity", "gamma": 0.0, "eta": 0.0},
                {
                    "name": "consensus:g=0.65:eta=0.25",
                    "gamma": 0.65,
                    "eta": 0.25,
                },
            ),
        )
        self.assertEqual(M.candidate_spec()["candidate_count"], 2)

    def test_identity_blend_is_bit_exact(self) -> None:
        rng = np.random.default_rng(1)
        baseline = rng.normal(size=(8, 12)).astype(np.float32)
        candidate = rng.normal(size=(8, 12)).astype(np.float32)
        output = M.row_norm_preserving_blend(baseline, candidate, 0.0)
        self.assertTrue(np.array_equal(output, baseline))

    def test_candidate_blend_preserves_row_norms(self) -> None:
        rng = np.random.default_rng(2)
        baseline = rng.normal(size=(11, 17)).astype(np.float32)
        candidate = rng.normal(size=(11, 17)).astype(np.float32)
        output = M.row_norm_preserving_blend(baseline, candidate, M.BLEND_ETA)
        np.testing.assert_allclose(
            np.linalg.norm(output, axis=1),
            np.linalg.norm(baseline, axis=1),
            rtol=2e-6,
            atol=2e-6,
        )

    def test_excluded_row_cannot_affect_any_fitted_basis_statistic(self) -> None:
        targets = ("t1", "t2", "t3", "t4", "t5", "held")
        genes = ("g1", "g2", "g3")
        donor_one = np.asarray(
            [
                [1.0, 0.0, 1.0],
                [0.0, 1.0, -1.0],
                [-1.0, 0.0, 0.5],
                [0.0, -1.0, -0.5],
                [1.0, 1.0, 0.25],
                [1e6, -1e6, 5e5],
            ]
        )
        donor_two = np.asarray(
            [
                [2.0, 0.0, -1.0],
                [0.0, 2.0, 1.0],
                [-2.0, 0.0, -0.25],
                [0.0, -2.0, 0.25],
                [2.0, 2.0, -0.5],
                [-7e5, 9e5, -4e5],
            ]
        )

        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            roots_a = {"d1": base / "a-d1", "d2": base / "a-d2"}
            roots_b = {"d1": base / "b-d1", "d2": base / "b-d2"}
            write_synthetic_donor(roots_a["d1"], targets, genes, donor_one)
            write_synthetic_donor(roots_a["d2"], targets, genes, donor_two)

            mutated_one = donor_one.copy()
            mutated_two = donor_two.copy()
            mutated_one[-1] = [-3e12, 4e12, 8e12]
            mutated_two[-1] = [9e12, -2e12, -6e12]
            write_synthetic_donor(roots_b["d1"], targets, genes, mutated_one)
            write_synthetic_donor(roots_b["d2"], targets, genes, mutated_two)

            guard_a = RecordingGuard()
            guard_b = RecordingGuard()
            basis_a, report_a = M.build_consensus_basis(
                ("d1", "d2"),
                roots_a,
                genes,
                {"held"},
                guard_a,
                chunk=2,
            )
            basis_b, report_b = M.build_consensus_basis(
                ("d1", "d2"),
                roots_b,
                genes,
                {"held"},
                guard_b,
                chunk=2,
            )

        np.testing.assert_array_equal(basis_a, basis_b)
        np.testing.assert_array_equal(
            basis_a.T @ basis_a,
            basis_b.T @ basis_b,
        )
        self.assertEqual(report_a["targets"], report_b["targets"])
        self.assertNotIn("held", report_a["targets"])
        self.assertEqual(
            report_a["consensus_mean_sha256"],
            report_b["consensus_mean_sha256"],
        )
        self.assertEqual(
            report_a["role_reference_median_norm"],
            report_b["role_reference_median_norm"],
        )
        self.assertEqual(
            report_a["donor_statistics"],
            report_b["donor_statistics"],
        )
        self.assertGreater(len(guard_a.labels), 0)
        self.assertGreater(len(guard_b.labels), 0)

    def test_numerical_rank_rejects_rows_greater_than_columns_singularity(self) -> None:
        with self.assertRaises(M.ContractError):
            M.regularize_spectrum(
                np.asarray([0.0, 1.0, 2.0]),
                n_rows=100,
                n_columns=3,
            )

    def test_numerical_rank_accepts_strictly_positive_spectrum(self) -> None:
        regularized, report = M.regularize_spectrum(
            np.asarray([1e-3, 1.0, 2.0]),
            n_rows=100,
            n_columns=3,
        )
        self.assertEqual(report["numerical_rank"], 3)
        self.assertTrue(report["regularized_positive_spectrum"])
        self.assertTrue(np.all(regularized > 0.0))

    def test_rss_and_free_memory_limits_cannot_be_bypassed(self) -> None:
        invalid = (
            (0.0, 25.0),
            (-1.0, 25.0),
            (18.1, 25.0),
            (18.0, 24.9),
            (18.0, 100.0),
        )
        for max_rss, min_free in invalid:
            with self.subTest(max_rss=max_rss, min_free=min_free):
                with self.assertRaises(M.ContractError):
                    M.MemoryGuard(max_rss, min_free)

        for function in (
            M.sha256_file,
            M.fit_allowed_donor_statistics,
            M.build_consensus_basis,
            M.fit_consensus_operator,
            M.atomic_json,
            M.atomic_npz,
            M.atomic_stage_directory,
        ):
            parameter = inspect.signature(function).parameters["guard"]
            self.assertIs(parameter.default, inspect.Parameter.empty)

        helpers = M.memory_guard_helpers()
        self.assertTrue(callable(helpers["process_rss"]))
        self.assertTrue(callable(helpers["free_percent"]))
        self.assertIn(M.MEMORY_GUARD_SCRIPT, M.source_paths())
        self.assertEqual(M.MEMORY_GUARD_SCRIPT, MEMORY_GUARD_PATH)

    def test_external_guard_is_mandatory_and_cannot_be_looser(self) -> None:
        keys = (
            M.EXTERNAL_GUARD_ACTIVE_ENV,
            M.EXTERNAL_GUARD_PATH_ENV,
            M.EXTERNAL_GUARD_MAX_RSS_ENV,
            M.EXTERNAL_GUARD_MIN_FREE_ENV,
        )
        with mock.patch.dict(os.environ, {}, clear=False):
            for key in keys:
                os.environ.pop(key, None)
            with self.assertRaises(M.ContractError):
                M.validate_external_guard(18.0, 25.0)

        token = "a" * 64
        request_read, request_write = os.pipe()
        response_read, response_write = os.pipe()
        valid = {
            M.EXTERNAL_GUARD_ACTIVE_ENV: "1",
            M.EXTERNAL_GUARD_PATH_ENV: str(MEMORY_GUARD_PATH.resolve()),
            M.EXTERNAL_GUARD_MAX_RSS_ENV: "18.0",
            M.EXTERNAL_GUARD_MIN_FREE_ENV: "25.0",
            M.EXTERNAL_GUARD_REQUEST_FD_ENV: str(request_write),
            M.EXTERNAL_GUARD_RESPONSE_FD_ENV: str(response_read),
            M.EXTERNAL_GUARD_TOKEN_ENV: token,
        }
        supervisor_errors: list[BaseException] = []

        def supervisor() -> None:
            try:
                request = os.read(request_read, len(token) + 1)
                self.assertEqual(request, f"{token}\n".encode())
                os.write(response_write, request)
            except BaseException as error:
                supervisor_errors.append(error)

        thread = threading.Thread(target=supervisor)
        thread.start()
        try:
            with mock.patch.dict(os.environ, valid, clear=False):
                report = M.validate_external_guard(18.0, 25.0)
        finally:
            thread.join(timeout=5)
            os.close(request_read)
            os.close(response_write)
        self.assertFalse(supervisor_errors)
        self.assertFalse(thread.is_alive())
        self.assertEqual(report["max_rss_gib"], 18.0)
        self.assertEqual(report["min_free_percent"], 25.0)
        self.assertEqual(
            report["live_handshake_verified"],
            "request-response-pipe",
        )

        for key, value in (
            (M.EXTERNAL_GUARD_MAX_RSS_ENV, "18.1"),
            (M.EXTERNAL_GUARD_MIN_FREE_ENV, "24.9"),
        ):
            invalid = dict(valid)
            invalid[key] = value
            with self.subTest(key=key):
                with mock.patch.dict(os.environ, invalid, clear=False):
                    with self.assertRaises(M.ContractError):
                        M.validate_external_guard(18.0, 25.0)

    def test_external_monitor_failures_terminate_the_process_group(self) -> None:
        class FakeProcess:
            pid = 424242

            @staticmethod
            def poll():
                return None

        process = FakeProcess()
        with (
            mock.patch.object(
                G,
                "process_group_members",
                side_effect=RuntimeError("ps failed"),
            ),
            mock.patch.object(G, "terminate_group") as terminate,
        ):
            status = G.monitor_process(
                process,
                max_rss_gib=18.0,
                min_free_percent=25.0,
                poll_seconds=0.001,
                report_file=None,
                quiet=True,
            )
        self.assertEqual(status, 2)
        terminate.assert_called_once()

        with (
            mock.patch.object(
                G,
                "process_group_members",
                return_value=[(424242, 1)],
            ),
            mock.patch.object(G, "free_percent", return_value=None),
            mock.patch.object(G, "terminate_group") as terminate,
        ):
            status = G.monitor_process(
                process,
                max_rss_gib=18.0,
                min_free_percent=25.0,
                poll_seconds=0.001,
                report_file=None,
                quiet=True,
            )
        self.assertEqual(status, 2)
        terminate.assert_called_once()

    def test_termination_escalates_and_verifies_group_extinction(self) -> None:
        class FakeProcess:
            pid = 424242

        process = FakeProcess()
        with (
            mock.patch.object(G.os, "killpg") as killpg,
            mock.patch.object(
                G,
                "process_group_members",
                side_effect=AssertionError("termination must not depend on ps"),
            ),
            mock.patch.object(
                G,
                "wait_for_group_exit",
                side_effect=(False, True),
            ) as wait,
        ):
            G.terminate_group(process, "test")
        self.assertEqual(
            [call.args[1] for call in killpg.call_args_list],
            [G.signal.SIGTERM, G.signal.SIGKILL],
        )
        self.assertEqual(wait.call_count, 2)

    def test_serializers_emit_memory_checkpoints(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            guard = RecordingGuard()
            M.atomic_json(root / "record.json", {"ok": True}, guard)
            M.atomic_npz(
                root / "arrays.npz",
                guard,
                value=np.asarray([1.0], dtype=np.float32),
            )
            with M.atomic_stage_directory(root / "final", guard) as staging:
                (staging / "item.txt").write_text("ok\n")
        self.assertTrue(
            any(label.startswith("before JSON serialization") for label in guard.labels)
        )
        self.assertTrue(
            any(label.startswith("after NPZ serialization") for label in guard.labels)
        )
        self.assertIn("before finalizing final", guard.labels)

    def test_authorization_identity_is_strict(self) -> None:
        authorization = {
            "schema_version": 2,
            "experiment": HERE.name,
            "stage": "source_only_prepare_v2",
            "confirmation": M.AUTHORIZATION_CONFIRMATION,
            "source_only_prepare_v2_authorized": True,
            "independent_no_truth_audit_passed": True,
            "development_outcomes_authorized": False,
            "held_outcomes_authorized": False,
            "sealed_613_410_values_authorized": False,
            "kolf_values_authorized": False,
            "reserved_expression_authorized": False,
            "challenge_controls_authorized": False,
            "cell_eval2_authorized": False,
            "counts_authorized": False,
            "build_authorized": False,
            "submission_authorized": False,
            "commit_authorized": False,
            "push_authorized": False,
        }
        M.validate_prepare_authorization_document(authorization)
        for key, value in (
            ("schema_version", 1),
            ("experiment", "other"),
            ("stage", "other"),
            ("confirmation", "other"),
        ):
            mutated = dict(authorization)
            mutated[key] = value
            with self.subTest(key=key):
                with self.assertRaises(M.ContractError):
                    M.validate_prepare_authorization_document(mutated)

    def test_v1_paths_are_quarantined_and_v2_writes_are_separate(self) -> None:
        protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
        frozen = protocol["v1_quarantine"]["sha256"]
        self.assertEqual(len(frozen), 11)
        self.assertTrue(
            all(
                relative == "SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json"
                or relative.startswith("work/preoutcome-v1/")
                for relative in frozen
            )
        )
        self.assertEqual(M.OUTPUT_DIR, HERE / "work/source-only-v2")
        self.assertNotEqual(M.OUTPUT_DIR, M.V1_OUTPUT_DIR)
        with tempfile.TemporaryDirectory() as temporary:
            existing = Path(temporary) / "already-exists"
            existing.mkdir()
            with self.assertRaises(M.ContractError):
                with M.atomic_stage_directory(existing, RecordingGuard()):
                    pass
        report = M.verify_v1_quarantine(RecordingGuard())
        self.assertTrue(report["directory_membership_verified"])

    def test_protocol_matches_code_and_remains_unauthorized(self) -> None:
        protocol = M.protocol()
        self.assertEqual(
            protocol["candidate_family"]["candidates"],
            list(M.CANDIDATES),
        )
        self.assertEqual(
            protocol["candidate_family"]["candidate_count_including_identity"], 2
        )
        self.assertFalse(protocol["candidate_family"]["new_grid_or_tuning"])
        self.assertEqual(protocol["status"], "SOURCE_ONLY_V2_UNAUTHORIZED")
        self.assertTrue(
            protocol["scope"][
                "independent_no_truth_audit_required_before_any_evaluation"
            ]
        )

    def test_forbidden_commands_and_imports_are_absent(self) -> None:
        report = M.forbidden_surface_audit()
        self.assertEqual(report["commands"], ["prepare", "audit"])
        self.assertEqual(report["forbidden_imports_found"], [])
        self.assertEqual(report["forbidden_functions_found"], [])
        for command in ("develop", "evaluate", "score", "counts", "build", "submit"):
            with self.subTest(command=command):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        M.build_parser().parse_args([command])

    def test_source_allowlist_contains_no_forbidden_value_path(self) -> None:
        M.ensure_source_allowlist()
        normalized = [str(path).lower() for path in M.source_paths()]
        self.assertFalse(any(path.endswith(".h5ad") for path in normalized))
        self.assertFalse(any("development-preoutcome" in path for path in normalized))
        self.assertFalse(any("development-result" in path for path in normalized))
        self.assertFalse(any("kolf" in path for path in normalized))
        self.assertFalse(any("/submissions/" in path for path in normalized))


if __name__ == "__main__":
    unittest.main(verbosity=2)
