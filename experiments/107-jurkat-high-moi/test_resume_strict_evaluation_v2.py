from __future__ import annotations

import inspect
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

import resume_strict_evaluation_v1 as V1
import resume_strict_evaluation_v2 as V2
import run_strict_evaluation as ORIGINAL
import strict_nadig_validation as STRICT


class EmptyDirectionalSetTests(unittest.TestCase):
    @staticmethod
    def control(genes: int) -> np.ndarray:
        proportions = np.full(genes, 1.0 / genes, dtype=np.float64)
        return np.log1p(ORIGINAL.TS * proportions)

    def test_empty_sets_are_zero_without_target_exclusion_or_reweighting(self) -> None:
        reference = ORIGINAL.ReferenceDE(
            genes=np.asarray(["PCCB", "A", "B"]),
            gate=np.asarray([True, True, True]),
            gate_indices=np.asarray([0, 1, 2]),
            control_mean_cpm=np.asarray([100.0, 100.0, 100.0]),
            targets=np.asarray(["PCCB", "OTHER", "NONEMPTY"]),
            lfc=np.asarray(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 1.0, -2.0],
                ]
            ),
            significant=(
                np.asarray([0]),
                np.asarray([], dtype=np.int64),
                np.asarray([1, 2]),
            ),
        )
        incumbent = np.zeros((3, 3), dtype=np.float64)
        candidate = np.asarray(
            [
                [0.0, 9.0, -9.0],
                [4.0, -4.0, 2.0],
                [0.0, 0.2, -0.1],
            ]
        )
        incumbent_values, incumbent_empty = (
            V2._directional_credit_with_empty_mask(  # noqa: SLF001
                incumbent,
                self.control(3),
                ("PCCB", "OTHER", "NONEMPTY"),
                reference,
            )
        )
        candidate_values, candidate_empty = (
            V2._directional_credit_with_empty_mask(  # noqa: SLF001
                candidate,
                self.control(3),
                ("PCCB", "OTHER", "NONEMPTY"),
                reference,
            )
        )
        np.testing.assert_array_equal(incumbent_empty, [True, True, False])
        np.testing.assert_array_equal(candidate_empty, incumbent_empty)
        np.testing.assert_array_equal(incumbent_values[:2], [0.0, 0.0])
        np.testing.assert_array_equal(candidate_values[:2], [0.0, 0.0])
        self.assertEqual(len(candidate_values), 3)
        self.assertEqual(
            float((candidate_values - incumbent_values)[:2].sum()),
            0.0,
        )

    def test_nonempty_targets_are_bit_identical_to_original(self) -> None:
        reference = ORIGINAL.ReferenceDE(
            genes=np.asarray(["A", "B", "C", "D"]),
            gate=np.asarray([True, True, True, True]),
            gate_indices=np.asarray([0, 1, 2, 3]),
            control_mean_cpm=np.asarray([20.0, 30.0, 40.0, 50.0]),
            targets=np.asarray(["T1", "T2"]),
            lfc=np.asarray(
                [
                    [1.0, -2.0, 0.0, 0.5],
                    [0.0, 3.0, -1.0, 0.0],
                ]
            ),
            significant=(np.asarray([0, 1, 3]), np.asarray([1, 2])),
        )
        effect = np.asarray(
            [[0.1, -0.2, 0.3, -0.1], [-0.1, 0.4, -0.2, 0.2]]
        )
        control = self.control(4)
        expected = ORIGINAL.directional_credit_per_target(
            effect, control, ("T1", "T2"), reference
        )
        observed, empty = V2._directional_credit_with_empty_mask(  # noqa: SLF001
            effect, control, ("T1", "T2"), reference
        )
        np.testing.assert_array_equal(observed, expected)
        np.testing.assert_array_equal(empty, [False, False])

    def test_implementation_is_generic_not_pccb_specific(self) -> None:
        source = inspect.getsource(V2._directional_credit_with_empty_mask)  # noqa: SLF001
        self.assertNotIn("PCCB", source)
        self.assertNotIn("EXPECTED_ERROR", source)

    def test_held_evaluator_retains_empty_target_and_uses_v1_pds(self) -> None:
        genes = np.asarray(["PCCB", "A", "B"])
        targets = np.asarray(["PCCB", "T2"])
        incumbent = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        candidate = np.asarray([[0.0, 0.8, 0.2], [0.0, 0.1, 0.9]])
        arrays = {
            "targets": targets,
            "genes": genes,
            "score_genes": np.asarray(["A", "B"]),
            "score_indices": np.asarray([1, 2]),
            "incumbent_effects": incumbent,
            "candidate_f100": candidate,
        }
        held_truth = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        reference = ORIGINAL.ReferenceDE(
            genes=genes,
            gate=np.asarray([True, True, True]),
            gate_indices=np.asarray([0, 1, 2]),
            control_mean_cpm=np.asarray([100.0, 100.0, 100.0]),
            targets=targets,
            lfc=np.asarray([[1.0, 0.0, 0.0], [0.0, 1.0, -1.0]]),
            significant=(np.asarray([0]), np.asarray([1, 2])),
        )
        control = self.control(3)
        expected_pds = V1.pds_on_frozen_score_axis(
            incumbent_full=incumbent,
            candidate_full=candidate,
            held_truth=held_truth,
            official_genes=genes,
            held_genes=genes,
            score_genes=arrays["score_genes"],
            score_indices=arrays["score_indices"],
        )
        with (
            mock.patch.object(ORIGINAL, "_load_integration", return_value=arrays),
            mock.patch.object(ORIGINAL, "_validate_integration_candidates"),
            mock.patch.object(
                ORIGINAL,
                "_load_held_effects",
                return_value=(genes, control, held_truth),
            ),
            mock.patch.object(ORIGINAL, "_load_reference_de", return_value=reference),
        ):
            observed = V2.evaluate_held_surface_v2(
                "synthetic",
                0.1,
                harness_path=Path("/synthetic.h5ad"),
                reference_de_path=Path("/synthetic-ref.npz"),
            )
        np.testing.assert_array_equal(observed["incumbent_pds"], expected_pds[0])
        np.testing.assert_array_equal(observed["candidate_pds"], expected_pds[1])
        np.testing.assert_array_equal(
            observed["empty_directional_credit_targets"], ["PCCB"]
        )
        np.testing.assert_array_equal(observed["targets"], targets)
        self.assertEqual(observed["incumbent_directional_credit"][0], 0.0)
        self.assertEqual(observed["candidate_directional_credit"][0], 0.0)

    def test_missing_reference_target_still_fails_closed(self) -> None:
        reference = ORIGINAL.ReferenceDE(
            genes=np.asarray(["A"]),
            gate=np.asarray([True]),
            gate_indices=np.asarray([0]),
            control_mean_cpm=np.asarray([100.0]),
            targets=np.asarray(["KNOWN"]),
            lfc=np.asarray([[1.0]]),
            significant=(np.asarray([0]),),
        )
        with self.assertRaisesRegex(V2.RecoveryV2Error, "missing target"):
            V2.directional_credit_per_target_v2(
                np.zeros((1, 1)),
                self.control(1),
                ("UNKNOWN",),
                reference,
            )


class V2AuthorizationTests(unittest.TestCase):
    def test_failure_amendment_binds_exact_v1_attempt(self) -> None:
        payload = json.loads(V2.V2_AMENDMENT_PATH.read_text())
        self.assertEqual(payload["failure"]["exact_error"], V2.EXPECTED_ERROR)
        bindings = payload["predecessor_bindings"]
        self.assertEqual(
            bindings["v1_resume_runner"]["sha256"],
            "aa8b0a8ced62bd9731776c01cbe1c410689af6c409dd605017de0f6a5d6f4c2d",
        )
        self.assertEqual(
            bindings["v1_failure_recovery_amendment"]["sha256"],
            "73d0f30b5e9e9fa941624e99fa067d58c90f30777a14688abadd38fd06198c66",
        )
        self.assertEqual(
            bindings["original_outcome_access_receipt"]["sha256"],
            "3ec3e08ab87d0c094457f8a3e76e1c0fac8cf5273ee675de0d9c2e1784a7078b",
        )
        self.assertEqual(
            bindings["v1_failed_recovery_lock"]["sha256"],
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        )

    def test_policy_preserves_every_frozen_degree_of_freedom(self) -> None:
        payload = json.loads(V2.V2_AMENDMENT_PATH.read_text())
        policy = payload["empty_set_policy"]
        self.assertEqual(policy["candidate_and_incumbent_absolute_credit"], 0.0)
        self.assertEqual(policy["candidate_minus_incumbent_delta"], 0.0)
        self.assertIn("retain_every_target", policy["target_handling"])
        self.assertIn("no_reweighting", policy["target_handling"])
        self.assertEqual(
            payload["recovery_contract"]["frozen_components_unchanged"],
            [
                "source_model",
                "integration_candidates",
                "fraction_grid",
                "nested_selection",
                "permutation_null",
                "cluster_bootstrap",
                "promotion_gates",
                "pds_score_axes",
                "exact_symbol_axis_projection",
                "target_order",
                "target_weighting",
            ],
        )

    def test_authentication_failure_precedes_outcome_loader(self) -> None:
        touched = False

        def outcome_loader(*args: object, **kwargs: object) -> np.ndarray:
            nonlocal touched
            touched = True
            raise AssertionError("outcome loader must not run")

        lease = V1.RecoveryLease(
            path=Path("/unused"),
            descriptor=-1,
            device=-1,
            inode=-1,
        )
        with mock.patch.object(
            V2,
            "authenticate_v2_amendment",
            side_effect=V2.RecoveryV2Error("bad v2 amendment"),
        ):
            with self.assertRaisesRegex(V2.RecoveryV2Error, "bad v2 amendment"):
                V2.authorize_then_reload_nadig_v2(
                    lease=lease,
                    outcome_loader=outcome_loader,
                )
        self.assertFalse(touched)

    def test_original_receipt_is_verified_before_outcome_reload(self) -> None:
        events: list[str] = []
        plan = {
            "nadig_target_count": 2,
            "nadig_gene_count": 3,
            "nadig_artifacts": {"binding": "frozen"},
        }
        receipt = {
            "source_manifest_sha256": "manifest",
            "source_model_sha256": "model",
            "nadig_artifacts": {"binding": "frozen"},
        }

        def authorize(*args: object, **kwargs: object) -> tuple[object, object]:
            events.append("authorize")
            return (
                SimpleNamespace(
                    manifest_sha256="manifest",
                    model_sha256="model",
                ),
                plan,
            )

        def verify(*args: object, **kwargs: object) -> dict[str, str]:
            events.append("verify-receipt")
            return {"binding": "frozen"}

        def load(*args: object, **kwargs: object) -> np.ndarray:
            events.append("load-outcome")
            return np.zeros((2, 3), dtype=np.float64)

        with (
            tempfile.TemporaryDirectory() as temporary,
            mock.patch.object(STRICT, "authorize_outcome_phase", side_effect=authorize),
            mock.patch.object(STRICT, "verify_nadig_artifacts", side_effect=verify),
        ):
            lock = Path(temporary) / "v2.lock"
            with V2._one_shot_v2_lock(lock) as lease:  # noqa: SLF001
                authenticated = {
                    "v2_lock_path": lock.resolve(),
                    "v1_auth": {
                        "original_auth": {"plan": plan},
                        "receipt": receipt,
                    },
                }
                V2._reload_nadig_after_v2_authentication(  # noqa: SLF001
                    authenticated,
                    lease=lease,
                    outcome_loader=load,
                )
        self.assertEqual(events, ["authorize", "verify-receipt", "load-outcome"])


class V2OutputInvariantTests(unittest.TestCase):
    def test_v2_lock_is_permanent_and_one_shot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            lock = Path(temporary) / ".strict-evaluation-recovery-v2.lock"
            with V2._one_shot_v2_lock(lock) as lease:  # noqa: SLF001
                self.assertEqual(lease.path, lock.resolve())
                self.assertEqual(lock.stat().st_size, 0)
            self.assertTrue(lock.is_file())
            with self.assertRaisesRegex(V2.RecoveryV2Error, "already claimed"):
                with V2._one_shot_v2_lock(lock):  # noqa: SLF001
                    pass

    def test_v1_lock_and_original_receipt_remain_immutable(self) -> None:
        payload = json.loads(V2.V2_AMENDMENT_PATH.read_text())
        bindings = payload["predecessor_bindings"]
        for name in (
            "v1_failed_recovery_lock",
            "original_outcome_access_receipt",
        ):
            record = bindings[name]
            path = V2.ROOT / record["path"]
            self.assertEqual(path.stat().st_size, record["bytes"])
            self.assertEqual(STRICT.sha256_file(path), record["sha256"])
            self.assertEqual(path.stat().st_mode & 0o222, 0)

    def test_v2_output_is_strictly_versioned(self) -> None:
        payload = json.loads(V2.V2_AMENDMENT_PATH.read_text())
        output = payload["recovery_contract"]["output_binding"]
        self.assertEqual(
            output["directory"],
            "experiments/107-jurkat-high-moi/work/strict-evaluation-results-v2",
        )
        self.assertEqual(V2.V2_RESULT_DIRECTORY.name, "strict-evaluation-results-v2")
        self.assertNotEqual(V2.V2_RESULT_DIRECTORY, ORIGINAL.RESULT_DIRECTORY)
        self.assertNotEqual(V2.V2_LOCK_PATH, V1.RECOVERY_LOCK_PATH)
        self.assertEqual(output["arrays"], ORIGINAL.RESULT_ARRAYS_NAME)
        self.assertEqual(output["summary"], ORIGINAL.RESULT_SUMMARY_NAME)
        self.assertEqual(output["manifest"], ORIGINAL.RESULT_MANIFEST_NAME)
        self.assertIn(
            "ORIGINAL.write_immutable_results",
            inspect.getsource(V2.run_v2_recovery_evaluation),
        )


if __name__ == "__main__":
    unittest.main()
