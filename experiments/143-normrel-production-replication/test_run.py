from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("experiment143_run", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


class Experiment143Tests(unittest.TestCase):
    def test_exact_normrel_formula_and_group_floor(self) -> None:
        keys = np.asarray([[1.0, 2.0], [4.0, 8.0]])
        available = np.ones_like(keys, dtype=bool)
        weights, floor = RUN.normrel_weights(keys, available)
        self.assertAlmostEqual(floor, 1.3)
        raw = np.maximum(keys, floor) ** -1.5
        expected = raw / raw.sum(axis=1, keepdims=True)
        np.testing.assert_allclose(weights, expected, rtol=0.0, atol=1e-15)

    def test_floor_uses_positive_available_values_only(self) -> None:
        keys = np.asarray([[0.0, 1.0, 1000.0], [2.0, 4.0, 2000.0]])
        available = np.asarray([[True, True, False], [True, True, False]])
        weights, floor = RUN.normrel_weights(keys, available)
        self.assertAlmostEqual(floor, 1.2)
        np.testing.assert_array_equal(weights[:, 2], 0.0)
        np.testing.assert_allclose(weights.sum(axis=1), 1.0)

    def test_one_donor_target_is_identity(self) -> None:
        keys = np.asarray([[2.0, 0.0], [0.0, 4.0]])
        available = np.asarray([[True, False], [False, True]])
        weights, _ = RUN.normrel_weights(keys, available)
        np.testing.assert_array_equal(weights, available.astype(np.float64))
        slots = np.arange(12, dtype=np.float64).reshape(2, 2, 3)
        expected = np.vstack([slots[0, 0], slots[1, 1]])
        np.testing.assert_array_equal(
            RUN.mix_slots_exact(slots, weights, available), expected
        )

    def test_k562_composite_key_uses_original_k562_response(self) -> None:
        centered = {
            "k562_gwps": np.asarray([[3.0, 4.0], [6.0, 8.0], [12.0, 16.0]]),
            "orion_hct116": np.asarray([[0.0, 10.0], [0.0, 20.0], [0.0, 40.0]]),
            "orion_hek293t": np.asarray([[8.0, 15.0], [16.0, 30.0], [32.0, 60.0]]),
        }
        with (
            patch.object(RUN, "RESPONSE_COUNT", 2),
            patch.object(RUN, "DEV_COUNT", len(centered["k562_gwps"])),
        ):
            keys, available, denominators = RUN.group_normrel_keys("non_k562", centered)
        np.testing.assert_allclose(keys[:, 0], [0.5, 1.0, 2.0])
        np.testing.assert_allclose(keys[:, 1], [0.5, 1.0, 2.0])
        np.testing.assert_allclose(keys[:, 2], [0.5, 1.0, 2.0])
        self.assertEqual(denominators["k562_composite"], 10.0)
        np.testing.assert_array_equal(available, True)

    def test_equal_weight_path_returns_frozen_baseline_bit_exactly(self) -> None:
        slots = np.arange(60, dtype=np.float32).reshape(5, 3, 4)
        available = np.ones((5, 3), dtype=bool)
        weights = RUN.equal_weights(available)
        baseline = np.mean(slots.astype(np.float64), axis=1).astype(np.float32)
        actual = RUN.mix_slots_exact(
            slots,
            weights,
            available,
            exact_equal_baseline=baseline,
        )
        self.assertTrue(np.array_equal(actual, baseline))
        self.assertEqual(actual.dtype, baseline.dtype)

    def test_nonidentity_weights_do_not_take_identity_shortcut(self) -> None:
        slots = np.asarray([[[1.0], [3.0]]])
        available = np.ones((1, 2), dtype=bool)
        weights = np.asarray([[0.75, 0.25]])
        frozen = np.asarray([[99.0]])
        actual = RUN.mix_slots_exact(
            slots,
            weights,
            available,
            exact_equal_baseline=frozen,
        )
        np.testing.assert_allclose(actual, [[1.5]])

    def test_hash_bound_dynamic_import_rejects_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            module_path = Path(directory) / "upstream.py"
            module_path.write_text("VALUE = 1\n")
            with (
                patch.object(RUN, "EXP142_RUN", module_path),
                patch.object(
                    RUN,
                    "frozen_bindings",
                    return_value={module_path.resolve(): "0" * 64},
                ),
                self.assertRaises(RUN.ContractError),
            ):
                RUN.load_experiment142()

    def test_authorization_hash_validation_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bound.txt"
            path.write_text("frozen")
            authorization = {"bindings": {str(path): RUN.sha256_file(path)}}
            RUN.validate_bound_paths(authorization, (path,), hash_files=True)
            path.write_text("changed")
            with self.assertRaises(RUN.ContractError):
                RUN.validate_bound_paths(authorization, (path,), hash_files=True)

    def test_failed_lock_acquisition_preserves_other_owner(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            lock = Path(directory) / "owned.lock"
            lock.write_text("other-process\n")
            with self.assertRaises(FileExistsError):
                with RUN.exclusive_lock(lock):
                    self.fail("existing lock must not be acquired")
            self.assertEqual(lock.read_text(), "other-process\n")

    def test_access_lock_precedes_truth_validation_and_loading(self) -> None:
        source = MODULE_PATH.read_text()
        body = source[source.index("def develop(") : source.index("def audit(")]
        lock_position = body.index("create_outcome_access_lock(authorization)")
        verification_position = body.index(
            "validate_development_authorization(verify_truth_files=True)"
        )
        load_position = body.index("truth = parent.load_development_truth(0.0)")
        self.assertLess(lock_position, verification_position)
        self.assertLess(verification_position, load_position)

    def test_prepare_and_audit_do_not_load_truth(self) -> None:
        source = MODULE_PATH.read_text()
        prepare_body = source[
            source.index("def prepare()") : source.index(
                "def validate_preoutcome_control_plane"
            )
        ]
        audit_body = source[source.index("def audit()") : source.index("def main()")]
        self.assertNotIn("load_development_truth", prepare_body)
        self.assertNotIn("load_development_truth", audit_body)
        self.assertNotIn("test_targets.txt", source)
        self.assertNotIn("source-test-access", source)
        self.assertNotIn("load_test_truth", source)

    def test_cli_has_only_prepare_develop_and_audit(self) -> None:
        source = MODULE_PATH.read_text()
        main_body = source[source.index("def main()") :]
        self.assertIn('choices=("prepare", "develop", "audit")', main_body)
        for forbidden in ("evaluate-test", "generate-counts", "build-challenge"):
            self.assertNotIn(forbidden, main_body)

    def test_raw_gate_set_is_complete_and_fail_closed(self) -> None:
        summary = {
            "macro": {"pds": 0.003, "signed_credit": 0.003},
            "families": {
                metric: {family: 0.001 for family in RUN.FAMILIES}
                for metric in ("pds", "signed_credit")
            },
        }
        cells = {
            metric: {family: [0.001] * 5 for family in RUN.FAMILIES}
            for metric in ("pds", "signed_credit")
        }
        bootstrap = {
            metric: {"lower_95": 0.0001, "upper_95": 0.01}
            for metric in ("pds", "signed_credit")
        }
        concentration = {
            metric: {"family": 0.5, "target_block": 0.5}
            for metric in ("pds", "signed_credit")
        }
        checks = RUN.promotion_checks(summary, cells, bootstrap, concentration)
        self.assertEqual(len(checks), 12)
        self.assertTrue(all(checks.values()))
        cells["pds"]["jurkat"][3] = -1e-3
        checks = RUN.promotion_checks(summary, cells, bootstrap, concentration)
        self.assertFalse(checks["pds_every_family_fold_cell_nonnegative"])

    def test_protocol_freezes_one_candidate_and_prohibitions(self) -> None:
        value = json.loads(Path(RUN.PROTOCOL_PATH).read_text())
        self.assertEqual(value["candidate"]["count"], 1)
        self.assertEqual(value["candidate"]["name"], "normrel613^-1.5")
        self.assertFalse(value["evaluation"]["adaptive_selection"])
        self.assertIn("positive_scalar_invariance", value["geometry"])
        self.assertTrue(all(value["prohibited"].values()))


if __name__ == "__main__":
    unittest.main()
