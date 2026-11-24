from __future__ import annotations

import copy
import importlib.util
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SPEC = importlib.util.spec_from_file_location(
    "experiment_119_source_audit",
    HERE / "source_audit.py",
)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


class ManifestTests(unittest.TestCase):
    def test_manifest_is_strictly_source_only(self) -> None:
        manifest = AUDIT.load_manifest()
        self.assertEqual(
            manifest["scope"],
            "source_acquisition_and_compact_source_only",
        )
        self.assertTrue(manifest["authorization"]["source_file_already_acquired"])
        for key, value in manifest["authorization"].items():
            if key != "source_file_already_acquired":
                self.assertFalse(value, key)

    def test_forbidden_authorization_fails_closed(self) -> None:
        manifest = AUDIT.load_manifest()
        drifted = copy.deepcopy(manifest)
        drifted["authorization"]["evaluation"] = True
        with self.assertRaisesRegex(
            AUDIT.SourceAuditError,
            "forbidden authorization",
        ):
            with mock.patch.object(AUDIT, "load_json", return_value=drifted):
                AUDIT.load_manifest(Path("/synthetic/manifest.json"))


class MappingTests(unittest.TestCase):
    def test_projection_is_exact_symbol_only(self) -> None:
        projection = AUDIT.exact_symbol_projection(
            ["ABC", "ABC-AS1", "OLD1"],
            ["ABC-AS1", "ABC", "NEW1"],
        )
        np.testing.assert_array_equal(projection, [1, 0, -1])

    def test_duplicate_source_symbol_is_rejected_not_merged(self) -> None:
        with self.assertRaisesRegex(AUDIT.SourceAuditError, "duplicates"):
            AUDIT.exact_symbol_projection(
                ["ABC", "ABC"],
                ["ABC", "DEF"],
            )

    def test_api_rounding_and_nan_comparison(self) -> None:
        maximum = AUDIT.compare_api_vector(
            [1.23457, None, -0.25],
            np.asarray([1.2345678, np.nan, -0.250004], dtype=np.float32),
            5.1e-6,
            "synthetic",
        )
        self.assertLessEqual(maximum, 5.1e-6)

    def test_api_missingness_disagreement_fails(self) -> None:
        with self.assertRaisesRegex(AUDIT.SourceAuditError, "missing-value"):
            AUDIT.compare_api_vector(
                [1.0, None],
                np.asarray([1.0, 0.0]),
                5.1e-6,
                "synthetic",
            )


class RealSourceContractTests(unittest.TestCase):
    def test_real_source_audit_matches_frozen_contract(self) -> None:
        report = AUDIT.build_report()
        self.assertEqual(
            report["status"],
            "official_source_authenticated_and_semantics_verified",
        )
        self.assertEqual(report["summary"]["target_rows"], 1000)
        self.assertEqual(report["summary"]["supported_target_rows"], 997)
        self.assertEqual(report["summary"]["canonical_exact_overlap"], 7843)
        self.assertEqual(report["summary"]["vcc_panel_target_overlap"], 21)
        self.assertEqual(report["summary"]["vcc_panel_supported_overlap"], 20)
        self.assertEqual(
            report["summary"]["all_nan_rows"],
            ["COQ10B", "HSD17B12", "NEK1", "Nontargeting"],
        )
        self.assertEqual(
            report["semantics_evidence"]["official_api_field"],
            "log2_fold_change",
        )
        self.assertTrue(report["semantics_evidence"]["core_nan_pattern_consistent"])
        self.assertFalse(report["safety"]["held_truth_opened"])
        self.assertFalse(report["safety"]["existing_outcome_opened"])


if __name__ == "__main__":
    unittest.main()
