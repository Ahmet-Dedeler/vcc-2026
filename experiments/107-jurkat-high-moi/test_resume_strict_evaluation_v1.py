from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

import resume_strict_evaluation_v1 as RECOVERY
import run_strict_evaluation as ORIGINAL
import strict_nadig_validation as STRICT


class AxisRecoveryTests(unittest.TestCase):
    def test_only_three_held_aliases_are_excluded_and_hspa14_remains(self) -> None:
        official = ("TBCE", "HSPA14", "TMSB15B", "X", "OFFICIAL_ONLY")
        held = (
            "TBCE",
            "TBCE-1",
            "HSPA14",
            "HSPA14-1",
            "TMSB15B",
            "TMSB15B-1",
            "X",
        )
        evidence = RECOVERY.build_axis_evidence(
            official_genes=official,
            held_genes=held,
            reference_genes=held,
            integration_genes=official,
            score_genes=("HSPA14", "X"),
        )
        self.assertEqual(
            evidence["held_only_genes"],
            ["TBCE-1", "HSPA14-1", "TMSB15B-1"],
        )
        self.assertEqual(evidence["held_only_gene_count"], 3)
        self.assertEqual(evidence["official_only_gene_count"], 1)
        self.assertEqual(evidence["representable_gene_count"], 4)
        self.assertEqual(
            evidence["hspa14_exact_symbol"],
            {
                "held_alias_index": 3,
                "held_index": 2,
                "official_index": 1,
                "representable_index": 1,
            },
        )

        official_indices, held_indices, representable, excluded = (
            RECOVERY.exact_symbol_projection(official, held)
        )
        np.testing.assert_array_equal(official_indices, [0, 1, 2, 3])
        np.testing.assert_array_equal(held_indices, [0, 2, 4, 6])
        self.assertEqual(representable, ("TBCE", "HSPA14", "TMSB15B", "X"))
        self.assertEqual(excluded, RECOVERY.EXPECTED_HELD_ONLY)

    def test_reference_projection_keeps_hspa14_without_merging_alias(self) -> None:
        held = np.asarray(
            [
                "TBCE",
                "TBCE-1",
                "HSPA14",
                "HSPA14-1",
                "TMSB15B",
                "TMSB15B-1",
                "X",
            ]
        )
        reference = ORIGINAL.ReferenceDE(
            genes=held,
            gate=np.ones(len(held), dtype=bool),
            gate_indices=np.arange(len(held), dtype=np.int64),
            control_mean_cpm=np.arange(1, len(held) + 1, dtype=np.float64),
            targets=np.asarray(["T"]),
            lfc=np.asarray([[10.0, 100.0, 20.0, 200.0, 30.0, 300.0, 40.0]]),
            significant=(np.arange(len(held), dtype=np.int64),),
        )
        _, held_indices, representable, _ = RECOVERY.exact_symbol_projection(
            ("TBCE", "HSPA14", "TMSB15B", "X", "OFFICIAL_ONLY"),
            held,
        )
        projected = RECOVERY.project_reference_de(reference, held_indices)
        self.assertEqual(tuple(projected.genes), representable)
        self.assertIn("HSPA14", projected.genes)
        self.assertNotIn("HSPA14-1", projected.genes)
        np.testing.assert_array_equal(projected.lfc, [[10.0, 20.0, 30.0, 40.0]])
        np.testing.assert_array_equal(projected.significant[0], [0, 1, 2, 3])

    def test_pds_is_identical_to_the_original_frozen_score_axis(self) -> None:
        official = np.asarray(["A", "HSPA14", "B", "C"])
        held = np.asarray(["A", "HSPA14", "HSPA14-1", "B", "C"])
        score = np.asarray(["A", "HSPA14", "B"])
        score_indices = np.asarray([0, 1, 2])
        incumbent = np.asarray(
            [[1.0, 0.0, 0.0, 2.0], [0.0, 1.0, 0.0, 2.0], [0.0, 0.0, 1.0, 2.0]]
        )
        candidate = np.asarray(
            [[1.0, 0.2, 0.0, 2.0], [0.1, 1.0, 0.0, 2.0], [0.0, 0.1, 1.0, 2.0]]
        )
        held_truth = np.asarray(
            [
                [1.0, 0.1, 999.0, 0.0, 5.0],
                [0.0, 1.0, -999.0, 0.1, 5.0],
                [0.1, 0.0, 777.0, 1.0, 5.0],
            ]
        )
        incumbent_pds, candidate_pds = RECOVERY.pds_on_frozen_score_axis(
            incumbent_full=incumbent,
            candidate_full=candidate,
            held_truth=held_truth,
            official_genes=official,
            held_genes=held,
            score_genes=score,
            score_indices=score_indices,
        )
        held_lookup = {gene: index for index, gene in enumerate(held)}
        truth_indices = np.asarray([held_lookup[gene] for gene in score])
        expected_truth = held_truth[:, truth_indices]
        np.testing.assert_array_equal(
            incumbent_pds,
            STRICT.pds_per_target(incumbent[:, score_indices], expected_truth),
        )
        np.testing.assert_array_equal(
            candidate_pds,
            STRICT.pds_per_target(candidate[:, score_indices], expected_truth),
        )


class AuthorizationTests(unittest.TestCase):
    def test_authentication_failure_precedes_outcome_loader(self) -> None:
        touched = False

        def outcome_loader(*args: object, **kwargs: object) -> np.ndarray:
            nonlocal touched
            touched = True
            raise AssertionError("outcome loader must not run")

        with mock.patch.object(
            RECOVERY,
            "authenticate_recovery_amendment",
            side_effect=RECOVERY.RecoveryError("bad recovery"),
        ):
            with self.assertRaisesRegex(RECOVERY.RecoveryError, "bad recovery"):
                RECOVERY.authorize_then_reload_nadig(
                    lease=RECOVERY.RecoveryLease(
                        path=Path("/unused"),
                        descriptor=-1,
                        device=-1,
                        inode=-1,
                    ),
                    outcome_loader=outcome_loader
                )
        self.assertFalse(touched)

    def test_existing_receipt_is_checked_before_reload(self) -> None:
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
            mock.patch.object(STRICT, "authorize_outcome_phase", side_effect=authorize),
            mock.patch.object(STRICT, "verify_nadig_artifacts", side_effect=verify),
            tempfile.TemporaryDirectory() as temporary,
        ):
            lock = Path(temporary) / "recovery.lock"
            with RECOVERY._one_shot_recovery_lock(lock) as lease:  # noqa: SLF001
                authenticated = {
                    "original_auth": {"plan": plan},
                    "receipt": receipt,
                    "recovery_lock_path": lock.resolve(),
                }
                RECOVERY._reload_nadig_after_existing_receipt(  # noqa: SLF001
                    authenticated,
                    lease=lease,
                    outcome_loader=load,
                )
        self.assertEqual(events, ["authorize", "verify-receipt", "load-outcome"])

    def test_closed_or_forged_recovery_lease_cannot_reload(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            lock = Path(temporary) / "recovery.lock"
            with RECOVERY._one_shot_recovery_lock(lock) as lease:  # noqa: SLF001
                pass
            authenticated = {"recovery_lock_path": lock.resolve()}
            with (
                mock.patch.object(STRICT, "authorize_outcome_phase") as authorize,
                self.assertRaisesRegex(RECOVERY.RecoveryError, "no longer valid"),
            ):
                RECOVERY._reload_nadig_after_existing_receipt(  # noqa: SLF001
                    authenticated,
                    lease=lease,
                )
            authorize.assert_not_called()

    def test_recovery_contract_forbids_post_outcome_tuning(self) -> None:
        payload = json.loads(RECOVERY.RECOVERY_AMENDMENT_PATH.read_text())
        contract = payload["recovery_contract"]
        self.assertEqual(contract["replay_limit"], 1)
        self.assertFalse(contract["alias_merging"])
        self.assertEqual(
            contract["frozen_components_unchanged"],
            [
                "source_model",
                "integration_candidates",
                "fraction_grid",
                "nested_selection",
                "permutation_null",
                "cluster_bootstrap",
                "promotion_gates",
                "pds_score_axes",
            ],
        )
        original = json.loads(ORIGINAL.AMENDMENT_PATH.read_text())
        self.assertEqual(
            original["execution_contract"]["fractions"],
            list(STRICT.FRACTIONS),
        )
        self.assertEqual(
            original["execution_contract"]["permutation_replicates"],
            STRICT.PERMUTATION_REPLICATES,
        )
        self.assertEqual(
            original["execution_contract"]["cluster_bootstrap_replicates"],
            STRICT.BOOTSTRAP_REPLICATES,
        )


class OutputInvariantTests(unittest.TestCase):
    def test_new_lock_allows_exactly_one_attempt_and_preserves_old_lock(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old_lock = root / ".strict-evaluation-run.lock"
            new_lock = root / ".strict-evaluation-recovery-v1.lock"
            old_lock.write_bytes(b"")
            old_hash = STRICT.sha256_file(old_lock)
            with RECOVERY._one_shot_recovery_lock(new_lock):  # noqa: SLF001
                self.assertTrue(new_lock.is_file())
                self.assertEqual(new_lock.stat().st_size, 0)
            self.assertEqual(STRICT.sha256_file(old_lock), old_hash)
            self.assertTrue(new_lock.is_file())
            with self.assertRaisesRegex(RECOVERY.RecoveryError, "already claimed"):
                with RECOVERY._one_shot_recovery_lock(new_lock):  # noqa: SLF001
                    pass

    def test_frozen_production_receipt_and_failed_lock_are_unchanged(self) -> None:
        payload = json.loads(RECOVERY.RECOVERY_AMENDMENT_PATH.read_text())
        predecessors = payload["predecessor_bindings"]
        for name in ("outcome_access_receipt", "failed_run_lock"):
            binding = predecessors[name]
            path = RECOVERY.ROOT / binding["path"]
            self.assertEqual(path.stat().st_size, binding["bytes"])
            self.assertEqual(STRICT.sha256_file(path), binding["sha256"])
        receipt_path = RECOVERY.ROOT / predecessors["outcome_access_receipt"]["path"]
        self.assertEqual(stat_mode(receipt_path) & 0o222, 0)

    def test_frozen_production_axis_evidence_is_exact(self) -> None:
        payload = json.loads(RECOVERY.RECOVERY_AMENDMENT_PATH.read_text())
        h1 = payload["axis_evidence"]["h1"]
        hek = payload["axis_evidence"]["orion_hek293t"]
        self.assertEqual(h1["held_only_genes"], list(RECOVERY.EXPECTED_HELD_ONLY))
        self.assertEqual(h1["held_only_gene_count"], 3)
        self.assertEqual(h1["representable_gene_count"], 18_077)
        self.assertEqual(h1["hspa14_exact_symbol"]["held_index"], 8676)
        self.assertEqual(h1["hspa14_exact_symbol"]["held_alias_index"], 8677)
        self.assertTrue(h1["score_axis_all_representable"])
        self.assertEqual(hek["held_only_genes"], [])
        self.assertEqual(hek["representable_gene_count"], 18_533)


def stat_mode(path: Path) -> int:
    return path.stat().st_mode


if __name__ == "__main__":
    unittest.main()
