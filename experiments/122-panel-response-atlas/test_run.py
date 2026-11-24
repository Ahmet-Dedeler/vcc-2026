from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


freeze_plans = load_module("experiment_122_freeze_plans", HERE / "freeze_plans.py")
sys.modules["freeze_plans"] = freeze_plans
run = load_module("experiment_122_run", HERE / "run.py")


def matrix(name: str, effect: np.ndarray, valid: np.ndarray | None = None):
    effect = np.asarray(effect, dtype=np.float64)
    if valid is None:
        valid = np.ones_like(effect, dtype=bool)
    else:
        valid = np.asarray(valid, dtype=bool)
    diagonal = np.arange(effect.shape[0])
    effect = effect.copy()
    valid = valid.copy()
    effect[diagonal, diagonal] = 0.0
    valid[diagonal, diagonal] = False
    return run.Matrix(name, effect, valid)


class FreezeTests(unittest.TestCase):
    def test_derangements_are_deterministic_and_valid(self) -> None:
        first = freeze_plans.derangements(8, 20, 122)
        second = freeze_plans.derangements(8, 20, 122)
        self.assertTrue(np.array_equal(first, second))
        freeze_plans.validate_bank(first, 8, 20)

    def test_synthetic_freeze_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE) as directory:
            output = Path(directory) / "freeze"
            freeze_plans.build(output)
            freeze_plans.verify(output)
            with self.assertRaises(freeze_plans.FreezeError):
                freeze_plans.build(output)


class AggregationTests(unittest.TestCase):
    def test_correlated_members_collapse_to_one_family(self) -> None:
        a = matrix("a", np.asarray([[0, 1, 3], [2, 0, 4], [5, 7, 0]]))
        b = matrix("b", np.asarray([[0, 3, 1], [4, 0, 6], [7, 5, 0]]))
        family = run.aggregate_family("group", (a, b), minimum_member_support=2)
        self.assertEqual(family.name, "group")
        self.assertTrue(np.array_equal(family.effect[0], [0.0, 2.0, 2.0]))
        self.assertTrue(np.all(np.diag(family.effect) == 0.0))

    def test_minimum_member_support_is_enforced(self) -> None:
        valid_a = np.ones((3, 3), dtype=bool)
        valid_b = np.ones((3, 3), dtype=bool)
        valid_b[0, 1] = False
        a = matrix("a", np.ones((3, 3)), valid_a)
        b = matrix("b", np.ones((3, 3)), valid_b)
        family = run.aggregate_family("group", (a, b), minimum_member_support=2)
        self.assertFalse(family.valid[0, 1])
        self.assertEqual(family.effect[0, 1], 0.0)


class CandidateTests(unittest.TestCase):
    def test_self_is_zero_and_confidence_is_separate(self) -> None:
        a = matrix("a", np.asarray([[0, 1, -2], [3, 0, -1], [2, 4, 0]]))
        b = matrix("b", np.asarray([[0, 2, -3], [1, 0, -4], [3, 5, 0]]))
        candidate = next(c for c in run.CANDIDATES if c.name == "equal_mean_agreement")
        residual, confidence, support = run.build_candidate((a, b), candidate)
        self.assertTrue(np.all(np.diag(residual) == 0.0))
        self.assertTrue(np.all(np.diag(confidence) == 0.0))
        self.assertTrue(np.all(confidence >= 0.0))
        self.assertTrue(np.all(np.diag(support) == 0))
        self.assertFalse(np.array_equal(np.abs(residual), confidence))

    def test_sign_consensus_requires_two_agreeing_families(self) -> None:
        a = matrix("a", np.asarray([[0, 1, 2], [1, 0, 2], [1, 2, 0]]))
        b = matrix("b", np.asarray([[0, -1, 3], [1, 0, -2], [1, 2, 0]]))
        candidate = next(
            c for c in run.CANDIDATES if c.name == "sign_consensus_median_min2"
        )
        residual, _, _ = run.build_candidate((a, b), candidate)
        self.assertEqual(residual[0, 1], 0.0)
        self.assertGreater(residual[0, 2], 0.0)

    def test_expansion_is_exactly_pds_null(self) -> None:
        panel = ("A", "B", "C")
        genes = ("X", "A", "Y", "B", "C", "Z")
        residual = np.asarray([[0, 1, -1], [2, 0, -2], [3, -3, 0]], dtype=float)
        confidence = np.abs(residual)
        full_residual = run.expand_panel_matrix(residual, panel, genes)
        full_confidence = run.expand_panel_matrix(confidence, panel, genes)
        report = run.invariant_report(
            residual,
            confidence,
            full_residual,
            full_confidence,
            panel,
            genes,
        )
        self.assertTrue(report["all_pass"])
        self.assertTrue(np.all(full_residual[:, [0, 2, 5]] == 0.0))


class SelectionTests(unittest.TestCase):
    def test_positive_shared_signal_selects_nonidentity(self) -> None:
        rng = np.random.default_rng(122)
        n = 30
        base = rng.normal(size=(n, n))
        np.fill_diagonal(base, 0.0)
        families = (
            matrix("replogle_k562", base),
            matrix("orion", base + rng.normal(scale=0.01, size=(n, n))),
            matrix("cd4", base + rng.normal(scale=0.01, size=(n, n))),
        )
        selected, summaries = run.evaluate_grid(
            families,
            tuple(f"T{i}" for i in range(n)),
        )
        self.assertNotEqual(selected.name, "identity")
        selected_summary = next(row for row in summaries if row["name"] == selected.name)
        self.assertTrue(selected_summary["eligible"])

    def test_target_permutation_breaks_target_alignment(self) -> None:
        rng = np.random.default_rng(7)
        n = 30
        base = rng.normal(size=(n, n))
        np.fill_diagonal(base, 0.0)
        families = tuple(
            matrix(name, base.copy())
            for name in ("replogle_k562", "orion", "cd4")
        )
        selected, summaries = run.evaluate_grid(
            families,
            tuple(f"T{i}" for i in range(n)),
        )
        observed = next(row for row in summaries if row["name"] == selected.name)
        permutation = np.roll(np.arange(n), 1)
        null_selected, null_summaries = run.evaluate_grid(
            families,
            tuple(f"T{i}" for i in range(n)),
            target_permutation=permutation,
        )
        null = next(row for row in null_summaries if row["name"] == null_selected.name)
        self.assertGreater(
            observed["worst_family_signed_credit"],
            null["worst_family_signed_credit"],
        )


class ContractTests(unittest.TestCase):
    def test_duplicate_nonpanel_gene_labels_are_allowed_but_panel_duplicates_fail(self) -> None:
        run.validate_axis(("A", "X", "X"), "genes", require_unique=False)
        panel_positions, source_positions = run.panel_indices(
            ("A", "X", "X"),
            ("A", "B"),
        )
        self.assertTrue(np.array_equal(panel_positions, [0]))
        self.assertTrue(np.array_equal(source_positions, [0]))
        with self.assertRaises(run.ExperimentError):
            run.panel_indices(("A", "A", "X"), ("A", "B"))

    def test_protocol_and_inventory_are_narrow(self) -> None:
        protocol = json.loads(
            (HERE / "preregistered_protocol.json").read_text(encoding="utf-8")
        )
        inventory = json.loads(
            (HERE / "SOURCE_INVENTORY.json").read_text(encoding="utf-8")
        )
        authorization = json.loads(
            (HERE / "SOURCE_ONLY_AUTHORIZATION.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            set(protocol["authorized_biological_families"]),
            {"replogle_k562", "orion", "cd4"},
        )
        self.assertTrue(inventory["coverage_summary"]["all_required_compact_panel_coordinates_recoverable"])
        self.assertFalse(authorization["authorization"]["read_a549"])
        self.assertFalse(
            authorization["authorization"]["read_dld1_experiment_116_matrix_values"]
        )
        self.assertFalse(
            authorization["authorization"]["read_vcc_held_harness_or_outcomes"]
        )


if __name__ == "__main__":
    unittest.main()
