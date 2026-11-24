from __future__ import annotations

import importlib.util
import math
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("experiment_118_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def synthetic_condition(
    name: str, lfc: np.ndarray, se: np.ndarray | None = None
) -> object:
    lfc = np.asarray(lfc, dtype=np.float64)
    if se is None:
        se = np.ones_like(lfc)
    valid = np.isfinite(lfc) & np.isfinite(se) & (se > 0)
    np.fill_diagonal(valid, False)
    return RUN.Condition(
        name=name,
        lfc=lfc.copy(),
        se=np.asarray(se, dtype=np.float64),
        fdr=np.where(valid, 0.01, 1.0),
        source_valid=valid,
        truth_valid=valid & (lfc != 0),
        target_qc=np.ones(lfc.shape[0], dtype=bool),
    )


class Experiment118Tests(unittest.TestCase):
    def test_candidate_grid_is_exact(self) -> None:
        self.assertEqual(len(RUN.candidates()), 13)
        self.assertEqual(RUN.candidates()[0].name, "identity")
        self.assertEqual(
            {(c.estimator, c.eta) for c in RUN.candidates()[1:]},
            {
                (estimator, eta)
                for estimator in RUN.ESTIMATORS
                for eta in RUN.ETAS
            },
        )

    def test_estimator_requires_two_conditions_and_excludes_self(self) -> None:
        a = np.array(
            [[9.0, 2.0, 4.0], [3.0, 9.0, -2.0], [1.0, 5.0, 9.0]]
        )
        b = np.array(
            [[9.0, 4.0, 8.0], [6.0, 9.0, -4.0], [2.0, 10.0, 9.0]]
        )
        c1 = synthetic_condition("a", a)
        c2 = synthetic_condition("b", b)
        c2.source_valid[0, 1] = False
        result = RUN.estimate_base([c1, c2], "unweighted")
        self.assertEqual(result[0, 1], 0.0)
        self.assertTrue(np.all(np.diag(result) == 0))
        self.assertTrue(
            all(math.fsum(map(float, row)) == 0.0 for row in result)
        )

    def test_panel_residual_sum_is_exact_zero(self) -> None:
        rng = np.random.default_rng(118)
        conditions = []
        for index in range(3):
            lfc = rng.normal(size=(11, 11))
            np.fill_diagonal(lfc, 1000.0 + index)
            conditions.append(synthetic_condition(str(index), lfc))
        for estimator in RUN.ESTIMATORS:
            result = RUN.estimate_base(conditions, estimator)
            self.assertTrue(
                all(math.fsum(map(float, row)) == 0.0 for row in result)
            )
            self.assertTrue(np.all(np.diag(result) == 0))

    def test_pds_visible_coordinates_are_exactly_untouched(self) -> None:
        panel = ["A", "B", "C"]
        genes = ["X", "A", "Y", "B", "C", "Z"]
        candidate_panel = np.zeros((13, 3, 3), dtype=np.float64)
        candidate_panel[1, 0] = np.array([0.0, 2.0, -2.0])
        candidate_panel[1, 1] = np.array([3.0, 0.0, -3.0])
        candidate_panel[1, 2] = np.array([4.0, -4.0, 0.0])
        full = RUN.full_residual(
            RUN.candidates()[1].name, candidate_panel, panel, genes
        )
        visible = [0, 2, 5]
        self.assertTrue(np.all(full[:, visible] == 0))
        panel_cols = [1, 3, 4]
        self.assertTrue(
            all(math.fsum(map(float, row)) == 0.0 for row in full[:, panel_cols])
        )
        self.assertTrue(np.all(full[np.arange(3), panel_cols] == 0))

    def test_identity_wins_exact_selection_tie(self) -> None:
        zero = np.zeros((4, 4), dtype=np.float64)
        conditions = [
            synthetic_condition("a", zero),
            synthetic_condition("b", zero),
            synthetic_condition("c", zero),
        ]
        # Give held truth one valid pair per target while training values stay
        # zero, so every candidate ties identity on both metrics.
        for condition in conditions:
            for target in range(4):
                response = (target + 1) % 4
                condition.lfc[target, response] = 1.0
                condition.truth_valid[target, response] = True
                condition.source_valid[target, response] = False
        selected, _ = RUN.evaluate_grid(conditions)
        self.assertEqual(selected.name, "identity")

    def test_permutation_banks_are_deranged_and_bin_preserving(self) -> None:
        original_draws = RUN.NULL_DRAWS
        try:
            RUN.NULL_DRAWS = 4
            bins = np.repeat(np.arange(10, dtype=np.int16), 4)
            target, response = RUN.freeze_permutations(bins)
            RUN.validate_permutations(target, response, bins)
        finally:
            RUN.NULL_DRAWS = original_draws

    def test_real_authorization_is_self_consistent(self) -> None:
        authorization = RUN.validate_authorization()
        self.assertTrue(authorization["source_only_authorized"])
        self.assertTrue(all(value is False for value in authorization["external_actions"].values()))


if __name__ == "__main__":
    unittest.main()
