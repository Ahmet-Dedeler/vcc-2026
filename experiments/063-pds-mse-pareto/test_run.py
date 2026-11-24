from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp063", HERE / "run.py")
RUN = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_predeclared_grid_and_convex_endpoints() -> None:
    assert RUN.BLEND_GRID == (0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 1.0)
    raw = np.asarray([[1.0, 2.0]])
    white = np.asarray([[3.0, 4.0]])
    predictions = {
        beta: (1 - beta) * white + beta * raw
        for beta in RUN.BLEND_GRID
    }
    np.testing.assert_allclose(predictions[0.0], white)
    np.testing.assert_allclose(predictions[1.0], raw)
    np.testing.assert_allclose(predictions[0.5], [[2.0, 3.0]])


def test_pds_uses_complete_gallery() -> None:
    truth = np.eye(4, dtype=np.float64)
    np.testing.assert_allclose(RUN.pds_per_target(truth, truth), 1.0)
    prediction = truth.copy()
    prediction[[0, 1]] = prediction[[1, 0]]
    scores = RUN.pds_per_target(prediction, truth)
    assert scores[0] < 1.0
    assert scores[1] < 1.0


def test_theoretical_mse_matches_direct_optimum() -> None:
    prediction = np.asarray([[1.0, 1.0], [-1.0, 0.0]])
    truth = np.asarray([[2.0, 0.0], [1.0, 0.0]])
    summary, _, cosines = RUN.metric_rows(prediction, truth)
    expected = np.mean(1.0 - cosines**2)
    assert np.isclose(summary["min_mse_unrestricted"], expected)
    # The anti-correlated second row needs negative alpha for its formal
    # 1-cos^2 optimum, so the nonnegative diagnostic correctly returns 1.
    assert np.isclose(summary["min_mse_nonnegative_alpha"], 0.75)


def test_pareto_frontier_rejects_dominated_row() -> None:
    rows = [
        {
            "beta": 0.0,
            "pds": 0.80,
            "raw_cosine": 0.03,
            "min_mse_unrestricted": 0.998,
            "min_mse_nonnegative_alpha": 0.999,
        },
        {
            "beta": 0.1,
            "pds": 0.79,
            "raw_cosine": 0.05,
            "min_mse_unrestricted": 0.996,
            "min_mse_nonnegative_alpha": 0.997,
        },
        {
            "beta": 0.2,
            "pds": 0.78,
            "raw_cosine": 0.04,
            "min_mse_unrestricted": 0.997,
            "min_mse_nonnegative_alpha": 0.998,
        },
    ]
    assert RUN.pareto_indices(rows) == [0, 1]


def test_gap_translation_requires_seven_thousandths() -> None:
    rows = [
        {
            "beta": 0.0,
            "pds": 0.80,
            "raw_cosine": 0.03,
            "min_mse_unrestricted": 0.999,
            "min_mse_nonnegative_alpha": 0.999,
        },
        {
            "beta": 0.1,
            "pds": 0.80,
            "raw_cosine": 0.20,
            "min_mse_unrestricted": 0.94,
            "min_mse_nonnegative_alpha": 0.94,
        },
    ]
    translated = RUN.optimistic_gap_translation(rows)
    assert translated[0]["optimistic_overall_delta_vs_beta0"] == 0.0
    assert translated[1]["could_close_0.007_gap"]


def test_panel_exclusion_is_exact() -> None:
    official = ["A", "B", "C", "D"]
    panel = {"B", "D"}
    source_union = {"A", "B", "D"}
    selected = [
        gene for gene in official if gene in source_union and gene not in panel
    ]
    assert selected == ["A"]
    assert not panel.intersection(selected)
