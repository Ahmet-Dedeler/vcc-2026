#!/usr/bin/env python3
"""Synthetic invariants for experiment 104."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp104_run", HERE / "run.py")
if SPEC is None or SPEC.loader is None:
    raise ImportError(HERE / "run.py")
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def test_masks() -> None:
    panel = RUN.PANEL[:2]
    safe = next(gene for gene in RUN.OFFICIAL if gene not in RUN.PANEL_SET)
    genes = np.asarray([panel[0], safe, panel[1], "NOT_ON_OFFICIAL_AXIS"])
    dmask = RUN.directional_mask(genes, panel[0])
    pmask = RUN.pds_mask(genes)
    assert dmask.tolist() == [False, True, True, False]
    assert pmask.tolist() == [False, True, False, False]


def test_transform_is_bounded_sign_preserving_and_norm_preserving() -> None:
    effect = np.asarray([[2.0, -1.0, 0.0, 3.0], [0.0, 0.0, 0.0, 0.0]])
    score = np.asarray([[2.0, -1.0, 9.0, 0.2], [4.0, 3.0, 2.0, 1.0]])
    active = np.asarray([True, True, True, True])
    supported = np.asarray([[True, True, False, True], [False, False, False, False]])
    result = RUN.transform_effect(effect, score, active, supported, beta=0.1)
    np.testing.assert_allclose(
        np.linalg.norm(result[0]), np.linalg.norm(effect[0]), rtol=1e-12, atol=1e-12
    )
    assert np.array_equal(np.sign(result), np.sign(effect))
    assert np.array_equal(result[effect == 0], effect[effect == 0])
    assert np.array_equal(result[1], effect[1])


def test_identity_beta_zero() -> None:
    rng = np.random.default_rng(4)
    effect = rng.normal(size=(3, 11))
    effect[:, 2] = 0
    score = rng.normal(size=effect.shape)
    active = np.ones(effect.shape[1], dtype=bool)
    supported = np.ones_like(effect, dtype=bool)
    result = RUN.transform_effect(effect, score, active, supported, beta=0.0)
    np.testing.assert_array_equal(result, effect)


def test_target_folds_are_whole_and_deterministic() -> None:
    targets = ["A", "B", "A", "C", "B"]
    observed = [RUN.target_fold(target) for target in targets]
    assert observed[0] == observed[2]
    assert observed[1] == observed[4]
    assert observed == [RUN.target_fold(target) for target in targets]
    assert all(0 <= value < RUN.N_FOLDS for value in observed)


def _blocks() -> list[object]:
    rows = []
    for context in ("left", "right"):
        for target_index in range(8):
            n = 10 + target_index
            x = np.column_stack(
                [
                    np.linspace(-1, 1, n),
                    np.full(n, target_index),
                ]
            )
            y = np.where(x[:, 0] >= 0, 1.0, -1.0)
            weight = np.linspace(0.1, 1.0, n)
            rows.append(
                RUN.TrainingBlock(
                    context=context,
                    target=f"T{target_index}",
                    x_agree=x,
                    x_full=np.column_stack([x, x[:, 0] ** 2]),
                    y=y,
                    weight=weight,
                )
            )
    return rows


def test_credit_model_and_ablation_dimensions() -> None:
    blocks = _blocks()
    models = RUN.fit_models(blocks)
    assert models["agreement"].coef.shape == (2,)
    assert models["full"].coef.shape == (3,)
    probe = np.asarray([[-0.8, 1.0], [0.8, 1.0]])
    pred = models["agreement"].predict(probe)
    assert pred[1] > pred[0]


def test_permutation_is_deterministic_and_changes_labels() -> None:
    blocks = _blocks()
    x1, y1, w1 = RUN.collect_blocks(blocks, "full", null_draw=3)
    x2, y2, w2 = RUN.collect_blocks(blocks, "full", null_draw=3)
    xr, yr, wr = RUN.collect_blocks(blocks, "full")
    np.testing.assert_array_equal(x1, x2)
    np.testing.assert_array_equal(y1, y2)
    np.testing.assert_array_equal(w1, w2)
    np.testing.assert_array_equal(x1, xr)
    np.testing.assert_array_equal(w1, wr)
    assert np.any(y1 != yr)
    assert sorted(y1.tolist()) == sorted(yr.tolist())


def test_write_scope_fails_closed() -> None:
    with tempfile.TemporaryDirectory() as directory:
        outside = Path(directory) / "bad.json"
        try:
            RUN.atomic_json(outside, {"bad": True})
        except ValueError:
            pass
        else:
            raise AssertionError("outside write was accepted")


def main() -> None:
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"{len(tests)} tests passed")


if __name__ == "__main__":
    main()
