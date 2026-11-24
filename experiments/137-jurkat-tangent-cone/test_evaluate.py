from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp137_evaluate", HERE / "evaluate.py")
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)

TANGENT_SPEC = importlib.util.spec_from_file_location(
    "exp137_eval_tangent", HERE / "tangent.py"
)
assert TANGENT_SPEC is not None and TANGENT_SPEC.loader is not None
TANGENT = importlib.util.module_from_spec(TANGENT_SPEC)
sys.modules[TANGENT_SPEC.name] = TANGENT
TANGENT_SPEC.loader.exec_module(TANGENT)


def test_derangement_has_no_fixed_points_and_is_deterministic() -> None:
    left = M.derangement(np.random.default_rng(137), 79)
    right = M.derangement(np.random.default_rng(137), 79)
    assert np.array_equal(left, right)
    assert np.array_equal(np.sort(left), np.arange(79))
    assert np.all(left != np.arange(79))


def test_reciprocal_clusters_are_deterministic() -> None:
    values = np.asarray(
        [[1.0, 0.0], [0.99, 0.01], [0.0, 1.0], [0.01, 0.99]],
        dtype=np.float64,
    )
    left = M.reciprocal_knn_clusters(
        values, fold_one_based=1, k=1, minimum_cosine=0.8
    )
    right = M.reciprocal_knn_clusters(
        values, fold_one_based=1, k=1, minimum_cosine=0.8
    )
    assert np.array_equal(left, right)
    assert left[0] == left[1]
    assert left[2] == left[3]
    assert left[0] != left[2]


def test_nested_statistics_selects_without_held_fold() -> None:
    pds = {}
    credit = {}
    for cap in TANGENT.CAP_DEGREES:
        gain = 0.01 if cap == 3.0 else 0.0
        pds[cap] = tuple(np.full(63, 0.5 + gain) for _ in range(5))
        credit[cap] = tuple(np.full(63, 0.2 + gain) for _ in range(5))
    summary = M.nested_statistics(pds, credit, TANGENT)
    assert summary["selected_caps_degrees"] == [3.0] * 5
    assert np.isclose(summary["macro_pds_delta"], 0.01)
    assert np.isclose(summary["macro_signed_credit_delta"], 0.01)
    assert summary["final_cap_degrees"] == 3.0


def test_development_scores_fix_unsupported_rows_at_half() -> None:
    rng = np.random.default_rng(137)
    incumbent = rng.normal(size=(394, 8))
    transferable = np.ones(394, dtype=bool)
    unsupported_rows = np.asarray([79, 142, 205])
    transferable[unsupported_rows] = False
    candidate = {
        "development_fold_id": np.repeat(np.arange(1, 6, dtype=np.int8), 63),
        "incumbent_score": incumbent,
        "transferable": transferable,
    }
    for cap in TANGENT.CAP_DEGREES:
        values = incumbent.copy()
        if cap:
            values[transferable] += cap * 0.01
        candidate[f"candidate_{M.cap_key(cap)}"] = values
    truth = rng.normal(size=(315, 8))
    pds, _, _ = M.development_scores(candidate, truth, TANGENT)
    for cap in TANGENT.CAP_DEGREES:
        assert pds[cap][0][0] == 0.5
        assert pds[cap][1][0] == 0.5
        assert pds[cap][2][0] == 0.5


def test_fix_unsupported_pds_is_shared_and_nonmutating() -> None:
    scores = np.asarray([0.1, 0.2, 0.3])
    transferable = np.asarray([True, False, True])
    fixed = M.fix_unsupported_pds(scores, transferable)
    assert np.array_equal(fixed, np.asarray([0.1, 0.5, 0.3]))
    assert np.array_equal(scores, np.asarray([0.1, 0.2, 0.3]))


def test_validated_test_cap_rejects_forged_or_mismatched_cap() -> None:
    candidate = {
        "incumbent_score": np.zeros((79, 4)),
        "candidate_cap_0300": np.zeros((79, 4)),
    }
    development = {"nested": {"final_cap_degrees": 3.0}}
    assert (
        M.validated_test_cap(
            {"selected_cap_degrees": 3.0}, development, candidate, TANGENT
        )
        == 3.0
    )
    for selected in (4.5, 6.0):
        try:
            M.validated_test_cap(
                {"selected_cap_degrees": selected}, development, candidate, TANGENT
            )
        except M.EvaluationError:
            pass
        else:
            raise AssertionError(f"forged cap {selected} was accepted")


def test_gain_concentration() -> None:
    deltas = (np.asarray([1.0, 1.0]), np.asarray([1.0, -1.0]))
    clusters = (
        np.asarray(["a", "b"]),
        np.asarray(["c", "d"]),
    )
    assert np.isclose(M.gain_concentration(deltas, clusters, 0.5), 1.0 / 3.0)
    assert np.isinf(M.gain_concentration(deltas, clusters, 0.0))


def test_access_lock_allows_only_exact_resume() -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "access.lock"
        payload = {"stage": "development", "binding": "abc"}
        M.create_or_validate_lock(path, payload)
        M.create_or_validate_lock(path, payload)
        try:
            M.create_or_validate_lock(path, {"stage": "development", "binding": "changed"})
        except M.EvaluationError:
            pass
        else:
            raise AssertionError("changed access lock was accepted")


if __name__ == "__main__":
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} tests passed")
