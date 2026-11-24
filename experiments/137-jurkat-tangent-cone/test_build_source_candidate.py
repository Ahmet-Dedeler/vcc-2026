from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "exp137_build_source_candidate", HERE / "build_source_candidate.py"
)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)

TANGENT_SPEC = importlib.util.spec_from_file_location(
    "exp137_test_tangent_runtime", HERE / "tangent.py"
)
assert TANGENT_SPEC is not None and TANGENT_SPEC.loader is not None
TANGENT = importlib.util.module_from_spec(TANGENT_SPEC)
sys.modules[TANGENT_SPEC.name] = TANGENT
TANGENT_SPEC.loader.exec_module(TANGENT)


def test_cap_key_is_stable() -> None:
    assert M.cap_key(0.0) == "cap_0000"
    assert M.cap_key(1.5) == "cap_0150"
    assert M.cap_key(6.0) == "cap_0600"


def test_build_candidates_preserves_identity_and_fallback() -> None:
    arrays = {
        "incumbent_score": np.asarray(
            [[1.0, 0.0, 9.0], [0.0, 1.0, 8.0]], dtype=np.float64
        ),
        "source_score": np.asarray(
            [[0.0, 1.0, 7.0], [1.0, 0.0, 6.0]], dtype=np.float64
        ),
        "transferable": np.asarray([True, False]),
        "coordinate_supported": np.asarray([True, True, False]),
    }
    built = M.build_candidates(arrays, TANGENT)
    zero = built["candidate_cap_0000"]
    assert np.array_equal(
        zero.view(np.uint8), arrays["incumbent_score"].view(np.uint8)
    )
    for cap in TANGENT.CAP_DEGREES[1:]:
        candidate = built[f"candidate_{M.cap_key(cap)}"]
        assert np.array_equal(
            candidate[1].view(np.uint8),
            arrays["incumbent_score"][1].view(np.uint8),
        )
        assert candidate[0, 2] == arrays["incumbent_score"][0, 2]


def test_build_plan_uses_all_394_targets() -> None:
    metadata, _ = M.validate_frozen_inputs()

    class SurfacePlan(SimpleNamespace):
        pass

    plan = M.build_plan(metadata, SimpleNamespace(SurfacePlan=SurfacePlan))
    assert len(plan.targets) == 394
    assert plan.targets[:2] == ("TFB2M", "RHPN1")
    assert plan.targets[78] == "UQCC2"
    assert len(plan.score_genes) == 223
    assert len(plan.operator_columns) == 17_585


if __name__ == "__main__":
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} tests passed")
