from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "exp137_freeze_metadata", HERE / "freeze_metadata.py"
)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def test_metadata_freeze_contract() -> None:
    arrays, report = M.build_freeze()
    assert report["split"]["fresh"] == 394
    assert report["split"]["test"] == 79
    assert report["split"]["development"] == 315
    assert report["split"]["test_transferable_targets"] == 79
    assert report["split"]["unsupported_fresh_targets"] == [
        "ZNF658",
        "LPA",
        "NOMO3",
    ]
    assert arrays["test_targets"][0] == "TFB2M"
    assert arrays["test_targets"][-1] == "UQCC2"
    assert arrays["development_fold_id"].shape == (315,)
    assert np.array_equal(
        np.bincount(arrays["development_fold_id"], minlength=6)[1:],
        np.repeat(63, 5),
    )
    assert arrays["score_genes"].shape == (223,)
    assert arrays["operator_columns"].shape == (17_585,)


if __name__ == "__main__":
    test_metadata_freeze_contract()
    print("1 test passed")
