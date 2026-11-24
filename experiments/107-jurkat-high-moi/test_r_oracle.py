#!/usr/bin/env python3
"""Compare the Python preprocessing port with pinned R package behavior."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import preprocess  # noqa: E402


EXPECTED_VERSIONS = {"DropletUtils": "1.31.2", "edgeR": "4.99.0"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def synthetic_matrix() -> sparse.csc_matrix:
    """Make a deterministic matrix with ambient and cell-like barcodes."""

    rng = np.random.default_rng(107_2026)
    ambient_totals = np.arange(21, 101, dtype=np.int64)
    cell_totals = np.r_[
        np.arange(120, 170, dtype=np.int64),
        np.arange(200, 230, dtype=np.int64),
    ]
    totals = np.r_[ambient_totals, cell_totals]
    ambient_probability = np.array(
        [
            0.28,
            0.19,
            0.14,
            0.10,
            0.08,
            0.06,
            0.05,
            0.035,
            0.025,
            0.015,
            0.01,
            0.095,
            0.0,
        ],
        dtype=np.float64,
    )
    ambient_probability /= ambient_probability.sum()

    columns: list[np.ndarray] = []
    for index, total in enumerate(totals):
        probability = ambient_probability.copy()
        if index >= ambient_totals.size:
            enriched = index % probability.size
            probability[enriched] += 0.18
            probability /= probability.sum()
        columns.append(rng.multinomial(int(total), probability))
    return sparse.csc_matrix(np.stack(columns, axis=1), dtype=np.int64)


def run() -> dict[str, object]:
    matrix = synthetic_matrix()
    with tempfile.TemporaryDirectory(prefix="vcc107-r-oracle-") as directory:
        temp = Path(directory)
        matrix_path = temp / "synthetic.mtx"
        output_path = temp / "oracle.json"
        mmwrite(matrix_path, matrix)
        subprocess.run(
            ["Rscript", str(HERE / "r_oracle.R"), str(matrix_path), str(output_path)],
            check=True,
        )
        oracle = json.loads(output_path.read_text())

    if oracle["versions"] != EXPECTED_VERSIONS:
        raise AssertionError(
            f"R oracle package versions drifted: {oracle['versions']!r}"
        )

    dense = np.asarray(matrix.toarray(), dtype=np.int64)
    totals = dense.sum(axis=0)
    rank = preprocess.barcode_ranks_thresholds(totals)
    feature_nonzero = dense.sum(axis=1) > 0
    ambient_counts = dense[:, totals <= 100].sum(axis=1)
    ambient = np.zeros(dense.shape[0], dtype=np.float64)
    ambient[feature_nonzero] = preprocess.safe_good_turing_proportions(
        ambient_counts[feature_nonzero]
    )
    tested = totals > 100
    observed = np.array(
        [
            preprocess.multinomial_data_log_probability(
                np.flatnonzero(dense[:, column]),
                dense[dense[:, column] > 0, column],
                ambient,
            )
            for column in np.flatnonzero(tested)
        ],
        dtype=np.float64,
    )
    bh_input = np.asarray(oracle["bh_input"], dtype=np.float64)
    bh = preprocess.benjamini_hochberg(bh_input)

    np.testing.assert_array_equal(totals, oracle["totals"])
    np.testing.assert_array_equal(np.flatnonzero(tested), oracle["tested"])
    np.testing.assert_allclose(rank.knee, oracle["knee"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(
        rank.inflection,
        oracle["inflection"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(ambient, oracle["ambient"], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(
        observed,
        oracle["observed_data_log_probability"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(bh, oracle["bh_adjusted"], rtol=0, atol=1e-15)

    return {
        "status": "pass",
        "synthetic_shape": list(matrix.shape),
        "synthetic_nnz": int(matrix.nnz),
        "versions": oracle["versions"],
        "r_oracle_sha256": sha256(HERE / "r_oracle.R"),
        "python_port_sha256": sha256(HERE / "preprocess.py"),
        "checks": [
            "column_totals",
            "barcode_ranks_knee",
            "barcode_ranks_inflection",
            "safe_good_turing_ambient_profile",
            "multinomial_observed_data_log_probability",
            "benjamini_hochberg",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
