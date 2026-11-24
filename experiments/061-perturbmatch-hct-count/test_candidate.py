"""Lightweight synthetic checks for experiment 061."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BUILD = load_module("exp061_build", "build_candidate.py")
COMPARE = load_module("exp061_compare", "compare_scores.py")


def test_blend_preserves_missing_gene_coordinates() -> None:
    baseline = np.asarray([10.0, 20.0, 30.0, 40.0])
    incumbent = np.asarray([1.0, 2.0, 3.0, 4.0])
    pm = np.asarray([8.0, 10.0])
    columns = np.asarray([0, 2])
    output = BUILD.blend_k562_slot(
        baseline,
        incumbent,
        pm,
        columns,
        effect_scale=2.0,
        blend=0.5,
        source_weight=0.5,
    )
    np.testing.assert_array_equal(output[[1, 3]], baseline[[1, 3]])
    np.testing.assert_allclose(output[[0, 2]], [13.5, 33.5])


def test_axis_mapping_uses_last_source_duplicate() -> None:
    targets, sources = BUILD.axis_mapping(
        ["A", "B", "A", "X"], ["B", "A", "C"]
    )
    np.testing.assert_array_equal(targets, [0, 1])
    np.testing.assert_array_equal(sources, [1, 2])


def write_agg(path: Path, values: dict[str, float]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["statistic", *COMPARE.METRICS])
        writer.writeheader()
        writer.writerow({"statistic": "mean", **values})


def test_score_comparison_has_all_six_raw_and_scaled_deltas() -> None:
    baseline = {
        "pds_cosine": 0.7,
        "expr_mse_unbiased_capped_norm": 1.1,
        "de_wilcoxon_direction_fidelity_yield_raw": 0.5,
        "de_wilcoxon_direction_reach_raw": 0.1,
        "de_wilcoxon_sig_jaccard": 0.03,
        "de_wilcoxon_lfc_nmae": 1.0,
    }
    candidate = dict(baseline)
    candidate["pds_cosine"] += 0.01
    candidate["expr_mse_unbiased_capped_norm"] -= 0.02
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        baseline_path = root / "baseline.csv"
        candidate_path = root / "candidate.csv"
        write_agg(baseline_path, baseline)
        write_agg(candidate_path, candidate)
        result = COMPARE.compare(baseline_path, candidate_path)
        assert set(result["metrics"]) == set(COMPARE.METRICS)
        assert np.isclose(
            result["metrics"]["pds_cosine"]["raw_candidate_minus_baseline"], 0.01
        )
        assert (
            result["metrics"]["expr_mse_unbiased_capped_norm"][
                "raw_favorable_delta"
            ]
            > 0
        )
        assert "candidate_minus_baseline" in result["avg6"]
        json.dumps(result)


def test_all_panel_coordinates_can_be_removed() -> None:
    official = ["P0", "G0", "P1", "G1", "G2"]
    panel = {"P0", "P1"}
    donor_union = {"P0", "P1", "G0", "G1", "G2"}
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in donor_union and gene not in panel
        ]
    )
    assert [official[index] for index in wcols] == ["G0", "G1", "G2"]


if __name__ == "__main__":
    test_blend_preserves_missing_gene_coordinates()
    test_axis_mapping_uses_last_source_duplicate()
    test_score_comparison_has_all_six_raw_and_scaled_deltas()
    test_all_panel_coordinates_can_be_removed()
    print(f"4 synthetic tests passed; peak RSS={BUILD.peak_rss_gib():.3f} GiB")
