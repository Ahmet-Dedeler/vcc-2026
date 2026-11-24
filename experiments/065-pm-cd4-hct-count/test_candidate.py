"""Lightweight synthetic checks for experiment 065."""

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


BUILD = load_module("exp065_build", "build_joint_candidate.py")
COMPARE = load_module("exp065_compare", "compare_scores.py")


def test_frozen_constants_match_sources() -> None:
    BUILD.assert_frozen_contract()
    assert BUILD.CD4_MODE == "gene_norm"
    assert BUILD.CD4_FRACTION == 0.05


def test_post_whiten_cd4_uses_exact_global_preserve_and_excludes_panel() -> None:
    pm = np.asarray(
        [
            [4.0, 3.0, 0.0, 2.0],
            [1.0, 0.0, 6.0, 2.0],
        ]
    )
    cd4 = BUILD.EXP062.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray(
            [
                [0.0, 5.0, 0.0, 7.0],
                [0.0, 2.0, 3.0, 0.0],
            ]
        ),
        gene_votes=np.asarray(
            [
                [0, 3, 0, 2],
                [0, 1, 3, 0],
            ],
            dtype=np.uint8,
        ),
        target_condition_count=np.asarray([3, 3], dtype=np.uint8),
        condition_details=[],
    )
    wcols = np.asarray([1, 2, 3], dtype=np.int64)
    panel_positions = np.asarray([0], dtype=np.int64)
    actual, metadata = BUILD.add_frozen_cd4_residual(
        pm, cd4, wcols, panel_positions
    )
    expected, expected_metadata = BUILD.EXP062.add_post_whiten_residual(
        pm, cd4, 0.05, wcols
    )
    np.testing.assert_array_equal(actual, expected)
    assert metadata == expected_metadata
    scale = metadata["post_add_global_preserve_scale"]
    np.testing.assert_allclose(actual[:, 0], pm[:, 0] * scale, rtol=0, atol=0)
    np.testing.assert_allclose(actual[0, 2], pm[0, 2] * scale, rtol=0, atol=0)


def test_panel_residual_fails_closed() -> None:
    pm = np.ones((1, 3))
    cd4 = BUILD.EXP062.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray([[1.0, 0.0, 0.0]]),
        gene_votes=np.asarray([[1, 0, 0]], dtype=np.uint8),
        target_condition_count=np.asarray([1], dtype=np.uint8),
        condition_details=[],
    )
    try:
        BUILD.add_frozen_cd4_residual(
            pm,
            cd4,
            np.asarray([1, 2], dtype=np.int64),
            np.asarray([0], dtype=np.int64),
        )
    except AssertionError as error:
        assert "VCC target coordinate" in str(error)
    else:
        raise AssertionError("panel-supported CD4 residual was accepted")


def test_incomplete_cd4_targets_are_not_extended_beyond_mean3() -> None:
    prepared = BUILD.EXP062.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
        gene_votes=np.asarray([[3, 3], [2, 2]], dtype=np.uint8),
        target_condition_count=np.asarray([3, 2], dtype=np.uint8),
        condition_details=[],
    )
    frozen = BUILD.freeze_complete_mean3(prepared)
    np.testing.assert_array_equal(frozen.mean_effect[0], [1.0, 2.0])
    np.testing.assert_array_equal(frozen.mean_effect[1], [0.0, 0.0])
    np.testing.assert_array_equal(frozen.gene_votes[1], [0, 0])


class FakeModel:
    names = [BUILD.EXP061.K562, BUILD.EXP061.HEK]
    target_genes = np.asarray(["G0", "G1", "G2", "G3"])
    effect_scale = 2.0
    _scale = {BUILD.EXP061.K562: 1.0}
    tcols = {BUILD.EXP061.K562: np.arange(4, dtype=np.int64)}
    dcols = {BUILD.EXP061.K562: np.arange(4, dtype=np.int64)}
    pert_index = {
        BUILD.EXP061.K562: {"T0": 0, "T1": 1},
        BUILD.EXP061.HEK: {"T0": 0, "T1": 1},
    }
    s = {
        BUILD.EXP061.K562: np.asarray(
            [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]
        )
    }

    def effect(self, target, _alpha, _shrink, *, weights):
        assert weights == {BUILD.EXP061.K562: 1.0, BUILD.EXP061.HEK: 1.0}
        return np.asarray([10.0, 20.0, 30.0, 40.0]) + (
            1.0 if target == "T1" else 0.0
        )


def test_exp061_missing_support_behavior_is_used_verbatim() -> None:
    pm = {
        "scale": 1.0,
        "target_columns": np.asarray([0, 2], dtype=np.int64),
        "source_columns": np.asarray([0, 1], dtype=np.int64),
        "mean": np.asarray([0.0, 0.0]),
        "delta": np.asarray([[8.0, 10.0]]),
        "pert_index": {"T0": 0},
    }
    baseline, candidate, audit = BUILD.EXP061.compute_effects(
        FakeModel(), pm, ("T0", "T1"), max_rss_gib=0.95
    )
    np.testing.assert_array_equal(candidate[1], baseline[1])
    np.testing.assert_array_equal(candidate[0, [1, 3]], baseline[0, [1, 3]])
    assert audit["missing_perturbmatch_targets"] == ["T1"]
    assert audit["genes_preserved_from_incumbent_slot"] == 2


def write_agg(path: Path, values: dict[str, float]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["statistic", *COMPARE.METRICS])
        writer.writeheader()
        writer.writerow({"statistic": "mean", **values})


def test_score_comparison_reports_all_six_metrics() -> None:
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
        assert len(result["metrics"]) == 6
        assert np.isclose(
            result["metrics"]["pds_cosine"]["raw_candidate_minus_baseline"], 0.01
        )
        assert (
            result["metrics"]["expr_mse_unbiased_capped_norm"][
                "raw_favorable_delta"
            ]
            > 0
        )
        json.dumps(result)


def test_cloud_runner_is_serial_and_exp061_gated() -> None:
    launcher = (HERE / "aws/launch_reuse_exp061.sh").read_text()
    runner = (HERE / "aws/run.sh").read_text()

    # The local launcher has permission to read output/* and refuses to dispatch
    # unless experiment 061 uploaded its result, aggregate, and zero exit code.
    assert "EXP061_OUTPUT_PREFIX" in launcher
    assert "${EXP061_OUTPUT_PREFIX}/results.json" in launcher
    assert "${EXP061_OUTPUT_PREFIX}/candidate-run/agg_results.csv" in launcher
    assert "${EXP061_OUTPUT_PREFIX}/exit-code.txt" in launcher
    assert 'if [[ "$exp061_exit" != "0" ]]' in launcher

    # The worker role cannot read output/*, so the cloud runner verifies the
    # already-local experiment-061 exit code and artifacts before taking locks.
    assert "EXP061_EXIT_FILE" in runner
    assert 'if [[ "$exp061_exit" != "0" ]]' in runner
    assert "cloud-output/results.json" in runner
    assert "cloud-output/candidate-run/agg_results.csv" in runner
    assert "flock -n 8" in runner
    assert "build_joint_candidate.py" in runner
    assert runner.index("build_joint_candidate.py") < runner.index("cell-eval2 run")
    assert runner.count("cell-eval2 run") == 1
    assert "build-memory-guard.json" in runner
    assert "score-memory-guard.json" in runner
    assert "shutdown -h" not in runner


if __name__ == "__main__":
    test_frozen_constants_match_sources()
    test_post_whiten_cd4_uses_exact_global_preserve_and_excludes_panel()
    test_panel_residual_fails_closed()
    test_incomplete_cd4_targets_are_not_extended_beyond_mean3()
    test_exp061_missing_support_behavior_is_used_verbatim()
    test_score_comparison_reports_all_six_metrics()
    test_cloud_runner_is_serial_and_exp061_gated()
    print(
        "7 synthetic tests passed; "
        f"peak RSS={BUILD.EXP061.peak_rss_gib():.3f} GiB"
    )
