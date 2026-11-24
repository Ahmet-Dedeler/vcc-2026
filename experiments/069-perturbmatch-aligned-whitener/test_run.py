#!/usr/bin/env python3
"""Low-memory synthetic tests for experiment 069."""

from __future__ import annotations

import csv
import hashlib
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


RUN = load_module("exp069_run", "run.py")
GATE = load_module("exp069_gate", "validate_exp061_gate.py")


def prepared(
    name: str,
    effects: np.ndarray,
    available: np.ndarray,
) -> object:
    return RUN.BASE.PreparedSource(
        name=name,
        effects=np.asarray(effects, dtype=np.float64),
        target_available=np.asarray(available, dtype=bool),
        mean_on_axis=np.zeros(effects.shape[1], dtype=np.float64),
        raw_median_norm=1.0,
        mapped_gene_count=effects.shape[1],
        duplicate_symbols=(),
    )


def write_bank(
    root: Path,
    name: str,
    genes: tuple[str, ...],
    perts: tuple[str, ...],
    values: np.ndarray,
) -> object:
    path = root / name
    path.mkdir()
    delta_path = path / "delta.npy"
    np.save(delta_path, np.asarray(values, dtype=np.float32))
    return RUN.BASE.Bank(
        name=name,
        path=path,
        genes=genes,
        perts=perts,
        pert_pos={perturbation: index for index, perturbation in enumerate(perts)},
        delta_path=delta_path,
        b_ctrl_path=path / "b_ctrl.npy",
        shape=values.shape,
    )


def test_frozen_constants() -> None:
    assert RUN.FOLDS == ("orion_hct116", "orion_hek293t")
    assert RUN.PM_FRACTION == 0.50
    assert RUN.CD4_MODE == "gene_norm"
    assert RUN.CD4_FRACTION == 0.05
    assert RUN.WHITEN_GAMMA == 0.65
    assert RUN.MACRO_GATE == 0.003


def test_pm_slot_target_and_gene_fallback() -> None:
    incumbent = prepared(
        "k562",
        np.asarray([[2.0, 4.0, 6.0, 8.0], [1.0, 3.0, 5.0, 7.0]]),
        np.asarray([True, True]),
    )
    pm = prepared(
        "pm",
        np.asarray([[10.0, 20.0, 30.0, 40.0], [9.0, 9.0, 9.0, 9.0]]),
        np.asarray([True, False]),
    )
    slot, audit = RUN.build_pm_slot_with_fallback(
        incumbent, pm, np.asarray([0, 2], dtype=np.int64)
    )
    np.testing.assert_array_equal(slot.effects[0], [6.0, 4.0, 18.0, 8.0])
    np.testing.assert_array_equal(slot.effects[1], incumbent.effects[1])
    np.testing.assert_array_equal(slot.target_available, [True, True])
    assert audit["targets_blended"] == 1
    assert audit["targets_incumbent_fallback"] == 1


def test_cd4_requires_all_three_conditions() -> None:
    cd4 = RUN.BASE.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray([[1.0, 2.0], [3.0, 4.0]]),
        gene_votes=np.asarray([[3, 3], [2, 2]], dtype=np.uint8),
        target_condition_count=np.asarray([3, 2], dtype=np.uint8),
        condition_details=[],
    )
    frozen = RUN.freeze_complete_cd4_mean3(cd4)
    np.testing.assert_array_equal(frozen.mean_effect[0], [1.0, 2.0])
    np.testing.assert_array_equal(frozen.mean_effect[1], [0.0, 0.0])
    np.testing.assert_array_equal(frozen.gene_votes[1], [0, 0])


def test_equal_source_mass_weights() -> None:
    k_rows, companion_rows = 7, 19
    k_weight = RUN.source_mass_row_weight(k_rows)
    companion_weight = RUN.source_mass_row_weight(companion_rows)
    assert np.isclose(k_rows * k_weight**2, 0.5)
    assert np.isclose(companion_rows * companion_weight**2, 0.5)


def test_aligned_whitener_excludes_exact_panel_and_uses_equal_mass() -> None:
    panel = tuple(f"P{index:03d}" for index in range(300))
    official = (*panel, "G0", "G1")
    genes = official
    incumbent_perts = (*panel, "X", "Y", "Z")
    perturbmatch_perts = (*panel, "X", "Y")
    companion_perts = (*panel, "U", "V", "W")

    def values(perts: tuple[str, ...], background: list[list[float]]) -> np.ndarray:
        matrix = np.zeros((len(perts), len(genes)), dtype=np.float32)
        matrix[: len(panel), -2:] = 1_000_000.0
        matrix[len(panel) :, -2:] = np.asarray(background, dtype=np.float32)
        matrix[:, : len(panel)] = 2_000_000.0
        return matrix

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        incumbent = write_bank(
            root,
            "k562",
            genes,
            incumbent_perts,
            values(incumbent_perts, [[2.0, 0.0], [-1.0, 1.0], [-1.0, -1.0]]),
        )
        perturbmatch = write_bank(
            root,
            "perturbmatch",
            genes,
            perturbmatch_perts,
            values(perturbmatch_perts, [[0.0, 2.0], [0.0, -2.0]]),
        )
        companion = write_bank(
            root,
            "companion",
            genes,
            companion_perts,
            values(companion_perts, [[1.0, 1.0], [-1.0, 1.0], [0.0, -2.0]]),
        )
        output = root / "whitener.npz"
        operator, wcols = RUN.fit_aligned_whitener(
            incumbent,
            perturbmatch,
            companion,
            official,
            panel,
            output,
            root / "scratch",
            chunk_rows=1,
            max_rss_gib=1.0,
        )

        np.testing.assert_array_equal(wcols, [300, 301])
        assert operator.evecs.shape == (2, 2)
        metadata = operator.metadata
        assert metadata["excluded_panel_rows_by_source"] == {
            "k562": 300,
            "perturbmatch": 300,
            "companion": 300,
        }
        assert metadata["basis_rows_by_source"] == {
            "k562_gwps_pm50_fallback": 3,
            "companion": 3,
        }
        assert metadata["k562_background_rows_blended"] == 2
        assert metadata["k562_background_rows_incumbent_fallback"] == 1
        assert metadata["source_mass"] == {
            "k562_gwps_pm50_fallback": 0.5,
            "companion": 0.5,
        }
        assert metadata["panel_rows_excluded"] is True
        assert metadata["panel_coordinates_excluded"] is True

        panel_set = set(panel)
        incumbent_stats = RUN.background_stats(
            incumbent,
            official,
            wcols,
            panel_set,
            chunk_rows=1,
            max_rss_gib=1.0,
        )
        perturbmatch_stats = RUN.background_stats(
            perturbmatch,
            official,
            wcols,
            panel_set,
            chunk_rows=1,
            max_rss_gib=1.0,
        )
        companion_stats = RUN.background_stats(
            companion,
            official,
            wcols,
            panel_set,
            chunk_rows=1,
            max_rss_gib=1.0,
        )
        ensemble_stats = RUN.ensemble_background_stats(
            incumbent,
            perturbmatch,
            incumbent_stats,
            perturbmatch_stats,
            official,
            wcols,
            panel_set,
            chunk_rows=1,
            max_rss_gib=1.0,
        )
        k562_block, matched = RUN.blend_background_block(
            incumbent,
            perturbmatch,
            incumbent_stats,
            perturbmatch_stats,
            incumbent_stats.rows,
            official,
            wcols,
            panel_set,
        )
        assert matched == 2
        k562_block -= ensemble_stats.mean
        k562_block *= RUN.BASE.EFFECT_SCALE / ensemble_stats.median_norm
        companion_block = RUN.normalized_source_block(
            companion,
            companion_stats,
            companion_stats.rows,
            len(wcols),
        )
        np.testing.assert_allclose(k562_block.mean(axis=0), 0.0, atol=1e-12)
        np.testing.assert_allclose(companion_block.mean(axis=0), 0.0, atol=1e-12)
        assert np.isclose(
            np.median(np.linalg.norm(k562_block, axis=1)),
            RUN.BASE.EFFECT_SCALE,
        )
        assert np.isclose(
            np.median(np.linalg.norm(companion_block, axis=1)),
            RUN.BASE.EFFECT_SCALE,
        )
        expected_covariance = (
            k562_block * RUN.source_mass_row_weight(len(incumbent_stats.rows))
        ).T @ (k562_block * RUN.source_mass_row_weight(len(incumbent_stats.rows))) + (
            companion_block * RUN.source_mass_row_weight(len(companion_stats.rows))
        ).T @ (companion_block * RUN.source_mass_row_weight(len(companion_stats.rows)))
        reconstructed_covariance = (
            operator.evecs.astype(np.float64)
            @ np.diag(operator.evals)
            @ operator.evecs.astype(np.float64).T
        )
        np.testing.assert_allclose(
            reconstructed_covariance,
            expected_covariance,
            rtol=1e-6,
            atol=1e-6,
        )


def test_panel_rows_are_excluded() -> None:
    class Bank:
        name = "tiny"
        perts = ("A", "PANEL", "B")

    rows = RUN.panel_background_rows(Bank(), {"PANEL"})
    np.testing.assert_array_equal(rows, [0, 2])


def test_stopping_gate_is_conjunctive() -> None:
    passed = RUN.stopping_gate(
        {"orion_hct116": 0.004, "orion_hek293t": 0.003},
        0.0035,
        {"ci95_percentile": [0.0001, 0.007]},
    )
    assert passed["passes"] is True
    failed = RUN.stopping_gate(
        {"orion_hct116": 0.004, "orion_hek293t": -0.0001},
        0.0035,
        {"ci95_percentile": [0.0001, 0.007]},
    )
    assert failed["passes"] is False


def test_bootstrap_matches_macro_delta() -> None:
    deltas = {
        "orion_hct116": np.asarray([0.01, 0.02, -0.01]),
        "orion_hek293t": np.asarray([0.03, 0.00, 0.01]),
    }
    result = RUN.paired_target_cluster_bootstrap(deltas, reps=2_000, seed=69)
    expected = np.mean(
        [
            np.mean(deltas["orion_hct116"]),
            np.mean(deltas["orion_hek293t"]),
        ]
    )
    assert np.isclose(result["estimate"], expected)
    assert result["targets"] == 3


def write_score(path: Path, marker: str) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["statistic", "marker"])
        writer.writerow(["mean", marker])


def test_exp061_gate_checks_avg6_and_hashes() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        baseline = root / "baseline.csv"
        candidate = root / "candidate.csv"
        write_score(baseline, "base")
        write_score(candidate, "candidate")
        build = root / "build.json"
        candidate_hash = "a" * 64
        build.write_text(
            json.dumps(
                {
                    "experiment": "061-perturbmatch-hct-count",
                    "status": "candidate_built",
                    "candidate_h5ad": {
                        "sha256": candidate_hash,
                        "bytes": 123,
                    },
                }
            )
        )
        checksum = root / "candidate.h5ad.sha256"
        checksum.write_text(f"{candidate_hash}  candidate.h5ad\n")
        exit_code = root / "exit-code.txt"
        exit_code.write_text("0\n")
        results = root / "results.json"
        payload = {
            "experiment": "061-perturbmatch-hct-count",
            "status": "scored",
            "avg6": {"candidate_minus_baseline": 0.001},
            "score_files": {
                "baseline": {"sha256": RUN.sha256_file(baseline)},
                "candidate": {"sha256": RUN.sha256_file(candidate)},
            },
        }
        results.write_text(json.dumps(payload))
        accepted = GATE.validate_gate(
            results_path=results,
            baseline_agg_path=baseline,
            candidate_agg_path=candidate,
            build_report_path=build,
            candidate_checksum_path=checksum,
            exit_code_path=exit_code,
            candidate_object_bytes=123,
        )
        assert accepted["status"] == "pass"
        payload["avg6"]["candidate_minus_baseline"] = -1e-9
        results.write_text(json.dumps(payload))
        try:
            GATE.validate_gate(
                results_path=results,
                baseline_agg_path=baseline,
                candidate_agg_path=candidate,
                build_report_path=build,
                candidate_checksum_path=checksum,
                exit_code_path=exit_code,
                candidate_object_bytes=123,
            )
        except ValueError as error:
            assert "avg6" in str(error)
        else:
            raise AssertionError("negative experiment-061 avg6 passed the gate")


def test_exp061_gate_rejects_corrupt_score_hash() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        baseline = root / "baseline.csv"
        candidate = root / "candidate.csv"
        write_score(baseline, "base")
        write_score(candidate, "candidate")
        candidate_hash = "b" * 64
        (root / "build.json").write_text(
            json.dumps(
                {
                    "experiment": "061-perturbmatch-hct-count",
                    "status": "candidate_built",
                    "candidate_h5ad": {
                        "sha256": candidate_hash,
                        "bytes": 321,
                    },
                }
            )
        )
        (root / "candidate.h5ad.sha256").write_text(
            f"{candidate_hash}  candidate.h5ad\n"
        )
        (root / "exit-code.txt").write_text("0\n")
        (root / "results.json").write_text(
            json.dumps(
                {
                    "experiment": "061-perturbmatch-hct-count",
                    "status": "scored",
                    "avg6": {"candidate_minus_baseline": 0.0},
                    "score_files": {
                        "baseline": {"sha256": "0" * 64},
                        "candidate": {
                            "sha256": hashlib.sha256(candidate.read_bytes()).hexdigest()
                        },
                    },
                }
            )
        )
        try:
            GATE.validate_gate(
                results_path=root / "results.json",
                baseline_agg_path=baseline,
                candidate_agg_path=candidate,
                build_report_path=root / "build.json",
                candidate_checksum_path=root / "candidate.h5ad.sha256",
                exit_code_path=root / "exit-code.txt",
                candidate_object_bytes=321,
            )
        except ValueError as error:
            assert "baseline aggregate hash mismatch" in str(error)
        else:
            raise AssertionError("corrupt score hash passed the gate")


def test_aws_launcher_and_runner_are_guarded() -> None:
    launcher = (HERE / "aws/launch.sh").read_text()
    runner = (HERE / "aws/run.sh").read_text()
    assert 'I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES"' in launcher
    assert "validate_exp061_gate.py" in launcher
    assert launcher.index("validate_exp061_gate.py") < launcher.index(
        "aws ssm send-command"
    )
    assert "candidate-object-bytes" in launcher
    assert "experiment 061 prerequisite failed" in launcher
    assert "validate_exp061_gate.py" in runner
    assert "memory_guard.py" in runner
    assert "--max-rss-gib 1.0" in runner
    assert "--max-rss-gib 55" in runner
    assert "cell-eval2" not in runner
    assert "vcc submit" not in launcher
    assert "vcc submit" not in runner
    assert "run-instances" not in launcher
    assert "start-instances" not in launcher


if __name__ == "__main__":
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
