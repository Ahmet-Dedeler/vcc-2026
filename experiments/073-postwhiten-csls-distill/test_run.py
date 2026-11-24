"""Synthetic tests for experiment 073.

No real donor matrix or 1.2 GiB experiment-062 operator is opened.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp073_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def write_csv(path: Path, column: str, values: list[str]) -> None:
    path.write_text(column + "\n" + "\n".join(values) + "\n")


def test_density_and_selection_primitives() -> None:
    scores = np.asarray(
        [
            [0.8, 0.7, 0.1],
            [0.6, 0.9, 0.2],
            [0.3, 0.4, 0.7],
        ],
        dtype=np.float64,
    )
    density = RUN.candidate_density(scores, 2)
    np.testing.assert_allclose(density, [0.7, 0.8, 0.45])
    assert RUN.pds_from_scores(scores) == 1.0

    candidates = [
        RUN.Candidate("bias", 5, 8, 0.1, 1.0),
        RUN.Candidate("bands", 5, 16, 0.1, 1.0),
    ]
    selected = RUN.choose_candidate(
        {
            candidates[0].key: 0.7014,
            candidates[1].key: 0.7018,
        },
        baseline=0.7000,
        candidates=candidates,
        min_gain=0.001,
    )
    # The two are within the 5e-4 tie window, so lower rank wins.
    assert selected == candidates[0].key
    assert (
        RUN.choose_candidate(
            {candidate.key: 0.7005 for candidate in candidates},
            baseline=0.7000,
            candidates=candidates,
            min_gain=0.001,
        )
        == "baseline"
    )


def test_frozen_operator_masks_both_sides() -> None:
    query = np.asarray(
        [
            [100.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            [-100.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        ],
        dtype=np.float32,
    )
    query[:, 0] = 0.0
    operator = RUN.Operator(
        evals=np.ones(6, dtype=np.float64),
        evecs=np.eye(6, dtype=np.float32),
        metadata={},
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        output_path = Path(directory) / "query.npy"
        metadata = RUN.apply_frozen_operator(
            query,
            operator,
            blocked_operator_positions=np.asarray([0], dtype=np.int64),
            selected_operator_positions=np.asarray([2, 4, 5], dtype=np.int64),
            output_path=output_path,
            transform_chunk_rows=1,
            max_rss_gib=0.95,
        )
        output = np.load(output_path)
    np.testing.assert_allclose(output, query[:, [2, 4, 5]])
    np.testing.assert_allclose(
        metadata["global_median_norm_preserve_scale"], 1.0, rtol=1e-6
    )


def test_dual_axis_scores_full_geometry_and_embeds_correction() -> None:
    query = np.eye(3, dtype=np.float64)
    truth = np.eye(3, dtype=np.float64)
    correction_in_score = np.asarray([0], dtype=np.int64)
    fitted = RUN.Fitted(
        basis=np.ones((1, 1), dtype=np.float64),
        bands=np.asarray([0.4], dtype=np.float64),
    )
    scores = RUN.fitted_score_matrix(
        query,
        truth,
        fitted,
        eta=0.5,
        correction_in_score=correction_in_score,
    )
    normalized_query = RUN.row_normalize(query)
    normalized_truth = RUN.row_normalize(truth)
    expected = normalized_query @ normalized_truth.T
    expected += 0.5 * 0.4 * (normalized_query[:, [0]] @ normalized_truth[:, [0]].T)
    np.testing.assert_allclose(scores, expected)
    assert RUN.pds_from_scores(scores) == 1.0
    correction_only = (
        normalized_query[:, correction_in_score]
        @ normalized_truth[:, correction_in_score].T
    )
    assert RUN.pds_from_scores(correction_only) < RUN.pds_from_scores(scores)


def make_distill_fixture(root: Path) -> argparse.Namespace:
    train_perts = [f"T{i:03d}" for i in range(600)]
    vcc_panel = [f"P{i:03d}" for i in range(300)]
    outer_perts = vcc_panel[:240]
    score_only_genes = [f"X{i:02d}" for i in range(3)]
    correction_genes = [f"G{i:02d}" for i in range(12)]
    score_genes = [*score_only_genes, *correction_genes]
    master_genes = [*train_perts, *score_genes]
    official = [*vcc_panel, *master_genes]
    train_mask = np.asarray(
        [False] * (len(train_perts) + len(score_only_genes))
        + [True] * len(correction_genes),
        dtype=bool,
    )

    truth_cache = root / "truth-cache"
    truth_cache.mkdir()
    (truth_cache / "master_genes.txt").write_text("\n".join(master_genes) + "\n")
    (truth_cache / "train_perts.txt").write_text("\n".join(train_perts) + "\n")
    (truth_cache / "outer_perts.txt").write_text("\n".join(outer_perts) + "\n")
    np.save(truth_cache / "train_mask.npy", train_mask)
    (truth_cache / "meta.json").write_text("{}\n")

    rng = np.random.default_rng(730)
    train_signal = rng.normal(size=(len(train_perts), len(correction_genes)))
    outer_signal = rng.normal(size=(len(outer_perts), len(score_genes)))
    train_full = np.zeros((len(train_perts), len(master_genes)), dtype=np.float32)
    outer_full = np.zeros((len(outer_perts), len(master_genes)), dtype=np.float32)
    for index, name in enumerate(RUN.AUX_DONORS):
        train_full[:] = 0.0
        train_full[:, : len(train_perts)] = 50.0
        train_full[:, -len(correction_genes) :] = train_signal + rng.normal(
            scale=0.08 + 0.01 * index, size=train_signal.shape
        )
        np.save(truth_cache / f"{name}_raw_train.npy", train_full)
    for index, name in enumerate(RUN.MAIN_DONORS):
        outer_full[:] = 0.0
        outer_full[:, : len(train_perts)] = 50.0
        outer_full[:, -len(score_genes) :] = outer_signal + rng.normal(
            scale=0.08 + 0.01 * index, size=outer_signal.shape
        )
        np.save(truth_cache / f"{name}_raw_outer.npy", outer_full)

    genes_csv = root / "gene_names.csv"
    panel_csv = root / "pert_counts.csv"
    write_csv(genes_csv, "gene_name", official)
    write_csv(panel_csv, "target_gene", vcc_panel)

    query_cache = root / "query-cache"
    query_cache.mkdir()
    folds: dict[str, object] = {}
    for index, name in enumerate(RUN.MAIN_DONORS):
        fold_dir = query_cache / name
        fold_dir.mkdir()
        train_query = np.column_stack(
            [
                rng.normal(size=(len(train_perts), len(score_only_genes))),
                train_signal,
            ]
        )
        train_query += rng.normal(scale=0.12 + 0.01 * index, size=train_query.shape)
        train_query = train_query.astype(np.float32)
        outer_query = (
            outer_signal
            + rng.normal(scale=0.12 + 0.01 * index, size=outer_signal.shape)
        ).astype(np.float32)
        train_path = fold_dir / "train_query.npy"
        outer_path = fold_dir / "outer_query.npy"
        np.save(train_path, train_query)
        np.save(outer_path, outer_query)
        folds[name] = {
            "panels": {
                "train": {"query_sha256": RUN.sha256_file(train_path)},
                "outer": {"query_sha256": RUN.sha256_file(outer_path)},
            }
        }
    metadata = {
        "schema_version": RUN.CACHE_SCHEMA,
        "training_targets": len(train_perts),
        "outer_targets": len(outer_perts),
        "outer_score_genes": len(score_genes),
        "correction_genes": len(correction_genes),
        "score_only_genes": len(score_only_genes),
        "train_outer_overlap": 0,
        "vcc_target_coordinate_overlap": 0,
        "training_target_coordinate_overlap": 0,
        "whitening_or_full_covariance_eigendecompositions": 0,
        "train_order_sha256": RUN.sha256_lines(train_perts),
        "outer_order_sha256": RUN.sha256_lines(outer_perts),
        "score_gene_order_sha256": RUN.sha256_lines(score_genes),
        "correction_gene_order_sha256": RUN.sha256_lines(correction_genes),
        "correction_in_score_sha256": RUN.sha256_array(
            np.arange(
                len(score_only_genes),
                len(score_genes),
                dtype=np.int64,
            )
        ),
        "folds": folds,
    }
    manifest_path = root / "artifact_manifest.json"
    operators: dict[str, object] = {}
    for name in RUN.MAIN_DONORS:
        fold = folds[name]
        fold["source_donors"] = [donor for donor in RUN.MAIN_DONORS if donor != name]
        fold["operator_npz_sha256"] = f"{name}-npz"
        fold["operator_json_sha256"] = f"{name}-json"
        operators[name] = {
            "source_donors": fold["source_donors"],
            "npz_sha256": fold["operator_npz_sha256"],
            "json_sha256": fold["operator_json_sha256"],
        }
    (query_cache / "metadata.json").write_text(json.dumps(metadata))
    manifest_path.write_text(json.dumps({"schema_version": 1, "operators": operators}))
    return argparse.Namespace(
        donor_root=root / "unused-donors",
        truth_cache=truth_cache,
        operator_dir=root / "unused-operators",
        artifact_manifest=manifest_path,
        genes_csv=genes_csv,
        panel_csv=panel_csv,
        query_cache=query_cache,
        output=root / "results.json",
        per_target_output=root / "per_target.csv",
        phase="distill",
        source_chunk_rows=16,
        transform_chunk_rows=4,
        max_rss_gib=0.95,
        methods=["bias", "bands", "combo"],
        ks=[5],
        ranks=[2],
        ridges=[0.1],
        etas=[0.5, 1.0],
        min_inner_gain=0.001,
        bootstrap_replicates=500,
        bootstrap_seed=7302026,
        guarded=True,
    )


def test_end_to_end_nested_distillation() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = make_distill_fixture(Path(directory))
        result = RUN.run_distillation(
            args, json.loads(args.artifact_manifest.read_text())
        )
        assert set(result["folds"]) == set(RUN.MAIN_DONORS)
        assert result["protocol"]["operator_refits"] == 0
        assert (
            result["protocol"]["whitening_or_full_covariance_eigendecompositions"] == 0
        )
        assert result["protocol"]["randomized_low_rank_distillation_svd"] is True
        assert result["protocol"]["train_outer_overlap"] == 0
        assert result["protocol"]["outer_score_genes"] == 15
        assert result["protocol"]["correction_genes"] == 12
        assert result["protocol"]["score_only_genes"] == 3
        audit = result["leakage_audit"]
        assert audit["training_and_outer_perturbations_disjoint"] is True
        assert audit["outer_donor_used_in_query"] is False
        assert audit["outer_donor_used_in_operator_fit"] is False
        assert audit["outer_donor_used_in_distillation_fit"] is False
        assert audit["outer_donor_used_in_hyperparameter_selection"] is False
        assert audit["inner_validation_teacher_used_in_inner_fit"] is False
        assert audit["vcc_target_coordinates_in_effective_map"] is False
        assert audit["training_target_coordinates_in_effective_map"] is False
        assert audit["outer_truth_transformed"] is False
        assert audit["outer_scores_used_to_select_candidate"] is False
        assert all(
            set(fold["inner"]["detail"]) == set(RUN.AUX_DONORS)
            for fold in result["folds"].values()
        )
        assert result["nested_macro"]["bootstrap"]["replicates"] == 500
        assert args.output.exists()
        assert args.per_target_output.exists()
        assert len(args.per_target_output.read_text().splitlines()) == 1 + 240 * len(
            RUN.MAIN_DONORS
        )
        assert RUN.peak_rss_gib() < 0.95


if __name__ == "__main__":
    test_density_and_selection_primitives()
    test_frozen_operator_masks_both_sides()
    test_dual_axis_scores_full_geometry_and_embeds_correction()
    test_end_to_end_nested_distillation()
    print(f"4 synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
