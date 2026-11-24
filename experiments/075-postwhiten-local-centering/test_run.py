"""Bounded synthetic tests for experiment 075.

These tests never open real donor matrices or frozen whitening operators.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp075_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def write_csv(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[column])
        writer.writeheader()
        writer.writerows([{column: value} for value in values])


def test_transform_is_query_only_permutation_equivariant_and_norm_preserving() -> None:
    rng = np.random.default_rng(75)
    query = rng.normal(size=(24, 17))
    query[0] = 0.0
    transformed, diagnostics = RUN.apply_reciprocal_local_centering(
        query, k=6, beta=0.4
    )
    np.testing.assert_allclose(
        np.linalg.norm(transformed, axis=1),
        np.linalg.norm(query, axis=1),
        atol=1e-10,
    )
    assert diagnostics["max_row_norm_error"] < 1e-10
    assert np.all(transformed[0] == 0.0)

    permutation = rng.permutation(len(query))
    inverse = np.argsort(permutation)
    permuted, _ = RUN.apply_reciprocal_local_centering(
        query[permutation], k=6, beta=0.4
    )
    np.testing.assert_allclose(transformed, permuted[inverse], atol=1e-10)


def test_fast_score_formula_matches_materialized_transform() -> None:
    rng = np.random.default_rng(751)
    query = rng.normal(size=(30, 19))
    truth = rng.normal(size=(30, 19))
    candidate = RUN.Candidate(k=8, beta=0.5)
    _, fast, _ = RUN.candidate_score_matrices(
        query,
        truth,
        [candidate],
        minimum_mutual=2,
    )
    transformed, _ = RUN.apply_reciprocal_local_centering(
        query, k=candidate.k, beta=candidate.beta
    )
    explicit = RUN.row_normalize(transformed) @ RUN.row_normalize(truth).T
    np.testing.assert_allclose(fast[candidate.key], explicit, atol=1e-10)


def test_robust_selector_requires_every_validation_cell_nonnegative() -> None:
    candidates = [
        RUN.Candidate(k=8, beta=0.25),
        RUN.Candidate(k=16, beta=0.25),
    ]
    cells = [
        {
            "deltas": {
                candidates[0].key: 0.004,
                candidates[1].key: 0.003,
            }
        },
        {
            "deltas": {
                candidates[0].key: -0.0001,
                candidates[1].key: 0.002,
            }
        },
    ]
    selected, audit = RUN.select_robust_candidate(
        cells, candidates, minimum_mean_gain=0.001
    )
    assert selected == candidates[1].key
    assert audit["eligible_count"] == 1


def make_fixture(root: Path) -> argparse.Namespace:
    train_perts = [f"T{i:03d}" for i in range(600)]
    vcc_panel = [f"P{i:03d}" for i in range(300)]
    outer_perts = vcc_panel[:240]
    active = [f"G{i:02d}" for i in range(18)]
    score_only = [f"S{i:02d}" for i in range(4)]
    score_genes = [*active, *score_only]
    master = [*train_perts, *score_genes]
    official = [*vcc_panel, *master]
    train_mask = np.asarray(
        [False] * len(train_perts) + [True] * len(active) + [False] * len(score_only),
        dtype=bool,
    )

    genes_csv = root / "gene_names.csv"
    panel_csv = root / "pert_counts.csv"
    write_csv(genes_csv, "gene_name", official)
    write_csv(panel_csv, "target_gene", vcc_panel)

    truth_cache = root / "truth-cache"
    truth_cache.mkdir()
    (truth_cache / "master_genes.txt").write_text("\n".join(master) + "\n")
    (truth_cache / "train_perts.txt").write_text("\n".join(train_perts) + "\n")
    (truth_cache / "outer_perts.txt").write_text("\n".join(outer_perts) + "\n")
    np.save(truth_cache / "train_mask.npy", train_mask)
    (truth_cache / "meta.json").write_text('{"fixture": true}\n')

    rng = np.random.default_rng(752)
    train_signal = rng.normal(size=(len(train_perts), len(active)))
    outer_signal = rng.normal(size=(len(outer_perts), len(score_genes)))
    for index, donor in enumerate(RUN.AUX_DONORS):
        values = np.zeros((len(train_perts), len(master)), dtype=np.float32)
        active_start = len(master) - len(score_genes)
        values[:, active_start : active_start + len(active)] = (
            train_signal
            + rng.normal(
                scale=0.10 + index * 0.02,
                size=train_signal.shape,
            )
        )
        np.save(truth_cache / f"{donor}_raw_train.npy", values)
    for index, donor in enumerate(RUN.MAIN_DONORS):
        values = np.zeros((len(outer_perts), len(master)), dtype=np.float32)
        values[:, -len(score_genes) :] = outer_signal + rng.normal(
            scale=0.11 + index * 0.02,
            size=outer_signal.shape,
        )
        np.save(truth_cache / f"{donor}_raw_outer.npy", values)

    query_cache = root / "query-cache"
    query_cache.mkdir()
    folds: dict[str, object] = {}
    operators: dict[str, object] = {}
    for index, held in enumerate(RUN.MAIN_DONORS):
        fold_dir = query_cache / held
        fold_dir.mkdir()
        train_query = (
            np.column_stack(
                [
                    train_signal
                    + rng.normal(
                        scale=0.15 + index * 0.01,
                        size=train_signal.shape,
                    ),
                    rng.normal(size=(len(train_perts), len(score_only))),
                ]
            )
        ).astype(np.float32)
        outer_query = (
            outer_signal
            + rng.normal(scale=0.15 + index * 0.01, size=outer_signal.shape)
        ).astype(np.float32)
        train_path = fold_dir / "train_query.npy"
        outer_path = fold_dir / "outer_query.npy"
        np.save(train_path, train_query)
        np.save(outer_path, outer_query)
        sources = [donor for donor in RUN.MAIN_DONORS if donor != held]
        operators[held] = {
            "source_donors": sources,
            "npz_sha256": f"{held}-npz",
            "json_sha256": f"{held}-json",
        }
        folds[held] = {
            "source_donors": sources,
            "operator_npz_sha256": f"{held}-npz",
            "operator_json_sha256": f"{held}-json",
            "panels": {
                "train": {"query_sha256": RUN.sha256_file(train_path)},
                "outer": {"query_sha256": RUN.sha256_file(outer_path)},
            },
        }
    metadata = {
        "schema_version": RUN.QUERY_CACHE_SCHEMA,
        "experiment": RUN.EXPECTED_UPSTREAM_EXPERIMENT,
        "gamma": RUN.EXPECTED_GAMMA,
        "training_targets": len(train_perts),
        "outer_targets": len(outer_perts),
        "train_outer_overlap": 0,
        "outer_score_genes": len(score_genes),
        "correction_genes": len(active),
        "vcc_target_coordinate_overlap": 0,
        "training_target_coordinate_overlap": 0,
        "whitening_or_full_covariance_eigendecompositions": 0,
        "train_order_sha256": RUN.sha256_lines(train_perts),
        "outer_order_sha256": RUN.sha256_lines(outer_perts),
        "score_gene_order_sha256": RUN.sha256_lines(score_genes),
        "correction_gene_order_sha256": RUN.sha256_lines(active),
        "correction_in_score_sha256": RUN.sha256_array(
            np.arange(len(active), dtype=np.int64)
        ),
        "folds": folds,
    }
    (query_cache / "metadata.json").write_text(json.dumps(metadata))

    truth_files = {
        path.name: RUN.sha256_file(path)
        for path in truth_cache.iterdir()
        if path.is_file()
    }
    upstream_manifest = {
        "schema_version": 1,
        "challenge": {
            "gene_names.csv": RUN.sha256_file(genes_csv),
            "pert_counts.csv": RUN.sha256_file(panel_csv),
        },
        "truth_cache": {"files": truth_files},
        "operators": operators,
    }
    manifest_path = root / "upstream_manifest.json"
    manifest_path.write_text(json.dumps(upstream_manifest))
    return argparse.Namespace(
        truth_cache=truth_cache,
        query_cache=query_cache,
        upstream_manifest=manifest_path,
        genes_csv=genes_csv,
        panel_csv=panel_csv,
        output=root / "results.json",
        per_target_output=root / "per_target.csv",
        ks=[4, 8],
        betas=[0.10, 0.25],
        minimum_mutual=2,
        minimum_inner_gain=0.0,
        bootstrap_replicates=200,
        bootstrap_seed=752,
        max_rss_gib=0.95,
        guarded=True,
    )


def test_end_to_end_protocol_and_leakage_audit() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = make_fixture(Path(directory))
        manifest = json.loads(args.upstream_manifest.read_text())
        RUN.verify_inputs(args, manifest)
        result = RUN.run_experiment(args, manifest)
        assert result["protocol"]["train_outer_overlap"] == 0
        assert result["method"]["query_only"] is True
        assert result["method"]["hidden_gallery_used"] is False
        assert result["leakage_audit"]["outer_truth_used_in_selection"] is False
        assert set(result["outer"]["folds"]) == set(RUN.MAIN_DONORS)
        assert len(result["inner"]["role_leave_one_out"]) == 3
        assert result["promotion_gate"]["predeclared"] is True
        assert args.output.exists()
        assert args.per_target_output.exists()


def test_cache_rejects_train_outer_overlap() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = make_fixture(Path(directory))
        outer_path = args.truth_cache / "outer_perts.txt"
        outer = list(RUN.read_lines(outer_path))
        outer[0] = "T000"
        outer_path.write_text("\n".join(outer) + "\n")
        try:
            RUN.validate_truth_cache(
                args.truth_cache,
                RUN.read_csv_column(args.genes_csv, "gene_name"),
                RUN.read_csv_column(args.panel_csv, "target_gene"),
            )
        except ValueError as error:
            assert "VCC-panel subset" in str(error) or "overlap" in str(error)
        else:
            raise AssertionError("train/outer overlap was not rejected")


def main() -> int:
    tests = [
        test_transform_is_query_only_permutation_equivariant_and_norm_preserving,
        test_fast_score_formula_matches_materialized_transform,
        test_robust_selector_requires_every_validation_cell_nonnegative,
        test_end_to_end_protocol_and_leakage_audit,
        test_cache_rejects_train_outer_overlap,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
