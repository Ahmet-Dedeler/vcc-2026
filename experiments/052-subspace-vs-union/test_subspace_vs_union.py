from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).with_name("run.py")
SPEC = importlib.util.spec_from_file_location("subspace_vs_union", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_pds_identity_is_perfect() -> None:
    rng = np.random.default_rng(1)
    truth = rng.normal(size=(30, 50))
    assert MODULE.pds(truth, truth) == 1.0


def test_union_operator_matches_explicit_formula() -> None:
    rng = np.random.default_rng(2)
    matrix = rng.normal(size=(12, 20)).astype(np.float32)
    vectors, _ = np.linalg.qr(rng.normal(size=(8, 8)))
    evals = np.linspace(0.5, 3.0, 8)
    positions = np.arange(4, 12)
    operator = MODULE.UnionOperator(
        positions,
        evals,
        vectors.astype(np.float32),
        gamma=0.65,
    )
    got = operator.apply(matrix)
    block = matrix[:, positions]
    weights = evals ** (-0.65 / 2.0)
    with np.errstate(all="ignore"):
        expected_block = ((block @ vectors) * weights) @ vectors.T
    assert np.isfinite(expected_block).all()
    expected_block *= np.median(np.linalg.norm(block, axis=1)) / np.median(
        np.linalg.norm(expected_block, axis=1)
    )
    expected = matrix.copy()
    expected[:, positions] = expected_block
    np.testing.assert_allclose(got, expected, rtol=2e-6, atol=2e-6)


def test_subspace_operator_matches_explicit_formula() -> None:
    rng = np.random.default_rng(3)
    matrix = rng.normal(size=(15, 30)).astype(np.float32)
    basis, _ = np.linalg.qr(rng.normal(size=(10, 4)))
    weights = np.asarray([0.2, 0.4, 0.7, 1.0], dtype=np.float32)
    positions = np.arange(5, 15)
    operator = MODULE.SubspaceOperator(
        positions,
        basis.astype(np.float32),
        weights,
    )
    got = operator.apply(matrix)
    block = matrix[:, positions]
    with np.errstate(all="ignore"):
        expected_block = block + ((block @ basis) * (weights - 1.0)) @ basis.T
    assert np.isfinite(expected_block).all()
    expected_block *= np.median(np.linalg.norm(block, axis=1)) / np.median(
        np.linalg.norm(expected_block, axis=1)
    )
    expected = matrix.copy()
    expected[:, positions] = expected_block
    np.testing.assert_allclose(got, expected, rtol=2e-6, atol=2e-6)


def test_composition_order_is_implemented_distinctly() -> None:
    rng = np.random.default_rng(4)
    matrix = rng.normal(size=(18, 16)).astype(np.float32)
    union_vectors, _ = np.linalg.qr(rng.normal(size=(12, 12)))
    sub_vectors, _ = np.linalg.qr(rng.normal(size=(7, 3)))
    union = MODULE.UnionOperator(
        np.arange(2, 14),
        np.linspace(0.3, 2.0, 12),
        union_vectors.astype(np.float32),
        0.65,
    )
    subspace = MODULE.SubspaceOperator(
        np.arange(5, 12),
        sub_vectors.astype(np.float32),
        np.asarray([0.2, 0.5, 1.0], dtype=np.float32),
    )
    variants = MODULE.transform_primary_variants(matrix, union, subspace)
    assert not np.allclose(
        variants["union_then_subspace"],
        variants["subspace_then_union"],
    )


def test_selection_falls_back_to_union() -> None:
    summary = {
        name: {
            "delta_vs_union": -0.01,
            "positive_folds_vs_union": 0,
            "direction_delta_vs_union": {"a->b": -0.01, "b->a": -0.01},
        }
        for name in MODULE.ELIGIBLE_OVER_UNION
    }
    assert MODULE.select_over_union(summary, 0.001, 6) == "union_only"


def test_fold_axes_remove_all_target_coordinates_from_both_operators() -> None:
    official = ["vcc", "train", "safe", "held_only", "source_only"]
    axes = MODULE.build_fold_axes(
        official=official,
        vcc=["vcc"],
        train=["train"],
        subspace_genes=["vcc", "train", "safe"],
        held_supported_official=np.asarray([0, 1, 2, 3]),
        source_supported_official=[np.asarray([0, 1, 2, 4])],
    )

    assert axes.score_genes == ["train", "safe", "held_only"]
    assert axes.union_genes == ["safe"]
    assert axes.subspace_genes == ["safe"]
    assert [axes.score_genes[index] for index in axes.inner_score_positions] == [
        "safe",
        "held_only",
    ]


def test_inner_folds_rank_queries_against_the_full_catalog() -> None:
    rng = np.random.default_rng(5)
    truth = rng.normal(size=(20, 12)).astype(np.float32)
    prediction = (truth + 0.7 * rng.normal(size=truth.shape)).astype(np.float32)
    variants = {
        name: prediction for name in (*MODULE.PRIMARY_VARIANTS, *MODULE.NULL_VARIANTS)
    }
    positions = np.arange(truth.shape[1])

    summary = MODULE.score_inner_seed(
        {"a->b": variants},
        {"a->b": truth},
        positions,
        seed=17,
        folds=5,
    )

    per_query = MODULE.pds_per_query(prediction, truth)
    fold_rows = MODULE.split_folds(
        len(truth),
        5,
        MODULE.stable_seed(("inner-folds", "a->b"), 17),
    )
    expected = [float(np.mean(per_query[rows])) for rows in fold_rows]
    np.testing.assert_allclose(summary["identity"]["fold_scores"], expected)
    assert summary["identity"]["mean"] == np.mean(expected)


def test_held_truth_can_be_extracted_without_background_centering(
    tmp_path: Path,
) -> None:
    np.save(
        tmp_path / "delta.npy",
        np.asarray([[11.0, 23.0]], dtype=np.float32),
    )
    donor = MODULE.Donor(
        name="held",
        path=tmp_path,
        genes=["g1", "g2"],
        perts=["p1"],
        gene_pos={"g1": 0, "g2": 1},
        pert_pos={"p1": 0},
        supported_official=np.asarray([0, 1]),
        supported_donor=np.asarray([0, 1]),
        background_rows=np.asarray([], dtype=np.int64),
        center=np.asarray([10.0, 20.0], dtype=np.float32),
    )

    centered = MODULE.extract_panel(
        donor,
        ["p1"],
        np.asarray([0, 1]),
        chunk_rows=1,
    )
    raw = MODULE.extract_panel(
        donor,
        ["p1"],
        np.asarray([0, 1]),
        chunk_rows=1,
        center_on_background=False,
    )

    np.testing.assert_array_equal(centered, [[1.0, 3.0]])
    np.testing.assert_array_equal(raw, [[11.0, 23.0]])


def test_covariance_path_is_memmapped_and_removed_after_diagonalization(
    tmp_path: Path,
) -> None:
    rng = np.random.default_rng(6)
    basis_path = tmp_path / "basis.npy"
    covariance_path = tmp_path / "covariance.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float32,
        shape=(32, 12),
    )
    basis[:] = rng.normal(size=basis.shape)
    expected = np.asarray(basis).T @ np.asarray(basis) / len(basis)

    covariance = MODULE.build_covariance(
        basis,
        covariance_path,
        soft_limit=0.95,
    )
    assert isinstance(covariance, np.memmap)
    np.testing.assert_allclose(covariance, expected, rtol=2e-6, atol=2e-6)
    evals, evecs = MODULE.diagonalize_union(
        covariance,
        covariance_path,
        soft_limit=0.95,
    )

    with np.errstate(all="ignore"):
        reconstructed = (evecs * evals) @ evecs.T
    assert np.isfinite(reconstructed).all()
    np.testing.assert_allclose(reconstructed, expected, rtol=2e-5, atol=2e-5)
    assert not covariance_path.exists()
    MODULE.close_memmap(basis)


def test_summary_never_selects_a_composition_from_outer_scores() -> None:
    outer_names = ["a", "b"]
    folds = {}
    for outer in outer_names:
        folds[outer] = {
            "outer_scores": {
                "identity": 0.50,
                "union_only": 0.70,
                "subspace_only": 0.80,
                "union_then_subspace": 0.99,
                "subspace_then_union": 0.60,
            },
            "seeds": {
                "1": {
                    "selected": "subspace_only",
                    "outer_selected_score": 0.80,
                    "outer_selected_delta_vs_union": 0.10,
                    "outer_matched_null_delta_vs_union": -0.01,
                }
            },
        }

    _, seed_macro, head_to_head = MODULE.summarize_folds(
        folds,
        outer_names,
        [1],
    )

    assert seed_macro["1"]["selected_variants_by_outer"] == {
        "a": "subspace_only",
        "b": "subspace_only",
    }
    assert head_to_head["nested_selection_gate"]
    assert "best_composition" not in head_to_head

    # A macro-average win is insufficient when the matched null wins one
    # outer fold.
    folds["b"]["seeds"]["1"]["outer_matched_null_delta_vs_union"] = 0.20
    _, seed_macro, head_to_head = MODULE.summarize_folds(
        folds,
        outer_names,
        [1],
    )
    assert (
        seed_macro["1"]["selected_delta_vs_union"]
        > seed_macro["1"]["matched_random_delta_vs_union"]
    )
    assert seed_macro["1"]["positive_outer_folds_vs_matched_random"] == 1
    assert not head_to_head["nested_selection_gate"]
