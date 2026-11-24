#!/usr/bin/env python3
"""Synthetic, memory-bounded tests for experiment 080."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
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


RUN = load_module("exp080_run_test", "run.py")
STAGE = load_module("exp080_stage_test", "stage_cache.py")


def write_csv(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[column])
        writer.writeheader()
        writer.writerows([{column: value} for value in values])


def cross_entropy(scores: np.ndarray, labels: np.ndarray, temperature: float) -> float:
    logits = scores / temperature
    logits -= logits.max(axis=1, keepdims=True)
    return float(
        np.mean(
            np.log(np.exp(logits).sum(axis=1)) - logits[np.arange(len(labels)), labels]
        )
    )


def test_analytic_listwise_gradient_matches_finite_difference() -> None:
    rng = np.random.default_rng(800)
    n, rank = 7, 3
    q = RUN.row_normalize(rng.normal(size=(n, 11)))
    t = RUN.row_normalize(rng.normal(size=(n, 11)))
    basis = RUN.randomized_right_basis(np.vstack([q, t]), rank=5, seed=800)
    x = q @ basis
    tp = t @ basis
    pair = RUN.TrainingPair(
        base_scores=q @ t.T,
        query_projection=x,
        truth_projection=tp,
        labels=np.arange(n, dtype=np.int64),
    )
    u = rng.normal(scale=0.1, size=(5, rank))
    v = rng.normal(scale=0.1, size=(5, rank))
    beta, temperature = 0.2, 0.15
    loss, gradient_u, gradient_v = RUN._softmax_cross_entropy_gradient(
        pair, u.copy(), v.copy(), beta=beta, temperature=temperature
    )

    def objective(local_u: np.ndarray, local_v: np.ndarray) -> float:
        latent = x @ local_v
        residual = latent @ local_u.T
        denominator = np.sqrt(
            np.maximum(
                1
                + 2 * beta * np.sum(x * residual, axis=1)
                + beta**2 * np.sum(residual * residual, axis=1),
                1e-18,
            )
        )
        scores = (pair.base_scores + beta * residual @ tp.T) / denominator[:, None]
        return cross_entropy(scores, pair.labels, temperature)

    epsilon = 1e-6
    for matrix, gradient, other, is_u in (
        (u, gradient_u, v, True),
        (v, gradient_v, u, False),
    ):
        for index in ((0, 0), (2, 1), (4, 2)):
            plus = matrix.copy()
            minus = matrix.copy()
            plus[index] += epsilon
            minus[index] -= epsilon
            if is_u:
                numerical = (objective(plus, other) - objective(minus, other)) / (
                    2 * epsilon
                )
            else:
                numerical = (objective(other, plus) - objective(other, minus)) / (
                    2 * epsilon
                )
            np.testing.assert_allclose(gradient[index], numerical, rtol=2e-5, atol=2e-6)
    assert np.isfinite(loss)


def test_target_split_is_hash_stable_disjoint_and_exhaustive() -> None:
    targets = [f"T{i:03d}" for i in range(101)]
    first = RUN.deterministic_target_split(targets, seed=80)
    second = RUN.deterministic_target_split(targets, seed=80)
    assert {name: values.tolist() for name, values in first.items()} == {
        name: values.tolist() for name, values in second.items()
    }
    sets = [set(values.tolist()) for values in first.values()]
    assert not any(
        left & right for i, left in enumerate(sets) for right in sets[i + 1 :]
    )
    assert set().union(*sets) == set(range(len(targets)))
    changed = RUN.deterministic_target_split(targets, seed=81)
    assert first["fit"].tolist() != changed["fit"].tolist()


def test_norm_decile_permutation_is_stratified_and_bijective() -> None:
    truth = np.diag(np.linspace(1.0, 100.0, 100))
    labels = RUN.permute_labels_within_norm_deciles(truth, seed=802)
    assert sorted(labels.tolist()) == list(range(100))
    for start in range(0, 100, 10):
        assert set(labels[start : start + 10].tolist()) == set(range(start, start + 10))


def test_coherent_permutation_plan_is_frozen_across_all_fit_stages() -> None:
    rng = np.random.default_rng(812)
    targets = [f"T{i:03d}" for i in range(60)]
    splits = RUN.deterministic_target_split(targets, seed=812)
    truths = {
        donor: rng.normal(size=(len(targets), 12)) for donor in RUN.TEACHER_DONORS
    }
    first = RUN.build_permutation_plan(truths, splits, replicate=3, seed=813)
    second = RUN.build_permutation_plan(truths, splits, replicate=3, seed=813)
    changed = RUN.build_permutation_plan(truths, splits, replicate=4, seed=813)
    assert first.plan_sha256 == second.plan_sha256
    assert first.plan_sha256 != changed.plan_sha256
    assert set(first.labels) == {"fit", "fit_tune", "all"}
    for stage in first.labels:
        for teacher in RUN.TEACHER_DONORS:
            observed = first.labels[stage][teacher]
            assert observed.flags.writeable is False
            assert np.array_equal(observed, second.labels[stage][teacher])
            assert sorted(observed.tolist()) == list(range(len(observed)))
            rows = (
                splits["fit"]
                if stage == "fit"
                else (
                    np.sort(np.concatenate([splits["fit"], splits["tune"]]))
                    if stage == "fit_tune"
                    else np.arange(len(targets), dtype=np.int64)
                )
            )
            expected = RUN.restrict_permutation(
                first.global_labels[teacher],
                rows,
            )
            np.testing.assert_array_equal(observed, expected)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "plans.npz"
        record = RUN.save_permutation_plans(path, [first, changed])
        assert record["replicates"] == 2
        assert record["sha256"] == RUN.sha256_file(path)
        with np.load(path, allow_pickle=False) as archive:
            for stage, by_teacher in first.labels.items():
                for teacher, labels in by_teacher.items():
                    key = f"r0003__{stage}__{teacher}"
                    np.testing.assert_array_equal(archive[key], labels)


def test_adapter_training_reduces_listwise_loss_and_preserves_rank_limit() -> None:
    rng = np.random.default_rng(803)
    n, d = 40, 24
    query = rng.normal(size=(n, d))
    transform = np.eye(d)
    transform[:4, 4:8] += 0.8 * np.eye(4)
    truth = query @ transform + rng.normal(scale=0.03, size=(n, d))
    candidate = RUN.Candidate(4, 0.20, 0.15, 1e-4)
    model = RUN.fit_adapter(
        query,
        {"teacher": truth},
        np.arange(d, dtype=np.int64),
        candidate,
        np.arange(n, dtype=np.int64),
        basis_dim=12,
        epochs=40,
        learning_rate=0.03,
        gradient_clip=5.0,
        seed=803,
    )
    assert model.u.shape == (12, 4)
    assert model.v.shape == (12, 4)
    assert model.fit_audit["last_loss"] < model.fit_audit["first_loss"]
    transformed = RUN.transform_query(query, np.arange(d, dtype=np.int64), model)
    np.testing.assert_allclose(np.linalg.norm(transformed, axis=1), 1.0, atol=1e-10)
    fast_scores, _ = RUN.adapter_scores_and_angles(
        query, truth, np.arange(d, dtype=np.int64), model
    )
    explicit_scores = transformed @ RUN.row_normalize(truth).T
    np.testing.assert_allclose(fast_scores, explicit_scores, atol=1e-10)


def test_fit_cache_is_semantics_preserving_for_every_adapter_rank() -> None:
    rng = np.random.default_rng(814)
    query = rng.normal(size=(36, 28))
    truths = {"teacher": query @ (np.eye(28) + 0.02 * rng.normal(size=(28, 28)))}
    rows = np.arange(30, dtype=np.int64)
    correction = np.arange(24, dtype=np.int64)
    cache = RUN.prepare_fit_cache(
        query,
        truths,
        correction,
        rows,
        basis_dim=20,
        basis_seed=815,
    )
    for rank in (4, 8, 16):
        candidate = RUN.Candidate(rank, 0.05, 0.07, 1e-3)
        kwargs = {
            "basis_dim": 20,
            "epochs": 4,
            "learning_rate": 0.01,
            "gradient_clip": 5.0,
            "seed": 816 + rank,
            "basis_seed": 815,
        }
        serial = RUN.fit_adapter(
            query,
            truths,
            correction,
            candidate,
            rows,
            **kwargs,
        )
        cached = RUN.fit_adapter(
            query,
            truths,
            correction,
            candidate,
            rows,
            prepared_cache=cache,
            **kwargs,
        )
        np.testing.assert_array_equal(cached.basis, serial.basis)
        np.testing.assert_array_equal(cached.u, serial.u)
        np.testing.assert_array_equal(cached.v, serial.v)
        assert cached.fit_audit["last_loss"] == serial.fit_audit["last_loss"]


def test_batch_equivalence_all_ranks_two_teachers_permuted_and_clipped() -> None:
    rng = np.random.default_rng(817)
    query = rng.normal(size=(28, 24))
    truths = {
        "teacher_a": rng.normal(size=(28, 24)),
        "teacher_b": rng.normal(size=(28, 24)),
    }
    rows = np.arange(24, dtype=np.int64)
    correction = np.arange(20, dtype=np.int64)
    cache = RUN.prepare_fit_cache(
        query,
        truths,
        correction,
        rows,
        basis_dim=20,
        basis_seed=818,
    )
    permutations = {
        "teacher_a": np.arange(len(rows) - 1, -1, -1, dtype=np.int64),
        "teacher_b": np.roll(np.arange(len(rows), dtype=np.int64), 7),
    }
    epochs = 6
    gradient_clip = 1e-8
    for rank in (4, 8, 16):
        candidates = [
            RUN.Candidate(rank, beta, temperature, regularization)
            for beta, temperature, regularization in (
                (0.05, 0.03, 1e-4),
                (0.20, 0.15, 1e-2),
            )
        ]
        batched = RUN.fit_adapter_batch(
            candidates,
            cache,
            epochs=epochs,
            learning_rate=0.02,
            gradient_clip=gradient_clip,
            seed=819 + rank,
            label_permutations=permutations,
        )
        for candidate, observed in zip(candidates, batched, strict=True):
            expected = RUN.fit_adapter(
                query,
                truths,
                correction,
                candidate,
                rows,
                basis_dim=20,
                epochs=epochs,
                learning_rate=0.02,
                gradient_clip=gradient_clip,
                seed=819 + rank,
                basis_seed=818,
                label_permutations=permutations,
                prepared_cache=cache,
            )
            np.testing.assert_allclose(observed.u, expected.u, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(observed.v, expected.v, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(
                observed.fit_audit["last_loss"],
                expected.fit_audit["last_loss"],
                rtol=1e-12,
                atol=1e-12,
            )
            assert observed.fit_audit["gradient_clipped_steps"] == epochs
            assert expected.fit_audit["gradient_clipped_steps"] == epochs
            assert observed.fit_audit["labels_permuted_within_norm_deciles"] is True
            assert observed.fit_audit["factor_penalty"] is False
            assert (
                observed.fit_audit["regularization_objective"]
                == "0.5 * lambda * ||U V^T||_F^2"
            )


def test_frozen_adapter_roundtrip_preserves_exact_weights() -> None:
    rng = np.random.default_rng(811)
    candidate = RUN.Candidate(4, 0.05, 0.07, 1e-3)
    adapter = RUN.Adapter(
        candidate=candidate,
        basis=rng.normal(size=(9, 6)),
        u=rng.normal(size=(6, 4)),
        v=rng.normal(size=(6, 4)),
        fit_audit={"rows": 12, "held_truth_opened": False},
    )
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        record = RUN.save_adapter(Path(directory) / "frozen.npz", adapter)
        loaded = RUN.load_frozen_adapter(record)
        assert loaded.candidate == candidate
        assert np.array_equal(loaded.basis, adapter.basis)
        assert np.array_equal(loaded.u, adapter.u)
        assert np.array_equal(loaded.v, adapter.v)


def test_candidate_selector_enforces_every_cell_and_angle_guard() -> None:
    good = RUN.Candidate(4, 0.02, 0.07, 1e-3)
    exactly_better = RUN.Candidate(16, 0.20, 0.15, 1e-4)
    bad_sign = RUN.Candidate(8, 0.05, 0.07, 1e-3)
    bad_angle = RUN.Candidate(16, 0.20, 0.07, 1e-3)
    cells = {
        good.key: [
            {"delta": 0.002, "median_angle_degrees": 2.0},
            {"delta": 0.003, "median_angle_degrees": 3.0},
        ],
        exactly_better.key: [
            {"delta": 0.0021, "median_angle_degrees": 2.0},
            {"delta": 0.0031, "median_angle_degrees": 3.0},
        ],
        bad_sign.key: [
            {"delta": 0.010, "median_angle_degrees": 2.0},
            {"delta": -0.001, "median_angle_degrees": 2.0},
        ],
        bad_angle.key: [
            {"delta": 0.02, "median_angle_degrees": 11.0},
            {"delta": 0.02, "median_angle_degrees": 9.0},
        ],
    }
    selected, audit = RUN.choose_candidate(
        cells,
        [good, exactly_better, bad_sign, bad_angle],
    )
    assert selected == exactly_better.key
    assert audit["eligible_count"] == 2
    assert audit["exact_tie_count"] == 1
    assert "tie_window" not in audit


def test_angular_guard_is_enforced_separately_in_every_fold() -> None:
    assert RUN.every_fold_angle_below_limit(
        {"orion_hct116": 9.9, "orion_hek293t": 8.0, "kolf2": 7.0}
    )
    assert not RUN.every_fold_angle_below_limit(
        {"orion_hct116": 2.0, "orion_hek293t": 10.0, "kolf2": 1.0}
    )


def test_full_gallery_pds_uses_all_600_candidates_for_subset_queries() -> None:
    scores = np.asarray(
        [
            [0.9, 0.8, 0.95, 0.1],
            [0.2, 0.7, 0.1, 0.8],
        ]
    )
    observed = RUN.pds_per_query_from_scores(scores, np.asarray([0, 3], dtype=np.int64))
    np.testing.assert_allclose(observed, [2 / 3, 1.0])


def test_source_clusters_and_fold_specific_bootstrap_do_not_mix_folds() -> None:
    rng = np.random.default_rng(804)
    first = RUN.row_normalize(rng.normal(size=(30, 12)))
    second = RUN.row_normalize(rng.normal(size=(17, 12)))
    labels_a, audit_a = RUN.source_effect_clusters(
        [first],
        reciprocal_k=4,
        cosine_threshold=0.4,
    )
    labels_b, audit_b = RUN.source_effect_clusters(
        [second],
        reciprocal_k=3,
        cosine_threshold=0.4,
    )
    assert len(labels_a) == 30 and len(labels_b) == 17
    assert audit_a["clusters"] >= 1 and audit_b["clusters"] >= 1
    deltas = {
        "fold_a": np.linspace(-0.01, 0.03, 30),
        "fold_b": np.linspace(0.00, 0.04, 17),
    }
    labels = {"fold_a": labels_a, "fold_b": labels_b}
    bootstrap = RUN.fold_specific_cluster_bootstrap(
        deltas, labels, replicates=200, seed=804
    )
    expected = np.mean([np.mean(deltas["fold_a"]), np.mean(deltas["fold_b"])])
    assert np.isclose(bootstrap["point_estimate"], expected)
    assert bootstrap["clusters_by_fold"] == {
        "fold_a": audit_a["clusters"],
        "fold_b": audit_b["clusters"],
    }


def test_generalized_pm_slot_preserves_missing_target_rows() -> None:
    incumbent = STAGE.EXP062.PreparedSource(
        name="k562",
        effects=np.asarray([[2.0, 4.0], [6.0, 8.0]]),
        target_available=np.asarray([True, True]),
        mean_on_axis=np.zeros(2),
        raw_median_norm=1.0,
        mapped_gene_count=2,
        duplicate_symbols=(),
    )
    perturbmatch = STAGE.EXP062.PreparedSource(
        name="pm",
        effects=np.asarray([[10.0, 20.0], [30.0, 40.0]]),
        target_available=np.asarray([True, False]),
        mean_on_axis=np.zeros(2),
        raw_median_norm=1.0,
        mapped_gene_count=2,
        duplicate_symbols=(),
    )
    slot, audit = STAGE.generalized_pm_slot(
        incumbent, perturbmatch, np.asarray([0, 1], dtype=np.int64)
    )
    np.testing.assert_allclose(slot.effects[0], [6.0, 12.0])
    np.testing.assert_allclose(slot.effects[1], incumbent.effects[1])
    assert audit["targets_pm_modified"] == 1
    assert audit["targets_incumbent_preserved"] == 1


def test_principal_operator_restriction_equals_dense_zeroed_application() -> None:
    rng = np.random.default_rng(810)
    width = 9
    q, _ = np.linalg.qr(rng.normal(size=(width, width)))
    evals = np.linspace(0.4, 2.0, width)
    operator = STAGE.EXP062.Operator(
        evals=evals,
        evecs=q,
        metadata={},
        path=Path("synthetic.npz"),
        metadata_path=Path("synthetic.json"),
    )
    full_columns = np.arange(width, dtype=np.int64)
    active = np.asarray([0, 2, 3, 6, 8], dtype=np.int64)
    matrix = rng.normal(size=(13, width))
    matrix[:, np.setdiff1d(full_columns, active)] = 0.0
    preserve = np.asarray([0, 2, 5, 9], dtype=np.int64)
    observed, audit = STAGE.apply_operator_restriction(
        matrix,
        operator,
        full_columns,
        active,
        preserve_rows=preserve,
    )

    weights = evals ** (-STAGE.EXP064.WHITEN_GAMMA / 2.0)
    dense_operator = (q * weights) @ q.T
    dense_retained = (matrix @ dense_operator)[:, active]
    before = np.median(np.linalg.norm(matrix[np.ix_(preserve, active)], axis=1))
    after = np.median(np.linalg.norm(dense_retained[preserve], axis=1))
    expected = dense_retained * (before / after)
    np.testing.assert_allclose(observed[:, active], expected, atol=1e-12)
    assert np.count_nonzero(observed[:, np.setdiff1d(full_columns, active)]) == 0
    assert audit["dense_zero_equivalence_proved"] is True

    changed_forbidden = matrix.copy()
    changed_forbidden[:, 1] = 1e9
    try:
        STAGE.apply_operator_restriction(
            changed_forbidden,
            operator,
            full_columns,
            active,
            preserve_rows=preserve,
        )
    except AssertionError:
        pass
    else:
        raise AssertionError("nonzero forbidden input was not rejected")


def write_compact_bank(
    path: Path,
    genes: list[str],
    perts: list[str],
    delta: np.ndarray,
    *,
    perturbmatch: bool = False,
    cd4: bool = False,
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", np.asarray(delta, dtype=np.float32))
    if not perturbmatch:
        np.save(path / "b_ctrl.npy", np.ones(len(genes), dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes) + "\n")
    (path / "perts.txt").write_text("\n".join(perts) + "\n")
    metadata: dict[str, object] = {}
    if perturbmatch:
        metadata.update(
            {
                "source_md5": "0017076a5cc7b9979febf3bd19bd2157",
                "vcc_target_coverage": 251,
            }
        )
    if cd4:
        metadata["effect_definition"] = "delta = ln(2) * source DESeq2 log2 fold change"
    (path / "meta.json").write_text(json.dumps(metadata))


def test_real_shape_contract_has_49_coordinate_less_targets_and_17806_genes() -> None:
    manifest = json.loads((HERE / "artifact_manifest.json").read_text())
    root = HERE.parents[1]
    donor_root = Path(os.environ.get("EXP080_REAL_DONOR_ROOT", root / "data/donors"))
    kolf_root = Path(
        os.environ.get("EXP080_REAL_KOLF_ROOT", root / "data/external/kolf2")
    )
    truth_cache = Path(
        os.environ.get(
            "EXP080_REAL_TRUTH_CACHE",
            root / "experiments/044-csls-distill/cache",
        )
    )
    official = RUN.read_csv_column(
        root / "data/challenge/gene_names.csv",
        "gene_name",
    )
    panel = RUN.read_csv_column(
        root / "data/challenge/pert_counts.csv",
        "target_gene",
    )
    train_targets = RUN.read_lines(truth_cache / "train_perts.txt")
    assert len(train_targets) == 600
    assert len(set(train_targets) - set(official)) == 49

    score_genes, score_positions, audit = STAGE.established_final_score_axis(
        official,
        panel,
        donor_root,
        manifest,
    )
    assert len(score_genes) == 17_806
    assert (
        audit["gene_order_sha256"] == manifest["final_score_axis"]["gene_order_sha256"]
    )
    assert (
        audit["official_positions_sha256"]
        == manifest["final_score_axis"]["official_positions_sha256"]
    )
    assert len(score_positions) == 17_806
    kolf_genes = set(RUN.read_gene_lines(kolf_root / "genes.txt"))
    assert (
        tuple(
            gene for gene in official if gene in kolf_genes and gene not in set(panel)
        )
        == score_genes
    )
    _, validated_targets, _, axes = RUN.validate_truth_cache(
        truth_cache,
        official,
        panel,
        donor_root,
        manifest,
    )
    assert len(validated_targets) == 600
    assert axes.panel_score_genes == score_genes


def test_compact_panel_truth_uses_complete_requested_axis_and_background_center() -> (
    None
):
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        donor_root = root / "donors"
        values = np.asarray(
            [
                [1.0, 3.0, 8.0],
                [2.0, 5.0, 4.0],
                [9.0, 7.0, 6.0],
            ],
            dtype=np.float32,
        )
        write_compact_bank(
            donor_root / "held",
            ["G0", "G1", "G0"],
            ["T0", "T1", "BG"],
            values,
        )
        observed, audit = RUN.load_compact_panel_truth(
            donor_root,
            "held",
            ["T1", "T0"],
            ["G0", "G1"],
            chunk_rows=2,
        )
        expected = values[[1, 0]][:, [2, 1]] - values[:, [2, 1]].mean(axis=0)
        np.testing.assert_allclose(observed, expected)
        assert audit["targets"] == 2
        assert audit["score_genes"] == 2
        assert audit["score_gene_order_sha256"] == RUN.sha256_lines(["G0", "G1"])


def test_blocked_coordinates_cannot_change_safe_centering_or_magnitude() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        genes = ["BLOCKED", "SAFE1", "SAFE2"]
        perts = ["T0", "T1", "BG"]
        base = np.asarray(
            [[1.0, 2.0, 7.0], [2.0, 8.0, 3.0], [3.0, 5.0, 9.0]],
            dtype=np.float32,
        )
        changed = base.copy()
        changed[:, 0] += np.asarray([1e6, -2e6, 3e6], dtype=np.float32)
        write_compact_bank(root / "base", genes, perts, base)
        write_compact_bank(root / "changed", genes, perts, changed)
        blocked = np.asarray([0], dtype=np.int64)
        prepared = []
        legacy = []
        for name in ("base", "changed"):
            bank = STAGE.EXP062.Bank.load(name, root / name)
            prepared.append(
                STAGE.prepare_main_source_safe(
                    bank,
                    ["T0", "T1"],
                    genes,
                    blocked,
                    chunk_rows=2,
                    max_rss_gib=0.95,
                )
            )
            legacy.append(
                STAGE.EXP062.prepare_main_source(
                    bank,
                    ["T0", "T1"],
                    genes,
                    blocked,
                    chunk_rows=2,
                    max_rss_gib=0.95,
                )
            )
        np.testing.assert_array_equal(prepared[0].effects, prepared[1].effects)
        assert prepared[0].raw_median_norm == prepared[1].raw_median_norm
        legacy_change = float(
            np.max(np.abs(legacy[0].effects[:, 1:] - legacy[1].effects[:, 1:]))
        )
        assert legacy_change > 5.548


def test_kolf_only_gain_cannot_satisfy_primary_macro_gate() -> None:
    deltas = {
        "orion_hct116": 0.0,
        "orion_hek293t": 0.0,
        "kolf2": 1.0,
    }
    assert RUN.primary_fold_mean(deltas) == 0.0
    assert RUN.primary_fold_mean(deltas) < 0.003


def test_arm_d_cache_stage_masks_both_target_panels_and_holds_out_donor() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        train_targets = [f"T{i:02d}" for i in range(20)]
        coordinate_train_targets = train_targets[:-2]
        outer_targets = [f"P{i:02d}" for i in range(15)]
        features = [f"G{i:02d}" for i in range(16)]
        official = [*outer_targets, *coordinate_train_targets, *features]
        all_targets = [*train_targets, *outer_targets]
        perts = [*all_targets, "BG0", "BG1", "BG2"]
        rng = np.random.default_rng(809)
        shared = rng.normal(size=(len(perts), len(official)))
        # Make target-gene coordinates huge; the cache must remove them.
        shared[:, : len(outer_targets) + len(train_targets)] += 1000.0

        donor_root = root / "donors"
        main_names = ("k562_gwps", "orion_hct116", "orion_hek293t")
        for offset, name in enumerate(main_names):
            write_compact_bank(
                donor_root / name,
                official,
                perts,
                shared + rng.normal(scale=0.02 + 0.01 * offset, size=shared.shape),
            )
        pm_root = root / "perturbmatch"
        write_compact_bank(
            pm_root,
            official,
            perts,
            shared + rng.normal(scale=0.025, size=shared.shape),
            perturbmatch=True,
        )
        cd4_root = root / "cd4"
        for offset, name in enumerate(STAGE.CD4_DONORS):
            write_compact_bank(
                cd4_root / name,
                official,
                outer_targets,
                shared[len(train_targets) : len(all_targets)]
                + rng.normal(
                    scale=0.03 + 0.01 * offset,
                    size=(len(outer_targets), len(official)),
                ),
                cd4=True,
            )

        genes_csv = root / "gene_names.csv"
        panel_csv = root / "pert_counts.csv"
        write_csv(genes_csv, "gene_name", official)
        write_csv(panel_csv, "target_gene", outer_targets)
        truth_cache = root / "truth"
        truth_cache.mkdir()
        master = [*coordinate_train_targets, *features]
        (truth_cache / "master_genes.txt").write_text("\n".join(master) + "\n")
        (truth_cache / "train_perts.txt").write_text("\n".join(train_targets) + "\n")
        (truth_cache / "outer_perts.txt").write_text("\n".join(outer_targets) + "\n")
        np.save(
            truth_cache / "train_mask.npy",
            np.asarray(
                [False] * len(coordinate_train_targets) + [True] * len(features)
            ),
        )
        (truth_cache / "meta.json").write_text('{"fixture": true}\n')

        operator_dir = root / "operators"
        operator_dir.mkdir()
        operators: dict[str, object] = {}
        panel_set = set(outer_targets)
        for held in STAGE.FINAL_FOLDS:
            if held == STAGE.KOLF_FOLD:
                source_donors = list(STAGE.MAIN_DONORS)
            else:
                companion = next(
                    name for name in STAGE.DEVELOPMENT_FOLDS if name != held
                )
                source_donors = ["k562_gwps", companion]
            columns = np.asarray(
                [index for index, gene in enumerate(official) if gene not in panel_set],
                dtype=np.int64,
            )
            npz_path = operator_dir / f"held_{held}.npz"
            json_path = operator_dir / f"held_{held}.json"
            np.savez(
                npz_path,
                evals=np.ones(len(columns), dtype=np.float64),
                evecs=np.eye(len(columns), dtype=np.float32),
            )
            json_path.write_text(
                json.dumps(
                    {
                        "axis": "union",
                        "donors": source_donors,
                        "n_cols": len(columns),
                        "wcols_sha256": STAGE.hashlib.sha256(
                            columns.tobytes()
                        ).hexdigest(),
                        "gamma": 0.65,
                    }
                )
            )
            operators[held] = {
                "source_donors": source_donors,
                "npz_sha256": STAGE.sha256_file(npz_path),
                "json_sha256": STAGE.sha256_file(json_path),
                "columns_sha256": STAGE.sha256_array(columns),
            }

        exp064_results = root / "exp064-results.json"
        exp064_results.write_text(
            json.dumps(
                {
                    "axis": {
                        "target_order": outer_targets,
                    }
                }
            )
        )

        def hashes(path: Path, filenames: list[str]) -> dict[str, str]:
            return {
                filename: STAGE.sha256_file(path / filename) for filename in filenames
            }

        manifest = {
            "schema_version": 2,
            "challenge": {
                "gene_names.csv": STAGE.sha256_file(genes_csv),
                "pert_counts.csv": STAGE.sha256_file(panel_csv),
            },
            "final_score_axis": {
                "gene_count": len(master),
                "gene_order_sha256": STAGE.sha256_lines(master),
                "official_positions_sha256": STAGE.sha256_array(
                    np.asarray(
                        [official.index(gene) for gene in master],
                        dtype=np.int64,
                    )
                ),
            },
            "donors": {
                name: hashes(
                    donor_root / name,
                    ["delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json"],
                )
                for name in main_names
            },
            "perturbmatch": hashes(
                pm_root, ["delta.npy", "genes.txt", "perts.txt", "meta.json"]
            ),
            "cd4": {
                name: hashes(
                    cd4_root / name,
                    ["delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json"],
                )
                for name in STAGE.CD4_DONORS
            },
            "truth_cache": {
                "files": {
                    path.name: STAGE.sha256_file(path)
                    for path in truth_cache.iterdir()
                    if path.is_file()
                }
            },
            "exp064_reference": {
                "results_sha256": STAGE.sha256_file(exp064_results),
                "gallery_targets": len(outer_targets),
                "target_order_sha256": STAGE.sha256_lines(outer_targets),
                "required_max_abs_error": 1e-6,
                "normalization_contract": "all_mapped_then_block_vcc",
                "independent_raw_bank_reconstruction_required": True,
                "learned_statistics_use_separate_retained_coordinate_features": True,
            },
            "kolf_confirmation": {
                "supported": True,
                "files": {},
            },
            "operators": operators,
        }
        manifest_path = root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))
        args = argparse.Namespace(
            donor_root=donor_root,
            perturbmatch=pm_root,
            cd4_root=cd4_root,
            truth_cache=truth_cache,
            operator_dir=operator_dir,
            exp064_results=exp064_results,
            artifact_manifest=manifest_path,
            genes_csv=genes_csv,
            panel_csv=panel_csv,
            output=root / "query-cache",
            chunk_rows=8,
            max_rss_gib=0.95,
            expected_training_targets=len(train_targets),
            minimum_outer_targets=len(outer_targets),
            exp064_atol=1e-6,
            guarded=True,
        )
        STAGE.verify_inputs(args, manifest)
        metadata = STAGE.stage(args, manifest)
        assert metadata["leakage_audit"]["held_truth_opened"] is False
        assert metadata["training_blocked_targets"] == len(train_targets) + len(
            outer_targets
        )
        assert metadata["training_blocked_target_coordinates"] == len(
            coordinate_train_targets
        ) + len(outer_targets)
        assert metadata["training_targets_without_official_coordinate"] == 2
        assert metadata["panel_blocked_target_coordinates"] == len(outer_targets)
        assert metadata["historical_target_coordinates_restored_for_panel"] == len(
            coordinate_train_targets
        )
        for held in STAGE.FINAL_FOLDS:
            assert held not in metadata["folds"][held]["source_donors"]
            train_query = np.load(args.output / held / "train_query.npy")
            cluster_query = np.load(args.output / held / "cluster_query.npy")
            panel_query = np.load(args.output / held / "panel_query.npy")
            assert train_query.shape == (len(train_targets), len(features))
            assert cluster_query.shape == (len(outer_targets), len(features))
            assert panel_query.shape == (
                len(outer_targets),
                len(master),
            )
            assert np.isfinite(train_query).all()
            assert np.isfinite(cluster_query).all()
            assert np.isfinite(panel_query).all()
            assert np.any(panel_query[:, : len(coordinate_train_targets)] != 0.0)
            restriction = metadata["folds"][held]["panels"]["train"][
                "operator_restriction"
            ]
            assert restriction["all_nonpanel_target_coordinates_excluded"] is True
            cluster_metadata = metadata["folds"][held]["panels"]["cluster"]
            assert (
                cluster_metadata["operator_restriction"][
                    "all_vcc_and_historical_target_coordinates_excluded"
                ]
                is True
            )
            assert cluster_metadata["support"]["used_for_source_clustering"] is True
        for held in STAGE.DEVELOPMENT_FOLDS:
            reproduction = metadata["folds"][held]["exp064_reproduction"]
            assert reproduction["passed"] is True
            assert reproduction["saved_query_score_axis_exact"] is True
            assert reproduction["retrieval_genes"] == len(master)
            assert reproduction["independent_from_staged_panel_preparation"] is True
            assert reproduction["max_abs_error"] <= 1e-6
            train_support = metadata["folds"][held]["panels"]["train"]["support"]
            panel_support = metadata["folds"][held]["panels"]["panel"]["support"]
            assert train_support["used_for_adapter_fit_or_selection"] is True
            assert panel_support["used_for_adapter_fit_or_selection"] is False
            assert "all mapped source" in panel_support["normalization_contract"]


def make_identity_fixture(
    root: Path, *, include_panel_truth: bool = False
) -> tuple[argparse.Namespace, dict[str, object]]:
    train_targets = [f"T{i:03d}" for i in range(600)]
    outer_targets = [f"P{i:03d}" for i in range(20)]
    score_genes = [f"G{i:03d}" for i in range(18)]
    master = [*train_targets, *score_genes]
    official = [*outer_targets, *master]
    panel = [*outer_targets, *[f"PX{i:03d}" for i in range(280)]]
    genes_csv = root / "gene_names.csv"
    panel_csv = root / "pert_counts.csv"
    write_csv(genes_csv, "gene_name", official + panel[len(outer_targets) :])
    write_csv(panel_csv, "target_gene", panel)

    truth_cache = root / "truth"
    truth_cache.mkdir()
    (truth_cache / "master_genes.txt").write_text("\n".join(master) + "\n")
    (truth_cache / "train_perts.txt").write_text("\n".join(train_targets) + "\n")
    (truth_cache / "outer_perts.txt").write_text("\n".join(outer_targets) + "\n")
    train_mask = np.asarray(
        [False] * len(train_targets) + [True] * len(score_genes), dtype=bool
    )
    np.save(truth_cache / "train_mask.npy", train_mask)
    (truth_cache / "meta.json").write_text('{"fixture": true}\n')
    rng = np.random.default_rng(805)
    signal = rng.normal(size=(len(train_targets), len(score_genes)))
    for donor in RUN.TEACHER_DONORS:
        values = np.zeros((len(train_targets), len(master)), dtype=np.float32)
        values[:, len(train_targets) :] = signal
        np.save(truth_cache / f"{donor}_raw_train.npy", values)
    donor_root = root / "donors"
    donor_perts = [*outer_targets, "BG0", "BG1", "BG2"]
    for donor in RUN.DEVELOPMENT_FOLDS:
        values = rng.normal(size=(len(donor_perts), len(master))).astype(np.float32)
        write_compact_bank(
            donor_root / donor,
            master,
            donor_perts,
            values,
        )

    query_cache = root / "query"
    query_cache.mkdir()
    folds: dict[str, object] = {}
    for held in RUN.FINAL_FOLDS:
        fold_dir = query_cache / held
        fold_dir.mkdir()
        train_query = np.asarray(signal, dtype=np.float32)
        panel_rows = len(panel) if held == RUN.KOLF_FOLD else len(outer_targets)
        cluster_query = rng.normal(size=(panel_rows, len(score_genes))).astype(
            np.float32
        )
        panel_query = rng.normal(size=(panel_rows, len(master))).astype(np.float32)
        np.save(fold_dir / "train_query.npy", train_query)
        np.save(fold_dir / "cluster_query.npy", cluster_query)
        np.save(fold_dir / "panel_query.npy", panel_query)
        if held == RUN.KOLF_FOLD:
            sources = ["k562_gwps", *RUN.DEVELOPMENT_FOLDS]
        else:
            source = next(name for name in RUN.DEVELOPMENT_FOLDS if name != held)
            sources = ["k562_gwps", source]
        folds[held] = {
            "held_donor": held,
            "source_donors": sources,
            "held_donor_used_in_prediction": False,
            "held_donor_used_in_operator": False,
            "operator_npz_sha256": f"{held}-npz",
            "operator_json_sha256": f"{held}-json",
            "panels": {
                "train": {
                    "query_sha256": RUN.sha256_file(fold_dir / "train_query.npy"),
                    "support": {
                        "used_for_adapter_fit_or_selection": True,
                    },
                    "operator_restriction": {
                        "all_nonpanel_target_coordinates_excluded": True,
                        "nonpanel_target_overlap_with_active_operator": 0,
                        "dense_zero_equivalence_proved": True,
                    },
                },
                "cluster": {
                    "query_sha256": RUN.sha256_file(fold_dir / "cluster_query.npy"),
                    "support": {
                        "used_for_adapter_fit_or_selection": False,
                        "used_for_source_clustering": True,
                    },
                    "operator_restriction": {
                        "all_vcc_and_historical_target_coordinates_excluded": True,
                    },
                    "score_gene_count": len(score_genes),
                    "score_gene_order_sha256": RUN.sha256_lines(score_genes),
                    "saved_query_is_complete_cluster_axis": True,
                },
                "panel": {
                    "query_sha256": RUN.sha256_file(fold_dir / "panel_query.npy"),
                    "support": {
                        "inference_only_fixed_baseline": True,
                        "used_for_adapter_fit_or_selection": False,
                    },
                    "score_gene_count": len(master),
                    "score_gene_order_sha256": RUN.sha256_lines(master),
                    "score_official_positions_sha256": RUN.sha256_array(
                        np.asarray(
                            [official.index(gene) for gene in master],
                            dtype=np.int64,
                        )
                    ),
                    "saved_query_is_complete_scored_axis": True,
                },
            },
        }
        if held in RUN.DEVELOPMENT_FOLDS:
            folds[held]["exp064_reproduction"] = {
                "passed": True,
                "independent_from_staged_panel_preparation": True,
                "max_abs_error": 0.0,
                "saved_query_score_axis_exact": True,
                "retrieval_genes": len(master),
                "retrieval_gene_order_sha256": RUN.sha256_lines(master),
            }
    correction_in_train_score = np.arange(len(score_genes), dtype=np.int64)
    correction_in_panel_score = np.arange(
        len(train_targets), len(master), dtype=np.int64
    )
    metadata = {
        "schema_version": RUN.CACHE_SCHEMA,
        "experiment": "080-listwise-metric-adapter",
        "arm": "D",
        "development_folds": list(RUN.DEVELOPMENT_FOLDS),
        "final_folds": list(RUN.FINAL_FOLDS),
        "training_targets": len(train_targets),
        "outer_targets": len(outer_targets),
        "training_outer_overlap": 0,
        "vcc_targets": len(panel),
        "train_score_genes": len(score_genes),
        "cluster_score_genes": len(score_genes),
        "panel_score_genes": len(master),
        "correction_genes": len(score_genes),
        "train_score_gene_order_sha256": RUN.sha256_lines(score_genes),
        "cluster_score_gene_order_sha256": RUN.sha256_lines(score_genes),
        "cluster_score_official_positions_sha256": RUN.sha256_array(
            np.asarray(
                [official.index(gene) for gene in score_genes],
                dtype=np.int64,
            )
        ),
        "panel_score_gene_order_sha256": RUN.sha256_lines(master),
        "panel_score_official_positions_sha256": RUN.sha256_array(
            np.asarray(
                [official.index(gene) for gene in master],
                dtype=np.int64,
            )
        ),
        "correction_gene_order_sha256": RUN.sha256_lines(score_genes),
        "correction_in_train_score_sha256": RUN.sha256_array(correction_in_train_score),
        "correction_in_panel_score_sha256": RUN.sha256_array(correction_in_panel_score),
        "train_order_sha256": RUN.sha256_lines(train_targets),
        "outer_order_sha256": RUN.sha256_lines(outer_targets),
        "panel_order_sha256": RUN.sha256_lines(panel),
        "final_score_axis": {
            "gene_count": len(master),
            "gene_order_sha256": RUN.sha256_lines(master),
            "official_positions_sha256": RUN.sha256_array(
                np.asarray(
                    [official.index(gene) for gene in master],
                    dtype=np.int64,
                )
            ),
            "development_donor_axes_equal": True,
            "excludes_exactly_vcc_coordinates_within_held_support": True,
        },
        "folds": folds,
    }
    (query_cache / "metadata.json").write_text(json.dumps(metadata))

    kolf_root = root / "kolf"
    kolf_files: dict[str, str] = {}
    if include_panel_truth:
        kolf_root.mkdir()
        np.save(
            kolf_root / "delta.npy",
            rng.normal(size=(len(panel), len(master))).astype(np.float32),
        )
        np.save(kolf_root / "b_ctrl.npy", np.ones(len(master), dtype=np.float32))
        (kolf_root / "genes.txt").write_text("\n".join(master) + "\n")
        (kolf_root / "perts.txt").write_text("\n".join(panel) + "\n")
        (kolf_root / "meta.json").write_text('{"fixture": true}\n')
        kolf_files = {
            path.name: RUN.sha256_file(path)
            for path in kolf_root.iterdir()
            if path.is_file()
        }

    manifest = {
        "schema_version": 2,
        "protocol_contracts": {
            "optimizer_regularization": "0.5 * lambda * ||U V^T||_F^2",
            "factor_penalty": False,
            "selection_objective": "exact_macro_tune_pds",
            "tie_window": None,
            "exact_tie_break": [
                "lower_rank",
                "lower_beta",
                "stronger_regularization",
                "temperature_nearest_0.07",
                "candidate_key",
            ],
            "batch_equivalence_regression": {
                "ranks": [4, 8, 16],
                "teachers": 2,
                "permuted_labels": True,
                "forced_gradient_clipping": True,
            },
            "source_clustering_query": "panel_rows_retained_coordinate_safe",
        },
        "challenge": {
            "gene_names.csv": RUN.sha256_file(genes_csv),
            "pert_counts.csv": RUN.sha256_file(panel_csv),
        },
        "final_score_axis": {
            "gene_count": len(master),
            "gene_order_sha256": RUN.sha256_lines(master),
            "official_positions_sha256": RUN.sha256_array(
                np.asarray(
                    [official.index(gene) for gene in master],
                    dtype=np.int64,
                )
            ),
        },
        "donors": {
            donor: {
                path.name: RUN.sha256_file(path)
                for path in (donor_root / donor).iterdir()
                if path.is_file()
            }
            for donor in RUN.DEVELOPMENT_FOLDS
        },
        "truth_cache": {
            "files": {
                path.name: RUN.sha256_file(path)
                for path in truth_cache.iterdir()
                if path.is_file()
            }
        },
        "exp064_reference": {
            "required_max_abs_error": 1e-6,
        },
        "operators": {
            held: {
                "source_donors": folds[held]["source_donors"],
                "npz_sha256": folds[held]["operator_npz_sha256"],
                "json_sha256": folds[held]["operator_json_sha256"],
            }
            for held in RUN.FINAL_FOLDS
        },
        "kolf_confirmation": {
            "supported": True,
            "expected_targets": len(panel),
            "files": kolf_files,
        },
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    args = argparse.Namespace(
        donor_root=donor_root,
        truth_cache=truth_cache,
        kolf_root=(kolf_root if include_panel_truth else root / "missing-kolf-truth"),
        query_cache=query_cache,
        artifact_manifest=manifest_path,
        genes_csv=genes_csv,
        panel_csv=panel_csv,
        output_dir=root / "output",
        ranks=[4],
        betas=[0.02],
        temperatures=[0.15],
        regularizations=[1e-2],
        basis_dim=8,
        epochs=3,
        learning_rate=0.01,
        gradient_clip=5.0,
        split_seed=805,
        optimizer_seed=806,
        null_seed=807,
        null_replicates=2,
        bootstrap_replicates=100,
        bootstrap_seed=808,
        cluster_reciprocal_k=3,
        cluster_cosine_threshold=0.8,
        truth_chunk_rows=8,
        max_rss_gib=0.95,
        guarded=True,
    )
    return args, manifest


def test_failed_pre_gate_never_opens_missing_outer_truth() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_identity_fixture(Path(directory))
        official = RUN.read_csv_column(args.genes_csv, "gene_name")
        panel = RUN.read_csv_column(args.panel_csv, "target_gene")
        _, train_targets, _, axes = RUN.validate_truth_cache(
            args.truth_cache,
            official,
            panel,
            args.donor_root,
            manifest,
        )
        assert not (set(axes.train_score_genes) & set(train_targets))
        assert set(axes.panel_score_genes) & set(train_targets)
        assert not (set(axes.panel_score_genes) & set(panel))
        result = RUN.run_experiment(args, manifest)
        assert result["selection"]["selected"] == RUN.IDENTITY
        assert result["pretest"]["gate"]["passed"] is False
        assert result["pretest"]["null_complete_grid_repetitions"] == 2
        assert result["outer"]["evaluated"] is False
        assert not any(args.truth_cache.glob("*_raw_outer.npy"))
        assert (args.output_dir / "selection_freeze.json").exists()
        assert (args.output_dir / "results.json").exists()
        frozen = json.loads((args.output_dir / "selection_freeze.json").read_text())
        assert len(frozen["complete_selection_nulls"]) == 2
        for null in frozen["complete_selection_nulls"]:
            assert null["selection"]["complete_grid_repeated"] is True
            assert null["selection"]["candidate_count_including_identity"] == 2


def test_missing_kolf_protocol_fails_closed_before_selection() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_identity_fixture(Path(directory))
        manifest["kolf_confirmation"] = {
            "supported": False,
            "reason": "synthetic missing artifact",
            "files": {},
        }
        try:
            RUN.run_experiment(args, manifest)
        except ValueError as error:
            assert "promotion must fail closed" in str(error)
        else:
            raise AssertionError("missing KOLF support did not fail closed")


def test_passing_path_freezes_models_before_truth_and_evaluates_kolf() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_identity_fixture(
            Path(directory), include_panel_truth=True
        )
        args.null_replicates = 1
        candidate = RUN.Candidate(4, 0.02, 0.15, 1e-2)
        original_selection = RUN.complete_candidate_selection
        original_pretest = RUN.evaluate_selected_pretest
        original_verify = RUN.verify_manifest_inputs
        original_clusters = RUN.source_effect_clusters
        authentication_observations: list[bool] = []
        cluster_widths: list[int] = []

        def selected_result(**kwargs):
            replicate = kwargs["permutation_replicate"]
            return {
                "selected": candidate.key,
                "candidate": candidate,
                "selection": {
                    "selected": candidate.key,
                    "candidate_summaries": {
                        candidate.key: {
                            "cell_deltas": [0.01, 0.01, 0.01, 0.01],
                            "cell_median_angles_degrees": [1.0, 1.0, 1.0, 1.0],
                            "minimum_delta": 0.01,
                            "macro_delta": 0.01,
                            "maximum_cell_median_angle_degrees": 1.0,
                        }
                    },
                },
                "cells": [],
                "permutation_replicate": replicate,
                "permutation_plan_sha256": (
                    None
                    if kwargs["permutation_plan"] is None
                    else kwargs["permutation_plan"].plan_sha256
                ),
                "complete_grid_repeated": True,
                "candidate_invariant_fit_cache_reused": True,
                "same_rank_candidates_batched": True,
                "candidate_count_including_identity": 2,
            }

        def pretest_result(**kwargs):
            replicate = kwargs["permutation_replicate"]
            delta = 0.01 if replicate is None else 0.0
            cells = [
                {
                    "outer_donor": fold,
                    "teacher": teacher,
                    "delta": delta,
                }
                for fold in RUN.DEVELOPMENT_FOLDS
                for teacher in RUN.TEACHER_DONORS
            ]
            return {
                "cells": cells,
                "macro_delta": delta,
                "fold_median_angles_degrees": {
                    fold: 1.0 for fold in RUN.DEVELOPMENT_FOLDS
                },
                "all_fold_angles_below_10_degrees": True,
            }

        def verify_with_order_check(local_args, local_manifest, *, include_outer_truth):
            if include_outer_truth:
                freeze = local_args.output_dir / "selection_freeze.json"
                authentication_observations.append(freeze.exists())
                assert freeze.exists()
                frozen = json.loads(freeze.read_text())
                assert frozen["final_models_fitted"] is True
                assert set(frozen["final_model_files"]["real"]) == set(RUN.FINAL_FOLDS)
                assert set(frozen["fold_cluster_files"]) == set(RUN.FINAL_FOLDS)
            return original_verify(
                local_args,
                local_manifest,
                include_outer_truth=include_outer_truth,
            )

        def leakage_safe_clusters(queries, *, reciprocal_k, cosine_threshold):
            assert len(queries) == 1
            cluster_widths.append(int(queries[0].shape[1]))
            assert queries[0].shape[1] == 18
            return original_clusters(
                queries,
                reciprocal_k=reciprocal_k,
                cosine_threshold=cosine_threshold,
            )

        RUN.complete_candidate_selection = selected_result
        RUN.evaluate_selected_pretest = pretest_result
        RUN.verify_manifest_inputs = verify_with_order_check
        RUN.source_effect_clusters = leakage_safe_clusters
        try:
            result = RUN.run_experiment(args, manifest)
        finally:
            RUN.complete_candidate_selection = original_selection
            RUN.evaluate_selected_pretest = original_pretest
            RUN.verify_manifest_inputs = original_verify
            RUN.source_effect_clusters = original_clusters

        assert authentication_observations == [True]
        assert cluster_widths == [18, 18, 18]
        assert result["outer"]["evaluated"] is True
        assert set(result["outer"]["folds"]) == set(RUN.FINAL_FOLDS)
        assert result["outer"]["folds"][RUN.KOLF_FOLD]["gallery"] == len(
            RUN.read_csv_column(args.panel_csv, "target_gene")
        )
        assert (
            result["outer"]["promotion_gate"]["kolf_independent_confirmation_included"]
            is True
        )


def test_no_adapter_fit_occurs_after_panel_truth_authentication_boundary() -> None:
    source = (HERE / "run.py").read_text()
    marker = "verify_manifest_inputs(args, manifest, include_outer_truth=True)"
    assert source.count(marker) == 1
    after_authentication = source.split(marker, 1)[1]
    assert "fit_adapter(" not in after_authentication
    assert "complete_candidate_selection(" not in after_authentication
    assert "fit_and_freeze_final_models(" not in after_authentication


def test_aws_scripts_pin_payloads_and_check_available_memory_and_collisions() -> None:
    launcher = (HERE / "aws" / "launch.sh").read_text()
    runner = (HERE / "aws" / "run.sh").read_text()
    assert "EXPECTED_SOURCE_SHA256" in runner
    assert "EXPECTED_RUNNER_SHA256" in runner
    assert "MemAvailable:" in runner
    assert "minimum_available_kib" in runner
    assert "source_archive_sha256" in launcher
    assert "runner_sha256" in launcher
    assert 'PYTHON_VERSION="3.12.12"' in runner
    assert 'NUMPY_VERSION="2.5.2"' in runner
    assert '--with "numpy==${NUMPY_VERSION}"' in runner
    assert "refusing to reuse existing work directory" in runner
    assert "aws s3api list-objects-v2" in launcher
    assert "refusing reused RUN_ID" in launcher
    assert "EXP080_REAL_DONOR_ROOT" in runner
    assert "EXP080_REAL_KOLF_ROOT" in runner
    assert "EXP080_REAL_TRUTH_CACHE" in runner
    assert "Read-only exp080 memory and process-collision preflight" in launcher
    for marker in (
        "[c]ell-eval",
        "[p]ython(3)? .*experiments/",
        "[u]v (run|tool|x)",
        "[t]orchrun",
        "[m]ake_submission",
        "[r]un_with_memory_guard",
    ):
        assert marker in runner
        assert marker in launcher


def main() -> int:
    tests = [
        test_analytic_listwise_gradient_matches_finite_difference,
        test_target_split_is_hash_stable_disjoint_and_exhaustive,
        test_norm_decile_permutation_is_stratified_and_bijective,
        test_coherent_permutation_plan_is_frozen_across_all_fit_stages,
        test_adapter_training_reduces_listwise_loss_and_preserves_rank_limit,
        test_fit_cache_is_semantics_preserving_for_every_adapter_rank,
        test_batch_equivalence_all_ranks_two_teachers_permuted_and_clipped,
        test_frozen_adapter_roundtrip_preserves_exact_weights,
        test_candidate_selector_enforces_every_cell_and_angle_guard,
        test_angular_guard_is_enforced_separately_in_every_fold,
        test_full_gallery_pds_uses_all_600_candidates_for_subset_queries,
        test_source_clusters_and_fold_specific_bootstrap_do_not_mix_folds,
        test_generalized_pm_slot_preserves_missing_target_rows,
        test_principal_operator_restriction_equals_dense_zeroed_application,
        test_real_shape_contract_has_49_coordinate_less_targets_and_17806_genes,
        test_compact_panel_truth_uses_complete_requested_axis_and_background_center,
        test_blocked_coordinates_cannot_change_safe_centering_or_magnitude,
        test_kolf_only_gain_cannot_satisfy_primary_macro_gate,
        test_arm_d_cache_stage_masks_both_target_panels_and_holds_out_donor,
        test_failed_pre_gate_never_opens_missing_outer_truth,
        test_missing_kolf_protocol_fails_closed_before_selection,
        test_passing_path_freezes_models_before_truth_and_evaluates_kolf,
        test_no_adapter_fit_occurs_after_panel_truth_authentication_boundary,
        test_aws_scripts_pin_payloads_and_check_available_memory_and_collisions,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
