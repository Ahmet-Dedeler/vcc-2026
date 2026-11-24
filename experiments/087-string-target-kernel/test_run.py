#!/usr/bin/env python3
"""Synthetic tests for experiment 087; no real matrices or STRING artifact are opened."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


RUN = load_module("exp087_run_test", "run.py")
STAGE = load_module("exp087_stage_test", "stage_cache.py")


def write_csv(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[column])
        writer.writeheader()
        writer.writerows([{column: value} for value in values])


def write_bank(
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
    (path / "meta.json").write_text(json.dumps(metadata) + "\n")


def with_realistic_duplicate_axis(
    genes: list[str],
    values: np.ndarray,
    *,
    start: int,
) -> tuple[list[str], np.ndarray]:
    """Add 21 duplicated symbols and 22 extra columns, matching Orion counts."""
    duplicate_indices = [(start + offset) % len(genes) for offset in range(21)]
    extra_indices = [*duplicate_indices, duplicate_indices[0]]
    extras = np.asarray(values[:, extra_indices], dtype=np.float64).copy()
    extras += np.arange(1, len(extra_indices) + 1, dtype=np.float64)[None, :]
    return (
        [*genes, *(genes[index] for index in extra_indices)],
        np.concatenate([values, extras], axis=1),
    )


def test_target_split_is_stable_disjoint_and_exhaustive() -> None:
    targets = [f"T{i:03d}" for i in range(101)]
    first = RUN.deterministic_target_split(targets, seed=82)
    second = RUN.deterministic_target_split(targets, seed=82)
    assert {key: value.tolist() for key, value in first.items()} == {
        key: value.tolist() for key, value in second.items()
    }
    groups = [set(value.tolist()) for value in first.values()]
    assert not any(
        left & right
        for index, left in enumerate(groups)
        for right in groups[index + 1 :]
    )
    assert set().union(*groups) == set(range(len(targets)))


def test_frozen_string_candidate_grid_includes_identity_and_rank_8() -> None:
    candidates = RUN.build_candidates(
        [8, 16, 32, 64],
        [1e-3, 1e-2, 1e-1, 1.0, 10.0],
        [0.02, 0.05, 0.10, 0.20],
    )
    assert len(candidates) == 80
    assert {candidate.rank for candidate in candidates} == {8, 16, 32, 64}
    assert all(candidate.key.startswith("string-ridge:") for candidate in candidates)
    assert len(candidates) + 1 == 81


def test_complete_gallery_pds_keeps_nonquery_competitors() -> None:
    scores = np.asarray(
        [
            [0.8, 0.1, 0.9, 0.0, 0.0, 0.0],
            [0.0, 0.7, 0.0, 0.8, 0.0, 0.0],
        ]
    )
    pds = RUN.pds_per_target_from_scores(scores, np.asarray([0, 1]))
    np.testing.assert_allclose(pds, [0.8, 0.8])
    assert not np.allclose(
        pds,
        RUN.pds_per_target_from_scores(scores[:, :2], np.asarray([0, 1])),
    )


def test_embedding_h5_loader_reads_index_and_requested_rows() -> None:
    import h5py

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        path = Path(directory) / "embedding.h5"
        with h5py.File(path, "w") as handle:
            handle.create_dataset("A", data=np.asarray([1.0, 4.0]))
            handle.create_dataset("B", data=np.asarray([2.0, 5.0]))
            handle.create_dataset("C", data=np.asarray([3.0, 6.0]))
        loaded = RUN.load_target_embeddings(
            path,
            ["C", "MISSING", "A"],
            expected_rows=3,
            expected_dim=2,
            read_chunk=1,
        )
        np.testing.assert_allclose(loaded.values[0], [3.0, 6.0])
        np.testing.assert_allclose(loaded.values[2], [1.0, 4.0])
        np.testing.assert_array_equal(loaded.values[1], np.zeros(2))
        assert loaded.available.tolist() == [True, False, True]
        assert loaded.audit["missing_rows_encoded_as_zero"] is True
        assert loaded.audit["missingness_flag"] == "available"


def test_production_manifest_pins_canonical_axes_and_control_hashes() -> None:
    manifest = json.loads((HERE / "artifact_manifest.json").read_text())
    assert manifest["final_score_axis"] == {
        "gene_count": 17_806,
        "gene_order_sha256": (
            "6eaba848c81b77dd1d4b83c0d61b7a8689c1e39762a34da7dcfe26ea1b5d765f"
        ),
        "official_positions_sha256": (
            "a000353b3dbeb8a55bdf7ae3bcc60abd13bbc52f714270413c10f7b1c23e1a02"
        ),
    }
    assert manifest["teacher_score_axis"] == {
        "source_donors": ["rpe1", "nadig_jurkat"],
        "gene_count": 5_582,
        "gene_order_sha256": (
            "015b59dc3d059f984bff2e5b2dd7a08b66978ca9e05cfb6a85b525670d9416d8"
        ),
        "train_mask_array_sha256": (
            "cb95ed14ebc5c9d9266e9dccdc04dd338c97798aea4e86756f6997d7d599da19"
        ),
        "teacher_to_dev_positions_sha256": (
            "7085202ace0b6bbeaf266566203421dc11d0227306c62f7bd8d51a3bed0b6c3d"
        ),
        "teacher_to_panel_positions_sha256": (
            "d3dd86516037d2edbdeec92e95799c5b7a233c5dae966ea7e8695106e63b4b43"
        ),
        "equals_ordered_common_support": True,
    }
    assert manifest["donors"]["k562_gwps"]["b_ctrl.npy"] == (
        "d06cf51d51cde1ae999e86bbd3be0bdc15471392bae153a1aa1e21dc5c81b957"
    )
    assert manifest["cd4"]["cd4_rest"]["b_ctrl.npy"] == (
        "5c47ef6b15a00a867b7b24034273bf3aba118d7cdf87ce58ec58427751313d81"
    )
    assert manifest["exp064_reference"] == {
        "results_sha256": (
            "723c7e94b3aee687b001be435bf7b7a7b55d6a54991404e0ce9dcdd3f59a6113"
        ),
        "gallery_targets": 221,
        "target_order_sha256": (
            "59e231eeae34ac02ac336bef4828bec6bba373224051af78ddb49914c3ec383f"
        ),
        "required_max_abs_error": 1e-6,
    }
    assert manifest["kolf_confirmation"]["expected_targets"] == 273
    assert manifest["embedding"] == {
        "repository": "RasmussenLab/SLIM",
        "repository_type": "git",
        "revision": "5a7e9ade5d0a",
        "filename": "data/gene_string_embeddings.v0.3.h5",
        "local_path": "data/external/slim/gene_string_embeddings.v0.3.h5",
        "citation_readme": "data/external/slim/README.md",
        "citation_readme_sha256": (
            "034e15175ef8472220ab918cac4cac511fad9bc9030350aa4d5309eacccdd6b4"
        ),
        "size_bytes": 19_502_608,
        "sha256": "789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d",
        "rows": 21_688,
        "dimensions": 64,
        "dtype": "float64",
        "layout": "one root-level dataset per gene symbol",
        "root_gene_order_sha256": (
            "881fe6f352f78cc857c6bff2e4b703ffecdff4566bb717cad0b025a6daeab17b"
        ),
        "feature_order_sha256": (
            "b21a85e099b34b6973584305ec84fe9087b9b3b33795acc37c68229a2c39e266"
        ),
        "vcc_panel_coverage": 298,
        "missing_vcc_panel_targets": ["EEF1A2", "KIF21B"],
    }
    assert manifest["development_split"] == {
        "seed": 8202026,
        "training_target_count": 600,
        "training_target_order_sha256": (
            "f97526dfe600acd21853ec217474e6a557d5c5eb8fdc6c81e0bfadfc9e5e788e"
        ),
        "fit": {
            "count": 360,
            "targets_sha256": (
                "17040e750280cf96a8241ec5b3cdb05e32176b0ca3073b386ddc852071b7739c"
            ),
        },
        "tune": {
            "count": 120,
            "targets_sha256": (
                "c7d68ccb71567daff155c4a17ae92c830e91532c7231112db48064aa58549d52"
            ),
        },
        "test": {
            "count": 120,
            "targets_sha256": (
                "b4fdee9f0e5c0ff3b4f9cefcce90261c4b7d27f4b7868c3a7c0eaeaa396ecaf8"
            ),
        },
    }
    assert manifest["held_panel_galleries"] == {
        "orion_hct116": {
            "targets": 290,
            "target_order_sha256": (
                "b4d62e504c487d5cea40548ecb8c1672fbd49c651cc039b182b14ab2325213c4"
            ),
        },
        "orion_hek293t": {
            "targets": 294,
            "target_order_sha256": (
                "b971c471a3533a436cb0f3ed9fcecfb9c4dac786a0a8182b28c1bb44b982a11b"
            ),
        },
        "kolf2": {
            "targets": 273,
            "target_order_sha256": (
                "27c77da3f1a59733468b2b9ee4ac42c7a058d1bf7b71e661fbad68448b893dfd"
            ),
        },
    }
    panel = RUN.read_csv_column(ROOT / "data/challenge/pert_counts.csv", "target_gene")
    assert len(panel) == 300
    for held in RUN.DEVELOPMENT_FOLDS:
        donor_targets = set(RUN.read_lines(ROOT / "data/donors" / held / "perts.txt"))
        complete = tuple(target for target in panel if target in donor_targets)
        assert len(complete) == manifest["held_panel_galleries"][held]["targets"]
        assert (
            RUN.sha256_lines(complete)
            == manifest["held_panel_galleries"][held]["target_order_sha256"]
        )


def test_real_axis_contract_has_109_master_gaps_and_17806_score_genes() -> None:
    manifest = json.loads((HERE / "artifact_manifest.json").read_text())
    official = RUN.read_csv_column(ROOT / "data/challenge/gene_names.csv", "gene_name")
    panel = RUN.read_csv_column(ROOT / "data/challenge/pert_counts.csv", "target_gene")
    cache = ROOT / "experiments/044-csls-distill/cache"
    master = RUN.read_lines(cache / "master_genes.txt")
    train_targets = RUN.read_lines(cache / "train_perts.txt")
    missing_master = tuple(sorted(set(train_targets) - set(master)))
    missing_official = tuple(sorted(set(train_targets) - set(official)))
    assert len(missing_master) == 109
    assert RUN.sha256_lines(missing_master) == (
        "b6820a8f5aa109e2fb7935bf67579e2dc9339d7830537994c40a201bc1a4268f"
    )
    assert len(missing_official) == 49
    assert RUN.sha256_lines(missing_official) == (
        "76fb8028ef4d38c4e880c05fbb21212406c5387bdb77b47fbbc4aebf5c224605"
    )
    score_genes, score_positions, audit = STAGE.established_final_score_axis(
        official,
        panel,
        ROOT / "data/donors",
        manifest,
    )
    assert len(score_genes) == 17_806
    assert len(score_positions) == 17_806
    assert (
        audit["gene_order_sha256"] == manifest["final_score_axis"]["gene_order_sha256"]
    )
    assert (
        audit["official_positions_sha256"]
        == manifest["final_score_axis"]["official_positions_sha256"]
    )
    assert len(set(score_genes) & set(train_targets)) == 551
    _, validated_targets, _, axes = RUN.validate_truth_cache(
        cache,
        official,
        panel,
        ROOT / "data/donors",
        manifest,
    )
    assert validated_targets == train_targets
    assert len(axes.panel_score_genes) == 17_806
    assert len(axes.dev_score_genes) == 6_821
    assert len(axes.teacher_score_genes) == 5_582
    assert RUN.sha256_array(axes.teacher_score_mask) == manifest[
        "teacher_score_axis"
    ]["train_mask_array_sha256"]
    assert RUN.sha256_lines(axes.teacher_score_genes) == manifest[
        "teacher_score_axis"
    ]["gene_order_sha256"]
    assert RUN.sha256_array(axes.teacher_to_dev_positions) == manifest[
        "teacher_score_axis"
    ]["teacher_to_dev_positions_sha256"]
    assert RUN.sha256_array(axes.teacher_to_panel_positions) == manifest[
        "teacher_score_axis"
    ]["teacher_to_panel_positions_sha256"]
    teacher_supports = []
    for teacher in RUN.TEACHER_DONORS:
        teacher_genes = set(
            RUN.read_lines(
                ROOT / "data/donors" / teacher / "genes.txt",
                require_unique=False,
            )
        )
        teacher_supports.append(teacher_genes)
    expected_teacher_genes = tuple(
        gene
        for gene in axes.dev_score_genes
        if all(gene in support for support in teacher_supports)
    )
    assert axes.teacher_score_genes == expected_teacher_genes
    for donor in RUN.DEVELOPMENT_FOLDS:
        rows = RUN.read_lines(
            ROOT / "data/donors" / donor / "genes.txt",
            require_unique=False,
        )
        _, duplicate_audit = RUN.last_occurrence_lookup(rows)
        assert duplicate_audit["axis_rows"] == 38_606
        assert duplicate_audit["unique_symbols"] == 38_584
        assert duplicate_audit["duplicate_symbols"] == 21
        assert duplicate_audit["duplicate_extra_rows"] == 22
        assert axes.donor_axis_audits[donor] == duplicate_audit


def test_blocked_coordinates_do_not_affect_source_statistics_or_safe_effects() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        genes = ["BLOCK", "SAFE1", "SAFE2", "SAFE3"]
        perts = [f"T{i}" for i in range(12)]
        rng = np.random.default_rng(820)
        base = rng.normal(size=(len(perts), len(genes))).astype(np.float32)
        changed = base.copy()
        changed[:, 0] += np.linspace(-10_000, 10_000, len(perts))
        write_bank(root / "a", genes, perts, base)
        write_bank(root / "b", genes, perts, changed)
        first = STAGE.prepare_main_source_safe(
            STAGE.EXP062.Bank.load("a", root / "a"),
            perts,
            genes,
            {"BLOCK"},
            chunk_rows=4,
            max_rss_gib=0.95,
        )
        second = STAGE.prepare_main_source_safe(
            STAGE.EXP062.Bank.load("b", root / "b"),
            perts,
            genes,
            {"BLOCK"},
            chunk_rows=4,
            max_rss_gib=0.95,
        )
        np.testing.assert_array_equal(first.effects, second.effects)
        assert first.raw_median_norm == second.raw_median_norm
        assert np.all(first.effects[:, 0] == 0)


def test_arm_d_formula_matches_direct_062_064_calls_exactly() -> None:
    rng = np.random.default_rng(821)
    rows, genes = 18, 9
    pm = STAGE.EXP062.PreparedSource(
        name="pm",
        effects=rng.normal(size=(rows, genes)),
        target_available=np.ones(rows, dtype=bool),
        mean_on_axis=np.zeros(genes),
        raw_median_norm=1.0,
        mapped_gene_count=genes,
        duplicate_symbols=(),
    )
    companion = STAGE.EXP062.PreparedSource(
        name="companion",
        effects=rng.normal(size=(rows, genes)),
        target_available=np.ones(rows, dtype=bool),
        mean_on_axis=np.zeros(genes),
        raw_median_norm=1.0,
        mapped_gene_count=genes,
        duplicate_symbols=(),
    )
    cd4 = STAGE.EXP062.CD4Prepared(
        mode=STAGE.EXP064.CD4_MODE,
        mean_effect=rng.normal(size=(rows, genes)),
        gene_votes=np.ones((rows, genes), dtype=np.uint8),
        target_condition_count=np.full(rows, len(STAGE.CD4_DONORS), dtype=np.uint8),
        condition_details=[],
    )
    q, _ = np.linalg.qr(rng.normal(size=(genes, genes)))
    operator = STAGE.EXP062.Operator(
        evals=np.linspace(0.5, 2.0, genes),
        evecs=q,
        metadata={},
        path=Path("synthetic.npz"),
        metadata_path=Path("synthetic.json"),
    )
    columns = np.arange(genes, dtype=np.int64)
    score = np.asarray([0, 2, 4, 6, 8], dtype=np.int64)
    observed, _, audit = STAGE.build_arm_d_panel(
        pm,
        [companion],
        cd4,
        operator,
        columns,
        columns,
        np.arange(rows, dtype=np.int64),
        score,
    )
    baseline, _ = STAGE.EXP062.build_baseline([pm, companion])
    whitened, _ = STAGE.EXP062.apply_operator(baseline, operator, columns)
    expected, _ = STAGE.EXP062.add_post_whiten_residual(
        whitened,
        cd4,
        STAGE.EXP064.CD4_FRACTION,
        columns,
    )
    np.testing.assert_array_equal(observed, expected[:, score].astype(np.float32))
    assert audit["operator_restriction"]["dense_zero_equivalence_proved"] is True


def test_principal_operator_restriction_equals_dense_zeroed_application() -> None:
    rng = np.random.default_rng(829)
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
    np.testing.assert_allclose(
        observed[:, active], dense_retained * (before / after), atol=1e-12
    )
    assert np.count_nonzero(observed[:, np.setdiff1d(full_columns, active)]) == 0
    assert audit["dense_zero_equivalence_proved"] is True


def test_real_stage_path_builds_all_folds_without_undefined_helpers() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        train_targets = [f"T{i:02d}" for i in range(20)]
        coordinate_targets = train_targets[:15]
        master_coordinate_targets = coordinate_targets[:12]
        outer_targets = [f"P{i:02d}" for i in range(15)]
        panel_targets = [*outer_targets, "PX00", "PX01", "PX02"]
        features = [f"G{i:02d}" for i in range(16)]
        official = [*panel_targets, *coordinate_targets, *features]
        all_targets = [*train_targets, *panel_targets]
        perts = [*all_targets, "BG0", "BG1", "BG2"]
        rng = np.random.default_rng(830)
        shared = rng.normal(size=(len(perts), len(official)))
        shared[:, : len(panel_targets) + len(train_targets)] += 1000.0

        donor_root = root / "donors"
        for offset, name in enumerate(STAGE.MAIN_DONORS):
            donor_genes = official
            donor_values = shared + rng.normal(
                scale=0.02 + 0.01 * offset,
                size=shared.shape,
            )
            if name in STAGE.DEVELOPMENT_FOLDS:
                donor_genes, donor_values = with_realistic_duplicate_axis(
                    official,
                    donor_values,
                    start=len(panel_targets),
                )
            write_bank(
                donor_root / name,
                donor_genes,
                perts,
                donor_values,
            )
        pm_root = root / "perturbmatch"
        write_bank(
            pm_root,
            official,
            perts,
            shared + rng.normal(scale=0.025, size=shared.shape),
            perturbmatch=True,
        )
        cd4_root = root / "cd4"
        for offset, name in enumerate(STAGE.CD4_DONORS):
            write_bank(
                cd4_root / name,
                official,
                panel_targets,
                shared[len(train_targets) : len(all_targets)]
                + rng.normal(
                    scale=0.03 + 0.01 * offset,
                    size=(len(panel_targets), len(official)),
                ),
                cd4=True,
            )

        genes_csv = root / "gene_names.csv"
        panel_csv = root / "pert_counts.csv"
        write_csv(genes_csv, "gene_name", official)
        write_csv(panel_csv, "target_gene", panel_targets)
        truth_cache = root / "truth"
        truth_cache.mkdir()
        master = [*master_coordinate_targets, *features]
        panel_score_genes = [*coordinate_targets, *features]
        (truth_cache / "master_genes.txt").write_text("\n".join(master) + "\n")
        (truth_cache / "train_perts.txt").write_text("\n".join(train_targets) + "\n")
        (truth_cache / "outer_perts.txt").write_text("\n".join(outer_targets) + "\n")
        np.save(
            truth_cache / "train_mask.npy",
            np.asarray(
                [False] * len(master_coordinate_targets) + [True] * len(features)
            ),
        )
        (truth_cache / "meta.json").write_text('{"fixture": true}\n')

        operator_dir = root / "operators"
        operator_dir.mkdir()
        operators: dict[str, object] = {}
        for held in STAGE.FINAL_FOLDS:
            source_donors = (
                list(STAGE.MAIN_DONORS)
                if held == STAGE.KOLF_FOLD
                else [
                    "k562_gwps",
                    next(name for name in STAGE.DEVELOPMENT_FOLDS if name != held),
                ]
            )
            columns = np.asarray(
                [
                    index
                    for index, gene in enumerate(official)
                    if gene not in set(panel_targets)
                ],
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
                        "wcols_sha256": STAGE.sha256_array(columns),
                        "gamma": STAGE.EXP064.WHITEN_GAMMA,
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
        exp064_results.write_text(json.dumps({"axis": {"target_order": outer_targets}}))

        def hashes(path: Path, filenames: list[str]) -> dict[str, str]:
            return {
                filename: STAGE.sha256_file(path / filename) for filename in filenames
            }

        manifest = {
            "schema_version": 1,
            "challenge": {
                "gene_names.csv": STAGE.sha256_file(genes_csv),
                "pert_counts.csv": STAGE.sha256_file(panel_csv),
            },
            "donors": {
                name: hashes(
                    donor_root / name,
                    ["delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json"],
                )
                for name in STAGE.MAIN_DONORS
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
            },
            "operators": operators,
            "final_score_axis": {
                "gene_count": len(panel_score_genes),
                "gene_order_sha256": STAGE.sha256_lines(panel_score_genes),
                "official_positions_sha256": STAGE.sha256_array(
                    np.asarray(
                        [official.index(gene) for gene in panel_score_genes],
                        dtype=np.int64,
                    )
                ),
            },
            "dependencies": {
                "exp062_run_sha256": STAGE.sha256_file(STAGE.EXP062_PATH),
                "exp064_run_sha256": STAGE.sha256_file(STAGE.EXP064_PATH),
            },
        }
        args = argparse.Namespace(
            donor_root=donor_root,
            perturbmatch=pm_root,
            cd4_root=cd4_root,
            truth_cache=truth_cache,
            operator_dir=operator_dir,
            exp064_results=exp064_results,
            artifact_manifest=root / "manifest.json",
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
        assert set(metadata["folds"]) == set(STAGE.FINAL_FOLDS)
        assert metadata["training_blocked_target_coordinates"] == len(
            coordinate_targets
        ) + len(panel_targets)
        assert metadata["training_targets_without_official_coordinate"] == (
            len(train_targets) - len(coordinate_targets)
        )
        assert metadata["training_targets_without_master_coordinate"] == (
            len(train_targets) - len(master_coordinate_targets)
        )
        assert metadata["final_score_axis"]["gene_count"] == len(panel_score_genes)
        for held in STAGE.DEVELOPMENT_FOLDS:
            duplicate_audit = metadata["final_score_axis"]["donor_gene_axes"][held]
            assert duplicate_audit["duplicate_symbols"] == 21
            assert duplicate_audit["duplicate_extra_rows"] == 22
        assert metadata["panel_blocked_target_coordinates"] == len(panel_targets)
        for held in STAGE.FINAL_FOLDS:
            train_query = np.load(args.output / held / "train_query.npy")
            panel_query = np.load(args.output / held / "panel_query.npy")
            assert train_query.shape == (len(train_targets), len(features))
            expected_rows = len(panel_targets)
            assert panel_query.shape == (expected_rows, len(panel_score_genes))
            saved_targets = RUN.read_lines(args.output / held / "panel_targets.txt")
            assert saved_targets == tuple(panel_targets)
            assert np.any(panel_query[:, : len(coordinate_targets)] != 0.0)
            restriction = metadata["folds"][held]["panels"]["train"][
                "operator_restriction"
            ]
            assert restriction["all_nonpanel_target_coordinates_excluded"] is True
        for held in STAGE.DEVELOPMENT_FOLDS:
            assert metadata["folds"][held]["exp064_reproduction"]["passed"] is True
            assert metadata["folds"][held]["exp064_reproduction"][
                "saved_query_score_axis_exact"
            ]


def test_held_truth_projection_is_uncentered_and_uses_last_duplicate_occurrence() -> (
    None
):
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        genes = ["A", "DUP", "B", "DUP"]
        perts = ["T0", "T1", "BG"]
        delta = np.asarray(
            [
                [1.0, 100.0, 10.0, 4.0],
                [2.0, 200.0, 20.0, 5.0],
                [3.0, 300.0, 30.0, 6.0],
            ],
            dtype=np.float32,
        )
        write_bank(root / "held", genes, perts, delta)
        values, audit = RUN.load_compact_panel_truth(
            root,
            "held",
            ["T0", "T1"],
            ["A", "DUP", "B"],
            chunk_rows=2,
        )
        expected = delta[:2, [0, 3, 2]]
        np.testing.assert_allclose(values, expected)
        assert audit["representation"] == "raw control-referenced donor delta, uncentered"
        assert audit["gene_axis"]["axis_rows"] == 4
        assert audit["gene_axis"]["unique_symbols"] == 3
        assert audit["gene_axis"]["duplicate_symbols"] == 1
        assert audit["gene_axis"]["duplicate_extra_rows"] == 1
        assert audit["selected_source_columns_sha256"] == RUN.sha256_array(
            np.asarray([0, 3, 2], dtype=np.int64)
        )


def test_held_truth_rejects_an_incomplete_panel_gallery() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        write_bank(
            root / "held",
            ["A", "B"],
            ["T0", "T1", "BG"],
            np.asarray(
                [
                    [1.0, 2.0],
                    [3.0, 4.0],
                    [5.0, 6.0],
                ],
                dtype=np.float32,
            ),
        )
        try:
            RUN.load_compact_panel_truth(
                root,
                "held",
                ["T0"],
                ["A", "B"],
                chunk_rows=2,
                complete_panel=["T0", "T1", "MISSING"],
            )
        except ValueError as error:
            assert "complete held-panel gallery" in str(error)
        else:
            raise AssertionError("incomplete held-panel gallery did not fail closed")


def test_teacher_truth_loader_rejects_a_missing_requested_gene() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        write_bank(
            root / "teacher",
            ["A", "B"],
            ["T0", "T1"],
            np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        )
        try:
            RUN.load_truth(
                root,
                "teacher",
                ["T0", "T1"],
                ["A", "MISSING"],
                chunk_rows=1,
            )
        except ValueError as error:
            assert "requested score genes are absent" in str(error)
        else:
            raise AssertionError("teacher truth missing a requested gene did not fail closed")


def make_model_fixture(
    *,
    n: int = 90,
    output_dim: int = 30,
    embedding_dim: int = 12,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray]:
    rng = np.random.default_rng(822)
    embedding = rng.normal(size=(n, embedding_dim)).astype(np.float32)
    biological = rng.normal(size=(embedding_dim, output_dim))
    residual = embedding @ biological * 0.08
    query = rng.normal(size=(n, output_dim))
    truths = {
        "rpe1": query + residual + rng.normal(scale=0.01, size=query.shape),
        "nadig_jurkat": query + residual + rng.normal(scale=0.01, size=query.shape),
    }
    available = np.ones(n, dtype=bool)
    available[7] = False
    embedding[~available] = 0.0
    return query, truths, embedding, available


def test_fit_statistics_use_only_fit_rows() -> None:
    query, truths, embedding, available = make_model_fixture()
    fit_rows = np.arange(60, dtype=np.int64)
    basis = RUN.fit_residual_basis(
        query, {"rpe1": truths["rpe1"]}, available, fit_rows, max_rank=8
    )
    first = RUN.fit_residual_family(embedding, available, basis, ridges=[0.1])
    changed = embedding.copy()
    changed[60:] = 1_000_000
    changed[~available] = 0.0
    second = RUN.fit_residual_family(changed, available, basis, ridges=[0.1])
    np.testing.assert_array_equal(first.embedding_mean, second.embedding_mean)
    np.testing.assert_array_equal(first.embedding_scale, second.embedding_scale)
    np.testing.assert_allclose(first.dual_by_ridge[0.1], second.dual_by_ridge[0.1])


def test_excluded_development_coordinates_cannot_affect_fit_pds_or_selection() -> (
    None
):
    query, truths, embedding, available = make_model_fixture(
        n=75,
        output_dim=8,
        embedding_dim=6,
    )
    rng = np.random.default_rng(8291)
    teacher_positions = np.asarray([0, 2, 3, 6, 8, 9, 11, 13], dtype=np.int64)
    full_queries: dict[str, np.ndarray] = {}
    for offset, outer in enumerate(RUN.OUTER_DONORS):
        full = rng.normal(size=(len(query), 14))
        full[:, teacher_positions] = query + offset * 0.001
        full_queries[outer] = full
    changed_queries = {outer: values.copy() for outer, values in full_queries.items()}
    excluded_positions = np.asarray(
        sorted(set(range(14)) - set(teacher_positions.tolist())),
        dtype=np.int64,
    )
    for offset, values in enumerate(changed_queries.values(), start=1):
        values[:, excluded_positions] = rng.normal(
            loc=offset * 1_000_000.0,
            scale=100_000.0,
            size=(len(query), len(excluded_positions)),
        )
    projected = {
        outer: RUN.project_query_axis(values, teacher_positions)
        for outer, values in full_queries.items()
    }
    projected_changed = {
        outer: RUN.project_query_axis(values, teacher_positions)
        for outer, values in changed_queries.items()
    }
    for outer in RUN.OUTER_DONORS:
        np.testing.assert_array_equal(projected[outer], projected_changed[outer])

    split = RUN.deterministic_target_split(
        [f"T{i:03d}" for i in range(len(query))],
        seed=8291,
    )
    candidates = RUN.build_candidates([2], [0.1, 1.0], [0.02, 0.05])
    first_bases = RUN.build_inner_residual_bases(
        projected,
        truths,
        available,
        split["fit"],
        max_rank=2,
    )
    second_bases = RUN.build_inner_residual_bases(
        projected_changed,
        truths,
        available,
        split["fit"],
        max_rank=2,
    )
    for key in first_bases:
        first = first_bases[key]
        second = second_bases[key]
        np.testing.assert_array_equal(first.residual_mean, second.residual_mean)
        np.testing.assert_array_equal(first.components, second.components)
        np.testing.assert_array_equal(first.coefficients, second.coefficients)

    held_key = (RUN.OUTER_DONORS[0], RUN.TEACHER_DONORS[0])
    first_family = RUN.fit_residual_family(
        embedding,
        available,
        first_bases[held_key],
        ridges=[0.1, 1.0],
    )
    second_family = RUN.fit_residual_family(
        embedding,
        available,
        second_bases[held_key],
        ridges=[0.1, 1.0],
    )
    first_grid = RUN.evaluate_candidate_grid(
        projected[RUN.OUTER_DONORS[0]],
        truths[RUN.TEACHER_DONORS[0]],
        embedding,
        available,
        first_family,
        candidates,
        split["tune"],
        split["tune"],
    )
    second_grid = RUN.evaluate_candidate_grid(
        projected_changed[RUN.OUTER_DONORS[0]],
        truths[RUN.TEACHER_DONORS[0]],
        embedding,
        available,
        second_family,
        candidates,
        split["tune"],
        split["tune"],
    )
    for candidate in candidates:
        np.testing.assert_array_equal(
            first_grid[candidate.key]["per_target_adapted"],
            second_grid[candidate.key]["per_target_adapted"],
        )
        np.testing.assert_array_equal(
            first_grid[candidate.key]["per_target_delta"],
            second_grid[candidate.key]["per_target_delta"],
        )

    first_selection = RUN.run_selection(
        projected,
        truths,
        embedding,
        available,
        split,
        candidates,
        max_rank=2,
        ridges=[0.1, 1.0],
        max_rss_gib=0.95,
        residual_bases=first_bases,
    )
    second_selection = RUN.run_selection(
        projected_changed,
        truths,
        embedding,
        available,
        split,
        candidates,
        max_rank=2,
        ridges=[0.1, 1.0],
        max_rss_gib=0.95,
        residual_bases=second_bases,
    )
    assert first_selection.selected_key == second_selection.selected_key
    assert first_selection.audit == second_selection.audit


def test_fast_scores_match_explicit_final_axis_and_angles() -> None:
    query, truths, embedding, available = make_model_fixture()
    fit_rows = np.arange(60, dtype=np.int64)
    basis = RUN.fit_residual_basis(
        query, {"rpe1": truths["rpe1"]}, available, fit_rows, max_rank=8
    )
    family = RUN.fit_residual_family(embedding, available, basis, ridges=[0.1])
    candidate = RUN.Candidate(rank=4, ridge=0.1, beta=0.1)
    frozen = RUN.freeze_predictor(family, candidate)
    extra = np.random.default_rng(823).normal(size=(len(query), 5))
    full_query = np.column_stack([extra[:, :2], query, extra[:, 2:]])
    full_truth = np.column_stack([extra[:, :2], truths["nadig_jurkat"], extra[:, 2:]])
    positions = np.arange(2, 2 + query.shape[1], dtype=np.int64)
    rows = np.asarray([61, 62, 63, 64], dtype=np.int64)
    fast_scores, _, fast_angles = RUN.candidate_scores(
        full_query[rows],
        full_truth,
        embedding[rows],
        available[rows],
        frozen,
        candidate,
        positions,
    )
    explicit = RUN.explicit_candidate_query(
        full_query[rows],
        embedding[rows],
        available[rows],
        frozen,
        candidate,
        positions,
    )
    expected_scores = RUN.row_normalize(explicit) @ RUN.row_normalize(full_truth).T
    np.testing.assert_allclose(fast_scores, expected_scores, atol=1e-9)
    expected_angles = np.degrees(
        np.arccos(
            np.clip(
                np.sum(
                    RUN.row_normalize(explicit) * RUN.row_normalize(full_query[rows]),
                    axis=1,
                ),
                -1,
                1,
            )
        )
    )
    np.testing.assert_allclose(fast_angles, expected_angles, atol=1e-8)


def test_final_residual_application_changes_only_teacher_axis_positions() -> None:
    query, truths, embedding, available = make_model_fixture(
        n=30,
        output_dim=6,
        embedding_dim=5,
    )
    basis = RUN.fit_residual_basis(
        query,
        truths,
        available,
        np.arange(20, dtype=np.int64),
        max_rank=3,
    )
    family = RUN.fit_residual_family(embedding, available, basis, ridges=[0.1])
    candidate = RUN.Candidate(rank=2, ridge=0.1, beta=0.05)
    predictor = RUN.freeze_predictor(family, candidate)
    rng = np.random.default_rng(8292)
    full_query = rng.normal(size=(5, 11))
    teacher_positions = np.asarray([0, 2, 4, 5, 8, 10], dtype=np.int64)
    full_query[:, teacher_positions] = query[:5]
    adapted = RUN.explicit_candidate_query(
        full_query,
        embedding[:5],
        available[:5],
        predictor,
        candidate,
        teacher_positions,
    )
    excluded_positions = np.asarray(
        sorted(set(range(full_query.shape[1])) - set(teacher_positions.tolist())),
        dtype=np.int64,
    )
    np.testing.assert_array_equal(
        adapted[:, excluded_positions],
        full_query[:, excluded_positions],
    )
    residual = predictor.residual_mean + RUN.frozen_predict_coefficients(
        predictor,
        embedding[:5],
        available[:5],
    ) @ predictor.components
    np.testing.assert_allclose(
        adapted[:, teacher_positions] - full_query[:, teacher_positions],
        candidate.beta * residual,
        atol=1e-12,
    )


def test_frozen_predictor_roundtrip_covers_model_and_identity() -> None:
    query, truths, embedding, available = make_model_fixture()
    basis = RUN.fit_residual_basis(
        query,
        truths,
        available,
        np.arange(60, dtype=np.int64),
        max_rank=4,
    )
    family = RUN.fit_residual_family(embedding, available, basis, ridges=[0.1])
    candidate = RUN.Candidate(rank=2, ridge=0.1, beta=0.05)
    predictor = RUN.freeze_predictor(family, candidate)
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        model_path = root / "model.npz"
        model_record = RUN.save_frozen_predictor(model_path, predictor)
        loaded = RUN.load_frozen_predictor(model_path, model_record)
        assert loaded is not None
        assert RUN.predictor_sha256(loaded) == RUN.predictor_sha256(predictor)

        identity_path = root / "identity.npz"
        identity_record = RUN.save_frozen_predictor(identity_path, None)
        assert RUN.load_frozen_predictor(identity_path, identity_record) is None


def test_nulls_preserve_strata_availability_and_are_deterministic() -> None:
    query, _, embedding, available = make_model_fixture()
    queries = {
        "orion_hct116": query,
        "orion_hek293t": query * 1.01,
    }
    counts = {
        name: np.asarray([1 + (index % 2) for index in range(len(query))])
        for name in queries
    }
    split = RUN.deterministic_target_split(
        [f"T{i:03d}" for i in range(len(query))], seed=824
    )
    strata = RUN.make_source_strata(
        queries,
        counts,
        RUN.partition_labels(len(query), split),
        norm_bins=3,
    )
    first, audit = RUN.permute_embeddings_within_strata(
        embedding, available, strata, seed=824
    )
    second, _ = RUN.permute_embeddings_within_strata(
        embedding, available, strata, seed=824
    )
    np.testing.assert_equal(first, second)
    assert np.all(first[~available] == 0.0)
    assert audit["strata"] == len(set(strata.tolist()))
    gaussian_a, _ = RUN.gaussian_embeddings_like(embedding, available, seed=825)
    gaussian_b, _ = RUN.gaussian_embeddings_like(embedding, available, seed=825)
    np.testing.assert_equal(gaussian_a, gaussian_b)
    assert gaussian_a.shape == embedding.shape
    assert np.all(gaussian_a[~available] == 0.0)


def test_null_pipeline_repeats_complete_selection() -> None:
    query, truths, embedding, available = make_model_fixture(n=75)
    queries = {name: query.copy() for name in RUN.OUTER_DONORS}
    split = RUN.deterministic_target_split(
        [f"T{i:03d}" for i in range(len(query))], seed=826
    )
    counts = {name: np.full(len(query), 2, dtype=np.uint8) for name in RUN.OUTER_DONORS}
    strata = RUN.make_source_strata(
        queries,
        counts,
        RUN.partition_labels(len(query), split),
        norm_bins=2,
    )
    candidates = RUN.build_candidates([2], [0.1], [0.05])
    inner = RUN.build_inner_residual_bases(
        queries, truths, available, split["fit"], max_rank=2
    )
    fit_tune = np.sort(np.concatenate([split["fit"], split["tune"]]))
    pretest = RUN.build_refit_residual_bases(
        queries, truths, available, fit_tune, max_rank=2
    )
    calls = 0
    original = RUN.run_selection

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    RUN.run_selection = counted
    try:
        records, selected = RUN.null_pretest_distribution(
            ("permuted", "gaussian"),
            2,
            train_embeddings=embedding,
            train_available=available,
            train_strata=strata,
            train_queries=queries,
            teacher_truth=truths,
            splits=split,
            candidates=candidates,
            max_rank=2,
            ridges=[0.1],
            null_seed=826,
            max_rss_gib=0.95,
            inner_residual_bases=inner,
            pretest_residual_bases=pretest,
        )
    finally:
        RUN.run_selection = original
    assert calls == 4
    assert set(records) == {"permuted", "gaussian"}
    assert len(selected) == 4
    assert all(
        record["embedding_audit"] for family in records.values() for record in family
    )
    assert all(
        record["complete_grid_repeated"]
        and record["candidate_count_including_identity"] == len(candidates) + 1
        for family in records.values()
        for record in family
    )


def test_fold_specific_cluster_bootstrap_and_strata() -> None:
    rng = np.random.default_rng(827)
    queries = {held: rng.normal(size=(24, 10)) for held in RUN.OUTER_DONORS}
    labels = {}
    audits = {}
    for held in RUN.OUTER_DONORS:
        labels[held], audits[held] = RUN.source_effect_clusters(
            queries[held], reciprocal_k=3, cosine_threshold=0.2
        )
        assert audits[held]["method"].startswith("outer-fold source-only")
    deltas = {
        RUN.OUTER_DONORS[0]: np.linspace(-0.01, 0.03, 24),
        RUN.OUTER_DONORS[1]: np.linspace(0.00, 0.04, 24),
    }
    bootstrap = RUN.outer_fold_cluster_bootstrap(
        deltas, labels, replicates=200, seed=827
    )
    expected = np.mean([np.mean(values) for values in deltas.values()])
    assert np.isclose(bootstrap["point_estimate"], expected)
    strata = RUN.stratum_summaries(
        deltas[RUN.OUTER_DONORS[0]],
        np.asarray([2] * 12 + [1] * 8 + [0] * 4),
        complete_count=2,
    )
    assert strata["complete_source"]["count"] == 12
    assert strata["partial_source"]["count"] == 12
    assert strata["direct_row"]["count"] == 20
    assert strata["missing_source"]["count"] == 4


def make_failed_pregate_fixture(
    root: Path,
    *,
    include_panel_truth: bool = False,
) -> tuple[argparse.Namespace, dict[str, object]]:
    import h5py

    train_targets = [f"T{i:03d}" for i in range(600)]
    outer_targets = [f"P{i:03d}" for i in range(20)]
    features = [f"G{i:03d}" for i in range(12)]
    common_teacher_features = features[:8]
    teacher_features = {
        "rpe1": features[:10],
        "nadig_jurkat": [*features[:8], *features[10:]],
    }
    extra_score_genes = [f"X{i:03d}" for i in range(5)]
    master = [*train_targets, *features]
    panel = [*outer_targets, *[f"PX{i:03d}" for i in range(280)]]
    panel_score_genes = [*master, *extra_score_genes]
    official = [*panel, *panel_score_genes]
    genes_csv = root / "gene_names.csv"
    panel_csv = root / "pert_counts.csv"
    write_csv(genes_csv, "gene_name", official)
    write_csv(panel_csv, "target_gene", panel)

    truth_cache = root / "truth"
    truth_cache.mkdir()
    (truth_cache / "master_genes.txt").write_text("\n".join(master) + "\n")
    (truth_cache / "train_perts.txt").write_text("\n".join(train_targets) + "\n")
    (truth_cache / "outer_perts.txt").write_text("\n".join(outer_targets) + "\n")
    np.save(
        truth_cache / "train_mask.npy",
        np.asarray(
            [False] * 600
            + [gene in set(common_teacher_features) for gene in features]
        ),
    )
    (truth_cache / "meta.json").write_text('{"fixture": true}\n')
    rng = np.random.default_rng(828)
    signal = rng.normal(size=(600, len(features))).astype(np.float32)
    embedding_values = rng.normal(size=(900, 8)).astype(np.float32)
    biological = embedding_values[:600] @ rng.normal(size=(8, len(features)))
    common = rng.normal(scale=0.01, size=signal.shape)
    teacher_values = {
        "rpe1": signal + 0.5 * biological + common,
        "nadig_jurkat": signal - 0.5 * biological + common,
    }
    donor_root = root / "donors"
    for donor in RUN.TEACHER_DONORS:
        selected = np.asarray(
            [features.index(gene) for gene in teacher_features[donor]],
            dtype=np.int64,
        )
        raw = np.zeros((600, len(master)), dtype=np.float32)
        raw[:, 600 + selected] = teacher_values[donor][:, selected]
        np.save(truth_cache / f"{donor}_raw_train.npy", raw)
        write_bank(
            donor_root / donor,
            teacher_features[donor],
            train_targets,
            teacher_values[donor][:, selected],
        )
    held_panel_targets = {
        RUN.DEVELOPMENT_FOLDS[0]: (*outer_targets, *panel[20:23]),
        RUN.DEVELOPMENT_FOLDS[1]: (*outer_targets, *panel[23:28]),
    }
    donor_axis_audits: dict[str, object] = {}
    for donor in RUN.DEVELOPMENT_FOLDS:
        donor_perts = [*held_panel_targets[donor], "BG0", "BG1", "BG2"]
        donor_values = rng.normal(size=(len(donor_perts), len(official))).astype(
            np.float32
        )
        donor_genes, donor_values = with_realistic_duplicate_axis(
            official,
            donor_values,
            start=len(panel),
        )
        write_bank(
            donor_root / donor,
            donor_genes,
            donor_perts,
            donor_values,
        )
        _, donor_axis_audits[donor] = RUN.last_occurrence_lookup(donor_genes)

    query_cache = root / "query"
    query_cache.mkdir()
    folds: dict[str, object] = {}
    dev_to_panel = np.arange(600, 612, dtype=np.int64)
    np.save(query_cache / "dev_to_panel_positions.npy", dev_to_panel)
    for held in RUN.FINAL_FOLDS:
        fold_dir = query_cache / held
        fold_dir.mkdir()
        train_query = signal.copy()
        saved_panel_targets = (
            panel if held == RUN.KOLF_FOLD else held_panel_targets[held]
        )
        panel_rows = len(saved_panel_targets)
        panel_query = rng.normal(size=(panel_rows, len(panel_score_genes))).astype(
            np.float32
        )
        source_count_train = np.full(600, 2, dtype=np.uint8)
        complete_count = 3 if held == RUN.KOLF_FOLD else 2
        source_count_panel = np.resize(
            np.asarray([complete_count, complete_count - 1], dtype=np.uint8),
            panel_rows,
        )
        np.save(fold_dir / "train_query.npy", train_query)
        np.save(fold_dir / "panel_query.npy", panel_query)
        np.save(fold_dir / "train_source_count.npy", source_count_train)
        np.save(fold_dir / "panel_source_count.npy", source_count_panel)
        panel_targets_path = fold_dir / "panel_targets.txt"
        panel_targets_path.write_text("\n".join(saved_panel_targets) + "\n")
        companion = (
            None
            if held == RUN.KOLF_FOLD
            else next(name for name in RUN.OUTER_DONORS if name != held)
        )
        source_donors = (
            ["k562_gwps", *RUN.DEVELOPMENT_FOLDS]
            if held == RUN.KOLF_FOLD
            else ["k562_gwps", companion]
        )
        folds[held] = {
            "held_donor": held,
            "source_donors": source_donors,
            "held_donor_used_in_prediction": False,
            "held_donor_used_in_operator": False,
            "operator_npz_sha256": f"{held}-npz",
            "operator_json_sha256": f"{held}-json",
            "panel_targets": len(saved_panel_targets),
            "panel_target_order_sha256": RUN.sha256_lines(saved_panel_targets),
            "panel_targets_file_sha256": RUN.sha256_file(panel_targets_path),
            "complete_held_panel": held in RUN.DEVELOPMENT_FOLDS,
            "panels": {
                "train": {
                    "query_sha256": RUN.sha256_file(fold_dir / "train_query.npy"),
                    "source_count_sha256": RUN.sha256_file(
                        fold_dir / "train_source_count.npy"
                    ),
                    "operator_restriction": {
                        "all_nonpanel_target_coordinates_excluded": True,
                        "nonpanel_target_overlap_with_active_operator": 0,
                        "dense_zero_equivalence_proved": True,
                    },
                },
                "panel": {
                    "query_sha256": RUN.sha256_file(fold_dir / "panel_query.npy"),
                    "source_count_sha256": RUN.sha256_file(
                        fold_dir / "panel_source_count.npy"
                    ),
                    "score_gene_count": len(panel_score_genes),
                    "score_gene_order_sha256": RUN.sha256_lines(panel_score_genes),
                    "score_official_positions_sha256": RUN.sha256_array(
                        np.asarray(
                            [official.index(gene) for gene in panel_score_genes],
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
                "saved_query_score_axis_exact": True,
                "retrieval_genes": len(panel_score_genes),
                "retrieval_gene_order_sha256": RUN.sha256_lines(panel_score_genes),
            }
    metadata = {
        "schema_version": RUN.CACHE_SCHEMA,
        "experiment": "087-string-target-kernel",
        "arm": "D",
        "development_folds": list(RUN.DEVELOPMENT_FOLDS),
        "final_folds": list(RUN.FINAL_FOLDS),
        "valid_outer_donors": list(RUN.OUTER_DONORS),
        "training_targets": 600,
        "outer_targets": 20,
        "training_outer_overlap": 0,
        "vcc_targets": 300,
        "master_genes": len(master),
        "dev_score_genes": len(features),
        "panel_score_genes": len(panel_score_genes),
        "dev_score_gene_order_sha256": RUN.sha256_lines(features),
        "panel_score_gene_order_sha256": RUN.sha256_lines(panel_score_genes),
        "panel_score_official_positions_sha256": RUN.sha256_array(
            np.asarray(
                [official.index(gene) for gene in panel_score_genes],
                dtype=np.int64,
            )
        ),
        "dev_to_panel_positions_sha256": RUN.sha256_array(dev_to_panel),
        "train_order_sha256": RUN.sha256_lines(train_targets),
        "outer_order_sha256": RUN.sha256_lines(outer_targets),
        "panel_order_sha256": RUN.sha256_lines(panel),
        "implementation": {
            "exp062_run_sha256": "exp062",
            "exp064_run_sha256": "exp064",
            "exact_arm_d_formula": True,
        },
        "coordinate_contract": {
            "development_blocked": {"before_dense_whitening": True},
            "final_panel_blocked": {"score_axis_excludes_exactly_vcc_targets": True},
        },
        "final_score_axis": {
            "gene_count": len(panel_score_genes),
            "gene_order_sha256": RUN.sha256_lines(panel_score_genes),
            "official_positions_sha256": RUN.sha256_array(
                np.asarray(
                    [official.index(gene) for gene in panel_score_genes],
                    dtype=np.int64,
                )
            ),
            "development_donor_axes_equal": True,
            "excludes_exactly_vcc_coordinates_within_held_support": True,
            "donor_gene_axes": donor_axis_audits,
        },
        "folds": folds,
        "panel_query_target_counts": {
            held: int(folds[held]["panel_targets"]) for held in RUN.FINAL_FOLDS
        },
        "panel_query_target_order_sha256": {
            held: str(folds[held]["panel_target_order_sha256"])
            for held in RUN.FINAL_FOLDS
        },
    }
    (query_cache / "metadata.json").write_text(json.dumps(metadata))

    embedding_path = root / "embeddings.h5"
    with h5py.File(embedding_path, "w") as handle:
        for target, values in zip(
            [*train_targets, *panel], embedding_values, strict=True
        ):
            handle.create_dataset(target, data=np.asarray(values, dtype=np.float64))

    kolf_root = root / "kolf2"
    kolf_files: dict[str, str] = {}
    if include_panel_truth:
        kolf_root.mkdir()
        np.save(
            kolf_root / "delta.npy",
            rng.normal(size=(len(panel), len(panel_score_genes))).astype(np.float32),
        )
        np.save(
            kolf_root / "b_ctrl.npy",
            np.ones(len(panel_score_genes), dtype=np.float32),
        )
        (kolf_root / "genes.txt").write_text("\n".join(panel_score_genes) + "\n")
        (kolf_root / "perts.txt").write_text("\n".join(panel) + "\n")
        (kolf_root / "meta.json").write_text('{"fixture": true}\n')
        kolf_files = {
            path.name: RUN.sha256_file(path)
            for path in kolf_root.iterdir()
            if path.is_file()
        }

    synthetic_split = RUN.deterministic_target_split(train_targets, seed=828)
    teacher_mask = np.asarray(
        [False] * 600 + [gene in set(common_teacher_features) for gene in features],
        dtype=bool,
    )
    teacher_to_dev = np.asarray(
        [features.index(gene) for gene in common_teacher_features],
        dtype=np.int64,
    )
    teacher_to_panel = np.asarray(
        [panel_score_genes.index(gene) for gene in common_teacher_features],
        dtype=np.int64,
    )
    manifest = {
        "schema_version": 1,
        "development_split": {
            "seed": 828,
            "training_target_count": len(train_targets),
            "training_target_order_sha256": RUN.sha256_lines(train_targets),
            **{
                name: {
                    "count": len(rows),
                    "targets_sha256": RUN.sha256_lines(
                        train_targets[index] for index in rows.tolist()
                    ),
                }
                for name, rows in synthetic_split.items()
            },
        },
        "challenge": {
            "gene_names.csv": RUN.sha256_file(genes_csv),
            "pert_counts.csv": RUN.sha256_file(panel_csv),
        },
        "truth_cache": {
            "files": {
                path.name: RUN.sha256_file(path)
                for path in truth_cache.iterdir()
                if path.is_file()
            }
        },
        "operators": {
            held: {
                "source_donors": folds[held]["source_donors"],
                "npz_sha256": folds[held]["operator_npz_sha256"],
                "json_sha256": folds[held]["operator_json_sha256"],
            }
            for held in RUN.FINAL_FOLDS
        },
        "donors": {
            donor: {
                path.name: RUN.sha256_file(path)
                for path in (donor_root / donor).iterdir()
                if path.is_file()
            }
            for donor in RUN.DEVELOPMENT_FOLDS
        },
        "teacher_truth_donors": {
            donor: {
                path.name: RUN.sha256_file(path)
                for path in (donor_root / donor).iterdir()
                if path.is_file()
            }
            for donor in RUN.TEACHER_DONORS
        },
        "final_score_axis": {
            "gene_count": len(panel_score_genes),
            "gene_order_sha256": RUN.sha256_lines(panel_score_genes),
            "official_positions_sha256": RUN.sha256_array(
                np.asarray(
                    [official.index(gene) for gene in panel_score_genes],
                    dtype=np.int64,
                )
            ),
        },
        "teacher_score_axis": {
            "source_donors": list(RUN.TEACHER_DONORS),
            "gene_count": len(common_teacher_features),
            "gene_order_sha256": RUN.sha256_lines(common_teacher_features),
            "train_mask_array_sha256": RUN.sha256_array(teacher_mask),
            "teacher_to_dev_positions_sha256": RUN.sha256_array(teacher_to_dev),
            "teacher_to_panel_positions_sha256": RUN.sha256_array(
                teacher_to_panel
            ),
            "equals_ordered_common_support": True,
        },
        "kolf_confirmation": {
            "supported": True,
            "expected_targets": len(panel),
            "files": kolf_files,
        },
        "held_panel_galleries": {
            held: {
                "targets": len(held_panel_targets[held]),
                "target_order_sha256": RUN.sha256_lines(held_panel_targets[held]),
            }
            for held in RUN.DEVELOPMENT_FOLDS
        }
        | {
            RUN.KOLF_FOLD: {
                "targets": len(panel),
                "target_order_sha256": RUN.sha256_lines(panel),
            }
        },
        "embedding": {
            "repository": "fixture",
            "revision": "fixture",
            "filename": embedding_path.name,
            "size_bytes": embedding_path.stat().st_size,
            "sha256": RUN.sha256_file(embedding_path),
        },
        "dependencies": {
            "exp062_run_sha256": "exp062",
            "exp064_run_sha256": "exp064",
        },
    }
    args = argparse.Namespace(
        truth_cache=truth_cache,
        donor_root=donor_root,
        kolf_root=kolf_root,
        query_cache=query_cache,
        embedding_h5=embedding_path,
        artifact_manifest=root / "manifest.json",
        genes_csv=genes_csv,
        panel_csv=panel_csv,
        output_dir=root / "output",
        ranks=[2],
        ridges=[0.1],
        betas=[0.05],
        split_seed=828,
        null_seed=829,
        null_replicates=2,
        null_norm_bins=2,
        bootstrap_replicates=100,
        bootstrap_seed=830,
        cluster_reciprocal_k=3,
        cluster_cosine_threshold=0.2,
        expected_embedding_rows=900,
        expected_embedding_dim=8,
        embedding_read_chunk=4,
        truth_chunk_rows=4,
        max_rss_gib=0.95,
        guarded=True,
    )
    args.artifact_manifest.write_text(json.dumps(manifest))
    return args, manifest


def test_failed_pregate_uses_complete_600_gallery_and_never_opens_outer_truth() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_failed_pregate_fixture(Path(directory))
        original = RUN.build_inner_residual_bases
        observed_partition_freeze: list[dict[str, object]] = []

        def assert_partition_frozen(*call_args, **call_kwargs):
            path = args.output_dir / "partition_freeze.json"
            assert path.exists()
            frozen = json.loads(path.read_text())
            assert frozen["frozen_before_any_fit"] is True
            assert frozen["split_seed"] == args.split_seed
            assert frozen["training_target_count"] == 600
            assert set(frozen["partitions"]) == {"fit", "tune", "test"}
            assert (
                sum(int(entry["count"]) for entry in frozen["partitions"].values())
                == 600
            )
            for entry in frozen["partitions"].values():
                assert entry["targets_sha256"] == RUN.sha256_lines(entry["targets"])
            observed_partition_freeze.append(frozen)
            return original(*call_args, **call_kwargs)

        RUN.build_inner_residual_bases = assert_partition_frozen
        try:
            result = RUN.run_experiment(args, manifest)
        finally:
            RUN.build_inner_residual_bases = original
        assert len(observed_partition_freeze) == 1
        assert result["selection"]["selected"] == RUN.IDENTITY
        assert result["outer"]["evaluated"] is False
        assert not any(args.truth_cache.glob("*_raw_outer.npy"))
        assert {cell["gallery_size"] for cell in result["pretest"]["cells"]} == {600}
        assert result["leakage_audit"][
            "development_queries_use_complete_600_target_gallery"
        ]


def test_missing_kolf_protocol_fails_closed_before_selection() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_failed_pregate_fixture(Path(directory))
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


def test_passing_path_persists_real_and_null_models_before_kolf_authentication() -> (
    None
):
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args, manifest = make_failed_pregate_fixture(
            Path(directory), include_panel_truth=True
        )
        args.null_replicates = 1
        candidate = RUN.Candidate(rank=2, ridge=0.1, beta=0.05)
        cells = tuple(
            {
                "outer_donor": outer,
                "held_teacher": teacher,
                "candidate": candidate.key,
                "delta": 0.01,
                "median_angle_degrees": 1.0,
            }
            for outer in RUN.DEVELOPMENT_FOLDS
            for teacher in RUN.TEACHER_DONORS
        )
        selection = RUN.SelectionResult(
            selected_key=candidate.key,
            selected=candidate,
            audit={"selected": candidate.key},
            cells=cells,
        )

        def fake_selection(*_args, **_kwargs):
            return selection

        def fake_refit(*_args, **_kwargs):
            evaluations = {}
            for outer in RUN.DEVELOPMENT_FOLDS:
                for teacher in RUN.TEACHER_DONORS:
                    evaluations[f"{outer}|{teacher}"] = {
                        "baseline": 0.5,
                        "adapted": 0.51,
                        "delta": 0.01,
                        "per_target_baseline": np.full(120, 0.5),
                        "per_target_adapted": np.full(120, 0.51),
                        "per_target_delta": np.full(120, 0.01),
                        "baseline_norms": np.ones(120),
                        "adapted_norms": np.ones(120),
                        "angles_degrees": np.ones(120),
                        "median_angle_degrees": 1.0,
                        "max_angle_degrees": 1.0,
                        "embedding_available": np.ones(120, dtype=bool),
                        "gallery_size": 600,
                    }
            return evaluations, {}

        def fake_nulls(*_args, **_kwargs):
            records = {
                kind: [
                    {
                        "replicate": 0,
                        "selected": candidate.key,
                        "macro_test_delta": 0.0,
                        "minimum_test_cell_delta": 0.0,
                        "embedding_audit": {"synthetic": True},
                    }
                ]
                for kind in ("permuted", "gaussian")
            }
            selected = {(kind, 0): candidate.key for kind in ("permuted", "gaussian")}
            return records, selected

        original_selection = RUN.run_selection
        original_refit = RUN.evaluate_refit
        original_nulls = RUN.null_pretest_distribution
        original_verify = RUN.verify_manifest_inputs
        authentication_observations: list[bool] = []

        def verify_with_freeze_check(
            local_args, local_manifest, *, include_outer_truth
        ):
            if include_outer_truth:
                freeze_path = local_args.output_dir / "final_model_freeze.json"
                assert freeze_path.exists()
                frozen = json.loads(freeze_path.read_text())
                assert set(frozen["real"]) == set(RUN.FINAL_FOLDS)
                for kind in ("permuted", "gaussian"):
                    assert len(frozen["nulls"][kind]) == 1
                    assert set(frozen["nulls"][kind][0]["models"]) == set(
                        RUN.FINAL_FOLDS
                    )
                    for record in frozen["nulls"][kind][0]["models"].values():
                        assert Path(str(record["path"])).exists()
                authentication_observations.append(True)
            return original_verify(
                local_args,
                local_manifest,
                include_outer_truth=include_outer_truth,
            )

        RUN.run_selection = fake_selection
        RUN.evaluate_refit = fake_refit
        RUN.null_pretest_distribution = fake_nulls
        RUN.verify_manifest_inputs = verify_with_freeze_check
        try:
            result = RUN.run_experiment(args, manifest)
        finally:
            RUN.run_selection = original_selection
            RUN.evaluate_refit = original_refit
            RUN.null_pretest_distribution = original_nulls
            RUN.verify_manifest_inputs = original_verify

        assert authentication_observations == [True]
        assert result["outer"]["evaluated"] is True
        assert set(result["outer"]["folds"]) == set(RUN.FINAL_FOLDS)
        assert result["outer"]["folds"][RUN.KOLF_FOLD]["gallery_size"] == len(
            RUN.read_csv_column(args.panel_csv, "target_gene")
        )
        assert result["outer"]["primary_folds"] == list(RUN.OUTER_DONORS)
        assert set(result["outer"]["primary_fold_deltas"]) == set(RUN.OUTER_DONORS)
        assert result["outer"]["folds"][RUN.DEVELOPMENT_FOLDS[0]]["targets"] == 23
        assert result["outer"]["folds"][RUN.DEVELOPMENT_FOLDS[1]]["targets"] == 25
        assert set(result["outer"]["cluster_bootstrap"]["clusters_by_fold"]) == set(
            RUN.OUTER_DONORS
        )
        expected_primary = np.mean(
            [result["outer"]["fold_deltas"][held] for held in RUN.OUTER_DONORS]
        )
        assert np.isclose(result["outer"]["macro_delta"], expected_primary)
        assert (
            result["outer"]["kolf_confirmation"]["included_in_primary_macro"] is False
        )
        assert (
            result["outer"]["kolf_confirmation"]["included_in_primary_bootstrap"]
            is False
        )
        assert (
            result["outer"]["kolf_confirmation"]["included_in_primary_null_p95"]
            is False
        )
        for held in RUN.DEVELOPMENT_FOLDS:
            gene_axis = result["outer"]["folds"][held]["held_truth"]["gene_axis"]
            assert gene_axis["duplicate_symbols"] == 21
            assert gene_axis["duplicate_extra_rows"] == 22
        for family in result["outer"]["nulls"].values():
            assert family["primary_folds"] == list(RUN.OUTER_DONORS)
            assert family["fold_specific_source_only_strata"] is True
            for record in family["records"]:
                expected_null_macro = np.mean(
                    [record["fold_deltas"][held] for held in RUN.OUTER_DONORS]
                )
                assert np.isclose(record["primary_macro_delta"], expected_null_macro)
                assert record["kolf_delta"] == record["fold_deltas"][RUN.KOLF_FOLD]
                for held in RUN.FINAL_FOLDS:
                    query = np.load(args.query_cache / held / "panel_query.npy")
                    counts = np.load(args.query_cache / held / "panel_source_count.npy")
                    length = len(query)
                    expected_strata = RUN.make_source_strata(
                        {held: query},
                        {held: counts},
                        np.asarray(["panel"] * length),
                        norm_bins=args.null_norm_bins,
                    )
                    audit = record["embedding_audit_by_fold"][held]
                    assert audit["source_fold"] == held
                    assert audit["source_only_strata"] is True
                    assert audit["strata_sha256"] == RUN.sha256_lines(
                        expected_strata.tolist()
                    )
        assert (
            result["outer"]["promotion_gate"]["kolf_independent_confirmation_included"]
            is True
        )


def test_no_model_fit_occurs_after_panel_truth_authentication_boundary() -> None:
    source = (HERE / "run.py").read_text()
    marker = "verify_manifest_inputs(args, manifest, include_outer_truth=True)"
    assert source.count(marker) == 1
    after_authentication = source.split(marker, 1)[1]
    for forbidden in (
        "fit_residual_basis(",
        "fit_residual_family(",
        "build_frozen_predictors(",
        "run_selection(",
    ):
        assert forbidden not in after_authentication


def test_aws_scripts_pin_payloads_and_check_memory_and_collisions() -> None:
    launcher = (HERE / "aws" / "launch.sh").read_text()
    runner = (HERE / "aws" / "run.sh").read_text()
    assert "EXPECTED_SOURCE_SHA256" in runner
    assert "EXPECTED_RUNNER_SHA256" in runner
    assert "source_archive_sha256" in launcher
    assert "runner_sha256" in launcher
    assert "Read-only exp087 memory and process-collision preflight" in launcher
    assert "HELD_OPERATOR_DIR=" in launcher
    assert '"${HELD_OPERATOR_DIR}/held_${held}.${suffix}"' in launcher
    assert (
        "s3://${BUCKET}/output/exp062-cd4-post-whiten/whiteners/held_${held}.${suffix}"
    ) not in launcher
    assert "data/external/slim/gene_string_embeddings.v0.3.h5" in launcher
    assert "data/external/slim/README.md" in launcher
    assert "gene_string_embeddings.v0.3.h5" in runner
    assert "789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d" in (
        runner
    )
    assert "huggingface.co" not in runner
    assert "STRING_pert_features.parquet" not in runner
    assert "--exp064-results" in runner
    assert "--exp064-atol 1e-6" in runner
    assert "--kolf-root" in runner
    assert '--donor-root "$INPUT/donors"' in runner
    assert "rpe1 nadig_jurkat" in runner
    for packaged_axis in (
        "data/donors/orion_hct116/genes.txt",
        "data/donors/orion_hct116/perts.txt",
        "data/donors/orion_hek293t/genes.txt",
        "data/donors/orion_hek293t/perts.txt",
        "data/donors/rpe1/genes.txt",
        "data/donors/rpe1/perts.txt",
        "data/donors/nadig_jurkat/genes.txt",
        "data/donors/nadig_jurkat/perts.txt",
        "experiments/044-csls-distill/cache/meta.json",
        "experiments/044-csls-distill/cache/master_genes.txt",
        "experiments/044-csls-distill/cache/train_perts.txt",
        "experiments/044-csls-distill/cache/outer_perts.txt",
        "experiments/044-csls-distill/cache/train_mask.npy",
    ):
        assert packaged_axis in launcher
    assert "code-sha256.txt" not in launcher
    assert "code-sha256.txt" not in runner
    assert launcher.count('assert_unused_s3_prefix "$SOURCE_PREFIX"') == 1
    assert launcher.count('assert_unused_s3_prefix "$OUTPUT_PREFIX"') == 1
    assert "refusing reused S3 prefix" in launcher
    assert launcher.index('assert_unused_s3_prefix "$SOURCE_PREFIX"') < launcher.index(
        'aws s3 cp "$ARCHIVE"'
    )
    assert "refusing reused or unverifiable output prefix" in runner
    assert runner.index("output_key_count=") < runner.index("upload_outputs()")
    assert "^[0-9]{8}T[0-9]{6}Z$" in launcher
    assert "^[0-9]{8}T[0-9]{6}Z$" in runner
    assert "os.path.realpath(sys.argv[1])" in runner
    assert '[[ -e "$WORK" || -L "$WORK" ]]' in runner
    assert 'rm -rf "$WORK"' not in runner

    base_environment = {
        **os.environ,
        "BUCKET": "fixture",
        "EXPECTED_SOURCE_SHA256": "0" * 64,
        "EXPECTED_RUNNER_SHA256": "0" * 64,
    }
    for unsafe in (".", "..", "20260823", "20260823T010203Z/escape"):
        completed = subprocess.run(
            ["bash", str(HERE / "aws" / "run.sh")],
            env={**base_environment, "RUN_ID": unsafe},
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 63
        assert "unsafe RUN_ID" in completed.stderr

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        work_root = Path(directory).resolve()
        run_id = "20260823T010203Z"
        existing = work_root / run_id
        existing.mkdir()
        sentinel = existing / "sentinel"
        sentinel.write_text("preserve")
        completed = subprocess.run(
            ["bash", str(HERE / "aws" / "run.sh")],
            env={
                **base_environment,
                "RUN_ID": run_id,
                "WORK_ROOT": str(work_root),
            },
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 70
        assert "refusing existing local run directory" in completed.stderr
        assert sentinel.read_text() == "preserve"

        completed = subprocess.run(
            ["bash", str(HERE / "aws" / "run.sh")],
            env={
                **base_environment,
                "RUN_ID": "20260823T010204Z",
                "WORK_ROOT": f"{work_root}/nested/..",
            },
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 70
        assert "absolute canonical non-root path" in completed.stderr
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
        test_target_split_is_stable_disjoint_and_exhaustive,
        test_frozen_string_candidate_grid_includes_identity_and_rank_8,
        test_complete_gallery_pds_keeps_nonquery_competitors,
        test_embedding_h5_loader_reads_index_and_requested_rows,
        test_production_manifest_pins_canonical_axes_and_control_hashes,
        test_real_axis_contract_has_109_master_gaps_and_17806_score_genes,
        test_blocked_coordinates_do_not_affect_source_statistics_or_safe_effects,
        test_arm_d_formula_matches_direct_062_064_calls_exactly,
        test_principal_operator_restriction_equals_dense_zeroed_application,
        test_real_stage_path_builds_all_folds_without_undefined_helpers,
        test_held_truth_projection_is_uncentered_and_uses_last_duplicate_occurrence,
        test_held_truth_rejects_an_incomplete_panel_gallery,
        test_teacher_truth_loader_rejects_a_missing_requested_gene,
        test_fit_statistics_use_only_fit_rows,
        test_excluded_development_coordinates_cannot_affect_fit_pds_or_selection,
        test_fast_scores_match_explicit_final_axis_and_angles,
        test_final_residual_application_changes_only_teacher_axis_positions,
        test_frozen_predictor_roundtrip_covers_model_and_identity,
        test_nulls_preserve_strata_availability_and_are_deterministic,
        test_null_pipeline_repeats_complete_selection,
        test_fold_specific_cluster_bootstrap_and_strata,
        test_failed_pregate_uses_complete_600_gallery_and_never_opens_outer_truth,
        test_missing_kolf_protocol_fails_closed_before_selection,
        test_passing_path_persists_real_and_null_models_before_kolf_authentication,
        test_no_model_fit_occurs_after_panel_truth_authentication_boundary,
        test_aws_scripts_pin_payloads_and_check_memory_and_collisions,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}", flush=True)
    print(f"{len(tests)} tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
