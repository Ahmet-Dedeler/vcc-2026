#!/usr/bin/env python3
"""Low-memory synthetic tests for experiment 070."""

from __future__ import annotations

import csv
import gzip
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


RUN = load_module("exp070_run", "run.py")
MANIFEST = load_module("exp070_manifest", "audit_manifest.py")


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
        raw_median_norm=RUN.BASE.EFFECT_SCALE,
        mapped_gene_count=effects.shape[1],
        duplicate_symbols=(),
    )


def write_bank(
    root: Path,
    name: str,
    genes: tuple[str, ...],
    perts: tuple[str, ...],
) -> object:
    path = root / name
    path.mkdir()
    delta = np.zeros((len(perts), len(genes)), dtype=np.float32)
    np.save(path / "delta.npy", delta)
    np.save(path / "b_ctrl.npy", np.zeros(len(genes), dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes) + "\n")
    (path / "perts.txt").write_text("\n".join(perts) + "\n")
    (path / "meta.json").write_text("{}\n")
    return RUN.BASE.Bank(
        name=name,
        path=path,
        genes=genes,
        perts=perts,
        pert_pos={target: index for index, target in enumerate(perts)},
        delta_path=path / "delta.npy",
        b_ctrl_path=path / "b_ctrl.npy",
        shape=delta.shape,
    )


def write_compact(
    path: Path,
    genes: tuple[str, ...],
    perts: tuple[str, ...],
    effects: np.ndarray,
) -> None:
    path.mkdir()
    source_perts = (*perts, "SOURCE_ONLY_1", "SOURCE_ONLY_2")
    source_center = np.asarray([0.5, -0.5], dtype=np.float64)
    source_variance = np.asarray([1.0, 4.0], dtype=np.float64)
    panel_center = effects.mean(axis=0)
    (path / "genes.txt").write_text("\n".join(genes) + "\n")
    (path / "perts.txt").write_text("\n".join(perts) + "\n")
    (path / "source_perts.txt").write_text("\n".join(source_perts) + "\n")
    with gzip.open(path / "effects.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target_gene", *genes])
        for target, row in zip(perts, effects, strict=True):
            writer.writerow([target, *map(str, row)])
    with gzip.open(path / "control.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "control_profile"])
        for gene in genes:
            writer.writerow([gene, "0"])
    for filename, column, values in (
        ("source_center.tsv.gz", "source_center", source_center),
        ("source_variance.tsv.gz", "source_variance", source_variance),
        ("panel_center.tsv.gz", "panel_center", panel_center),
    ):
        with gzip.open(path / filename, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["gene", column])
            for gene, value in zip(genes, values, strict=True):
                writer.writerow([gene, value])
    with (path / "target_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target_gene",
                "manifest_cell_count",
                "observed_cell_count",
                "count_difference",
                "low_cell_manifest",
                "low_cell_observed",
            ],
        )
        writer.writeheader()
        for index, target in enumerate(perts):
            count = 10 if index == 0 else 50
            writer.writerow(
                {
                    "target_gene": target,
                    "manifest_cell_count": count,
                    "observed_cell_count": count,
                    "count_difference": 0,
                    "low_cell_manifest": index == 0,
                    "low_cell_observed": index == 0,
                }
            )
    (path / "metadata_column_audit.csv").write_text(
        "column,target_coverage,target_cells,control_cells,preference\n"
        f"Gene_assignment,{len(perts)},{len(perts)},100,1\n"
    )
    (path / "object_audit.txt").write_text("class=Seurat\n")
    (path / "r_session_info.txt").write_text("R synthetic\n")
    artifact_names = {
        "genes.txt",
        "perts.txt",
        "source_perts.txt",
        "effects.tsv.gz",
        "control.tsv.gz",
        "source_center.tsv.gz",
        "source_variance.tsv.gz",
        "panel_center.tsv.gz",
        "target_audit.csv",
        "metadata_column_audit.csv",
        "object_audit.txt",
        "r_session_info.txt",
    }
    artifacts_md5 = {
        name: hashlib.md5((path / name).read_bytes()).hexdigest()
        for name in artifact_names
    }
    (path / "meta.json").write_text(
        json.dumps(
            {
                "experiment": "070-viperturb-k562",
                "zenodo_record": "18460279",
                "source_file": "genome_wide_filtered.rds",
                "source_md5": RUN.VIP_SOURCE_MD5,
                "manifest_file": "genome_wide_manifest.txt",
                "manifest_md5": RUN.VIP_MANIFEST_MD5,
                "manifest_filtered_vcc_coverage": len(perts),
                "filtered_targets": list(perts),
                "normalization_scope": "all_filtered_source_perturbations",
                "n_source_perturbations": len(source_perts),
                "source_perturbation_order_md5": hashlib.md5(
                    (path / "source_perts.txt").read_bytes()
                ).hexdigest(),
                "source_median_centered_norm": 2.0,
                "panel_median_centered_norm": float(
                    np.median(np.linalg.norm(effects - panel_center, axis=1))
                ),
                "n_perturbations": len(perts),
                "n_genes": len(genes),
                "expression_source": "layer:counts",
                "expression_kind": "counts",
                "effect_definition": "synthetic",
                "target_column": "Gene_assignment",
                "n_control_cells": 100,
                "manifest_low_cell_targets": [perts[0]],
                "observed_low_cell_targets": [perts[0]],
                "low_cell_threshold": RUN.LOW_CELL_THRESHOLD,
                "manifest_observed_count_mismatches": 0,
                "artifacts_md5": artifacts_md5,
            }
        )
        + "\n"
    )


def test_frozen_constants() -> None:
    assert RUN.FOLDS == ("orion_hct116", "orion_hek293t")
    assert RUN.PM_FRACTION == 0.50
    assert RUN.VIP_FRACTION == 0.50
    assert RUN.CD4_MODE == "gene_norm"
    assert RUN.CD4_FRACTION == 0.05
    assert RUN.WHITEN_GAMMA == 0.65
    assert RUN.MIN_MODIFIED_TARGETS == 70
    assert RUN.MIN_NON_LOW_TARGETS == 20
    assert RUN.LOW_CELL_THRESHOLD == 30
    assert RUN.NULL_REPS == 1_000
    assert RUN.MIN_NULL_REPS == 100
    assert RUN.EXPECTED_VIP_COVERAGE == 124


def test_manifest_exact_124_and_low_cell_audit() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        panel = tuple(f"P{index:03d}" for index in range(300))
        with (root / "panel.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["target_gene"])
            writer.writerows([[target] for target in panel])
        with (root / "manifest.txt").open("w") as handle:
            handle.write("gene cell_count file is_in_filtered_obj\n")
            for part in ("partA", "partB", "partC"):
                handle.write(f"NO-TARGET 100 {part} TRUE\n")
            for index, target in enumerate(panel):
                included = index < 124
                count = 8 if index == 0 else 21 if index == 1 else 50
                handle.write(
                    f"{target} {count} partA {'TRUE' if included else 'FALSE'}\n"
                )
        result = MANIFEST.audit_manifest(
            root / "manifest.txt",
            root / "panel.csv",
            expected_filtered_coverage=124,
            expected_source_perturbations=124,
            low_cell_threshold=30,
        )
        assert result["source_perturbations"] == 124
        assert result["manifest_control_labels"] == ["NO-TARGET"]
        assert result["manifest_control_rows"] == 3
        assert result["filtered_vcc_coverage"] == 124
        assert result["low_cell_count"] == 2
        assert result["cell_count"]["bins"]["lt_10"] == 1
        assert result["cell_count"]["bins"]["20_to_29"] == 1
        for name, values in {
            "k562.txt": panel[:100],
            "pm.txt": panel[:90],
            "hct.txt": panel[:200],
            "hek.txt": panel[:180],
        }.items():
            (root / name).write_text("\n".join(values) + "\n")
        with_overlap = MANIFEST.add_overlap_preflight(
            result,
            root / "panel.csv",
            k562_perts_path=root / "k562.txt",
            perturbmatch_perts_path=root / "pm.txt",
            hct116_perts_path=root / "hct.txt",
            hek293t_perts_path=root / "hek.txt",
            minimum_modified_targets=70,
        )
        assert with_overlap["overlap_preflight"]["full_held_gallery_targets"] == 180
        assert with_overlap["overlap_preflight"]["expected_modified_targets"] == 90


def test_source_catalogue_matches_observed_seurat_v5_schema() -> None:
    rows = [
        {
            "gene": "NO-TARGET",
            "cell_count": 7_949,
            "file": part,
            "is_in_filtered_obj": True,
        }
        for part in ("partA", "partB", "partC")
    ]
    rows.extend(
        [
            {
                "gene": "A",
                "cell_count": 11,
                "file": "partA",
                "is_in_filtered_obj": True,
            },
            {
                "gene": "B",
                "cell_count": 17,
                "file": "partB",
                "is_in_filtered_obj": True,
            },
            {
                "gene": "NOT_FILTERED",
                "cell_count": 20,
                "file": "partC",
                "is_in_filtered_obj": False,
            },
        ]
    )
    catalogue = MANIFEST.filtered_source_catalogue(
        rows,
        expected_source_perturbations=2,
    )
    assert catalogue == {
        "targets": ("A", "B"),
        "control_labels": ("NO-TARGET",),
        "control_rows": 3,
    }

    duplicated = [
        *rows,
        {
            "gene": "A",
            "cell_count": 9,
            "file": "partC",
            "is_in_filtered_obj": True,
        },
    ]
    try:
        MANIFEST.filtered_source_catalogue(
            duplicated,
            expected_source_perturbations=3,
        )
    except ValueError as error:
        assert "duplicated perturbations: A" in str(error)
    else:
        raise AssertionError("duplicated real perturbation was accepted")


def test_manifest_rejects_wrong_filtered_coverage() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        (root / "panel.csv").write_text("target_gene\nA\nB\n")
        (root / "manifest.txt").write_text(
            "gene cell_count file is_in_filtered_obj\n"
            "NO-TARGET 100 partA TRUE\n"
            "NO-TARGET 100 partB TRUE\n"
            "NO-TARGET 100 partC TRUE\n"
            "A 20 partA TRUE\n"
            "B 30 partA FALSE\n"
        )
        try:
            MANIFEST.audit_manifest(
                root / "manifest.txt",
                root / "panel.csv",
                expected_filtered_coverage=2,
                expected_source_perturbations=1,
            )
        except ValueError as error:
            assert "coverage" in str(error)
        else:
            raise AssertionError("wrong manifest coverage was accepted")


def test_compact_loader_authenticates_and_aligns() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        compact = root / "compact"
        genes = ("G1", "G2")
        perts = ("P1", "P2", "P3")
        write_compact(
            compact,
            genes,
            perts,
            np.asarray([[1.0, 0.0], [0.0, 2.0], [-1.0, -2.0]]),
        )
        prepared_vip, panel_vip, audit = RUN.load_viperturb_compact(
            compact,
            ("P1", "P2", "P3", "P4"),
            ("P1", "P2", "P3", "P4", "G1", "G2"),
            np.asarray([0, 1, 2, 3]),
            expected_coverage=3,
            expected_source_md5=RUN.VIP_SOURCE_MD5,
            expected_manifest_md5=RUN.VIP_MANIFEST_MD5,
        )
        np.testing.assert_array_equal(
            prepared_vip.target_available, [True, True, True, False]
        )
        np.testing.assert_array_equal(prepared_vip.effects[:, :4], 0.0)
        expected_source = (
            np.asarray([[1.0, 0.0], [0.0, 2.0], [-1.0, -2.0]]) - np.asarray([0.5, -0.5])
        ) * (RUN.BASE.EFFECT_SCALE / 2.0)
        np.testing.assert_allclose(prepared_vip.effects[:3, 4:], expected_source)
        assert not np.allclose(prepared_vip.effects, panel_vip.effects)
        assert audit.target_rows["P1"]["low_cell_manifest"] is True
        np.testing.assert_allclose(audit.source_variance[4:], [1.0, 4.0])
        assert len(audit.compact_hashes) == 12


def test_compact_loader_rejects_panel_only_normalization() -> None:
    with tempfile.TemporaryDirectory() as directory:
        compact = Path(directory) / "compact"
        genes = ("G1", "G2")
        perts = ("P1", "P2", "P3")
        write_compact(compact, genes, perts, np.eye(3, 2))
        metadata_path = compact / "meta.json"
        metadata = json.loads(metadata_path.read_text())
        metadata["normalization_scope"] = "retained_panel_rows"
        metadata_path.write_text(json.dumps(metadata) + "\n")
        try:
            RUN.load_viperturb_compact(
                compact,
                perts,
                (*perts, *genes),
                np.asarray([0, 1, 2]),
                expected_coverage=3,
                expected_source_md5=RUN.VIP_SOURCE_MD5,
                expected_manifest_md5=RUN.VIP_MANIFEST_MD5,
            )
        except ValueError as error:
            assert "normalization_scope" in str(error)
        else:
            raise AssertionError("panel-only normalization metadata was accepted")


def test_compact_loader_accepts_only_recognized_legacy_r_scale() -> None:
    with tempfile.TemporaryDirectory() as directory:
        compact = Path(directory) / "compact"
        genes = ("G1", "G2")
        perts = ("P1", "P2", "P3")
        effects = np.asarray(
            [
                [1.0, 10.0],
                [2.0, 20.0],
                [4.0, 40.0],
            ]
        )
        write_compact(compact, genes, perts, effects)
        metadata_path = compact / "meta.json"
        metadata = json.loads(metadata_path.read_text())
        panel_center = effects.mean(axis=0)
        canonical = RUN.panel_median_centered_norm(effects, panel_center)
        legacy = RUN.legacy_r_recycled_panel_norm(effects, panel_center)
        assert not np.isclose(canonical, legacy, rtol=1e-6, atol=1e-6)
        resolved, validation, legacy_value = RUN.resolve_panel_median_norm(
            effects,
            panel_center,
            round(canonical, 4),
        )
        assert resolved == canonical
        assert validation == "canonical_jsonlite_four_decimal"
        assert legacy_value is None
        metadata["panel_median_centered_norm"] = round(legacy, 4)
        metadata_path.write_text(json.dumps(metadata) + "\n")

        _, panel_vip, audit = RUN.load_viperturb_compact(
            compact,
            ("P1", "P2", "P3", "P4"),
            ("P1", "P2", "P3", "P4", "G1", "G2"),
            np.asarray([0, 1, 2, 3]),
            expected_coverage=3,
            expected_source_md5=RUN.VIP_SOURCE_MD5,
            expected_manifest_md5=RUN.VIP_MANIFEST_MD5,
        )
        assert panel_vip.raw_median_norm == canonical
        assert (
            audit.metadata["panel_median_centered_norm_validation"]
            == "legacy_r_recycling_jsonlite_four_decimal"
        )
        assert audit.metadata["panel_median_centered_norm_stored"] == round(legacy, 4)

        metadata["panel_median_centered_norm"] = legacy + 0.001
        metadata_path.write_text(json.dumps(metadata) + "\n")
        try:
            RUN.load_viperturb_compact(
                compact,
                ("P1", "P2", "P3", "P4"),
                ("P1", "P2", "P3", "P4", "G1", "G2"),
                np.asarray([0, 1, 2, 3]),
                expected_coverage=3,
                expected_source_md5=RUN.VIP_SOURCE_MD5,
                expected_manifest_md5=RUN.VIP_MANIFEST_MD5,
            )
        except ValueError as error:
            assert "panel-centered scale artifact is inconsistent" in str(error)
        else:
            raise AssertionError("unrecognized panel scale mismatch was accepted")


def test_assay_balanced_slot_and_fallback() -> None:
    incumbent = prepared(
        "incumbent",
        np.asarray([[2.0, 4.0], [10.0, 20.0]]),
        np.asarray([True, True]),
    )
    perturbmatch = prepared(
        "pm",
        np.asarray([[6.0, 8.0], [30.0, 40.0]]),
        np.asarray([True, True]),
    )
    vip = prepared(
        "vip",
        np.asarray([[10.0, 12.0], [50.0, 60.0]]),
        np.asarray([True, False]),
    )
    pm_slot, _ = RUN.build_pm_slot_with_fallback(
        incumbent, perturbmatch, np.asarray([0, 1])
    )
    candidate, modified, audit = RUN.build_vip_slot(
        incumbent,
        perturbmatch,
        pm_slot,
        vip,
        np.asarray([0, 1]),
        ("P1", "P2"),
    )
    np.testing.assert_allclose(pm_slot.effects[0], [4.0, 6.0])
    np.testing.assert_allclose(candidate.effects[0], [7.0, 9.0])
    np.testing.assert_allclose(candidate.effects[1], pm_slot.effects[1])
    np.testing.assert_array_equal(modified, [True, False])
    assert audit["modified_target_count"] == 1


def test_fixed_scales_preserve_unmodified_rows() -> None:
    operator = RUN.BASE.Operator(
        evals=np.ones(2),
        evecs=np.eye(2),
        metadata={},
        path=Path("synthetic.npz"),
        metadata_path=Path("synthetic.json"),
    )
    baseline = np.asarray([[1.0, 2.0], [3.0, 4.0]])
    candidate = baseline.copy()
    candidate[0] = [5.0, 6.0]
    baseline_post, whitening = RUN.BASE.apply_operator(
        baseline, operator, np.asarray([0, 1])
    )
    candidate_post = RUN.apply_operator_with_frozen_scale(
        candidate,
        operator,
        np.asarray([0, 1]),
        whitening["applied_scale"],
    )
    cd4 = RUN.BASE.CD4Prepared(
        mode="gene_norm",
        mean_effect=np.asarray([[0.2, -0.1], [0.4, 0.3]]),
        gene_votes=np.ones((2, 2), dtype=np.uint8),
        target_condition_count=np.asarray([3, 3], dtype=np.uint8),
        condition_details=[],
    )
    baseline_d, preserve = RUN.BASE.add_post_whiten_residual(
        baseline_post, cd4, RUN.CD4_FRACTION, np.asarray([0, 1])
    )
    candidate_d = RUN.add_cd4_with_frozen_scale(
        candidate_post, cd4, preserve["post_add_global_preserve_scale"]
    )
    np.testing.assert_array_equal(candidate_d[1], baseline_d[1])


def test_fold_pipeline_keeps_unmodified_target_exact() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        panel = ("P1", "P2", "P3")
        official = (*panel, "G1", "G2", "G3")
        banks = {
            name: write_bank(root, name, official, panel)
            for name in ("k562_gwps", "orion_hct116", "orion_hek293t")
        }
        held_truth = np.zeros((3, 6))
        held_truth[:, 3:] = np.eye(3)
        companion = np.zeros((3, 6))
        companion[:, 3:] = np.asarray(
            [[0.8, 0.1, 0.0], [0.0, 0.8, 0.1], [0.1, 0.0, 0.8]]
        )
        pm_effects = companion.copy()
        vip_effects = companion.copy()
        vip_effects[[0, 2], 3:] = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        prepared_sources = {
            "orion_hct116": prepared(
                "orion_hct116", held_truth, np.ones(3, dtype=bool)
            ),
            "orion_hek293t": prepared(
                "orion_hek293t", companion, np.ones(3, dtype=bool)
            ),
        }
        pm_slot = prepared("pm_slot", pm_effects, np.ones(3, dtype=bool))
        modified = np.asarray([True, False, True])
        vip_source = prepared("vip_source", vip_effects, modified)
        candidate_effects = pm_effects.copy()
        candidate_effects[modified] = (1.0 - RUN.VIP_FRACTION) * pm_effects[
            modified
        ] + RUN.VIP_FRACTION * vip_effects[modified]
        vip_slot = prepared("vip_slot", candidate_effects, np.ones(3, dtype=bool))
        cd4 = RUN.BASE.CD4Prepared(
            mode="gene_norm",
            mean_effect=np.asarray(
                [
                    [0, 0, 0, 0.02, -0.01, 0.03],
                    [0, 0, 0, -0.01, 0.01, 0.02],
                    [0, 0, 0, 0.01, 0.02, -0.02],
                ],
                dtype=np.float64,
            ),
            gene_votes=np.ones((3, 6), dtype=np.uint8),
            target_condition_count=np.asarray([3, 3, 3], dtype=np.uint8),
            condition_details=[],
        )
        wcols = np.asarray([3, 4, 5], dtype=np.int64)
        whitener = root / "held_orion_hct116.npz"
        evecs, _ = np.linalg.qr(
            np.asarray(
                [
                    [1.0, 2.0, 3.0],
                    [2.0, -1.0, 1.0],
                    [0.5, 1.5, -2.0],
                ]
            )
        )
        np.savez(
            whitener,
            evals=np.asarray([0.5, 1.0, 2.0]),
            evecs=evecs.astype(np.float32),
        )
        whitener.with_suffix(".json").write_text(
            json.dumps(
                {
                    "axis": "union",
                    "donors": ["k562_gwps", "orion_hek293t"],
                    "gamma": 0.65,
                    "n_cols": 3,
                    "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
                }
            )
        )
        vip_audit = RUN.VIPAudit(
            metadata={},
            target_rows={
                "P1": {
                    "manifest_cell_count": 10,
                    "observed_cell_count": 10,
                    "low_cell_manifest": True,
                },
                "P3": {
                    "manifest_cell_count": 50,
                    "observed_cell_count": 50,
                    "low_cell_manifest": False,
                },
            },
            compact_hashes={},
            source_variance=np.ones(6),
        )
        fold, rows, delta, cache = RUN.evaluate_fold(
            "orion_hct116",
            banks,
            prepared_sources,
            pm_slot,
            vip_source,
            vip_slot,
            modified,
            wcols,
            cd4,
            official,
            panel,
            panel,
            vip_audit,
            whitener,
            expected_whitener_sha256="",
            max_rss_gib=1.0,
            min_operator_genes=1,
            min_retrieval_genes=1,
        )
        assert fold["gallery_targets"] == 3
        assert rows[1]["modified"] is False
        assert rows[1]["delta"] == 0.0
        assert delta[1] == 0.0
        assert fold["low_cell_audit"]["modified_low_cell_targets"] == 1
        replay = RUN.candidate_pds_from_vip_rows(
            cache,
            vip_source.effects[np.ix_(modified, wcols)],
            pm_slot.effects[np.ix_(modified, wcols)],
        )
        np.testing.assert_array_equal(
            replay - cache.baseline_per_target,
            delta,
        )


def test_nulls_are_deterministic_stratified_and_bounded() -> None:
    baseline_logits = np.eye(4)
    cache = RUN.FoldCache(
        held="synthetic",
        baseline_logits=baseline_logits,
        vip_projection=np.asarray(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        ),
        row_factors=np.ones(4),
        baseline_per_target=RUN.pds_per_target_from_logits(baseline_logits),
    )
    caches = {held: cache for held in RUN.FOLDS}
    vip_rows = np.eye(4)
    pm_rows = np.zeros((4, 4))
    quartiles = np.zeros(4, dtype=np.int64)
    first = RUN.run_target_label_null(
        caches,
        vip_rows,
        pm_rows,
        quartiles,
        reps=20,
        seed=123,
        max_rss_gib=1.0,
    )
    second = RUN.run_target_label_null(
        caches,
        vip_rows,
        pm_rows,
        quartiles,
        reps=20,
        seed=123,
        max_rss_gib=1.0,
    )
    np.testing.assert_array_equal(first, second)
    assert np.unique(first).size > 1

    values = np.asarray([5, 5, 1, 9, 2, 8, 3, 7], dtype=float)
    labels = RUN.balanced_rank_strata(values, 4)
    permutation = RUN.stratified_permutation(labels, np.random.default_rng(7))
    np.testing.assert_array_equal(labels[permutation], labels)

    gene_first = RUN.run_gene_coordinate_null(
        caches,
        vip_rows,
        pm_rows,
        np.asarray([0, 0, 1, 1]),
        reps=20,
        seed=456,
        max_rss_gib=1.0,
    )
    gene_second = RUN.run_gene_coordinate_null(
        caches,
        vip_rows,
        pm_rows,
        np.asarray([0, 0, 1, 1]),
        reps=20,
        seed=456,
        max_rss_gib=1.0,
    )
    np.testing.assert_array_equal(gene_first, gene_second)


def test_bootstrap_sensitivity_and_stopping_gate_are_frozen() -> None:
    deltas = {
        "orion_hct116": np.full(70, 0.004),
        "orion_hek293t": np.full(70, 0.004),
    }
    bootstrap = RUN.paired_target_cluster_bootstrap(deltas, reps=2_000, seed=70)
    macro = float(
        np.mean([deltas["orion_hct116"].mean(), deltas["orion_hek293t"].mean()])
    )
    gate = RUN.stopping_gate(
        70,
        {
            "orion_hct116": float(deltas["orion_hct116"].mean()),
            "orion_hek293t": float(deltas["orion_hek293t"].mean()),
        },
        macro,
        bootstrap,
        {"p95": 0.001},
        {"p95": 0.002},
        RUN.low_cell_sensitivity(
            deltas,
            np.ones(70, dtype=bool),
            np.zeros(70, dtype=bool),
            minimum_non_low_targets=20,
        ),
    )
    assert gate["passes"] is True
    assert (
        RUN.stopping_gate(
            69,
            {"orion_hct116": 0.004, "orion_hek293t": 0.004},
            0.004,
            {"ci95_percentile": [0.001, 0.007]},
            {"p95": 0.001},
            {"p95": 0.002},
            {"gain_persists": True, "minimum_non_low_targets": 20},
        )["passes"]
        is False
    )
    assert (
        RUN.stopping_gate(
            70,
            {"orion_hct116": 0.004, "orion_hek293t": 0.004},
            0.004,
            {"ci95_percentile": [0.001, 0.007]},
            {"p95": 0.001},
            {"p95": 0.002},
            {"gain_persists": False, "minimum_non_low_targets": 20},
        )["passes"]
        is False
    )
    assert (
        RUN.stopping_gate(
            70,
            {"orion_hct116": 0.004, "orion_hek293t": 0.004},
            0.004,
            {"ci95_percentile": [0.001, 0.007]},
            {"p95": 0.004},
            {"p95": 0.002},
            {"gain_persists": True, "minimum_non_low_targets": 20},
        )["passes"]
        is False
    )


def test_r_extractor_contract_is_compact_and_audited() -> None:
    source = (HERE / "extract_viperturb.R").read_text()
    assert "readRDS(input)" in source
    assert "LayerData" in source
    assert "GetAssayData" in source
    assert "EXPECTED_FILTERED_COVERAGE <- 124L" in source
    assert "EXPECTED_SOURCE_PERTURBATIONS <- 6724L" in source
    assert "source_catalogue_from_manifest" in source
    assert "targets <- labels[!control_rows]" in source
    assert "duplicate_targets <- unique(targets[duplicated(targets)])" in source
    assert "unmapped_non_controls" in source
    assert "group_levels <- c(CONTROL_LABEL, source_targets)" in source
    assert 'normalization_scope = "all_filtered_source_perturbations"' in source
    assert "source_center.tsv.gz" in source
    assert "source_variance.tsv.gz" in source
    assert "panel_center.tsv.gz" in source
    assert "source_median_centered_norm" in source
    assert "panel_median_centered_norm" in source
    assert "panel_centered <- sweep(" in source
    assert "digits = NA" in source
    assert "group_chunk_size" in source
    assert "target_audit.csv" in source
    assert "low_cell_manifest" in source
    assert "selected_matrix %*% indicator" in source
    assert "writeMM" not in source
    assert "SaveH5Seurat" not in source
    assert "Convert(" not in source


def test_r_schema_regression_fixture_is_wired() -> None:
    fixture = (HERE / "test_extract_schema.R").read_text()
    assert 'c("NO-TARGET", "NO-TARGET", "NO-TARGET", "A", "B"' in fixture
    assert "nCount_RNA" in fixture
    assert "nFeature_RNA" in fixture
    assert "nCount_GDO" in fixture
    assert "guide =" in fixture
    assert "gene =" in fixture
    assert "sample =" in fixture
    assert "nFeature_GDO" in fixture
    assert 'identical(assignment$column, "gene")' in fixture
    assert "duplicated perturbations: A" in fixture


def test_binary_r_environment_is_pinned() -> None:
    environment = (HERE / "environment.yml").read_text()
    assert "conda-forge" in environment
    assert "r-base=4.4.3" in environment
    assert "r-seurat=5.3.1" in environment
    assert "r-seuratobject=5.1.0" in environment


def test_aws_launcher_and_runner_are_guarded() -> None:
    launcher = (HERE / "aws/launch.sh").read_text()
    runner = (HERE / "aws/run.sh").read_text()
    evaluate_launcher = (HERE / "aws/evaluate_only.sh").read_text()
    evaluate_runner = (HERE / "aws/evaluate_only_remote.sh").read_text()
    assert 'I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES"' in launcher
    assert "INSTANCE_ID must name an already-running worker" in launcher
    assert "run-instances" not in launcher
    assert "start-instances" not in launcher
    assert "vcc submit" not in launcher
    assert "vcc submit" not in runner
    assert "exec 9>/work/vcc-heavy.lock" in runner
    assert "require_free_disk" in runner
    assert 'MIN_FREE_GIB="${MIN_FREE_GIB:-24}"' in runner
    assert RUN.VIP_SOURCE_MD5 in runner
    assert "3611371127" in runner
    assert "md5sum --check" in runner
    assert "micromamba create" in runner
    assert "curl-minimal" in runner
    assert "dnf install -y bzip2 gzip tar time util-linux" not in runner
    assert "environment.yml" in runner
    assert "extract_viperturb.R" in runner
    assert "test_extract_schema.R" in runner
    assert "--minimum-modified-targets 70" in runner
    assert '--null-reps "$NULL_REPS"' in runner
    assert "--max-rss-gib 1.0" in runner
    assert "uv run cell-eval2" not in runner
    assert "README.md" in launcher
    assert "--setenv=NULL_REPS=${NULL_REPS}" in launcher
    assert 'I_HAVE_EXPLICIT_CLEARANCE:-}" != "YES"' in evaluate_launcher
    assert "run-instances" not in evaluate_launcher
    assert "start-instances" not in evaluate_launcher
    assert "vcc submit" not in evaluate_launcher
    assert "vcc submit" not in evaluate_runner
    assert "genome_wide_filtered.rds" not in evaluate_runner
    assert "extract_viperturb.R" not in evaluate_runner
    assert 'BASE_RUN_ID="${BASE_RUN_ID:-20260823T063300Z}"' in evaluate_launcher
    assert 'COMPACT="${BASE_WORK}/output/viperturb_compact"' in evaluate_runner
    assert "curl -fL" not in evaluate_runner
    assert "Rscript" not in evaluate_runner
    assert "exec 9>/work/vcc-heavy.lock" in evaluate_runner
    assert "--max-rss-gib 1.0" in evaluate_runner


if __name__ == "__main__":
    tests = [
        value
        for name, value in sorted(globals().items())
        if name.startswith("test_") and callable(value)
    ]
    for test in tests:
        test()
    print(f"{len(tests)} synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
