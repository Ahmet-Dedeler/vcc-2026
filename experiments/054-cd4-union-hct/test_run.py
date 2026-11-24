"""Strict synthetic tests for experiment 054.

These tests create only tiny temporary banks, a tiny H5AD, and a tiny audited
whitener. They never open any real VCC, Orion, or CD4 artifact.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent


def load_module():
    spec = importlib.util.spec_from_file_location("exp054_run", HERE / "run.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


RUN = load_module()


def write_bank(
    path: Path,
    genes: list[str],
    perts: list[str],
    delta: np.ndarray,
    *,
    cd4: bool,
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", np.asarray(delta, dtype=np.float32))
    np.save(path / "b_ctrl.npy", np.linspace(1.0, 2.0, len(genes), dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes))
    (path / "perts.txt").write_text("\n".join(perts))
    if cd4:
        (path / "meta.json").write_text(
            json.dumps(
                {
                    "effect_definition": (
                        "delta = ln(2) * source DESeq2 log2 fold change"
                    )
                }
            )
        )


def write_hct(
    path: Path,
    genes: list[str],
    targets: list[str],
    signatures: np.ndarray,
) -> None:
    rows: list[np.ndarray] = []
    labels: list[str] = []
    control = np.full(len(genes), 40, dtype=np.int32)
    for target_index, target in enumerate(targets):
        profile = np.maximum(
            1,
            np.rint(control * np.exp(signatures[target_index])).astype(np.int32),
        )
        rows.extend([profile, profile + 1, profile])
        labels.extend([target] * 3)
    rows.extend([control, control + 1, control, control + 2])
    labels.extend([RUN.CONTROL] * 4)
    ad.AnnData(
        X=sp.csr_matrix(np.vstack(rows), dtype=np.float32),
        obs=pd.DataFrame(
            {RUN.PERT_COL: labels},
            index=[f"cell-{index}" for index in range(len(labels))],
        ),
        var=pd.DataFrame(index=genes),
    ).write_h5ad(path)


def synthetic_fixture(root: Path) -> argparse.Namespace:
    panel = ["P0", "P1", "P2", "P3", "P4"]
    features = [f"G{i}" for i in range(12)]
    genes = [*panel, *features]
    hct_targets = ["P0", "P1", "P2", "P3"]
    backgrounds = ["BG0", "BG1", "BG2"]
    source_perts = [*hct_targets, *backgrounds]

    panel_csv = root / "pert_counts.csv"
    genes_csv = root / "gene_names.csv"
    pd.DataFrame({"target_gene": panel}).to_csv(panel_csv, index=False)
    pd.DataFrame({"gene_name": genes}).to_csv(genes_csv, index=False)

    rng = np.random.default_rng(54)
    signatures = rng.normal(scale=0.12, size=(len(hct_targets), len(genes)))
    signatures[:, : len(panel)] = 0.0
    hct_real = root / "orion-hct116.h5ad"
    write_hct(hct_real, genes, hct_targets, signatures)
    hct_meta = root / "orion-hct116.h5ad.meta.json"
    hct_meta.write_text(json.dumps({"perts": hct_targets}))

    k562 = root / "k562_gwps"
    hek = root / "orion_hek293t"
    for index, (name, path) in enumerate((("k562_gwps", k562), ("orion_hek293t", hek))):
        donor_genes = genes.copy()
        if name == "orion_hek293t":
            # Duplicate G2 checks experiment-031-compatible norm and last-write mapping.
            donor_genes.append("G2")
        matrix = rng.normal(
            scale=0.03,
            size=(len(source_perts), len(donor_genes)),
        )
        matrix[: len(hct_targets), : len(genes)] += signatures
        if len(donor_genes) > len(genes):
            matrix[:, -1] = matrix[:, genes.index("G2")] + 0.01
        write_bank(path, donor_genes, source_perts, matrix, cd4=False)

    cd4_root = root / "cd4_compact"
    coverage = {
        "cd4_rest": ["P0", "P1", "P2"],
        "cd4_stim8hr": ["P0", "P1", "P2"],
        "cd4_stim48hr": ["P0", "P1"],
    }
    cd4_genes = [*panel, *features[:7]]
    for condition_index, name in enumerate(RUN.CD4_NAMES):
        perts = [*coverage[name], *backgrounds]
        matrix = rng.normal(scale=0.02, size=(len(perts), len(cd4_genes)))
        for row, target in enumerate(coverage[name]):
            hct_index = hct_targets.index(target)
            matrix[row] += signatures[hct_index, : len(cd4_genes)] * (
                0.9 + 0.05 * condition_index
            )
        write_bank(cd4_root / name, cd4_genes, perts, matrix, cd4=True)

    # Build exactly the metadata run.py derives, then store a nontrivial
    # orthonormal operator so the absent-target post-whitening assertion matters.
    source_banks = (
        RUN.Bank.load("k562_gwps", k562),
        RUN.Bank.load("orion_hek293t", hek),
    )
    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    wcols = np.asarray(
        sorted(
            index
            for index, gene in enumerate(genes)
            if gene in donor_union and gene not in set(panel)
        ),
        dtype=np.int64,
    )
    q, _ = np.linalg.qr(rng.normal(size=(len(wcols), len(wcols))))
    evals = np.linspace(0.5, 2.0, len(wcols), dtype=np.float64)
    whitener_dir = root / "out_hct_clean"
    whitener_dir.mkdir()
    whitener = whitener_dir / "whitener_union.npz"
    np.savez(whitener, evals=evals, evecs=q.astype(np.float32))
    whitener.with_suffix(".json").write_text(
        json.dumps(
            {
                "axis": "union",
                "donors": list(RUN.SOURCE_NAMES),
                "n_cols": len(wcols),
                "wcols_sha256": RUN.hashlib.sha256(wcols.tobytes()).hexdigest(),
                "excluded_targets": sorted(panel),
            }
        )
    )

    return argparse.Namespace(
        hct_real=hct_real,
        hct_meta=hct_meta,
        k562=k562,
        hek293t=hek,
        cd4_root=cd4_root,
        panel_csv=panel_csv,
        genes_csv=genes_csv,
        whitener=whitener,
        output=root / "results.json",
        scratch_parent=root / "scratch",
        source_chunk_rows=2,
        hct_chunk_rows=2,
        whitening_feature_block=3,
        max_rss_gib=0.95,
        min_hct_targets=4,
        min_cd4_covered_targets=3,
        min_whitener_cols=10,
        guarded=True,
    )


def test_end_to_end_fixed_protocol_and_missingness_controls() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = synthetic_fixture(Path(directory))
        result = RUN.run(args)

        assert result["protocol"]["cd4_fraction"] == 0.5
        assert result["protocol"]["selection_on_hct"] is False
        assert result["protocol"]["count_emission"] is False
        assert result["protocol"]["cell_eval2"] is False
        assert result["axis_audit"]["wcols"] == 12
        assert result["cd4_coverage"]["any_condition_targets"] == 3
        assert result["cd4_coverage"]["all_three_condition_targets"] == 2
        assert result["cd4_coverage"]["absent_targets"] == 1
        assert result["cd4_coverage"]["absent_target_symbols"] == ["P3"]

        candidates = result["candidates"]
        assert len(candidates) == 8
        assert sum(bool(candidate["is_primary"]) for candidate in candidates) == 1
        assert {candidate["cd4_mode"] for candidate in candidates} == {
            "raw",
            "gene_norm",
        }
        assert {candidate["missing_gene_policy"] for candidate in candidates} == {
            "preserve_baseline",
            "zero_vote",
        }
        assert {
            candidate["whitening_renormalization_policy"] for candidate in candidates
        } == {"production_self", "baseline_fixed"}
        for candidate in candidates:
            assert candidate["pre_whiten_absent_target_max_abs_diff"] == 0.0
            if candidate["whitening_renormalization_policy"] == "baseline_fixed":
                assert candidate["post_whiten_absent_target_max_abs_diff"] == 0.0
                assert candidate["post_whiten_absent_target_exact_equal"] is True
                assert (
                    candidate["whitening"]["applied_scale"]
                    == result["baseline"]["whitening"]["applied_scale"]
                )
            else:
                assert (
                    candidate["whitening"]["applied_scale"]
                    == candidate["whitening"]["candidate_self_scale"]
                )
                assert candidate["post_whiten_absent_target_max_abs_diff"] >= 0.0
            assert candidate["pds"]["all_hct_targets"]["n_targets"] == 4
            assert candidate["pds"]["cd4_covered_targets"]["n_targets"] == 3
            assert (
                candidate["pds"]["shared_all_three_condition_targets"]["n_targets"] == 2
            )
        assert args.output.exists()
        assert args.output.with_name("summary.csv").exists()
        assert args.output.with_name("per_target_pds.csv").exists()
        assert RUN.peak_rss_gib() < 1.0
        primary = next(candidate for candidate in candidates if candidate["is_primary"])
        assert primary["whitening_renormalization_policy"] == "production_self"


def test_production_self_can_move_absent_targets_but_baseline_fixed_cannot() -> None:
    rng = np.random.default_rng(540)
    baseline = rng.normal(size=(4, 6))
    candidate = baseline.copy()
    candidate[:3] *= 4.0
    wcols = np.arange(6, dtype=np.int64)

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
        np.save(root / "evals.npy", np.linspace(0.4, 2.0, 6))
        np.save(root / "evecs.npy", q.astype(np.float32))
        operator = RUN.ExtractedWhitener(
            evals_path=root / "evals.npy",
            evecs_path=root / "evecs.npy",
            size=6,
        )
        baseline_post, baseline_info = RUN.apply_whitener(
            baseline,
            operator,
            wcols,
            renormalization_policy="production_self",
            baseline_scale=None,
            feature_block=2,
            max_rss_gib=0.95,
        )
        production_post, production_info = RUN.apply_whitener(
            candidate,
            operator,
            wcols,
            renormalization_policy="production_self",
            baseline_scale=float(baseline_info["applied_scale"]),
            feature_block=2,
            max_rss_gib=0.95,
        )
        diagnostic_post, diagnostic_info = RUN.apply_whitener(
            candidate,
            operator,
            wcols,
            renormalization_policy="baseline_fixed",
            baseline_scale=float(baseline_info["applied_scale"]),
            feature_block=2,
            max_rss_gib=0.95,
        )

        assert production_info["applied_scale"] != baseline_info["applied_scale"]
        assert diagnostic_info["applied_scale"] == baseline_info["applied_scale"]
        assert not np.array_equal(production_post[3], baseline_post[3])
        np.testing.assert_array_equal(diagnostic_post[3], baseline_post[3])


def test_missing_gene_policies_are_structurally_distinct() -> None:
    baseline = np.asarray([[2.0, 4.0, 8.0], [3.0, 6.0, 9.0]])
    cd4 = RUN.CD4Prepared(
        mode="raw",
        mean_effect=np.asarray([[10.0, 20.0, 0.0], [0.0, 0.0, 0.0]]),
        gene_votes=np.asarray([[3, 1, 0], [0, 0, 0]], dtype=np.uint8),
        target_condition_count=np.asarray([3, 0], dtype=np.uint8),
        condition_details=[],
    )
    preserve = RUN.blend_cd4(baseline, cd4, "preserve_baseline")
    zero_vote = RUN.blend_cd4(baseline, cd4, "zero_vote")

    np.testing.assert_allclose(preserve[0], [6.0, 12.0, 8.0])
    np.testing.assert_allclose(zero_vote[0], [6.0, 12.0, 4.0])
    np.testing.assert_array_equal(preserve[1], baseline[1])
    np.testing.assert_array_equal(zero_vote[1], baseline[1])


def test_baseline_matches_existing_multidonor_implementation() -> None:
    from vcc2026.donor_transfer import Donor
    from vcc2026.multi_donor import MultiDonorModel

    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = synthetic_fixture(Path(directory))
        genes = RUN.read_csv_column(args.genes_csv, "gene_name")
        targets = tuple(sorted(json.loads(args.hct_meta.read_text())["perts"]))
        banks = (
            RUN.Bank.load("k562_gwps", args.k562),
            RUN.Bank.load("orion_hek293t", args.hek293t),
        )
        prepared = [
            RUN.prepare_main_source(
                bank,
                targets,
                genes,
                chunk_rows=2,
                max_rss_gib=0.95,
            )
            for bank in banks
        ]
        observed, _ = RUN.build_current_baseline(prepared)

        model = MultiDonorModel(
            {
                "k562_gwps": Donor.load(args.k562),
                "orion_hek293t": Donor.load(args.hek293t),
            },
            np.asarray(genes),
            gene_norm=False,
            eb=False,
        )
        expected = np.vstack(
            [
                model.effect(
                    target,
                    1.0,
                    0.0,
                    weights={name: 1.0 for name in model.names},
                )
                for target in targets
            ]
        )
        # Streaming float64 accumulation differs from the old full-array
        # float32 reduction only at sub-micro-unit roundoff.
        np.testing.assert_allclose(observed, expected, rtol=2e-6, atol=1e-6)


def test_whitener_audit_rejects_any_metadata_drift() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = synthetic_fixture(Path(directory))
        metadata_path = args.whitener.with_suffix(".json")
        metadata = json.loads(metadata_path.read_text())
        metadata["excluded_targets"] = metadata["excluded_targets"][:-1]
        metadata_path.write_text(json.dumps(metadata))

        panel = RUN.read_csv_column(args.panel_csv, "target_gene")
        genes = RUN.read_csv_column(args.genes_csv, "gene_name")
        sources = (
            RUN.Bank.load("k562_gwps", args.k562),
            RUN.Bank.load("orion_hek293t", args.hek293t),
        )
        try:
            RUN.audit_whitener(args.whitener, genes, panel, sources)
        except ValueError as error:
            assert "excluded_targets" in str(error)
        else:
            raise AssertionError("drifted whitener metadata was accepted")


def test_cli_refuses_unguarded_real_run() -> None:
    args = RUN.build_parser().parse_args([])
    try:
        RUN.validate_args(args)
    except SystemExit as error:
        assert "Refusing an unguarded" in str(error)
    else:
        raise AssertionError("unguarded real-data invocation was accepted")


if __name__ == "__main__":
    test_end_to_end_fixed_protocol_and_missingness_controls()
    test_production_self_can_move_absent_targets_but_baseline_fixed_cannot()
    test_missing_gene_policies_are_structurally_distinct()
    test_baseline_matches_existing_multidonor_implementation()
    test_whitener_audit_rejects_any_metadata_drift()
    test_cli_refuses_unguarded_real_run()
    print(f"6 synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
