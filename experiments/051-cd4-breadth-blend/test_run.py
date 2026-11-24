from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def load_module():
    spec = importlib.util.spec_from_file_location("exp051_run", HERE / "run.py")
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
    breadth: np.ndarray | None = None,
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", np.asarray(delta, dtype=np.float32))
    np.save(path / "b_ctrl.npy", np.ones(len(genes), dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes))
    (path / "perts.txt").write_text("\n".join(perts))
    if not cd4:
        return
    (path / "meta.json").write_text(
        json.dumps(
            {"effect_definition": ("delta = ln(2) * source DESeq2 log2 fold change")}
        )
    )
    if breadth is None:
        breadth = np.arange(1, len(perts) + 1)
    pd.DataFrame(
        {
            "target_gene": perts,
            "n_total_de_genes": breadth,
            "n_downstream": breadth,
            "ontarget_significant": [True] * len(perts),
            "low_target_gex": [False] * len(perts),
            "distal_offtarget_flag": [False] * len(perts),
            "neighboring_gene_KD": [False] * len(perts),
            "n_guides": [2] * len(perts),
            "n_cells_target": np.arange(len(perts)) + 100,
        }
    ).to_csv(path / "obs.csv", index=False)


def synthetic_inputs(root: Path) -> dict[str, Path]:
    panel = [f"P{i}" for i in range(8)]
    features = [f"G{i}" for i in range(12)]
    official = [*panel, *features]
    background = ["BG0", "BG1"]
    perts = [*panel, *background]
    panel_path = root / "panel.csv"
    genes_path = root / "genes.csv"
    pd.DataFrame({"target_gene": panel}).to_csv(panel_path, index=False)
    pd.DataFrame({"gene_name": official}).to_csv(genes_path, index=False)

    rng = np.random.default_rng(51)
    signatures = rng.normal(size=(len(panel), len(features))).astype(np.float32)
    signatures /= np.linalg.norm(signatures, axis=1, keepdims=True)

    main_paths: dict[str, Path] = {}
    for offset, name in enumerate(("k562_gwps", "orion_hct116", "orion_hek293t")):
        matrix = np.zeros((len(perts), len(official)), dtype=np.float32)
        matrix[: len(panel), : len(panel)] = np.eye(len(panel)) * (100 + offset)
        matrix[: len(panel), len(panel) :] = signatures + rng.normal(
            scale=0.03 + offset * 0.01, size=signatures.shape
        )
        matrix[len(panel) :] = rng.normal(scale=0.02, size=(2, len(official)))
        path = root / name
        write_bank(path, official, perts, matrix, cd4=False)
        main_paths[name] = path

    cd4_root = root / "cd4"
    breadth = np.asarray([1, 2, 4, 8, 16, 32, 64, 128], dtype=float)
    for offset, name in enumerate(RUN.CD4_NAMES):
        matrix = np.zeros((len(panel), len(official)), dtype=np.float32)
        matrix[:, : len(panel)] = np.eye(len(panel)) * (1000 + offset)
        quality = np.linspace(0.35, 1.0, len(panel))[:, None]
        matrix[:, len(panel) :] = signatures * quality
        write_bank(
            cd4_root / name,
            official,
            panel,
            matrix,
            cd4=True,
            breadth=breadth + offset,
        )
    return {
        **main_paths,
        "cd4_root": cd4_root,
        "panel": panel_path,
        "genes": genes_path,
    }


def make_args(paths: dict[str, Path], output: Path) -> argparse.Namespace:
    return argparse.Namespace(
        cd4_root=paths["cd4_root"],
        k562=paths["k562_gwps"],
        hct116=paths["orion_hct116"],
        hek293t=paths["orion_hek293t"],
        panel_csv=paths["panel"],
        genes_csv=paths["genes"],
        output=output,
        scratch_parent=output.parent,
        fractions=[0.2, 1 / 3],
        slopes=[0.5, 1.0],
        qc_strengths=[0.0, 0.5],
        modes=["raw", "gene_norm"],
        cd4_sources=["cd4_stim8hr", "cd4_mean3"],
        weight_lower=0.02,
        weight_upper=0.8,
        tie_tolerance=0.001,
        shuffle_replicates=5,
        seed=51,
        chunk_rows=3,
        whitening_block_rows=2,
        whitener_dir=None,
        whitening_gamma=0.65,
        max_rss_gib=0.0,
        min_shared_targets=8,
        min_master_genes=10,
        guarded=True,
    )


def write_identity_whiteners(
    root: Path,
    panel: list[str],
    feature_genes: list[str],
) -> None:
    for held in RUN.MAIN_NAMES:
        fold = root / held
        fold.mkdir(parents=True)
        sources = [name for name in RUN.MAIN_NAMES if name != held]
        (fold / "genes.txt").write_text("\n".join(feature_genes))
        np.save(fold / "evals.npy", np.ones(len(feature_genes), dtype=np.float64))
        np.save(
            fold / "evecs.npy",
            np.eye(len(feature_genes), dtype=np.float32),
        )
        (fold / "meta.json").write_text(
            json.dumps(
                {
                    "schema": RUN.WHITENER_SCHEMA,
                    "held_out_donor": held,
                    "fit_uses_held_truth": False,
                    "fit_uses_challenge_outcomes": False,
                    "axis_kind": "production_union",
                    "source_donors": sources,
                    "excluded_targets_sha256": RUN.sha256_lines(panel),
                    "operator_genes_sha256": RUN.sha256_lines(feature_genes),
                }
            )
        )


def test_adaptive_weights_are_smooth_monotone_and_mean_matched() -> None:
    confidence = RUN.ConfidenceData(
        source="synthetic",
        breadth_log=np.log1p([1, 2, 4, 8, 16, 32]),
        quality=np.full(6, 0.5),
        provenance={},
    )
    weights = RUN.adaptive_weights(
        confidence,
        fraction=1 / 3,
        slope=1.0,
        qc_strength=0.0,
        lower=0.02,
        upper=0.8,
    )
    assert np.all(np.diff(weights) > 0)
    assert np.isclose(weights.mean(), 1 / 3, atol=1e-10)
    assert weights.min() >= 0.02
    assert weights.max() <= 0.8

    shuffled = confidence.breadth_log[[3, 0, 5, 1, 4, 2]]
    shuffled_weights = RUN.adaptive_weights(
        confidence,
        fraction=1 / 3,
        slope=1.0,
        qc_strength=0.0,
        lower=0.02,
        upper=0.8,
        breadth_override=shuffled,
    )
    np.testing.assert_allclose(np.sort(weights), np.sort(shuffled_weights))


def test_end_to_end_keeps_k562_out_of_selection_and_excludes_panel() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        paths = synthetic_inputs(root)
        output = root / "results.json"
        args = make_args(paths, output)
        panel = pd.read_csv(paths["panel"])["target_gene"].astype(str).tolist()
        features = [
            gene
            for gene in pd.read_csv(paths["genes"])["gene_name"].astype(str)
            if gene not in set(panel)
        ]
        whitener_root = root / "whiteners"
        write_identity_whiteners(whitener_root, panel, features)
        args.whitener_dir = whitener_root
        result = RUN.run(args)

        development_held = {row["held"] for row in result["development_rows"]}
        assert development_held == {"orion_hct116", "orion_hek293t"}
        assert result["protocol"]["selection_uses_k562_pds"] is False
        assert result["protocol"]["k562_holdout_evaluations"] == 1
        assert result["protocol"]["master_union_genes"] == 12
        assert result["protocol"]["whitening"]["enabled"] is True
        assert {row["geometry"] for row in result["development_rows"]} == {
            "none",
            "audited_union:g=0.65",
        }
        assert result["k562_holdout"]["held"] == "k562_gwps"
        assert result["k562_holdout"]["evaluation_count"] == 1
        roles = {row["role"] for row in result["k562_holdout"]["rows"]}
        assert roles == {"selected_fixed", "selected_adaptive"}
        assert (
            result["k562_holdout"]["adaptive_shuffled_breadth_control"]["replicates"]
            == 5
        )
        assert output.exists()
        assert output.with_name("development_rows.csv").exists()
        assert output.with_name("k562_holdout_rows.csv").exists()


def test_whitener_contract_accepts_clean_identity_and_rejects_held_truth() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        fold = root / "orion_hct116"
        fold.mkdir(parents=True)
        genes = ["G1", "G2"]
        (fold / "genes.txt").write_text("\n".join(genes))
        np.save(fold / "evals.npy", np.ones(2, dtype=np.float64))
        np.save(fold / "evecs.npy", np.eye(2, dtype=np.float32))
        metadata = {
            "schema": RUN.WHITENER_SCHEMA,
            "held_out_donor": "orion_hct116",
            "fit_uses_held_truth": False,
            "fit_uses_challenge_outcomes": False,
            "axis_kind": "production_union",
            "source_donors": ["k562_gwps", "orion_hek293t"],
            "excluded_targets_sha256": RUN.sha256_lines(["P1"]),
            "operator_genes_sha256": RUN.sha256_lines(genes),
        }
        (fold / "meta.json").write_text(json.dumps(metadata))
        whitener = RUN.AuditedWhitener.load(
            root,
            "orion_hct116",
            genes,
            genes,
            ["k562_gwps", "orion_hek293t"],
            RUN.sha256_lines(["P1"]),
            gamma=0.65,
            max_rss_gib=0.0,
            block_rows=2,
        )
        matrix = np.asarray([[1.0, 2.0], [-3.0, 4.0]], dtype=np.float32)
        np.testing.assert_allclose(whitener.apply(matrix), matrix)

        # Deliberate leak: the held donor appears among fit sources.
        metadata["source_donors"] = [
            "k562_gwps",
            "orion_hek293t",
            "orion_hct116",
        ]
        (fold / "meta.json").write_text(json.dumps(metadata))
        try:
            RUN.AuditedWhitener.load(
                root,
                "orion_hct116",
                genes,
                genes,
                ["k562_gwps", "orion_hek293t"],
                RUN.sha256_lines(["P1"]),
                gamma=0.65,
                max_rss_gib=0.0,
                block_rows=2,
            )
        except ValueError as error:
            assert "source_donors" in str(error)
        else:
            raise AssertionError("held-truth whitener was accepted")


def test_cli_refuses_unguarded_real_run() -> None:
    parser = RUN.build_parser()
    args = parser.parse_args(["--cd4-root", "/tmp/not-used"])
    try:
        RUN.validate_args(args)
    except SystemExit as error:
        assert "Refusing an unguarded" in str(error)
    else:
        raise AssertionError("unguarded CLI run was accepted")


if __name__ == "__main__":
    test_adaptive_weights_are_smooth_monotone_and_mean_matched()
    test_end_to_end_keeps_k562_out_of_selection_and_excludes_panel()
    test_whitener_contract_accepts_clean_identity_and_rejects_held_truth()
    test_cli_refuses_unguarded_real_run()
    print("4 synthetic checks passed")
