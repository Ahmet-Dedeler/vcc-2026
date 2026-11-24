from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


EXTRACT = load_module("extract_cd4", "extract_cd4.py")
EVALUATE = load_module("evaluate_cd4_lodo", "evaluate_lodo.py")


def test_extract_converts_scale_and_resolves_duplicate_gene(tmp_path: Path) -> None:
    conditions = list(EXTRACT.CONDITIONS)
    obs = pd.DataFrame(
        {
            "target_contrast_gene_name": ["A", "B"] * len(conditions),
            "target_contrast": ["ENSGA", "ENSGB"] * len(conditions),
            "culture_condition": np.repeat(conditions, 2),
            "n_cells_target": [100, 90] * len(conditions),
            "n_guides": [2, 2] * len(conditions),
            "ontarget_significant": [True, True] * len(conditions),
        },
        index=[f"row_{i}" for i in range(6)],
    )
    var = pd.DataFrame(
        {
            "gene_name": ["A", "B", "DUP", "DUP"],
            "gene_ids": ["ENSGA", "ENSGB", "ENSGD1", "ENSGD2"],
        },
        index=["va", "vb", "vd1", "vd2"],
    )
    log_fc = np.arange(24, dtype=np.float32).reshape(6, 4) / 10
    base_mean = np.ones((6, 4), dtype=np.float32)
    base_mean[:, 2] = 2
    base_mean[:, 3] = 20
    source = ad.AnnData(X=np.zeros_like(log_fc), obs=obs, var=var)
    source.layers["log_fc"] = log_fc
    source.layers["zscore"] = log_fc + 1
    source.layers["lfcSE"] = np.full_like(log_fc, 0.2)
    source.layers["adj_p_value"] = np.full_like(log_fc, 0.01)
    source.layers["p_value"] = np.full_like(log_fc, 0.001)
    source.layers["baseMean"] = base_mean
    input_path = tmp_path / "source.h5ad"
    source.write_h5ad(input_path)

    panel_path = tmp_path / "panel.csv"
    genes_path = tmp_path / "genes.csv"
    pd.DataFrame({"target_gene": ["B", "A"]}).to_csv(panel_path, index=False)
    pd.DataFrame({"gene_name": ["A", "B", "DUP", "MISS"]}).to_csv(
        genes_path, index=False
    )
    output = tmp_path / "out"
    args = argparse.Namespace(
        input=str(input_path),
        output=str(output),
        panel_csv=str(panel_path),
        genes_csv=str(genes_path),
        row_block=1,
        max_rss_gib=0.0,
        min_targets=2,
        min_union_targets=2,
        min_genes=3,
        min_finite_fraction=1.0,
        overwrite=False,
        dry_run=False,
    )
    manifest = EXTRACT.run_extraction(args)

    rest = output / "cd4_rest"
    delta = np.load(rest / "delta.npy")
    assert delta.shape == (2, 3)
    np.testing.assert_allclose(delta[:, :2], log_fc[[1, 0], :2] * np.log(2))
    np.testing.assert_allclose(delta[:, 2], log_fc[[1, 0], 3] * np.log(2))
    assert (rest / "perts.txt").read_text().splitlines() == ["B", "A"]
    assert "" not in (rest / "genes.txt").read_text().split("\n")
    assert "" not in (rest / "perts.txt").read_text().split("\n")
    gene_map = pd.read_csv(output / "gene_map.csv")
    assert gene_map.loc[gene_map["gene_name"] == "DUP", "source_var_index"].item() == 3
    assert manifest["panel_union_coverage"] == 2


def test_pds_and_candidate_blending() -> None:
    truth = np.eye(8)
    assert EVALUATE.pds(truth, truth) == 1.0
    baseline_parts = {"a": truth, "b": truth}
    cd4_parts = {
        "cd4_rest": truth,
        "cd4_stim8hr": truth,
        "cd4_stim48hr": truth,
    }
    candidates = EVALUATE.candidate_predictions(baseline_parts, cd4_parts, [1 / 3])
    assert "cd4_mean3:fraction=0.333333" in candidates
    np.testing.assert_allclose(candidates["baseline"], truth)


def test_bank_resolves_duplicate_gene_by_stronger_basal_column(
    tmp_path: Path,
) -> None:
    path = tmp_path / "duplicate_bank"
    path.mkdir()
    genes = ["G1", "DUP", "DUP", "G2"]
    perts = ["P1", "P2"]
    np.save(path / "delta.npy", np.zeros((2, 4), dtype=np.float32))
    np.save(path / "b_ctrl.npy", np.asarray([1.0, 0.1, 2.0, 1.0], dtype=np.float32))
    (path / "genes.txt").write_text("\n".join(genes))
    (path / "perts.txt").write_text("\n".join(perts))

    bank = EVALUATE.Bank.load("duplicate", path)
    assert bank.duplicate_gene_symbols == 1
    assert bank.gene_pos["DUP"] == 2


def test_lodo_runner_excludes_all_panel_coordinates(tmp_path: Path) -> None:
    panel = ["P1", "P2", "P3", "P4"]
    genes = [*panel, "G1", "G2", "G3"]
    panel_path = tmp_path / "panel.csv"
    genes_path = tmp_path / "genes.csv"
    pd.DataFrame({"target_gene": panel}).to_csv(panel_path, index=False)
    pd.DataFrame({"gene_name": genes}).to_csv(genes_path, index=False)

    base = np.zeros((4, 7), dtype=np.float32)
    base[:, 4:] = np.asarray(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2], [-1, -1, -1]], dtype=np.float32
    )

    def write_bank(path: Path, delta: np.ndarray, *, cd4: bool) -> None:
        path.mkdir(parents=True)
        np.save(path / "delta.npy", delta)
        np.save(path / "b_ctrl.npy", np.ones(len(genes), dtype=np.float32))
        (path / "genes.txt").write_text("\n".join(genes))
        (path / "perts.txt").write_text("\n".join(panel))
        if cd4:
            (path / "meta.json").write_text(
                '{"effect_definition":"delta = ln(2) * source DESeq2 log2 fold change"}'
            )

    main_paths = {}
    for i, name in enumerate(("k562", "hct116", "hek293t")):
        path = tmp_path / name
        write_bank(path, base + i * 0.01, cd4=False)
        main_paths[name] = path
    cd4_root = tmp_path / "cd4"
    for i, name in enumerate(EVALUATE.CD4_NAMES):
        write_bank(cd4_root / name, base + i * 0.005, cd4=True)

    args = argparse.Namespace(
        panel_csv=str(panel_path),
        genes_csv=str(genes_path),
        k562=str(main_paths["k562"]),
        hct116=str(main_paths["hct116"]),
        hek293t=str(main_paths["hek293t"]),
        cd4_root=str(cd4_root),
        cd4_fractions=[1 / 3],
        chunk_rows=2,
        max_rss_gib=0.0,
        min_shared_targets=4,
        min_shared_genes=3,
        output=str(tmp_path / "results.json"),
    )
    result = EVALUATE.run(args)
    assert result["protocol"]["genes"] == 3
    assert result["protocol"]["perturbations"] == 4
    assert result["somatic_selected_then_k562"] is not None


if __name__ == "__main__":
    # Keep even temporary test writes inside this experiment directory.
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        test_extract_converts_scale_and_resolves_duplicate_gene(Path(directory))
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        test_lodo_runner_excludes_all_panel_coordinates(Path(directory))
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        test_bank_resolves_duplicate_gene_by_stronger_basal_column(Path(directory))
    test_pds_and_candidate_blending()
    print("4 synthetic checks passed")
