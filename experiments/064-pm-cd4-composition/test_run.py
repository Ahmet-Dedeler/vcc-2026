"""Synthetic tests for experiment 064; no real donor arrays are read."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp064_run", HERE / "run.py")
assert SPEC is not None and SPEC.loader is not None
RUN = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUN
SPEC.loader.exec_module(RUN)


def write_bank(
    path: Path,
    genes: list[str],
    perts: list[str],
    delta: np.ndarray,
    *,
    cd4: bool = False,
    perturbmatch: bool = False,
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", np.asarray(delta, dtype=np.float32))
    if not perturbmatch:
        np.save(
            path / "b_ctrl.npy",
            np.linspace(1.0, 2.0, len(genes)).astype(np.float32),
        )
    (path / "genes.txt").write_text("\n".join(genes))
    (path / "perts.txt").write_text("\n".join(perts))
    metadata: dict[str, object] = {}
    if cd4:
        metadata["effect_definition"] = (
            "delta = ln(2) * source DESeq2 log2 fold change"
        )
    if perturbmatch:
        metadata.update(
            {
                "source_md5": "0017076a5cc7b9979febf3bd19bd2157",
                "vcc_target_coverage": 251,
            }
        )
    (path / "meta.json").write_text(json.dumps(metadata))


def write_csv(path: Path, column: str, values: list[str]) -> None:
    path.write_text(column + "\n" + "\n".join(values) + "\n")


def fixture(root: Path) -> argparse.Namespace:
    panel = ["P0", "P1", "P2", "P3"]
    features = [f"G{i}" for i in range(14)]
    official = [*panel, *features]
    perts = [*panel, "BG0", "BG1", "BG2", "BG3"]
    rng = np.random.default_rng(64)
    shared = rng.normal(scale=0.3, size=(len(perts), len(official)))
    shared[:, : len(panel)] += 100.0
    main_root = root / "main"
    for index, name in enumerate(
        ("k562_gwps", "orion_hct116", "orion_hek293t")
    ):
        write_bank(
            main_root / name,
            official,
            perts,
            shared + rng.normal(scale=0.02 + index * 0.004, size=shared.shape),
        )
    write_bank(
        root / "perturbmatch",
        official,
        perts,
        shared + rng.normal(scale=0.025, size=shared.shape),
        perturbmatch=True,
    )
    cd4_root = root / "cd4"
    for index, name in enumerate(RUN.CD4_NAMES):
        write_bank(
            cd4_root / name,
            official,
            perts,
            shared + rng.normal(scale=0.03 + index * 0.004, size=shared.shape),
            cd4=True,
        )
    panel_csv = root / "pert_counts.csv"
    genes_csv = root / "gene_names.csv"
    write_csv(panel_csv, "target_gene", panel)
    write_csv(genes_csv, "gene_name", official)

    whitener_dir = root / "whiteners"
    whitener_dir.mkdir()
    for held in RUN.FOLDS:
        companion = next(name for name in RUN.FOLDS if name != held)
        path = whitener_dir / f"held_{held}.npz"
        np.savez(
            path,
            evals=np.linspace(0.5, 1.5, len(features)),
            evecs=np.eye(len(features), dtype=np.float32),
        )
        (path.with_suffix(".json")).write_text(
            json.dumps(
                {
                    "axis": "union",
                    "donors": ["k562_gwps", companion],
                    "n_cols": len(features),
                    "wcols_sha256": RUN.hashlib.sha256(
                        np.arange(
                            len(panel), len(official), dtype=np.int64
                        ).tobytes()
                    ).hexdigest(),
                    "gamma": 0.65,
                }
            )
        )
    output = root / "output" / "results.json"
    return argparse.Namespace(
        k562=main_root / "k562_gwps",
        hct116=main_root / "orion_hct116",
        hek293t=main_root / "orion_hek293t",
        perturbmatch=root / "perturbmatch",
        cd4_root=cd4_root,
        panel_csv=panel_csv,
        genes_csv=genes_csv,
        whitener_dir=whitener_dir,
        output=output,
        chunk_rows=2,
        max_rss_gib=0.95,
        expected_panel_targets=4,
        expected_shared_targets=4,
        min_pm_shared_genes=10,
        min_operator_genes=10,
        min_retrieval_genes=10,
        bootstrap_reps=500,
        bootstrap_seed=6402026,
        expected_hct_whitener_sha256="",
        expected_hek_whitener_sha256="",
        instance_id="i-synthetic",
        instance_type="tiny",
        market="synthetic",
        s3_prefix="s3://synthetic/exp064/",
        guarded=True,
    )


def test_pm_slot_blends_only_shared_columns() -> None:
    incumbent = RUN.BASE.PreparedSource(
        name="old",
        effects=np.asarray([[2.0, 4.0, 6.0]]),
        target_available=np.asarray([True]),
        mean_on_axis=np.zeros(3),
        raw_median_norm=1.0,
        mapped_gene_count=3,
        duplicate_symbols=(),
    )
    pm = RUN.BASE.PreparedSource(
        name="pm",
        effects=np.asarray([[10.0, 20.0, 30.0]]),
        target_available=np.asarray([True]),
        mean_on_axis=np.zeros(3),
        raw_median_norm=1.0,
        mapped_gene_count=3,
        duplicate_symbols=(),
    )
    slot = RUN.build_pm_slot(incumbent, pm, np.asarray([0, 2]))
    np.testing.assert_allclose(slot.effects, [[6.0, 4.0, 18.0]])


def test_factorial_effects_are_exact() -> None:
    effects = RUN.arm_statistics({"A": 1.0, "B": 1.2, "C": 1.3, "D": 1.7})
    assert np.isclose(effects["pm_main_without_cd4_B_minus_A"], 0.2)
    assert np.isclose(effects["cd4_main_without_pm_C_minus_A"], 0.3)
    assert np.isclose(effects["interaction_D_minus_B_minus_C_plus_A"], 0.2)


def test_paired_bootstrap_point_estimate() -> None:
    values = {
        "orion_hct116": np.asarray([0.1, -0.1, 0.3]),
        "orion_hek293t": np.asarray([0.2, 0.0, 0.1]),
    }
    result = RUN.paired_target_cluster_bootstrap(values, reps=500, seed=64)
    expected = np.mean(np.stack(list(values.values()), axis=1))
    assert np.isclose(result["estimate"], expected)
    assert result["targets"] == 3


def test_end_to_end_frozen_factorial() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = fixture(Path(directory))
        result = RUN.run(args)
        assert result["protocol"]["design"] == "frozen 2x2 factorial; no tuning"
        assert result["protocol"]["perturbmatch_fraction"] == 0.5
        assert result["protocol"]["cd4_fraction"] == 0.05
        assert result["axis"]["gallery_targets"] == 4
        assert len(result["folds"]) == 2
        assert all(set(fold["scores"]) == set(RUN.ARMS) for fold in result["folds"])
        for fold in result["folds"]:
            expected = (
                fold["scores"]["D"]
                - fold["scores"]["B"]
                - fold["scores"]["C"]
                + fold["scores"]["A"]
            )
            assert np.isclose(
                fold["effects"]["interaction_D_minus_B_minus_C_plus_A"],
                expected,
            )
        assert result["leakage_audit"]["hyperparameters_retuned"] is False
        assert result["leakage_audit"]["k562_truth_evaluated"] is False
        assert (args.output.parent / "fold_summary.csv").exists()
        assert (args.output.parent / "per_target_pds.csv").exists()
        assert RUN.peak_rss_gib() < 1.0


if __name__ == "__main__":
    test_pm_slot_blends_only_shared_columns()
    test_factorial_effects_are_exact()
    test_paired_bootstrap_point_estimate()
    test_end_to_end_frozen_factorial()
    print(f"4 synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")

