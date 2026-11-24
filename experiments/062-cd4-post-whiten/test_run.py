"""Synthetic tests for experiment 062; no real donor arrays are read."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("exp062_run", HERE / "run.py")
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
) -> None:
    path.mkdir(parents=True)
    np.save(path / "delta.npy", np.asarray(delta, dtype=np.float32))
    np.save(path / "b_ctrl.npy", np.linspace(1.0, 2.0, len(genes)).astype(np.float32))
    (path / "genes.txt").write_text("\n".join(genes))
    (path / "perts.txt").write_text("\n".join(perts))
    (path / "meta.json").write_text(
        json.dumps(
            {
                "effect_definition": (
                    "delta = ln(2) * source DESeq2 log2 fold change"
                    if cd4
                    else "synthetic centered effect"
                )
            }
        )
    )


def write_csv(path: Path, column: str, values: list[str]) -> None:
    path.write_text(column + "\n" + "\n".join(values) + "\n")


def fixture(root: Path) -> argparse.Namespace:
    panel = ["P0", "P1", "P2", "P3"]
    features = [f"G{i}" for i in range(14)]
    official = [*panel, *features]
    perts = [*panel, "BG0", "BG1", "BG2", "BG3"]
    rng = np.random.default_rng(62)
    shared = rng.normal(scale=0.3, size=(len(perts), len(official)))
    shared[:, : len(panel)] += 100.0
    main_root = root / "main"
    for index, name in enumerate(RUN.MAIN_NAMES):
        write_bank(
            main_root / name,
            official,
            perts,
            shared + rng.normal(scale=0.025 + index * 0.004, size=shared.shape),
        )
    cd4_root = root / "cd4"
    for index, name in enumerate(RUN.CD4_NAMES):
        write_bank(
            cd4_root / name,
            official,
            perts,
            shared + rng.normal(scale=0.03 + index * 0.005, size=shared.shape),
            cd4=True,
        )
    panel_csv = root / "pert_counts.csv"
    genes_csv = root / "gene_names.csv"
    write_csv(panel_csv, "target_gene", panel)
    write_csv(genes_csv, "gene_name", official)
    output_dir = root / "output"
    return argparse.Namespace(
        k562=main_root / "k562_gwps",
        hct116=main_root / "orion_hct116",
        hek293t=main_root / "orion_hek293t",
        cd4_root=cd4_root,
        panel_csv=panel_csv,
        genes_csv=genes_csv,
        output=output_dir / "results.json",
        whitener_dir=output_dir / "whiteners",
        scratch=root / "scratch",
        chunk_rows=2,
        max_rss_gib=0.95,
        expected_panel_targets=4,
        min_shared_targets=4,
        min_operator_genes=10,
        min_retrieval_genes=10,
        instance_id="i-synthetic",
        instance_type="tiny",
        market="synthetic",
        s3_prefix="s3://synthetic/exp062/",
        guarded=True,
    )


def test_selection_rule_is_conservative_and_deterministic() -> None:
    rows = []
    values = {
        ("raw", 0.01): [0.003, 0.004],
        ("raw", 0.02): [0.001, 0.020],
        ("gene_norm", 0.01): [0.003, 0.004],
    }
    for (mode, fraction), deltas in values.items():
        for held, delta in zip(RUN.DEVELOPMENT_HELD, deltas, strict=True):
            rows.append(
                {
                    "held": held,
                    "mode": mode,
                    "fraction": fraction,
                    "delta_vs_baseline": delta,
                }
            )
    selected = RUN.choose_candidate(rows)
    assert selected["mode"] == "gene_norm"
    assert selected["fraction"] == 0.01
    assert selected["worst_fold_delta"] == 0.003


def test_post_whiten_addition_uses_original_axis_and_preserves_direction() -> None:
    baseline = np.asarray([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
    cd4 = RUN.CD4Prepared(
        mode="raw",
        mean_effect=np.asarray([[0.0, 4.0, 9.0], [4.0, 0.0, 9.0]]),
        gene_votes=np.asarray([[1, 1, 0], [1, 1, 0]], dtype=np.uint8),
        target_condition_count=np.asarray([3, 3], dtype=np.uint8),
        condition_details=[],
    )
    candidate, metadata = RUN.add_post_whiten_residual(
        baseline, cd4, 0.1, np.asarray([0, 1], dtype=np.int64)
    )
    raw_candidate = baseline.copy()
    raw_candidate[cd4.gene_votes > 0] += 0.1 * cd4.mean_effect[
        cd4.gene_votes > 0
    ]
    ratios = candidate[raw_candidate != 0] / raw_candidate[raw_candidate != 0]
    np.testing.assert_allclose(ratios, metadata["post_add_global_preserve_scale"])
    assert np.all(candidate[:, 2] == 0.0)


def test_end_to_end_nested_protocol() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = fixture(Path(directory))
        result = RUN.run(args)
        assert result["protocol"]["whitening_gamma"] == 0.65
        assert len(result["development_folds"]) == 2
        assert {
            fold["held"] for fold in result["development_folds"]
        } == set(RUN.DEVELOPMENT_HELD)
        assert all(
            len(fold["candidates"]) == len(RUN.MODES) * len(RUN.FRACTIONS)
            for fold in result["development_folds"]
        )
        final = result["k562_holdout"]
        assert final["held"] == RUN.FINAL_HELD
        assert final["number_of_candidates_evaluated"] == 1
        assert len(final["candidates"]) == 1
        assert final["candidates"][0]["mode"] == result["selection"]["mode"]
        assert (
            final["candidates"][0]["fraction"]
            == result["selection"]["fraction"]
        )
        assert all(
            fold["operator_genes"] == 14
            and fold["retrieval_genes"] == 14
            for fold in [*result["development_folds"], final]
        )
        assert result["leakage_audit"][
            "k562_target_role_scores_used_for_selection"
        ] is False
        assert result["leakage_audit"]["k562_candidate_grid_evaluated"] is False
        assert len(list(args.whitener_dir.glob("*.npz"))) == 3
        assert args.output.exists()
        assert RUN.peak_rss_gib() < 1.0


if __name__ == "__main__":
    test_selection_rule_is_conservative_and_deterministic()
    test_post_whiten_addition_uses_original_axis_and_preserves_direction()
    test_end_to_end_nested_protocol()
    print(f"3 synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")

