"""Small synthetic checks for experiment 060; no real donor arrays are read."""

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
SPEC = importlib.util.spec_from_file_location("exp060_run", HERE / "run.py")
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


def fixture(root: Path) -> argparse.Namespace:
    panel = ["P0", "P1", "P2", "P3"]
    features = [f"G{i}" for i in range(14)]
    official = [*panel, *features]
    perts = [*panel, "BG0", "BG1", "BG2"]
    rng = np.random.default_rng(60)
    signal = rng.normal(scale=0.2, size=(len(perts), len(official)))
    signal[:, : len(panel)] += 20.0

    hct = root / "hct"
    hek = root / "hek"
    k562 = root / "k562"
    write_bank(hct, official, perts, signal + rng.normal(scale=0.02, size=signal.shape))
    write_bank(hek, official, perts, signal + rng.normal(scale=0.02, size=signal.shape))
    write_bank(k562, official, perts, signal + rng.normal(scale=0.02, size=signal.shape))

    cd4_root = root / "cd4"
    for index, name in enumerate(RUN.CD4_NAMES):
        write_bank(
            cd4_root / name,
            official,
            perts,
            signal + rng.normal(scale=0.03 + index * 0.005, size=signal.shape),
            cd4=True,
        )

    panel_csv = root / "pert_counts.csv"
    genes_csv = root / "gene_names.csv"
    pd.DataFrame({"target_gene": panel}).to_csv(panel_csv, index=False)
    pd.DataFrame({"gene_name": official}).to_csv(genes_csv, index=False)
    return argparse.Namespace(
        hct116=hct,
        hek293t=hek,
        k562=k562,
        cd4_root=cd4_root,
        panel_csv=panel_csv,
        genes_csv=genes_csv,
        output=root / "results.json",
        whitener=root / "whitener_union.npz",
        scratch=root / "scratch",
        chunk_rows=2,
        max_rss_gib=0.95,
        expected_panel_targets=4,
        min_shared_targets=4,
        min_operator_genes=10,
        min_retrieval_genes=10,
        include_raw_diagnostic=True,
        instance_id="i-synthetic",
        instance_type="tiny",
        s3_prefix="s3://synthetic/exp060/",
    )


def test_end_to_end_protocol_and_leakage_audit() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        args = fixture(Path(directory))
        result = RUN.run(args)
        assert result["protocol"]["gallery_size"] == 4
        assert result["protocol"]["whitening_gamma"] == 0.65
        assert result["protocol"]["selection_on_k562"] is False
        assert result["leakage_audit"]["held_k562_in_fit_sources"] is False
        assert result["leakage_audit"]["truth_transformed"] is False
        assert result["leakage_audit"]["challenge_targets_in_operator_axis"] == 0
        assert result["leakage_audit"]["challenge_targets_in_retrieval_axis"] == 0
        assert result["axis_audit"]["operator_genes"] == 14
        assert result["axis_audit"]["retrieval_genes"] == 14
        assert result["whitener"]["metadata"]["donors"] == list(RUN.SOURCE_NAMES)
        assert len(result["candidates"]) == 2
        primary = next(item for item in result["candidates"] if item["is_primary"])
        assert primary["candidate"] == "cd4_mean3_f0.5_gene_norm"
        assert args.output.exists()
        assert args.output.with_name("summary.csv").exists()
        assert args.output.with_name("per_target_pds.csv").exists()
        assert not args.scratch.exists()
        assert RUN.peak_rss_gib() < 1.0


def test_target_coordinates_cannot_affect_retrieval() -> None:
    prediction = np.asarray([[1000.0, -1000.0, 1.0, 0.0], [-900.0, 900.0, 0.0, 1.0]])
    truth = np.asarray([[1.0, 0.0], [0.0, 1.0]])
    score = RUN.pds_per_target(prediction[:, 2:], truth)
    np.testing.assert_allclose(score, [1.0, 1.0])


if __name__ == "__main__":
    test_end_to_end_protocol_and_leakage_audit()
    test_target_coordinates_cannot_affect_retrieval()
    print(f"2 synthetic tests passed; peak RSS={RUN.peak_rss_gib():.3f} GiB")
