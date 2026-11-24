"""Small synthetic tests for the streaming Feng effect extractor."""

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
SPEC = importlib.util.spec_from_file_location(
    "extract_feng_effects", HERE / "extract_feng_effects.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_csv(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[column])
        writer.writeheader()
        writer.writerows({column: value} for value in values)


def make_source(path: Path) -> tuple[int, str]:
    rows = [
        ("A", "G1", "ENSG1", -0.5, 0.01),
        ("A", "G2", "ENSG2", 0.2, 0.5),
        ("A", "G4", "ENSG4", 0.1, 1.0),
        ("B", "G1", "ENSG1", -0.2, 0.2),
        ("B", "G2", "ENSG2", 0.4, 0.03),
        ("B", "G4", "ENSG4", -0.1, 0.9),
        ("C", "G1", "ENSG1", 0.3, 0.1),
        ("C", "G2", "ENSG2", -0.6, 0.001),
        ("C", "G4", "ENSG4", 0.0, 1.0),
    ]
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(value.decode() for value in MODULE.SOURCE_HEADER)
        writer.writerows(rows)
    payload = path.read_bytes()
    return len(payload), hashlib.md5(payload).hexdigest()


def test_extract_preserves_panel_and_gene_order() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        panel_csv = root / "panel.csv"
        genes_csv = root / "genes.csv"
        source = root / "source.tsv.gz"
        output = root / "out"
        write_csv(panel_csv, "target_gene", ["C", "MISSING", "A"])
        write_csv(genes_csv, "gene_name", ["G4", "G3", "G1", "G2"])
        size, md5 = make_source(source)

        metadata = MODULE.extract(
            source=str(source),
            panel_csv=panel_csv,
            genes_csv=genes_csv,
            output_dir=output,
            expected_size=size,
            expected_md5=md5,
            expected_rows=9,
            expected_targets=3,
            expected_genes=3,
            timeout=10,
            max_rss_gib=0.9,
            overwrite=False,
        )

        assert (output / "targets.txt").read_text().splitlines() == ["C", "A"]
        assert (output / "genes.txt").read_text().splitlines() == [
            "G4",
            "G3",
            "G1",
            "G2",
        ]
        lfc = np.load(output / "lfc.npy")
        centered = np.load(output / "centered_lfc.npy")
        source_mean = np.load(output / "source_gene_mean.npy")
        pval = np.load(output / "pval_adj.npy")
        observed = np.load(output / "observed.npy")
        assert lfc.shape == (2, 4)
        np.testing.assert_allclose(
            lfc,
            [[0.0, np.nan, 0.3, -0.6], [0.1, np.nan, -0.5, 0.2]],
            equal_nan=True,
        )
        np.testing.assert_allclose(
            pval,
            [[1.0, np.nan, 0.1, 0.001], [1.0, np.nan, 0.01, 0.5]],
            equal_nan=True,
        )
        np.testing.assert_allclose(
            source_mean,
            [0.0, np.nan, (-0.5 - 0.2 + 0.3) / 3, (0.2 + 0.4 - 0.6) / 3],
            equal_nan=True,
            atol=1e-7,
        )
        np.testing.assert_allclose(
            centered,
            [
                [0.0, np.nan, 0.3 - (-0.4 / 3), -0.6],
                [0.1, np.nan, -0.5 - (-0.4 / 3), 0.2],
            ],
            equal_nan=True,
            atol=1e-7,
        )
        assert observed[:, 1].sum() == 0
        assert metadata["selection"]["covered_targets"] == 2
        assert metadata["selection"]["missing_target_symbols"] == ["MISSING"]
        assert metadata["source_schema"] == {
            "rows": 9,
            "targets": 3,
            "features_in_file": 3,
            "biological_expressed_genes_reported_in_paper": 6_471,
            "synthetic_features": [],
            "rectangular": True,
            "header": [
                "Target",
                "Expressed_Gene_Symbol",
                "Expressed_Gene_Ens_ID",
                "lfc",
                "pval_adj",
            ],
        }
        on_disk = json.loads((output / "meta.json").read_text())
        assert on_disk["integrity_checks"]["md5"] is True
        assert MODULE.peak_rss_gib() < 0.9


def test_axis_mismatch_is_rejected() -> None:
    with tempfile.TemporaryDirectory(dir=HERE) as directory:
        root = Path(directory)
        panel_csv = root / "panel.csv"
        genes_csv = root / "genes.csv"
        source = root / "source.tsv.gz"
        output = root / "out"
        write_csv(panel_csv, "target_gene", ["A"])
        write_csv(genes_csv, "gene_name", ["G1", "G2"])
        with gzip.open(source, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(value.decode() for value in MODULE.SOURCE_HEADER)
            writer.writerows(
                [
                    ("A", "G1", "ENSG1", 0.0, 1.0),
                    ("A", "G2", "ENSG2", 0.0, 1.0),
                    ("B", "G2", "ENSG2", 0.0, 1.0),
                    ("B", "G1", "ENSG1", 0.0, 1.0),
                ]
            )

        try:
            MODULE.extract(
                source=str(source),
                panel_csv=panel_csv,
                genes_csv=genes_csv,
                output_dir=output,
                expected_size=None,
                expected_md5=None,
                expected_rows=None,
                expected_targets=None,
                expected_genes=None,
                timeout=10,
                max_rss_gib=0.9,
                overwrite=False,
            )
        except ValueError as exc:
            assert "source axis mismatch" in str(exc)
        else:
            raise AssertionError("malformed rectangular source was accepted")
