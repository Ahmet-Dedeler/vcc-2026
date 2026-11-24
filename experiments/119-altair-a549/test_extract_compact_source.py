from __future__ import annotations

import csv
import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SPEC = importlib.util.spec_from_file_location(
    "experiment_119_extract_compact_source",
    HERE / "extract_compact_source.py",
)
assert SPEC is not None and SPEC.loader is not None
EXTRACT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EXTRACT
SPEC.loader.exec_module(EXTRACT)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv_axis(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([column])
        writer.writerows([[value] for value in values])


def make_source(path: Path, duplicate_symbols: bool = False) -> None:
    source_symbols = (
        ["G2", "G2", "UNMAPPED", "G3"]
        if duplicate_symbols
        else ["G2", "G1", "UNMAPPED", "G3"]
    )
    x = np.asarray(
        [
            [2.0, np.nan, 4.0, 5.0],
            [np.nan, np.nan, np.nan, np.nan],
            [1.0, -1.0, 3.0, np.nan],
        ],
        dtype=np.float32,
    )
    fdr = np.asarray(
        [
            [3.0, np.nan, 2.0, np.nan],
            [np.nan, np.nan, np.nan, np.nan],
            [4.0, 1.0, np.nan, np.nan],
        ],
        dtype=np.float32,
    )
    obs = pd.DataFrame(
        {
            "role": pd.Categorical(["targeting", "control", "targeting"]),
            "gene_id": pd.Categorical(["ENSG_T2", "non-targeting", "ENSG_T1"]),
            "gene_symbol": ["T2", "Nontargeting", "T1"],
            "n_guides": [4, 2, 4],
            "n_cells": [20, 100, 30],
        },
        index=["T2", "Nontargeting", "T1"],
    )
    var = pd.DataFrame(
        {"gene_symbol": source_symbols},
        index=["ENSG2", "ENSG1", "ENSGU", "ENSG3"],
    )
    source = ad.AnnData(X=x, obs=obs, var=var)
    source.layers["neg_log10_fdr"] = fdr
    source.write_h5ad(path)


def make_contract(
    directory: Path,
    *,
    duplicate_symbols: bool = False,
) -> tuple[object, dict[str, object], dict[str, object]]:
    source = directory / "source.h5ad"
    manifest = directory / "SOURCE_MANIFEST.json"
    audit = directory / "SOURCE_AUDIT.json"
    receipt = directory / "SOURCE_ACQUISITION_RECEIPT.json"
    genes = directory / "genes.csv"
    panel = directory / "panel.csv"
    make_source(source, duplicate_symbols=duplicate_symbols)
    manifest.write_text('{"synthetic":"manifest"}\n', encoding="utf-8")
    audit.write_text('{"synthetic":"audit"}\n', encoding="utf-8")
    receipt.write_text('{"synthetic":"receipt"}\n', encoding="utf-8")
    write_csv_axis(genes, "gene_name", ["G1", "G2", "MISSING"])
    write_csv_axis(panel, "target_gene", ["T1", "MISSING_TARGET"])
    contract = EXTRACT.ExtractionContract(
        source_h5ad=source,
        source_manifest=manifest,
        source_audit=audit,
        acquisition_receipt=receipt,
        canonical_genes_csv=genes,
        panel_csv=panel,
        expected_source_size=source.stat().st_size,
        expected_source_sha256=sha256(source),
        expected_source_manifest_sha256=sha256(manifest),
        expected_canonical_genes_sha256=sha256(genes),
        expected_panel_sha256=sha256(panel),
        expected_source_rows=3,
        expected_target_rows=2,
        expected_source_genes=4,
        expected_canonical_genes=3,
        expected_panel_targets=2,
        expected_vcc_overlap=1,
        expected_vcc_supported_overlap=1,
    )
    full_axis = ["T2", "Nontargeting", "T1"]
    report: dict[str, object] = {
        "source": {
            "size_bytes": contract.expected_source_size,
            "sha256": contract.expected_source_sha256,
        },
        "summary": {
            "canonical_exact_overlap": 2,
        },
        "axis_bindings": {
            "source_perturbations_sha256": EXTRACT.sha256_lines(full_axis),
            "source_ensembl_ids_sha256": EXTRACT.sha256_lines(
                ["ENSG2", "ENSG1", "ENSGU", "ENSG3"]
            ),
            "source_gene_symbols_sha256": EXTRACT.sha256_lines(
                ["G2", "G1", "UNMAPPED", "G3"]
                if not duplicate_symbols
                else ["G2", "G2", "UNMAPPED", "G3"]
            ),
        },
    }
    acquisition: dict[str, object] = {
        "source_manifest_sha256": contract.expected_source_manifest_sha256,
    }
    return contract, report, acquisition


class ProductionContractTests(unittest.TestCase):
    def test_production_contract_is_source_only(self) -> None:
        contract = EXTRACT.production_contract()
        self.assertEqual(contract.expected_source_size, 281_474_102)
        self.assertEqual(contract.expected_target_rows, 1000)
        self.assertEqual(contract.expected_source_genes, 9999)
        self.assertEqual(contract.expected_canonical_genes, 18_533)
        self.assertEqual(contract.expected_vcc_overlap, 21)
        self.assertEqual(contract.expected_vcc_supported_overlap, 20)

    def test_output_scope_is_enforced(self) -> None:
        with self.assertRaisesRegex(EXTRACT.ExtractionError, "must be a child"):
            EXTRACT.safe_output(Path("/tmp/altair-forbidden-output"))


class SyntheticExtractionTests(unittest.TestCase):
    def test_exact_projection_preserves_nans_and_row_order(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE, prefix=".compact-test-") as name:
            root = Path(name)
            contract, report, receipt = make_contract(root)
            output = root / "compact"
            with mock.patch.object(
                EXTRACT.source_audit,
                "verify_frozen_reports",
                return_value=(report, receipt),
            ):
                EXTRACT.extract_compact_source(
                    contract,
                    output,
                    chunk_rows=1,
                )
            meta = EXTRACT.verify_compact_source(output, contract)

            self.assertEqual(
                (output / "perts.txt").read_text().splitlines(),
                ["T2", "T1"],
            )
            self.assertEqual(
                (output / "genes.txt").read_text().splitlines(),
                ["G1", "G2", "MISSING"],
            )
            lfc = np.load(output / "log2_fold_change.npy", allow_pickle=False)
            measured = np.load(output / "measured_mask.npy", allow_pickle=False)
            fdr = np.load(output / "neg_log10_fdr.npy", allow_pickle=False)
            fdr_mask = np.load(output / "fdr_mask.npy", allow_pickle=False)
            np.testing.assert_allclose(
                lfc,
                [[np.nan, 2.0, np.nan], [-1.0, 1.0, np.nan]],
                equal_nan=True,
            )
            np.testing.assert_array_equal(
                measured,
                [[False, True, False], [True, True, False]],
            )
            np.testing.assert_allclose(
                fdr,
                [[np.nan, 3.0, np.nan], [1.0, 4.0, np.nan]],
                equal_nan=True,
            )
            np.testing.assert_array_equal(
                fdr_mask,
                [[False, True, False], [True, True, False]],
            )
            self.assertEqual(meta["axes"]["source_exact_symbol_matches"], 2)
            self.assertEqual(meta["axes"]["vcc_target_rows"], 1)
            self.assertFalse(meta["mapping"]["alias_merge"])

    def test_duplicate_source_symbols_fail_instead_of_merging(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE, prefix=".compact-test-") as name:
            root = Path(name)
            contract, report, receipt = make_contract(
                root,
                duplicate_symbols=True,
            )
            output = root / "compact"
            with mock.patch.object(
                EXTRACT.source_audit,
                "verify_frozen_reports",
                return_value=(report, receipt),
            ):
                with self.assertRaisesRegex(EXTRACT.ExtractionError, "duplicates"):
                    EXTRACT.extract_compact_source(contract, output)
            self.assertFalse(output.exists())

    def test_existing_output_is_not_overwritten(self) -> None:
        with tempfile.TemporaryDirectory(dir=HERE, prefix=".compact-test-") as name:
            root = Path(name)
            contract, _, _ = make_contract(root)
            output = root / "compact"
            output.mkdir()
            with self.assertRaisesRegex(EXTRACT.ExtractionError, "overwrite"):
                EXTRACT.extract_compact_source(contract, output)


if __name__ == "__main__":
    unittest.main()
