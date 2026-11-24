"""Synthetic tests for the fail-closed GSE295214 compact donor extractor."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SPEC = importlib.util.spec_from_file_location(
    "extract_compact_donor",
    HERE / "extract_compact_donor.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_csv_axis(path: Path, column: str, values: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([column])
        writer.writerows([[value] for value in values])


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_h5ad(
    path: Path,
    *,
    bad_mask: bool = False,
    varying_control: bool = False,
) -> None:
    source_ids = ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]
    source_symbols = ["G1", "DUP", "DUP", "UNMAPPED"]
    control_cpm = np.asarray([1000.0, 2000.0, 1000.0, 6000.0])
    control_rows = np.broadcast_to(control_cpm, (4, 4)).copy()
    if varying_control:
        control_rows[1] = [1001.0, 1999.0, 1000.0, 6000.0]
    perturbed_cpm = np.asarray(
        [
            [500.0, 2000.0, 1500.0, 6000.0],  # C
            [2000.0, 1000.0, 1000.0, 6000.0],  # A
            [1000.0, 1000.0, 2000.0, 6000.0],  # B
            [1000.0, 2000.0, 1000.0, 6000.0],  # other dataset
        ],
        dtype=np.float64,
    )
    control_log = np.log1p(control_rows).astype(np.float32)
    effect = (np.log1p(perturbed_cpm).astype(np.float32) - control_log).astype(
        np.float32
    )
    measure_mask = np.ones((4, 4), dtype=np.float32)
    measure_mask[0, 2] = 0.0  # Partial duplicate measurement => DUP unsupported.
    if bad_mask:
        measure_mask[1, 0] = 0.5
    obs = pd.DataFrame(
        {
            "dataset": [
                "Pan2025_hESC",
                "Pan2025_hESC",
                "Pan2025_hESC",
                "Other",
            ],
            "perturbed_gene_name": ["C", "A", "B", "OTHER"],
        },
        index=["row-c", "row-a", "row-b", "row-other"],
    )
    var = pd.DataFrame({"gene_name": source_symbols}, index=source_ids)
    source = ad.AnnData(X=effect, obs=obs, var=var)
    source.layers["control"] = control_log
    source.layers["measure_mask"] = measure_mask
    source.write_h5ad(path)


def make_contract(
    root: Path,
    *,
    bad_mask: bool = False,
    varying_control: bool = False,
) -> tuple[object, object]:
    source = root / "source.h5ad"
    manifest = root / "SOURCE_MANIFEST.json"
    amendment = root / "PSEUDOBULK_ACQUISITION_AMENDMENT.json"
    receipt = root / "pseudobulk-download-receipt.json"
    targets = root / "targets.txt"
    genes = root / "genes.csv"
    panel = root / "panel.csv"
    make_h5ad(
        source,
        bad_mask=bad_mask,
        varying_control=varying_control,
    )
    manifest.write_text('{"synthetic":"manifest"}\n', encoding="utf-8")
    authorization = {
        "asset_name": "public_pseudobulk_database",
        "download_public_pseudobulk_database": True,
        "download_geo_raw_archive": False,
        "confirmation_token": "synthetic",
        "max_download_bytes": source.stat().st_size,
    }
    amendment.write_text(
        json.dumps({"authorization": authorization}, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    targets.write_text("A\nB\nC\n", encoding="utf-8")
    write_csv_axis(genes, "gene_name", ["DUP", "G1", "MISSING"])
    write_csv_axis(panel, "target_gene", ["A", "C", "NOT_PRESENT"])
    source_hash = sha256(source)
    contract = MODULE.ExtractionContract(
        source_h5ad=source,
        download_receipt=receipt,
        source_manifest=manifest,
        amendment=amendment,
        target_axis=targets,
        canonical_genes_csv=genes,
        panel_csv=panel,
        expected_source_size=source.stat().st_size,
        expected_source_sha256=source_hash,
        expected_source_manifest_sha256=sha256(manifest),
        expected_amendment_sha256=sha256(amendment),
        expected_target_axis_sha256=sha256(targets),
        expected_canonical_genes_sha256=sha256(genes),
        expected_panel_sha256=sha256(panel),
        expected_rows=3,
        expected_vcc_targets=2,
        expected_canonical_genes=3,
        expected_panel_targets=3,
    )
    receipt.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "experiment": "115-gse295214-h1",
                "dataset": "GSE295214",
                "status": "verified_public_pseudobulk_download",
                "completed_at_utc": "2026-08-30T12:00:00Z",
                "source_manifest_sha256": contract.expected_source_manifest_sha256,
                "acquisition_amendment_sha256": contract.expected_amendment_sha256,
                "asset": {
                    "path": str(source.relative_to(MODULE.ROOT)),
                    "size_bytes": contract.expected_source_size,
                    "sha256": source_hash,
                    "repository_revision": ("1d3f28d966772dcf75870beac052923efdc078c5"),
                    "remote_lfs_sha256": source_hash,
                },
                "authorization": authorization,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    authenticated = MODULE.authenticate_download(contract)
    return contract, authenticated


def temporary_root() -> tempfile.TemporaryDirectory[str]:
    return tempfile.TemporaryDirectory(dir=HERE, prefix=".compact-test-")


class CompactDonorTests(unittest.TestCase):
    def test_production_contract_is_bound_to_frozen_inputs(self) -> None:
        contract = MODULE.production_contract()
        self.assertEqual(contract.expected_source_size, 4_425_260_319)
        self.assertEqual(contract.expected_rows, 2_354)
        self.assertEqual(contract.expected_vcc_targets, 64)
        self.assertEqual(contract.expected_canonical_genes, 18_533)
        self.assertEqual(contract.expected_panel_targets, 300)

    def test_authentication_rejects_receipt_hash_drift(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, _ = make_contract(root)
            receipt = json.loads(contract.download_receipt.read_text())
            receipt["asset"]["sha256"] = "0" * 64
            contract.download_receipt.write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(MODULE.ExtractionError, "asset binding"):
                MODULE.authenticate_download(contract)

    def test_output_scope_is_enforced(self) -> None:
        with self.assertRaisesRegex(MODULE.ExtractionError, "must be a child"):
            MODULE._safe_output(Path("/tmp/gse295214-forbidden-output"))

    def test_conversion_duplicate_collapse_masks_and_vcc_flags(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root)
            output = root / "compact"
            MODULE.extract_compact_donor(
                authenticated,
                contract,
                output,
                chunk_rows=2,
            )
            MODULE.verify_compact_donor(output)

            self.assertEqual(
                (output / "perts.txt").read_text().splitlines(),
                ["A", "B", "C"],
            )
            self.assertEqual(
                (output / "genes.txt").read_text().splitlines(),
                ["DUP", "G1", "MISSING"],
            )
            self.assertEqual(
                (output / "vcc_targets.txt").read_text().splitlines(),
                ["A", "C"],
            )
            np.testing.assert_array_equal(
                np.load(output / "vcc_target_mask.npy", allow_pickle=False),
                [True, False, True],
            )
            np.testing.assert_array_equal(
                np.load(
                    output / "source_covered_gene_mask.npy",
                    allow_pickle=False,
                ),
                [True, True, False],
            )

            b_ctrl = np.load(output / "b_ctrl.npy", allow_pickle=False)
            expected_control = np.asarray(
                [np.log1p(5.0 * 3000.0), np.log1p(5.0 * 1000.0), 0.0],
                dtype=np.float32,
            )
            np.testing.assert_allclose(b_ctrl, expected_control, rtol=1e-6)

            delta = np.load(output / "delta.npy", allow_pickle=False)
            measured = np.load(output / "measured_mask.npy", allow_pickle=False)
            expected_a = [
                np.log1p(5.0 * 2000.0) - np.log1p(5.0 * 3000.0),
                np.log1p(5.0 * 2000.0) - np.log1p(5.0 * 1000.0),
                0.0,
            ]
            np.testing.assert_allclose(
                delta[0],
                expected_a,
                rtol=1e-6,
                atol=1e-7,
            )
            np.testing.assert_array_equal(measured[0], [True, True, False])
            np.testing.assert_array_equal(measured[2], [False, True, False])
            self.assertEqual(delta[2, 0], 0.0)
            self.assertTrue(np.all(delta[~measured] == 0.0))

            rows = (output / "rows.tsv").read_text().splitlines()
            self.assertEqual(rows[1].split("\t")[2:], ["row-a", "A", "1"])
            self.assertEqual(rows[3].split("\t")[2:], ["row-c", "C", "1"])
            meta = json.loads((output / "meta.json").read_text())
            self.assertEqual(meta["axes"]["rows"], 3)
            self.assertEqual(meta["axes"]["source_duplicate_feature_rows"], 1)
            self.assertFalse(meta["safety"]["held_truth_opened"])
            self.assertEqual(
                meta["implementation_sha256"],
                MODULE.sha256_file(HERE / "extract_compact_donor.py"),
            )
            receipt = json.loads(
                (output / MODULE.COMPACT_RECEIPT_NAME).read_text(encoding="utf-8")
            )
            self.assertEqual(
                receipt["implementation"]["sha256"],
                MODULE.sha256_file(HERE / "extract_compact_donor.py"),
            )

    def test_outputs_are_byte_deterministic(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root)
            first = root / "first"
            second = root / "second"
            MODULE.extract_compact_donor(
                authenticated,
                contract,
                first,
                chunk_rows=1,
            )
            MODULE.extract_compact_donor(
                authenticated,
                contract,
                second,
                chunk_rows=3,
            )
            names = sorted(path.name for path in first.iterdir())
            self.assertEqual(
                names,
                sorted(path.name for path in second.iterdir()),
            )
            for name in names:
                self.assertEqual(
                    (first / name).read_bytes(),
                    (second / name).read_bytes(),
                    name,
                )

    def test_verifier_detects_artifact_corruption(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root)
            output = root / "compact"
            MODULE.extract_compact_donor(authenticated, contract, output)
            with (output / "genes.txt").open("ab") as handle:
                handle.write(b"CORRUPTION\n")
            with self.assertRaisesRegex(MODULE.ExtractionError, "binding changed"):
                MODULE.verify_compact_donor(output)

    def test_nonbinary_measure_mask_fails_closed(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root, bad_mask=True)
            with self.assertRaisesRegex(
                MODULE.ExtractionError,
                "not exactly binary",
            ):
                MODULE.extract_compact_donor(
                    authenticated,
                    contract,
                    root / "compact",
                )

    def test_missing_required_layer_fails_closed(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, _ = make_contract(root)
            with h5py.File(contract.source_h5ad, "r+") as handle:
                del handle["layers/measure_mask"]
            source_hash = sha256(contract.source_h5ad)
            receipt = json.loads(contract.download_receipt.read_text())
            receipt["asset"]["size_bytes"] = contract.source_h5ad.stat().st_size
            receipt["asset"]["sha256"] = source_hash
            receipt["asset"]["remote_lfs_sha256"] = source_hash
            contract.download_receipt.write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            modified = MODULE.ExtractionContract(
                **{
                    **contract.__dict__,
                    "expected_source_size": contract.source_h5ad.stat().st_size,
                    "expected_source_sha256": source_hash,
                }
            )
            authenticated = MODULE.authenticate_download(modified)
            with self.assertRaisesRegex(
                MODULE.ExtractionError,
                "layers/measure_mask",
            ):
                MODULE.extract_compact_donor(
                    authenticated,
                    modified,
                    root / "compact",
                )

    def test_nonidentical_h1_control_rows_fail_closed(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root, varying_control=True)
            with self.assertRaisesRegex(MODULE.ExtractionError, "not identical"):
                MODULE.extract_compact_donor(
                    authenticated,
                    contract,
                    root / "compact",
                )

    def test_frozen_target_axis_drift_fails_closed(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, authenticated = make_contract(root)
            contract.target_axis.write_text("A\nC\nD\n", encoding="utf-8")
            drifted = MODULE.ExtractionContract(
                **{
                    **contract.__dict__,
                    "expected_target_axis_sha256": sha256(contract.target_axis),
                }
            )
            with self.assertRaisesRegex(MODULE.ExtractionError, "differ"):
                MODULE.extract_compact_donor(
                    authenticated,
                    drifted,
                    root / "compact",
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
