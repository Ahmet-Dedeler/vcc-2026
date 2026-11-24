"""Synthetic and fail-closed tests for Experiment 120."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "myllia_builder",
    HERE / "build_compact_donor.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def temporary_root() -> tempfile.TemporaryDirectory[str]:
    return tempfile.TemporaryDirectory(dir=HERE, prefix=".test-")


def write_strings(group: h5py.Group, name: str, values: list[str]) -> None:
    group.create_dataset(
        name,
        data=np.asarray(values, dtype=h5py.string_dtype("utf-8")),
    )


def write_csr(group: h5py.Group, matrix: sparse.csr_matrix) -> None:
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = np.asarray(matrix.shape, dtype=np.int32)
    group.create_dataset("data", data=matrix.data.astype(np.float64))
    group.create_dataset("indices", data=matrix.indices.astype(np.int32))
    group.create_dataset("indptr", data=matrix.indptr.astype(np.int32))


def write_frame(
    handle: h5py.File,
    name: str,
    index: list[str],
    columns: dict[str, np.ndarray | list[str]],
) -> h5py.Group:
    group = handle.create_group(name)
    group.attrs["_index"] = "_index"
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.1.0"
    write_strings(group, "_index", index)
    for key, values in columns.items():
        if np.asarray(values).dtype.kind in {"O", "S", "U"}:
            write_strings(group, key, list(map(str, values)))
        else:
            group.create_dataset(key, data=np.asarray(values))
    return group


def make_h5ad(
    path: Path,
    channel: str,
    genes: list[str],
    guide_ids: list[str],
    symbols: list[str],
    counts: np.ndarray,
) -> dict[str, object]:
    matrix = sparse.csr_matrix(counts.astype(np.float64))
    obs_ids = [f"{channel}-cell-{index}" for index in range(len(guide_ids))]
    ncount = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    nfeature = np.diff(matrix.indptr).astype(np.int32)
    with h5py.File(path, "w") as handle:
        write_csr(handle.create_group("X"), matrix)
        write_frame(
            handle,
            "obs",
            obs_ids,
            {
                "channel": [channel] * len(guide_ids),
                "cl": ["A549"] * len(guide_ids),
                "nCount_RNA": ncount,
                "nFeature_RNA": nfeature,
                "sgrna_id": guide_ids,
                "sgrna_symbol": symbols,
            },
        )
        handle.create_group("obsm")
        raw = handle.create_group("raw")
        write_csr(raw.create_group("X"), matrix)
        write_frame(raw, "var", genes, {})
        write_frame(handle, "var", genes, {"features": genes})
        handle.create_group("varm")
    return {
        "cells": len(guide_ids),
        "gene_desert_cells": symbols.count(MODULE.DESERT_SYMBOL),
        "genes": len(genes),
        "gene_axis_sha256": MODULE.axis_sha256(genes),
        "ncount_sum": int(np.sum(ncount)),
        "nnz": int(matrix.nnz),
        "ntc_cells": symbols.count(MODULE.CONTROL_SYMBOL),
        "obs_axis_sha256": MODULE.axis_sha256(obs_ids),
        "observed_guides": len(set(guide_ids)),
    }


def make_library(path: Path) -> None:
    rows: list[list[object]] = []
    for target_index in range(218):
        symbol = f"T{target_index}"
        for rank in range(1, 4):
            rows.append(
                [
                    f"{symbol}_{rank}",
                    symbol,
                    rank,
                    "gene",
                    10_000 + target_index,
                    f"G{target_index:03d}{rank}",
                    100,
                ]
            )
    for rank in range(1, 21):
        rows.append(
            [
                f"gene_desert_{rank}",
                MODULE.DESERT_SYMBOL,
                rank,
                "ctrl_T",
                "",
                f"D{rank:02d}",
                100,
            ]
        )
    for rank in range(1, 21):
        rows.append(
            [
                f"non-targeting_{rank}",
                MODULE.CONTROL_SYMBOL,
                rank,
                "ctrl_NT",
                "",
                f"N{rank:02d}",
                100,
            ]
        )
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["id", "symbol", "rank", "class", "entrez", "guide", "count"])
        writer.writerows(rows)


def make_synthetic_contract(root: Path) -> tuple[object, dict[str, object]]:
    genes = ["A", "B", "C", "T0", "T1"]
    guide_ids = [
        "non-targeting_1",
        "gene_desert_1",
        "T0_1",
        "T0_2",
        "T1_1",
        "T1_2",
    ]
    symbols = [
        MODULE.CONTROL_SYMBOL,
        MODULE.DESERT_SYMBOL,
        "T0",
        "T0",
        "T1",
        "T1",
    ]
    ch1_counts = np.asarray(
        [
            [10, 20, 30, 40, 50],
            [9, 21, 30, 40, 50],
            [10, 30, 20, 40, 10],
            [10, 40, 10, 40, 20],
            [10, 25, 35, 40, 5],
            [10, 35, 25, 40, 5],
        ]
    )
    ch2_counts = np.asarray(
        [
            [20, 20, 20, 40, 50],
            [19, 21, 20, 40, 50],
            [20, 30, 10, 40, 10],
            [20, 40, 10, 30, 20],
            [20, 25, 25, 40, 5],
            [20, 35, 15, 40, 5],
        ]
    )
    channel_paths = {
        "ch_1": root / "ch1.h5ad",
        "ch_2": root / "ch2.h5ad",
    }
    expected_channels = {
        "ch_1": make_h5ad(
            channel_paths["ch_1"],
            "ch_1",
            genes,
            guide_ids,
            symbols,
            ch1_counts,
        ),
        "ch_2": make_h5ad(
            channel_paths["ch_2"],
            "ch_2",
            genes,
            guide_ids,
            symbols,
            ch2_counts,
        ),
    }
    library = root / "library.csv.gz"
    make_library(library)
    canonical = (
        "A",
        "B",
        "C",
        "T0",
        "T1",
        *(f"P{index}" for index in range(299)),
    )
    panel = ("A", "T1", *(f"P{index}" for index in range(298)))
    manifest = {
        "expected_audit": {
            "canonical_common_covered_genes": 5,
            "canonical_transferable_nonpanel_genes": 3,
            "guide_channel_rows_min_10_cells": 8,
            "paired_guides_min_10_each_channel": 4,
            "qc_passing_targets": 2,
            "qc_passing_vcc_targets": 1,
            "vcc_exact_overlap": 1,
        },
        "inputs": {
            "channel_1_h5ad": {"sha256": MODULE.sha256_file(channel_paths["ch_1"])},
            "channel_2_h5ad": {"sha256": MODULE.sha256_file(channel_paths["ch_2"])},
            "description_pdf": {"sha256": "0" * 64},
            "sgrna_library": {"sha256": MODULE.sha256_file(library)},
        },
    }
    contract = MODULE.Contract(
        manifest=manifest,
        provenance_sha256=MODULE.EXPECTED_PROVENANCE_SHA256,
        canonical_genes=canonical,
        panel=panel,
        channel_paths=channel_paths,
        library_path=library,
        pdf_path=root / "description.pdf",
        min_cells=1,
        min_paired_guides=2,
    )
    _, library_by_id = MODULE.load_library(library)
    metadata = {
        channel: MODULE.inspect_channel(
            channel,
            channel_paths[channel],
            expected_channels[channel],
            library_by_id,
        )
        for channel in MODULE.CHANNEL_ORDER
    }
    return contract, metadata


class CompactDonorTests(unittest.TestCase):
    def test_production_contract_is_hash_bound_and_source_only(self) -> None:
        contract = MODULE.load_contract()
        self.assertEqual(len(contract.canonical_genes), 18_533)
        self.assertEqual(len(contract.panel), 300)
        self.assertEqual(contract.min_cells, 10)
        self.assertFalse(contract.manifest["authorization"]["evaluation"])
        self.assertFalse(contract.manifest["authorization"]["open_vcc_held_truth"])
        self.assertFalse(
            contract.manifest["authorization"]["open_existing_outcome_directories"]
        )

    def test_output_scope_is_enforced(self) -> None:
        with self.assertRaisesRegex(MODULE.SourceError, "must be a child"):
            MODULE.safe_child(Path("/tmp/myllia-forbidden"))

    def test_exact_symbol_projection_rejects_duplicate_destinations(self) -> None:
        with self.assertRaisesRegex(MODULE.SourceError, "map multiply"):
            MODULE.canonical_projection(["A", "A"], {"A": 0})

    def test_robust_pooling_and_panel_exclusion(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, metadata = make_synthetic_contract(root)
            output = root / "compact"
            MODULE.build_compact_donor(
                contract,
                metadata,
                output,
                chunk_rows=2,
            )
            MODULE.verify_compact_donor(output, contract=contract)

            self.assertEqual(
                (output / "perts.txt").read_text().splitlines(),
                ["T0", "T1"],
            )
            panel_mask = np.load(
                output / "panel_coordinate_mask.npy",
                allow_pickle=False,
            )
            self.assertEqual(int(panel_mask.sum()), 300)
            for name in (
                "delta.npy",
                "guide_delta.npy",
                "guide_channel_delta.npy",
            ):
                values = np.load(output / name, allow_pickle=False)
                self.assertTrue(np.all(values[:, panel_mask] == 0.0))

            with (output / "guide_rows.tsv").open() as handle:
                guide_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (output / "guide_channel_rows.tsv").open() as handle:
                channel_rows = list(csv.DictReader(handle, delimiter="\t"))
            guide_delta = np.load(output / "guide_delta.npy", allow_pickle=False)
            channel_delta = np.load(
                output / "guide_channel_delta.npy",
                allow_pickle=False,
            )
            for guide_index, guide_row in enumerate(guide_rows):
                matching = [
                    index
                    for index, row in enumerate(channel_rows)
                    if row["guide_id"] == guide_row["guide_id"]
                ]
                self.assertEqual(len(matching), 2)
                np.testing.assert_allclose(
                    guide_delta[guide_index],
                    np.mean(channel_delta[matching], axis=0),
                    rtol=1e-6,
                    atol=1e-7,
                )

            target_delta = np.load(output / "delta.npy", allow_pickle=False)
            for target_index, target in enumerate(["T0", "T1"]):
                matching = [
                    index
                    for index, row in enumerate(guide_rows)
                    if row["target"] == target
                ]
                np.testing.assert_allclose(
                    target_delta[target_index],
                    np.median(guide_delta[matching], axis=0),
                    rtol=1e-6,
                    atol=1e-7,
                )

    def test_build_is_byte_deterministic_across_chunk_sizes(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, metadata = make_synthetic_contract(root)
            first = root / "first"
            second = root / "second"
            MODULE.build_compact_donor(contract, metadata, first, chunk_rows=1)
            MODULE.build_compact_donor(contract, metadata, second, chunk_rows=4)
            self.assertEqual(
                sorted(path.name for path in first.iterdir()),
                sorted(path.name for path in second.iterdir()),
            )
            for name in sorted(path.name for path in first.iterdir()):
                self.assertEqual(
                    (first / name).read_bytes(),
                    (second / name).read_bytes(),
                    name,
                )

    def test_verifier_detects_corruption(self) -> None:
        with temporary_root() as directory:
            root = Path(directory)
            contract, metadata = make_synthetic_contract(root)
            output = root / "compact"
            MODULE.build_compact_donor(contract, metadata, output, chunk_rows=3)
            with (output / "genes.txt").open("ab") as handle:
                handle.write(b"CORRUPTION\n")
            with self.assertRaisesRegex(MODULE.SourceError, "binding changed"):
                MODULE.verify_compact_donor(output, contract=contract)


if __name__ == "__main__":
    unittest.main(verbosity=2)
