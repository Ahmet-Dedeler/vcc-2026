#!/usr/bin/env python3
"""Stream the public PerturbMatch K562 matrix into a compact donor bank.

The source is a gzip CSV with genes as rows and perturbations as columns.  This
extractor never materializes the CSV or a dense copy in memory: it writes one
gene column at a time into a Fortran-ordered ``(perturbation, gene)`` NPY
memmap.  Hashes are computed over the compressed source bytes as they arrive.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import resource
import sys
import tempfile
import urllib.request
from collections.abc import Iterable
from pathlib import Path
from typing import BinaryIO

import numpy as np
import pandas as pd

SOURCE_URL = (
    "https://zenodo.org/api/records/21232999/files/lfc_matrix_ALL.csv.gz/content"
)
SOURCE_BYTES = 748_423_850
SOURCE_MD5 = "0017076a5cc7b9979febf3bd19bd2157"
ZENODO_RECORD = "https://zenodo.org/records/21232999"
ZENODO_DOI = "10.5281/zenodo.21232999"
ANALYSIS_COMMIT = "41e1040ddefa381133b293911c55bc9218a2e7dc"


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


class HashingReader:
    """Minimal file-like wrapper that hashes every compressed byte read."""

    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.md5 = hashlib.md5()
        self.sha256 = hashlib.sha256()
        self.bytes_read = 0

    def read(self, size: int = -1) -> bytes:
        block = self.raw.read(size)
        if block:
            self.md5.update(block)
            self.sha256.update(block)
            self.bytes_read += len(block)
        return block

    def readable(self) -> bool:
        return True


def strip_k5(value: str) -> str:
    if not value.endswith("_K5"):
        raise ValueError(f"unexpected PerturbMatch column label {value!r}")
    return value[:-3]


def write_lines(path: Path, values: Iterable[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values))


def read_panel(path: Path) -> tuple[str, ...]:
    frame = pd.read_csv(path)
    if "target_gene" not in frame:
        raise ValueError(f"{path}: missing target_gene column")
    values = tuple(frame["target_gene"].dropna().astype(str).str.strip())
    if len(values) != 300 or len(set(values)) != 300:
        raise ValueError(f"{path}: expected 300 unique panel targets")
    return values


def extract_stream(
    compressed: BinaryIO,
    output: Path,
    panel: tuple[str, ...],
    *,
    expected_genes: int,
    expected_bytes: int | None,
    expected_md5: str | None,
) -> dict[str, object]:
    """Extract one compressed source stream into ``output``."""
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing {output}")

    hashing = HashingReader(compressed)
    with tempfile.TemporaryDirectory(
        prefix=".perturbmatch-extract-", dir=output.parent
    ) as temporary_name:
        temporary = Path(temporary_name)
        with gzip.GzipFile(fileobj=hashing, mode="rb") as gz:
            text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
            reader = csv.reader(text)
            header = next(reader)
            if not header or header[0] != "gene":
                raise ValueError("source header must begin with 'gene'")
            source_columns = tuple(header[1:])
            perturbations = tuple(strip_k5(value) for value in source_columns)
            if len(perturbations) != len(set(perturbations)):
                raise ValueError("source contains duplicate perturbation labels")

            raw_path = temporary / "gene_major.f32"
            genes: list[str] = []
            with raw_path.open("wb") as raw:
                for gene_index, row in enumerate(reader):
                    if len(row) != len(header):
                        raise ValueError(
                            f"gene row {gene_index + 2} has {len(row)} fields; "
                            f"expected {len(header)}"
                        )
                    gene = row[0].strip()
                    if not gene:
                        raise ValueError(f"empty gene label on row {gene_index + 2}")
                    values = np.asarray(row[1:], dtype=np.float32)
                    if not np.isfinite(values).all():
                        raise ValueError(f"non-finite values for gene {gene}")
                    values.tofile(raw)
                    genes.append(gene)
                    if (gene_index + 1) % 250 == 0:
                        print(
                            f"genes={gene_index + 1} peak_rss={peak_rss_gib():.3f} GiB",
                            flush=True,
                        )
            text.detach()

        if expected_genes > 0 and len(genes) != expected_genes:
            raise ValueError(
                f"source has {len(genes)} genes; expected {expected_genes}"
            )
        if len(genes) != len(set(genes)):
            raise ValueError("source contains duplicate gene labels")
        raw = np.memmap(
            raw_path,
            mode="r",
            dtype=np.float32,
            shape=(len(genes), len(perturbations)),
        )
        delta = np.lib.format.open_memmap(
            temporary / "delta.npy",
            mode="w+",
            dtype=np.float32,
            shape=(len(perturbations), len(genes)),
            fortran_order=True,
        )
        for lo in range(0, len(genes), 64):
            hi = min(lo + 64, len(genes))
            delta[:, lo:hi] = raw[lo:hi].T
        delta.flush()
        del delta
        del raw
        raw_path.unlink()

        source_md5 = hashing.md5.hexdigest()
        source_sha256 = hashing.sha256.hexdigest()
        if expected_bytes is not None and hashing.bytes_read != expected_bytes:
            raise ValueError(
                f"compressed byte count {hashing.bytes_read} != {expected_bytes}"
            )
        if expected_md5 is not None and source_md5 != expected_md5:
            raise ValueError(f"source md5 {source_md5} != {expected_md5}")

        panel_set = set(panel)
        panel_covered = tuple(target for target in panel if target in perturbations)
        panel_missing = tuple(target for target in panel if target not in perturbations)
        write_lines(temporary / "genes.txt", genes)
        write_lines(temporary / "perts.txt", perturbations)
        write_lines(temporary / "panel_perts.txt", panel_covered)
        metadata: dict[str, object] = {
            "dataset": (
                "Replogle et al. (2022) K562 genome-wide CRISPRi "
                "Perturb-seq reanalysis by Yeung et al. (2026)"
            ),
            "source_url": SOURCE_URL,
            "zenodo_record": ZENODO_RECORD,
            "doi": ZENODO_DOI,
            "source_file": "lfc_matrix_ALL.csv.gz",
            "source_bytes": hashing.bytes_read,
            "source_md5": source_md5,
            "source_sha256": source_sha256,
            "license": "CC BY 4.0",
            "analysis_repository": (
                "https://github.com/Genentech/perturbmatch_analysis"
            ),
            "analysis_commit": ANALYSIS_COMMIT,
            "raw_data": "GEO GSE146194",
            "effect_definition": (
                "limma-voom log-fold-change from PerturbMatch matched controls; "
                "ALL combines singlets and guide multiplets at K=5"
            ),
            "matrix_layout": {
                "delta.npy": "float32 Fortran-order (perturbation, gene)",
                "source_csv": "gzip CSV (gene rows, perturbation_K5 columns)",
            },
            "n_perturbations": len(perturbations),
            "n_genes": len(genes),
            "vcc_target_coverage": len(panel_covered),
            "vcc_targets_covered": list(panel_covered),
            "vcc_targets_missing": list(panel_missing),
            "all_perturbations_are_unique": True,
            "all_genes_are_unique": True,
            "peak_rss_gib": peak_rss_gib(),
        }
        if len(panel_set & set(perturbations)) != len(panel_covered):
            raise AssertionError("panel coverage accounting failed")
        (temporary / "meta.json").write_text(json.dumps(metadata, indent=2) + "\n")
        os.replace(temporary, output)
        return metadata


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = read_panel(Path(args.panel_csv))
    request = urllib.request.Request(
        args.url,
        headers={"User-Agent": "vcc2026-perturbmatch-extractor/1.0"},
    )
    with urllib.request.urlopen(request, timeout=args.timeout) as response:
        metadata = extract_stream(
            response,
            Path(args.output),
            panel,
            expected_genes=args.expected_genes,
            expected_bytes=args.expected_bytes,
            expected_md5=args.expected_md5,
        )
    print(json.dumps(metadata, indent=2), flush=True)
    return metadata


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=SOURCE_URL)
    parser.add_argument(
        "--panel-csv", default=str(root / "data/challenge/pert_counts.csv")
    )
    parser.add_argument("--output", default=str(root / "data/external/perturbmatch"))
    parser.add_argument(
        "--expected-genes",
        type=int,
        default=9761,
        help="strict expected row count; use 0 to derive it from the stream",
    )
    parser.add_argument("--expected-bytes", type=int, default=SOURCE_BYTES)
    parser.add_argument("--expected-md5", default=SOURCE_MD5)
    parser.add_argument("--timeout", type=int, default=180)
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
