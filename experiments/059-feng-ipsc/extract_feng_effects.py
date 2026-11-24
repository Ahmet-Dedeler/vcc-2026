#!/usr/bin/env python3
"""Stream the official Feng genome-wide LFC table into a VCC-panel slice.

The source is long-form gzip with 43,187,656 target-by-feature rows. The
published paper reports 6,471 expressed genes and 43,180,983 biological
target-gene effects; the deposited file adds one synthetic selection marker,
``PURO``, for every target.
This extractor never materializes or stores the source table. It:

1. hashes the compressed byte stream while decompressing it;
2. validates the published rectangular target/gene layout;
3. retains only VCC panel targets and the canonical VCC gene axis;
4. writes small float32 arrays plus provenance metadata.

Peak memory is dominated by two 300 x 18,533 float32 work arrays (~45 MiB).
The published ``lfc`` values are preserved exactly apart from float32 rounding;
no log-base conversion or significance threshold is applied.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import io
import json
import math
import os
import resource
import shutil
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any, BinaryIO
from urllib.parse import urlsplit, urlunsplit
from urllib.request import Request, urlopen

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

SOURCE_URL = "https://ndownloader.figshare.com/files/49291654"
SOURCE_ARTICLE = "https://doi.org/10.6084/m9.figshare.26819743"
SOURCE_FILE_ID = 49_291_654
SOURCE_FILE_NAME = "GenomeWideScreen_LFC_byGene.tsv.gz"
SOURCE_SIZE = 794_610_300
SOURCE_MD5 = "03cb8a3b4328de50d45746a9d784bbc5"
SOURCE_ROWS = 43_187_656
SOURCE_TARGETS = 6_673
SOURCE_GENES = 6_472
SOURCE_HEADER = (
    b"Target",
    b"Expressed_Gene_Symbol",
    b"Expressed_Gene_Ens_ID",
    b"lfc",
    b"pval_adj",
)


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    import csv

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing required column {column!r}")
        values = tuple(value for row in reader if (value := str(row[column]).strip()))
    if not values:
        raise ValueError(f"{path}: empty {column!r} axis")
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


class HashingReader(io.RawIOBase):
    """Count and hash compressed bytes while gzip consumes a stream."""

    def __init__(self, raw: BinaryIO) -> None:
        self.raw = raw
        self.md5 = hashlib.md5()
        self.bytes_read = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray | memoryview) -> int:
        payload = self.raw.read(len(buffer))
        if not payload:
            return 0
        size = len(payload)
        buffer[:size] = payload
        self.md5.update(payload)
        self.bytes_read += size
        return size


@contextlib.contextmanager
def open_source(
    source: str, timeout: float
) -> Iterator[tuple[BinaryIO, dict[str, Any]]]:
    """Open a local gzip file or a streamed HTTPS response."""
    local = Path(source).expanduser()
    if local.exists():
        handle = local.open("rb")
        try:
            yield (
                handle,
                {
                    "kind": "local",
                    "resolved_source": str(local.resolve()),
                    "content_length": local.stat().st_size,
                },
            )
        finally:
            handle.close()
        return

    if not source.startswith("https://"):
        raise ValueError("source must be an existing local file or an HTTPS URL")
    request = Request(
        source,
        headers={
            "Accept-Encoding": "identity",
            "User-Agent": "vcc2026-feng-stream-extractor/1",
        },
    )
    response = urlopen(request, timeout=timeout)
    try:
        length = response.headers.get("Content-Length")
        resolved = urlsplit(response.geturl())
        yield (
            response,
            {
                "kind": "https",
                "resolved_source": urlunsplit(
                    (resolved.scheme, resolved.netloc, resolved.path, "", "")
                ),
                "content_length": int(length) if length else None,
                "etag": response.headers.get("ETag"),
                "last_modified": response.headers.get("Last-Modified"),
            },
        )
    finally:
        response.close()


def parse_float(payload: bytes, *, row_number: int, field: str) -> float:
    try:
        value = float(payload)
    except ValueError as exc:
        raise ValueError(
            f"source row {row_number}: invalid {field} {payload!r}"
        ) from exc
    if not math.isfinite(value):
        raise ValueError(f"source row {row_number}: non-finite {field}")
    return value


def write_lines(path: Path, values: tuple[str, ...] | list[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values))


def extract(
    *,
    source: str,
    panel_csv: Path,
    genes_csv: Path,
    output_dir: Path,
    expected_size: int | None,
    expected_md5: str | None,
    expected_rows: int | None,
    expected_targets: int | None,
    expected_genes: int | None,
    timeout: float,
    max_rss_gib: float,
    overwrite: bool,
) -> dict[str, Any]:
    panel = read_csv_axis(panel_csv, "target_gene")
    genes = read_csv_axis(genes_csv, "gene_name")
    panel_pos = {target.encode(): index for index, target in enumerate(panel)}
    gene_pos = {gene.encode(): index for index, gene in enumerate(genes)}

    if output_dir.exists():
        if not overwrite:
            raise FileExistsError(f"{output_dir} exists; pass --overwrite")
        if output_dir.is_symlink() or not output_dir.is_dir():
            raise ValueError(f"refusing to replace non-directory {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    lfc = np.full((len(panel), len(genes)), np.nan, dtype=np.float32)
    pval = np.full((len(panel), len(genes)), np.nan, dtype=np.float32)
    observed = np.zeros((len(panel), len(genes)), dtype=np.bool_)
    target_seen = np.zeros(len(panel), dtype=np.bool_)

    source_gene_axis: list[bytes] = []
    source_ensembl_axis: list[bytes] = []
    seen_source_targets: set[bytes] = set()
    current_target: bytes | None = None
    current_position = 0
    row_count = 0
    source_lfc_sum: list[float] | np.ndarray = []
    duplicate_selected_rows = 0
    selected_source_rows = 0
    response_metadata: dict[str, Any]

    with open_source(source, timeout) as (raw_source, response_metadata):
        hashing = HashingReader(raw_source)
        buffered = io.BufferedReader(hashing, buffer_size=1 << 20)
        with gzip.GzipFile(fileobj=buffered, mode="rb") as decompressed:
            header_line = decompressed.readline()
            if not header_line:
                raise ValueError("empty gzip source")
            header = tuple(header_line.rstrip(b"\r\n").split(b"\t"))
            if header != SOURCE_HEADER:
                raise ValueError(f"unexpected source header: {header!r}")

            for row_number, line in enumerate(decompressed, start=2):
                row_count += 1
                parts = line.rstrip(b"\r\n").split(b"\t")
                if len(parts) != 5:
                    raise ValueError(
                        f"source row {row_number}: expected 5 fields, got {len(parts)}"
                    )
                target, gene, ensembl, lfc_raw, pval_raw = parts
                lfc_value = parse_float(lfc_raw, row_number=row_number, field="lfc")

                if target != current_target:
                    if current_target is not None:
                        if current_position != len(source_gene_axis):
                            raise ValueError(
                                f"target {current_target!r} has {current_position} "
                                f"genes; expected {len(source_gene_axis)}"
                            )
                        if target in seen_source_targets:
                            raise ValueError(
                                f"source target {target!r} is not contiguous"
                            )
                    seen_source_targets.add(target)
                    current_target = target
                    current_position = 0
                    if len(seen_source_targets) == 2:
                        source_lfc_sum = np.asarray(source_lfc_sum, dtype=np.float64)

                if len(seen_source_targets) == 1:
                    source_gene_axis.append(gene)
                    source_ensembl_axis.append(ensembl)
                    assert isinstance(source_lfc_sum, list)
                    source_lfc_sum.append(lfc_value)
                else:
                    assert isinstance(source_lfc_sum, np.ndarray)
                    if current_position >= len(source_gene_axis):
                        raise ValueError(
                            f"target {target!r} has more genes than first target"
                        )
                    if (
                        gene != source_gene_axis[current_position]
                        or ensembl != source_ensembl_axis[current_position]
                    ):
                        raise ValueError(
                            f"source axis mismatch at target {target!r}, "
                            f"position {current_position}"
                        )
                    source_lfc_sum[current_position] += lfc_value

                target_index = panel_pos.get(target)
                if target_index is not None:
                    target_seen[target_index] = True
                    destination = gene_pos.get(gene)
                    if destination is not None:
                        if observed[target_index, destination]:
                            duplicate_selected_rows += 1
                        else:
                            pval_value = parse_float(
                                pval_raw, row_number=row_number, field="pval_adj"
                            )
                            if not 0.0 <= pval_value <= 1.0:
                                raise ValueError(
                                    f"source row {row_number}: pval_adj outside [0,1]"
                                )
                            lfc[target_index, destination] = lfc_value
                            pval[target_index, destination] = pval_value
                            observed[target_index, destination] = True
                            selected_source_rows += 1

                current_position += 1
                if row_count % 1_000_000 == 0:
                    enforce_memory_limit(
                        max_rss_gib, f"after {row_count:,} source rows"
                    )

        # gzip may buffer compressed bytes past the logical end. Drain the
        # compressed stream so size/hash cover the complete source object.
        while buffered.read(1 << 20):
            pass

    if current_target is None:
        raise ValueError("source contains no data rows")
    if current_position != len(source_gene_axis):
        raise ValueError(
            f"final target {current_target!r} has {current_position} genes; "
            f"expected {len(source_gene_axis)}"
        )

    source_md5 = hashing.md5.hexdigest()
    source_bytes = hashing.bytes_read
    integrity_checks = {
        "size": expected_size is None or source_bytes == expected_size,
        "md5": expected_md5 is None or source_md5 == expected_md5.lower(),
        "rows": expected_rows is None or row_count == expected_rows,
        "targets": (
            expected_targets is None or len(seen_source_targets) == expected_targets
        ),
        "genes": expected_genes is None or len(source_gene_axis) == expected_genes,
        "duplicate_selected_rows": duplicate_selected_rows == 0,
    }
    failed = [name for name, passed in integrity_checks.items() if not passed]
    if failed:
        raise ValueError(
            "source integrity/schema checks failed: "
            + ", ".join(failed)
            + f"; bytes={source_bytes}, md5={source_md5}, rows={row_count}, "
            + f"targets={len(seen_source_targets)}, genes={len(source_gene_axis)}"
        )

    covered_indices = np.flatnonzero(target_seen)
    covered_targets = [panel[index] for index in covered_indices]
    compact_lfc = lfc[covered_indices]
    compact_pval = pval[covered_indices]
    compact_observed = observed[covered_indices]
    measured_gene_mask = compact_observed.any(axis=0)
    if compact_observed.size:
        row_observed_counts = compact_observed.sum(axis=1)
        if not np.all(row_observed_counts == row_observed_counts[0]):
            raise ValueError("selected targets do not share one measured-gene axis")

    source_symbols = tuple(value.decode() for value in source_gene_axis)
    source_ensembl = tuple(value.decode() for value in source_ensembl_axis)
    if isinstance(source_lfc_sum, list):
        source_lfc_sum = np.asarray(source_lfc_sum, dtype=np.float64)
    source_lfc_mean = source_lfc_sum / len(seen_source_targets)
    source_to_vcc = np.asarray(
        [gene_pos.get(value, -1) for value in source_gene_axis], dtype=np.int32
    )
    vcc_source_mean = np.full(len(genes), np.nan, dtype=np.float32)
    mapped_source = source_to_vcc >= 0
    vcc_source_mean[source_to_vcc[mapped_source]] = source_lfc_mean[
        mapped_source
    ].astype(np.float32)
    centered_lfc = compact_lfc - vcc_source_mean[None, :]

    with tempfile.TemporaryDirectory(
        dir=output_dir.parent, prefix=f".{output_dir.name}."
    ) as temporary:
        work = Path(temporary)
        np.save(work / "lfc.npy", compact_lfc)
        np.save(work / "centered_lfc.npy", centered_lfc)
        np.save(work / "pval_adj.npy", compact_pval)
        np.save(work / "observed.npy", compact_observed)
        np.save(work / "measured_gene_mask.npy", measured_gene_mask)
        np.save(work / "source_gene_mean.npy", vcc_source_mean)
        np.save(work / "source_to_vcc.npy", source_to_vcc)
        write_lines(work / "targets.txt", covered_targets)
        write_lines(work / "genes.txt", list(genes))
        with (work / "source_gene_axis.tsv").open("w") as handle:
            handle.write("source_index\tgene_symbol\tensembl_id\tvcc_index\n")
            for index, (symbol, ensembl, destination) in enumerate(
                zip(source_symbols, source_ensembl, source_to_vcc, strict=True)
            ):
                handle.write(f"{index}\t{symbol}\t{ensembl}\t{destination}\n")

        self_lfc: dict[str, float] = {}
        self_pval: dict[str, float] = {}
        for row, target in enumerate(covered_targets):
            destination = gene_pos.get(target.encode())
            if destination is None or not compact_observed[row, destination]:
                continue
            self_lfc[target] = float(compact_lfc[row, destination])
            self_pval[target] = float(compact_pval[row, destination])

        files = {}
        for path in sorted(work.iterdir()):
            if path.name == "meta.json":
                continue
            files[path.name] = {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        metadata: dict[str, Any] = {
            "format_version": 1,
            "created_by": "experiments/059-feng-ipsc/extract_feng_effects.py",
            "source": {
                "article_doi": SOURCE_ARTICLE,
                "file_id": SOURCE_FILE_ID,
                "file_name": SOURCE_FILE_NAME,
                "requested_source": source,
                "resolved_source": response_metadata,
                "compressed_bytes": source_bytes,
                "compressed_md5": source_md5,
                "license": "MIT",
                "published_value_semantics": (
                    "Authors' pooled genome-wide target-level lfc and "
                    "Benjamini-Hochberg-adjusted p-value. Values are preserved "
                    "without log-base conversion or significance filtering."
                ),
                "centered_value_semantics": (
                    "centered_lfc = published lfc minus the per-expressed-gene "
                    "mean across all 6,673 source targets"
                ),
            },
            "source_schema": {
                "rows": row_count,
                "targets": len(seen_source_targets),
                "features_in_file": len(source_gene_axis),
                "biological_expressed_genes_reported_in_paper": 6_471,
                "synthetic_features": [
                    symbol
                    for symbol, ensembl in zip(
                        source_symbols, source_ensembl, strict=True
                    )
                    if symbol == ensembl and not ensembl.startswith("ENSG")
                ],
                "rectangular": True,
                "header": [value.decode() for value in SOURCE_HEADER],
            },
            "selection": {
                "vcc_panel_targets": len(panel),
                "covered_targets": len(covered_targets),
                "covered_target_symbols": covered_targets,
                "missing_target_symbols": [
                    target
                    for index, target in enumerate(panel)
                    if not target_seen[index]
                ],
                "vcc_genes": len(genes),
                "measured_vcc_genes": int(measured_gene_mask.sum()),
                "selected_effect_rows": selected_source_rows,
                "observed_values_per_target": (
                    int(compact_observed[0].sum()) if len(covered_targets) else 0
                ),
            },
            "quality": {
                "self_lfc_available": len(self_lfc),
                "self_lfc_median": (
                    float(np.median(list(self_lfc.values()))) if self_lfc else None
                ),
                "self_lfc_negative_fraction": (
                    float(np.mean(np.asarray(list(self_lfc.values())) < 0))
                    if self_lfc
                    else None
                ),
                "self_pval_below_0_1_fraction": (
                    float(np.mean(np.asarray(list(self_pval.values())) < 0.1))
                    if self_pval
                    else None
                ),
            },
            "integrity_checks": integrity_checks,
            "peak_rss_gib": peak_rss_gib(),
            "files": files,
        }
        (work / "meta.json").write_text(json.dumps(metadata, indent=2) + "\n")

        if output_dir.exists():
            shutil.rmtree(output_dir)
        os.replace(work, output_dir)

    metadata["files"]["meta.json"] = {
        "bytes": (output_dir / "meta.json").stat().st_size,
        "sha256": sha256_file(output_dir / "meta.json"),
    }
    enforce_memory_limit(max_rss_gib, "completed extraction")
    return metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default=SOURCE_URL)
    parser.add_argument(
        "--panel-csv", type=Path, default=ROOT / "data/challenge/pert_counts.csv"
    )
    parser.add_argument(
        "--genes-csv", type=Path, default=ROOT / "data/challenge/gene_names.csv"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/external/feng-ipsc/genomewide_vcc_panel",
    )
    parser.add_argument("--expected-size", type=int, default=SOURCE_SIZE)
    parser.add_argument("--expected-md5", default=SOURCE_MD5)
    parser.add_argument("--expected-rows", type=int, default=SOURCE_ROWS)
    parser.add_argument("--expected-targets", type=int, default=SOURCE_TARGETS)
    parser.add_argument("--expected-genes", type=int, default=SOURCE_GENES)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--max-rss-gib", type=float, default=0.90)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    metadata = extract(
        source=args.source,
        panel_csv=args.panel_csv,
        genes_csv=args.genes_csv,
        output_dir=args.output_dir,
        expected_size=args.expected_size,
        expected_md5=args.expected_md5,
        expected_rows=args.expected_rows,
        expected_targets=args.expected_targets,
        expected_genes=args.expected_genes,
        timeout=args.timeout,
        max_rss_gib=args.max_rss_gib,
        overwrite=args.overwrite,
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
