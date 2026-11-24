#!/usr/bin/env python3
"""Range-extract VCC target rows from the official KOLF2.1J web matrix.

The source is commit-pinned and its metadata, builder, and browser loader are
hash-checked before any binary row is decoded. A response that does not honor
HTTP Range is rejected before its body is read.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import os
import resource
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

REPOSITORY = "https://github.com/Y-Doctor/KOLF2.1J_Perturbation_Cell_Atlas"
COMMIT = "8f84531a784f24b14d7e2c9b80c375d7e78fd46c"
RAW_ROOT = (
    "https://raw.githubusercontent.com/Y-Doctor/"
    f"KOLF2.1J_Perturbation_Cell_Atlas/{COMMIT}"
)
META_PATH = "docs/data/panel/all/meta.json"
BUILDER_PATH = "psp/notebooks/mde_webapp/webapp_build/build_panel_all_perts.py"
LOADER_PATH = "docs/js/main.js"
LICENSE_PATH = "LICENSE"

EXPECTED_SHA256 = {
    META_PATH: "ed4c569933cb62675db09e418290b5866632efcd5627ad8a629c6ca69da13cfe",
    BUILDER_PATH: "d5f37e34d20252012bb541b5ca25eec5f759ca137b28b65c4256f9eaee5110a0",
    LOADER_PATH: "55af64324b9b27065b2b4073df178582bcc7689100a7d50fdeff9d94b6ca0cc0",
    LICENSE_PATH: "55a4257f3c9354852a4534ea4bb94660160abcb870cb3fe37c0446ae56dd2629",
}

CHUNK_GIT_BLOBS = {
    0: "a605c0117ae665e68d146419771acf397d476819",
    1: "833fa5009cb2944f1174eea86f7a0d663de2e0e8",
    2: "cf366973b6eb2083ddd11c269f99baaf774e718f",
    3: "daa88999463b1c4b13fbc52460ca18749e5d7252",
    4: "2c0e26d899719179b883496a47a0b89d1a0fe313",
}

EXPECTED_META = {
    "n_perts": 10_167,
    "n_genes": 38_606,
    "chunk_size": 2_500,
    "n_chunks": 5,
    "row_bytes": 38_606,
    "scale": 0.04,
}


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024**3 if os.uname().sysname == "Darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float) -> None:
    if limit_gib > 0 and peak_rss_gib() > limit_gib:
        raise MemoryError(
            f"peak RSS {peak_rss_gib():.2f} GiB exceeded {limit_gib:.2f} GiB"
        )


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, block_bytes: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def fetch_small_asset(path: str, *, timeout: float) -> tuple[bytes, dict[str, str]]:
    url = f"{RAW_ROOT}/{path}"
    request = Request(
        url,
        headers={
            "Accept-Encoding": "identity",
            "User-Agent": "vcc2026-kolf-range-extractor/1",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
        headers = {key.lower(): value for key, value in response.headers.items()}
    expected = EXPECTED_SHA256[path]
    actual = sha256_bytes(payload)
    if actual != expected:
        raise ValueError(f"{path}: SHA-256 {actual} != pinned {expected}")
    return payload, headers


def validate_official_contract(
    meta_payload: bytes,
    builder_payload: bytes,
    loader_payload: bytes,
) -> dict[str, Any]:
    """Validate the schema and the exact producer/consumer byte contract."""
    meta = json.loads(meta_payload)
    required = {
        "perts",
        "genes",
        "scale",
        "chunk_size",
        "n_chunks",
        "n_perts",
        "n_genes",
        "row_bytes",
    }
    if set(meta) != required:
        raise ValueError(
            f"unexpected metadata keys: got {sorted(meta)}, expected {sorted(required)}"
        )
    for key, expected in EXPECTED_META.items():
        if meta[key] != expected:
            raise ValueError(f"metadata {key}={meta[key]!r}, expected {expected!r}")
    if len(meta["perts"]) != meta["n_perts"]:
        raise ValueError("metadata perturbation count does not match n_perts")
    if len(meta["genes"]) != meta["n_genes"]:
        raise ValueError("metadata gene count does not match n_genes")
    if len(set(meta["perts"])) != len(meta["perts"]):
        raise ValueError("official perturbation axis contains duplicates")
    if len(set(meta["genes"])) != len(meta["genes"]):
        raise ValueError("official gene axis contains duplicates")
    if meta["row_bytes"] != meta["n_genes"]:
        raise ValueError("one-byte int8 row contract is not satisfied")
    expected_chunks = (meta["n_perts"] + meta["chunk_size"] - 1) // meta["chunk_size"]
    if meta["n_chunks"] != expected_chunks:
        raise ValueError("n_chunks is inconsistent with n_perts and chunk_size")

    builder = builder_payload.decode("utf-8")
    loader = loader_payload.decode("utf-8")
    builder_markers = (
        "int8 row-major",
        "q = np.clip(np.round(M / SCALE), -127, 127).astype(np.int8)",
        "chunk.tobytes(order='C')",
        "'row_bytes':  int(n_genes)",
    )
    loader_markers = (
        "Rows are fetched on demand via HTTP Range",
        "const offset   = (pi % cs) * rb;",
        "Range: `bytes=${offset}-${offset + rb - 1}`",
        "row = new Int8Array(buf);",
        "row[gi] * scale",
    )
    missing_builder = [text for text in builder_markers if text not in builder]
    missing_loader = [text for text in loader_markers if text not in loader]
    if missing_builder or missing_loader:
        raise ValueError(
            "official source no longer proves the expected binary contract: "
            f"builder={missing_builder}, loader={missing_loader}"
        )
    return meta


def expected_chunk_size(meta: dict[str, Any], chunk_index: int) -> int:
    first = chunk_index * int(meta["chunk_size"])
    rows = min(int(meta["chunk_size"]), int(meta["n_perts"]) - first)
    if rows <= 0:
        raise IndexError(f"chunk {chunk_index} is outside the metadata range")
    return rows * int(meta["row_bytes"])


class StrictRangeClient:
    """Persistent HTTPS client that accepts only exact 206 responses."""

    def __init__(self, root_url: str, *, timeout: float, retries: int) -> None:
        parsed = urlsplit(root_url)
        if parsed.scheme != "https":
            raise ValueError("range source must use HTTPS")
        self.host = parsed.netloc
        self.root_path = parsed.path.rstrip("/")
        self.timeout = timeout
        self.retries = retries
        self.connection: http.client.HTTPSConnection | None = None

    def close(self) -> None:
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def _connect(self) -> http.client.HTTPSConnection:
        if self.connection is None:
            self.connection = http.client.HTTPSConnection(
                self.host, timeout=self.timeout
            )
        return self.connection

    def get(
        self,
        relative_path: str,
        start: int,
        end: int,
        total_size: int,
    ) -> tuple[bytes, dict[str, str]]:
        expected_length = end - start + 1
        expected_range = f"bytes {start}-{end}/{total_size}"
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                connection = self._connect()
                connection.request(
                    "GET",
                    f"{self.root_path}/{relative_path}",
                    headers={
                        "Range": f"bytes={start}-{end}",
                        "Accept-Encoding": "identity",
                        "User-Agent": "vcc2026-kolf-range-extractor/1",
                    },
                )
                response = connection.getresponse()
                headers = {key.lower(): value for key, value in response.getheaders()}
                if response.status != 206:
                    # Do not consume a possible full-chunk 200 response.
                    response.close()
                    self.close()
                    raise RuntimeError(
                        f"{relative_path}: server returned HTTP {response.status}; "
                        "strict range extraction refuses the response body"
                    )
                if headers.get("content-range") != expected_range:
                    response.close()
                    self.close()
                    raise RuntimeError(
                        f"{relative_path}: Content-Range "
                        f"{headers.get('content-range')!r} != {expected_range!r}"
                    )
                payload = response.read(expected_length + 1)
                response.close()
                if len(payload) != expected_length:
                    raise RuntimeError(
                        f"{relative_path}: got {len(payload)} bytes, "
                        f"expected {expected_length}"
                    )
                return payload, headers
            except (OSError, http.client.HTTPException, RuntimeError) as exc:
                last_error = exc
                self.close()
                if isinstance(exc, RuntimeError) and "HTTP 200" in str(exc):
                    raise
                if attempt == self.retries:
                    break
                time.sleep(min(2**attempt, 8))
        assert last_error is not None
        raise last_error


def read_axis(path: Path, column: str) -> list[str]:
    table = pd.read_csv(path)
    if column not in table:
        raise ValueError(f"{path}: missing required column {column!r}")
    values = table[column].dropna().astype(str).str.strip().tolist()
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: axis is empty or contains duplicate labels")
    return values


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, Any]:
    output = Path(args.output)
    if output.exists() and any(output.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{output} is non-empty; pass --overwrite")
    output.mkdir(parents=True, exist_ok=True)
    provenance_dir = output / "official"
    provenance_dir.mkdir(exist_ok=True)

    fetched_assets: dict[str, dict[str, Any]] = {}
    payloads: dict[str, bytes] = {}
    for path in (META_PATH, BUILDER_PATH, LOADER_PATH, LICENSE_PATH):
        payload, headers = fetch_small_asset(path, timeout=args.timeout_seconds)
        payloads[path] = payload
        local = provenance_dir / Path(path).name
        local.write_bytes(payload)
        fetched_assets[path] = {
            "url": f"{RAW_ROOT}/{path}",
            "sha256": sha256_bytes(payload),
            "bytes": len(payload),
            "etag": headers.get("etag"),
            "last_modified": headers.get("last-modified"),
            "local_path": str(local.relative_to(output)),
        }

    meta = validate_official_contract(
        payloads[META_PATH],
        payloads[BUILDER_PATH],
        payloads[LOADER_PATH],
    )
    panel = read_axis(Path(args.panel_csv), "target_gene")
    official_genes = read_axis(Path(args.genes_csv), "gene_name")
    source_pert_pos = {name: index for index, name in enumerate(meta["perts"])}
    source_gene_pos = {name: index for index, name in enumerate(meta["genes"])}
    vcc_gene_pos = {name: index for index, name in enumerate(official_genes)}
    selected_perts = [target for target in panel if target in source_pert_pos]
    missing_perts = [target for target in panel if target not in source_pert_pos]
    selected_genes = [gene for gene in official_genes if gene in source_gene_pos]
    missing_genes = [gene for gene in official_genes if gene not in source_gene_pos]
    if len(selected_perts) != args.expected_targets:
        raise ValueError(
            f"target overlap {len(selected_perts)} != expected {args.expected_targets}"
        )
    if len(selected_genes) < args.min_genes:
        raise ValueError(
            f"gene overlap {len(selected_genes)} < required {args.min_genes}"
        )

    row_bytes = int(meta["row_bytes"])
    total_range_bytes = len(selected_perts) * row_bytes
    if total_range_bytes > args.max_range_bytes:
        raise ValueError(
            f"planned range bytes {total_range_bytes:,} exceed cap "
            f"{args.max_range_bytes:,}"
        )

    source_rows_tmp = output / "source_rows_int8.npy.tmp"
    source_rows_path = output / "source_rows_int8.npy"
    mapped_tmp = output / "delta.npy.tmp"
    mapped_path = output / "delta.npy"
    for path in (source_rows_tmp, mapped_tmp):
        if path.exists():
            path.unlink()
    source_rows = np.lib.format.open_memmap(
        source_rows_tmp,
        mode="w+",
        dtype=np.int8,
        shape=(len(selected_perts), int(meta["n_genes"])),
    )
    source_base = f"{RAW_ROOT}/docs/data/panel/all"
    client = StrictRangeClient(
        source_base,
        timeout=args.timeout_seconds,
        retries=args.retries,
    )
    row_audit: list[dict[str, Any]] = []
    chunk_headers: dict[int, dict[str, Any]] = {}
    try:
        for output_row, target in enumerate(selected_perts):
            source_row = source_pert_pos[target]
            chunk = source_row // int(meta["chunk_size"])
            local_row = source_row % int(meta["chunk_size"])
            start = local_row * row_bytes
            end = start + row_bytes - 1
            chunk_size = expected_chunk_size(meta, chunk)
            payload, headers = client.get(
                f"chunk_{chunk}.bin",
                start,
                end,
                chunk_size,
            )
            source_rows[output_row] = np.frombuffer(payload, dtype=np.int8)
            chunk_headers.setdefault(
                chunk,
                {
                    "url": f"{source_base}/chunk_{chunk}.bin",
                    "expected_bytes": chunk_size,
                    "git_blob_oid": CHUNK_GIT_BLOBS[chunk],
                    "etag": headers.get("etag"),
                    "last_modified": headers.get("last-modified"),
                },
            )
            row_audit.append(
                {
                    "output_row": output_row,
                    "target_gene": target,
                    "source_row": source_row,
                    "chunk": chunk,
                    "byte_start": start,
                    "byte_end": end,
                    "content_range": headers.get("content-range"),
                    "etag": headers.get("etag"),
                    "sha256": sha256_bytes(payload),
                }
            )
            if (output_row + 1) % 25 == 0 or output_row + 1 == len(selected_perts):
                print(
                    f"fetched {output_row + 1}/{len(selected_perts)} rows; "
                    f"bytes={((output_row + 1) * row_bytes):,}; "
                    f"peak_rss={peak_rss_gib():.2f} GiB",
                    flush=True,
                )
            enforce_memory_limit(args.max_rss_gib)
    finally:
        client.close()
    source_rows.flush()
    del source_rows
    os.replace(source_rows_tmp, source_rows_path)

    gene_columns = np.asarray(
        [source_gene_pos[gene] for gene in selected_genes], dtype=np.int64
    )
    source_rows = np.load(source_rows_path, mmap_mode="r")
    mapped = np.lib.format.open_memmap(
        mapped_tmp,
        mode="w+",
        dtype=np.float32,
        shape=(len(selected_perts), len(selected_genes)),
    )
    scale = float(meta["scale"])
    saturated = 0
    for lo in range(0, len(selected_perts), args.map_chunk_rows):
        hi = min(lo + args.map_chunk_rows, len(selected_perts))
        block = np.asarray(source_rows[lo:hi, gene_columns], dtype=np.int8)
        saturated += int(np.count_nonzero(np.abs(block.astype(np.int16)) == 127))
        mapped[lo:hi] = block.astype(np.float32) * scale
        enforce_memory_limit(args.max_rss_gib)
    mapped.flush()
    del mapped, source_rows
    os.replace(mapped_tmp, mapped_path)

    np.save(output / "b_ctrl.npy", np.ones(len(selected_genes), dtype=np.float32))
    (output / "perts.txt").write_text("\n".join(selected_perts) + "\n")
    (output / "genes.txt").write_text("\n".join(selected_genes) + "\n")
    write_csv(output / "target_row_map.csv", row_audit)
    write_csv(
        output / "gene_map.csv",
        [
            {
                "output_column": index,
                "gene_name": gene,
                "vcc_gene_index": vcc_gene_pos[gene],
                "source_gene_index": source_gene_pos[gene],
            }
            for index, gene in enumerate(selected_genes)
        ],
    )

    contract = {
        "source_repository": REPOSITORY,
        "source_commit": COMMIT,
        "source_kind": "quantized per-perturbation z-score matrix",
        "effect_definition": (
            "delta = int8 source value * 0.04; official builder input "
            "zscore_filt_log1p.parquet"
        ),
        "binary_contract": {
            "dtype": "int8",
            "layout": "row-major C order",
            "row_bytes": row_bytes,
            "scale": scale,
            "clipping": "round(value / 0.04), clipped to [-127, 127]",
        },
        "vcc_target_coverage": len(selected_perts),
        "vcc_gene_coverage": len(selected_genes),
        "missing_vcc_targets": missing_perts,
        "missing_vcc_genes": missing_genes,
        "range_request_count": len(row_audit),
        "range_bytes_downloaded": total_range_bytes,
        "full_matrix_bytes_avoided": int(meta["n_perts"]) * row_bytes,
        "saturated_mapped_entries": saturated,
        "b_ctrl_definition": (
            "all-ones compatibility placeholder; no basal profile is present "
            "in the official web matrix and this experiment does not use it"
        ),
        "peak_rss_gib": peak_rss_gib(),
    }
    (output / "meta.json").write_text(json.dumps(contract, indent=2) + "\n")

    artifacts = {}
    for name in (
        "source_rows_int8.npy",
        "delta.npy",
        "b_ctrl.npy",
        "perts.txt",
        "genes.txt",
        "target_row_map.csv",
        "gene_map.csv",
        "meta.json",
    ):
        path = output / name
        artifacts[name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    provenance = {
        "source_repository": REPOSITORY,
        "source_commit": COMMIT,
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "official_assets": fetched_assets,
        "chunk_objects": chunk_headers,
        "artifacts": artifacts,
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(contract, indent=2), flush=True)
    return contract


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=str(ROOT / "data/external/kolf2"),
    )
    parser.add_argument(
        "--panel-csv",
        default=str(ROOT / "data/challenge/pert_counts.csv"),
    )
    parser.add_argument(
        "--genes-csv",
        default=str(ROOT / "data/challenge/gene_names.csv"),
    )
    parser.add_argument("--expected-targets", type=int, default=273)
    parser.add_argument("--min-genes", type=int, default=18_000)
    parser.add_argument("--max-range-bytes", type=int, default=12_000_000)
    parser.add_argument("--map-chunk-rows", type=int, default=32)
    parser.add_argument("--timeout-seconds", type=float, default=60)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--max-rss-gib", type=float, default=3.5)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.map_chunk_rows < 1:
        raise SystemExit("--map-chunk-rows must be positive")
    if not 0 < args.max_rss_gib <= 4:
        raise SystemExit("--max-rss-gib must be in (0, 4]")
    run(args)


if __name__ == "__main__":
    main()
