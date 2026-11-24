#!/usr/bin/env python3
"""Leakage-safe context-control variance transport under fixed union whitening.

The held donor contributes control-cell variability at inference, exactly as a
challenge context may. Its perturbation effects are not touched until the
source-only candidate has been selected.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import re
import resource
import struct
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

MAIN_NAMES = ("k562_gwps", "orion_hct116", "orion_hek293t")
BANK_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
K562_H5AD = (
    ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
)
HCT116_H5AD = ROOT / "data/harness/orion-hct116.h5ad"
ORION_DATA = ROOT / "data/external/orion/data"
ORION_GENE_META = ROOT / "data/external/orion/metadata/gene_metadata.parquet"
WHITENER_DIR = (
    ROOT
    / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
)
CACHE_NPZ = HERE / "control_stats.npz"
CACHE_JSON = HERE / "control_stats.json"
ORION_WORKER_DIR = HERE / "orion-control-workers"
RESULTS_JSON = HERE / "results.json"
PER_TARGET_CSV = HERE / "per_target_pds.csv"

TARGET_SUM = 50_000.0
WHITEN_GAMMA = 0.65
EFFECT_SCALE = 6.817
ETAS = (0.25, 0.5, 1.0)
METHODS = ("sd", "mad")
PLACEMENTS = ("pre", "post")
RATIO_CLIP = (0.5, 2.0)
CONTROL_REPLICATES = ("a", "b")
CONTROL_PATTERNS = (
    r"^non.?targeting",
    r"^ntc",
    r"^control",
    r"^safe.?harbou?r",
    r"^neg",
)


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached "
            f"the internal {limit_gib:.3f} GiB stop"
        )


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(x.strip() for x in path.read_text().splitlines() if x.strip())
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column!r}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate {column!r} values")
    return values


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def decode_strings(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind == "S":
        return values.astype("U")
    return values.astype(str)


def categorical_code(
    handle: Any,
    path: str,
    label: str,
) -> tuple[np.ndarray, int]:
    group = handle[path]
    categories = decode_strings(np.asarray(group["categories"]))
    matches = np.flatnonzero(categories == label)
    if len(matches) != 1:
        raise ValueError(f"{path}: expected one category {label!r}, got {categories}")
    return np.asarray(group["codes"]), int(matches[0])


def normalize_dense_rows(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float32)
    totals = counts.sum(axis=1, dtype=np.float64)
    if np.any(totals <= 0):
        raise ValueError("control sample contains a zero-library cell")
    counts *= (1_000_000.0 / totals).astype(np.float32)[:, None]
    np.log1p(counts, out=counts)
    return counts


def robust_scales(sample: np.ndarray, gene_block: int = 256) -> dict[str, np.ndarray]:
    """Return winsorized SD and scaled MAD of log1p(CPM) controls."""
    n_genes = sample.shape[1]
    sd = np.zeros(n_genes, dtype=np.float32)
    mad = np.zeros(n_genes, dtype=np.float32)
    for lo in range(0, n_genes, gene_block):
        hi = min(n_genes, lo + gene_block)
        block = np.asarray(sample[:, lo:hi], dtype=np.float32)
        quantiles = np.quantile(block, (0.01, 0.99), axis=0)
        clipped = np.clip(block, quantiles[0], quantiles[1])
        sd[lo:hi] = np.std(clipped, axis=0, ddof=1).astype(np.float32)
        median = np.median(block, axis=0)
        mad[lo:hi] = (
            1.4826 * np.median(np.abs(block - median), axis=0)
        ).astype(np.float32)
    return {"sd": sd, "mad": mad}


def select_dense_h5_strata(
    control_rows: np.ndarray,
    row_chunk: int,
    sample_per_rep: int,
    n_strata: int,
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    if sample_per_rep % n_strata:
        raise ValueError("sample_per_rep must be divisible by n_strata")
    quota = sample_per_rep // n_strata
    chunk_ids = np.unique(control_rows // row_chunk)
    selected_chunks = np.unique(
        chunk_ids[
            np.rint(np.linspace(0, len(chunk_ids) - 1, n_strata)).astype(int)
        ]
    )
    if len(selected_chunks) != n_strata:
        raise ValueError("could not select the requested number of dense strata")
    plan: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for chunk_id in selected_chunks:
        rows = control_rows[control_rows // row_chunk == chunk_id]
        needed = 2 * quota
        if len(rows) < needed:
            raise ValueError(
                f"dense chunk {chunk_id} has {len(rows)} controls, need {needed}"
            )
        positions = np.rint(np.linspace(0, len(rows) - 1, needed)).astype(int)
        picked = rows[positions]
        if len(np.unique(picked)) != needed:
            raise ValueError(f"dense chunk {chunk_id}: duplicate sampled rows")
        plan[int(chunk_id)] = (picked[0::2], picked[1::2])
    return plan


def dense_h5_control_samples(
    path: Path,
    label_col: str,
    control_label: str,
    n_genes: int,
    *,
    sample_per_rep: int,
    n_strata: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    import h5py

    samples = [
        np.empty((sample_per_rep, n_genes), dtype=np.float32),
        np.empty((sample_per_rep, n_genes), dtype=np.float32),
    ]
    offsets = [0, 0]
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        if not isinstance(x, h5py.Dataset):
            raise TypeError(f"{path}: expected dense X dataset")
        codes, control_code = categorical_code(
            handle, f"obs/{label_col}", control_label
        )
        control_rows = np.flatnonzero(codes == control_code)
        row_chunk = int(x.chunks[0] if x.chunks else 1024)
        plan = select_dense_h5_strata(
            control_rows, row_chunk, sample_per_rep, n_strata
        )
        for done, (chunk_id, rep_rows) in enumerate(sorted(plan.items()), start=1):
            lo = chunk_id * row_chunk
            hi = min(x.shape[0], lo + row_chunk)
            block = np.asarray(x[lo:hi, :], dtype=np.float32)
            for rep_index, rows in enumerate(rep_rows):
                selected = np.asarray(block[rows - lo], dtype=np.float32)
                selected = normalize_dense_rows(selected)
                end = offsets[rep_index] + len(rows)
                samples[rep_index][offsets[rep_index] : end] = selected
                offsets[rep_index] = end
            del block
            gc.collect()
            enforce_memory_limit(max_rss_gib, f"K562 control stratum {done}")
            print(
                f"  K562 control stratum {done}/{len(plan)} "
                f"samples={offsets}",
                flush=True,
            )
    if offsets != [sample_per_rep, sample_per_rep]:
        raise RuntimeError(f"K562 sample counts are {offsets}")
    return samples[0], samples[1], {
        "source": str(path.resolve()),
        "available_controls": int(len(control_rows)),
        "sample_per_replicate": sample_per_rep,
        "strata": n_strata,
        "hdf5_row_chunk": row_chunk,
    }


def read_csr_row(group: Any, row: int) -> tuple[np.ndarray, np.ndarray]:
    indptr = group["indptr"]
    start, stop = int(indptr[row]), int(indptr[row + 1])
    return (
        np.asarray(group["indices"][start:stop], dtype=np.int64),
        np.asarray(group["data"][start:stop], dtype=np.float32),
    )


def csr_h5_control_samples(
    path: Path,
    label_col: str,
    control_label: str,
    n_genes: int,
    *,
    sample_per_rep: int,
    seed: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    import h5py

    combined = np.zeros((2 * sample_per_rep, n_genes), dtype=np.float32)
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        if not isinstance(x, h5py.Group):
            raise TypeError(f"{path}: expected sparse X group")
        codes, control_code = categorical_code(
            handle, f"obs/{label_col}", control_label
        )
        control_rows = np.flatnonzero(codes == control_code)
        if len(control_rows) < len(combined):
            raise ValueError(f"{path}: not enough control rows")
        rng = np.random.default_rng(seed)
        selected = rng.choice(control_rows, len(combined), replace=False)
        for destination, row in enumerate(selected):
            indices, values = read_csr_row(x, int(row))
            total = float(values.sum(dtype=np.float64))
            if total <= 0:
                raise ValueError(f"{path}: zero library at row {row}")
            combined[destination, indices] = np.log1p(
                np.float32(1_000_000.0 / total) * values
            )
            if destination and destination % 256 == 0:
                enforce_memory_limit(
                    max_rss_gib, f"HCT116 controls {destination}"
                )
    return combined[:sample_per_rep], combined[sample_per_rep:], {
        "source": str(path.resolve()),
        "available_controls": int(len(control_rows)),
        "sample_per_replicate": sample_per_rep,
        "selection": "disjoint seeded sample without replacement",
        "seed": seed,
    }


def is_control(value: str) -> bool:
    lowered = str(value).lower()
    return any(re.match(pattern, lowered) for pattern in CONTROL_PATTERNS)


def orion_token_to_official(official: Sequence[str]) -> np.ndarray:
    import pyarrow.parquet as pq

    table = pq.read_table(
        ORION_GENE_META, columns=["gene_token_id", "gene_name"]
    ).to_pydict()
    max_token = max(map(int, table["gene_token_id"]))
    mapping = np.full(max_token + 1, -1, dtype=np.int32)
    official_pos = {gene: index for index, gene in enumerate(official)}
    for token, gene in zip(
        table["gene_token_id"], table["gene_name"], strict=True
    ):
        destination = official_pos.get(str(gene))
        if destination is not None:
            mapping[int(token)] = destination
    return mapping


def selected_orion_files(line: str, n_strata: int) -> tuple[list[Path], int]:
    files = sorted(
        ORION_DATA.glob(f"{line}_Batch*.parquet"),
        key=lambda path: int(re.search(r"Batch(\d+)", path.name)[1]),
    )
    selected_indices = np.unique(
        np.rint(np.linspace(0, len(files) - 1, n_strata)).astype(int)
    )
    if len(selected_indices) != n_strata:
        raise ValueError("could not select the requested Orion strata")
    return [files[index] for index in selected_indices], len(files)


def orion_control_samples(
    line: str,
    official: Sequence[str],
    *,
    sample_per_rep: int,
    n_strata: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    if sample_per_rep % n_strata:
        raise ValueError("sample_per_rep must be divisible by n_strata")
    quota = sample_per_rep // n_strata
    selected_files, available_files = selected_orion_files(line, n_strata)
    samples = [
        np.zeros((sample_per_rep, len(official)), dtype=np.float32),
        np.zeros((sample_per_rep, len(official)), dtype=np.float32),
    ]
    offsets = [0, 0]
    file_counts: dict[str, int] = {}
    worker_peaks: dict[str, float] = {}
    for file_index, path in enumerate(selected_files, start=1):
        wanted = 2 * quota
        worker_output = ORION_WORKER_DIR / f"{path.stem}.npy"
        worker_metadata = worker_output.with_suffix(".json")
        if not worker_output.exists() or not worker_metadata.exists():
            raise FileNotFoundError(
                f"missing {worker_output}; run --stage orion-workers first"
            )
        worker_detail = json.loads(worker_metadata.read_text())
        expected_fingerprint = file_fingerprint(path)
        if worker_detail.get("source_fingerprint") != expected_fingerprint:
            raise ValueError(f"{path}: stale worker source fingerprint")
        if int(worker_detail.get("wanted", -1)) != wanted:
            raise ValueError(f"{path}: stale worker sample count")
        extracted = np.load(worker_output, mmap_mode="r")
        if extracted.shape != (wanted, len(official)):
            raise ValueError(f"{path}: worker returned shape {extracted.shape}")
        samples[0][offsets[0] : offsets[0] + quota] = extracted[0::2]
        samples[1][offsets[1] : offsets[1] + quota] = extracted[1::2]
        offsets[0] += quota
        offsets[1] += quota
        worker_peaks[path.name] = float(worker_detail["peak_rss_gib"])
        del extracted
        file_counts[path.name] = wanted
        gc.collect()
        enforce_memory_limit(
            max_rss_gib, f"HEK293T control stratum {file_index}"
        )
        print(
            f"  HEK293T control stratum {file_index}/{len(selected_files)} "
            f"samples={offsets}",
            flush=True,
        )
    if offsets != [sample_per_rep, sample_per_rep]:
        raise RuntimeError(f"HEK293T sample counts are {offsets}")
    return samples[0], samples[1], {
        "source": str(ORION_DATA.resolve()),
        "line": line,
        "available_files": available_files,
        "sampled_files": [path.name for path in selected_files],
        "sample_per_replicate": sample_per_rep,
        "controls_read_per_file": file_counts,
        "worker_peak_rss_gib": worker_peaks,
    }


def extract_orion_control_file(
    path: Path,
    output: Path,
    *,
    wanted: int,
    max_rss_gib: float,
) -> dict[str, object]:
    """Extract a tiny control sample in a short-lived process."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    official = read_csv_column(
        ROOT / "data/challenge/gene_names.csv", "gene_name"
    )
    token_map = orion_token_to_official(official)
    parquet = pq.ParquetFile(path)
    metadata = parquet.read(
        columns=["gene_target", "pass_guide_filter"],
        use_threads=False,
    )
    targets = metadata.column("gene_target").combine_chunks()
    passed = metadata.column("pass_guide_filter").combine_chunks()
    selected_rows = [
        row
        for row in range(len(targets))
        if int(passed[row].as_py()) == 1
        and is_control(str(targets[row].as_py()))
    ][:wanted]
    if len(selected_rows) != wanted:
        raise ValueError(
            f"{path}: found {len(selected_rows)} filtered controls, need {wanted}"
        )
    del metadata, targets, passed
    pa.default_memory_pool().release_unused()
    gc.collect()
    enforce_memory_limit(max_rss_gib, f"{path.name} metadata")

    ids_by_selected_row: list[np.ndarray | None] = [None] * wanted
    selected_cursor = 0
    row_cursor = 0
    for batch in parquet.iter_batches(
        batch_size=16,
        columns=["gene_token_id"],
        use_threads=False,
    ):
        ids_column = batch.column(0)
        batch_stop = row_cursor + batch.num_rows
        while (
            selected_cursor < wanted
            and selected_rows[selected_cursor] < batch_stop
        ):
            local_row = selected_rows[selected_cursor] - row_cursor
            ids_by_selected_row[selected_cursor] = (
                ids_column[local_row]
                .values.to_numpy(zero_copy_only=False)
                .astype(np.int64, copy=True)
            )
            selected_cursor += 1
        row_cursor = batch_stop
        del ids_column, batch
        enforce_memory_limit(max_rss_gib, f"{path.name} token scan")
        if selected_cursor == wanted:
            break
    if selected_cursor != wanted or any(
        values is None for values in ids_by_selected_row
    ):
        raise RuntimeError(f"{path}: incomplete token extraction")
    pa.default_memory_pool().release_unused()
    gc.collect()
    enforce_memory_limit(max_rss_gib, f"{path.name} tokens")

    extracted = np.zeros((wanted, len(official)), dtype=np.float32)
    selected_cursor = 0
    row_cursor = 0
    for batch in parquet.iter_batches(
        batch_size=16,
        columns=["gene_expression"],
        use_threads=False,
    ):
        counts_column = batch.column(0)
        batch_stop = row_cursor + batch.num_rows
        while (
            selected_cursor < wanted
            and selected_rows[selected_cursor] < batch_stop
        ):
            local_row = selected_rows[selected_cursor] - row_cursor
            counts_array = (
                counts_column[local_row]
                .values.to_numpy(zero_copy_only=False)
                .astype(np.float32, copy=False)
            )
            ids_array = ids_by_selected_row[selected_cursor]
            if ids_array is None or len(ids_array) != len(counts_array):
                raise RuntimeError(f"{path}: token/expression mismatch")
            official_columns = token_map[ids_array]
            keep = official_columns >= 0
            total = float(counts_array.sum(dtype=np.float64))
            if total <= 0:
                raise ValueError(f"{path}: zero control library")
            extracted[selected_cursor, official_columns[keep]] = np.log1p(
                np.float32(1_000_000.0 / total) * counts_array[keep]
            )
            ids_by_selected_row[selected_cursor] = None
            selected_cursor += 1
        row_cursor = batch_stop
        del counts_column, batch
        enforce_memory_limit(max_rss_gib, f"{path.name} expression scan")
        if selected_cursor == wanted:
            break
    if selected_cursor != wanted:
        raise RuntimeError(f"{path}: incomplete expression extraction")
    np.save(output, extracted, allow_pickle=False)
    payload = {
        "source": str(path.resolve()),
        "source_fingerprint": file_fingerprint(path),
        "wanted": wanted,
        "selected_rows": selected_rows,
        "selected_rows_sha256": sha256_bytes(
            np.asarray(selected_rows, dtype=np.int64).tobytes()
        ),
        "peak_rss_gib": peak_rss_gib(),
    }
    output.with_suffix(".json").write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def build_orion_worker_cache(
    line: str,
    *,
    sample_per_rep: int,
    n_strata: int,
    max_rss_gib: float,
) -> dict[str, object]:
    if sample_per_rep % n_strata:
        raise ValueError("sample_per_rep must be divisible by n_strata")
    selected_files, available_files = selected_orion_files(line, n_strata)
    wanted = 2 * (sample_per_rep // n_strata)
    ORION_WORKER_DIR.mkdir(parents=True, exist_ok=True)
    worker_peaks = {}
    for index, path in enumerate(selected_files, start=1):
        output = ORION_WORKER_DIR / f"{path.stem}.npy"
        metadata_path = output.with_suffix(".json")
        output.unlink(missing_ok=True)
        metadata_path.unlink(missing_ok=True)
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--stage",
                "orion-worker",
                "--orion-file",
                str(path),
                "--orion-output",
                str(output),
                "--orion-wanted",
                str(wanted),
                "--max-rss-gib",
                str(max_rss_gib),
            ],
            check=True,
        )
        worker_detail = json.loads(metadata_path.read_text())
        worker_peaks[path.name] = float(worker_detail["peak_rss_gib"])
        print(
            f"  cached HEK293T stratum {index}/{len(selected_files)} "
            f"worker_peak={worker_peaks[path.name]:.3f} GiB",
            flush=True,
        )
    return {
        "line": line,
        "available_files": available_files,
        "sampled_files": [path.name for path in selected_files],
        "wanted_per_file": wanted,
        "worker_peak_rss_gib": worker_peaks,
        "peak_rss_gib": peak_rss_gib(),
    }


def embed_to_official(
    values: np.ndarray,
    source_genes: Sequence[str],
    official: Sequence[str],
) -> np.ndarray:
    output = np.full(len(official), np.nan, dtype=np.float32)
    official_pos = {gene: index for index, gene in enumerate(official)}
    for source_index, gene in enumerate(source_genes):
        destination = official_pos.get(gene)
        if destination is not None:
            output[destination] = values[source_index]
    return output


def stats_stability(
    first: np.ndarray,
    second: np.ndarray,
) -> dict[str, float | int]:
    valid = (
        np.isfinite(first)
        & np.isfinite(second)
        & (first > 0)
        & (second > 0)
    )
    if valid.sum() < 2:
        return {
            "valid_genes": int(valid.sum()),
            "pearson": math.nan,
            "median_abs_log_ratio": math.nan,
        }
    a = np.log(first[valid])
    b = np.log(second[valid])
    return {
        "valid_genes": int(valid.sum()),
        "pearson": float(np.corrcoef(a, b)[0, 1]),
        "median_abs_log_ratio": float(np.median(np.abs(a - b))),
    }


def control_cache_signature(
    sample_per_rep: int,
    n_strata: int,
) -> dict[str, object]:
    files = {
        "k562_h5ad": file_fingerprint(K562_H5AD),
        "hct116_h5ad": file_fingerprint(HCT116_H5AD),
        "orion_gene_metadata": file_fingerprint(ORION_GENE_META),
    }
    hek_files = sorted(
        ORION_DATA.glob("HEK293T_Batch*.parquet"),
        key=lambda path: int(re.search(r"Batch(\d+)", path.name)[1]),
    )
    files["hek_file_count"] = len(hek_files)
    files["hek_total_bytes"] = sum(path.stat().st_size for path in hek_files)
    return {
        "version": 1,
        "sample_per_replicate": sample_per_rep,
        "strata": n_strata,
        "target_sum": 1_000_000.0,
        "sd_definition": "1%-99% winsorized sample SD of log1p(CPM)",
        "mad_definition": "1.4826 * median absolute deviation of log1p(CPM)",
        "files": files,
    }


def build_control_cache(
    *,
    sample_per_rep: int,
    n_strata: int,
    max_rss_gib: float,
) -> dict[str, object]:
    signature = control_cache_signature(sample_per_rep, n_strata)
    official = read_csv_column(
        ROOT / "data/challenge/gene_names.csv", "gene_name"
    )
    arrays: dict[str, np.ndarray] = {}
    details: dict[str, object] = {}

    print("extracting K562 control variability", flush=True)
    k562_genes = read_lines(BANK_PATHS["k562_gwps"] / "genes.txt")
    sample_a, sample_b, detail = dense_h5_control_samples(
        K562_H5AD,
        "perturbation",
        "control",
        len(k562_genes),
        sample_per_rep=sample_per_rep,
        n_strata=n_strata,
        max_rss_gib=max_rss_gib,
    )
    for rep, sample in zip(CONTROL_REPLICATES, (sample_a, sample_b), strict=True):
        for method, values in robust_scales(sample).items():
            arrays[f"k562_gwps_{rep}_{method}"] = embed_to_official(
                values, k562_genes, official
            )
    del sample_a, sample_b
    details["k562_gwps"] = detail
    gc.collect()
    enforce_memory_limit(max_rss_gib, "K562 control statistics")

    print("extracting HCT116 control variability", flush=True)
    sample_a, sample_b, detail = csr_h5_control_samples(
        HCT116_H5AD,
        "target_gene",
        "non-targeting",
        len(official),
        sample_per_rep=sample_per_rep,
        seed=670116,
        max_rss_gib=max_rss_gib,
    )
    for rep, sample in zip(CONTROL_REPLICATES, (sample_a, sample_b), strict=True):
        for method, values in robust_scales(sample).items():
            arrays[f"orion_hct116_{rep}_{method}"] = values
    del sample_a, sample_b
    details["orion_hct116"] = detail
    gc.collect()
    enforce_memory_limit(max_rss_gib, "HCT116 control statistics")

    print("extracting HEK293T control variability", flush=True)
    sample_a, sample_b, detail = orion_control_samples(
        "HEK293T",
        official,
        sample_per_rep=sample_per_rep,
        n_strata=n_strata,
        max_rss_gib=max_rss_gib,
    )
    for rep, sample in zip(CONTROL_REPLICATES, (sample_a, sample_b), strict=True):
        for method, values in robust_scales(sample).items():
            arrays[f"orion_hek293t_{rep}_{method}"] = values
    del sample_a, sample_b
    details["orion_hek293t"] = detail
    gc.collect()
    enforce_memory_limit(max_rss_gib, "HEK293T control statistics")

    stability = {}
    for donor in MAIN_NAMES:
        stability[donor] = {
            method: stats_stability(
                arrays[f"{donor}_a_{method}"],
                arrays[f"{donor}_b_{method}"],
            )
            for method in METHODS
        }
    payload = {
        **signature,
        "details": details,
        "stability": stability,
        "peak_rss_gib": peak_rss_gib(),
    }
    np.savez_compressed(CACHE_NPZ, **arrays)
    for name in details["orion_hek293t"]["sampled_files"]:
        stem = Path(name).stem
        (ORION_WORKER_DIR / f"{stem}.npy").unlink(missing_ok=True)
        (ORION_WORKER_DIR / f"{stem}.json").unlink(missing_ok=True)
    try:
        ORION_WORKER_DIR.rmdir()
    except OSError:
        pass
    details["orion_hek293t"]["worker_cache_removed_after_use"] = True
    CACHE_JSON.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote {CACHE_NPZ.relative_to(ROOT)}; "
        f"peak RSS={peak_rss_gib():.3f} GiB",
        flush=True,
    )
    return payload


def load_control_cache(
    sample_per_rep: int,
    n_strata: int,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    if not CACHE_NPZ.exists() or not CACHE_JSON.exists():
        raise FileNotFoundError("control cache is missing; run --stage controls")
    metadata = json.loads(CACHE_JSON.read_text())
    expected = control_cache_signature(sample_per_rep, n_strata)
    for key in (
        "version",
        "sample_per_replicate",
        "strata",
        "target_sum",
        "sd_definition",
        "mad_definition",
        "files",
    ):
        if metadata.get(key) != expected.get(key):
            raise ValueError(f"stale control cache: mismatch at {key}")
    with np.load(CACHE_NPZ) as archive:
        arrays = {key: np.asarray(archive[key]) for key in archive.files}
    return arrays, metadata


@dataclass(frozen=True)
class NpyLayout:
    shape: tuple[int, int]
    dtype: np.dtype
    data_offset: int
    row_bytes: int


def read_npy_layout(path: Path) -> NpyLayout:
    with path.open("rb") as handle:
        version = np.lib.format.read_magic(handle)
        if version == (1, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(handle)
        elif version == (2, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_2_0(handle)
        else:
            raise ValueError(f"{path}: unsupported NPY version {version}")
        data_offset = handle.tell()
    if fortran or len(shape) != 2:
        raise ValueError(f"{path}: expected a C-contiguous two-dimensional array")
    dtype = np.dtype(dtype)
    return NpyLayout(
        shape=(int(shape[0]), int(shape[1])),
        dtype=dtype,
        data_offset=data_offset,
        row_bytes=int(shape[1]) * dtype.itemsize,
    )


def iter_npy_row_blocks(
    path: Path,
    chunk_rows: int,
) -> Iterator[tuple[int, np.ndarray]]:
    layout = read_npy_layout(path)
    with path.open("rb") as handle:
        handle.seek(layout.data_offset)
        for lo in range(0, layout.shape[0], chunk_rows):
            count = min(chunk_rows, layout.shape[0] - lo)
            payload = handle.read(count * layout.row_bytes)
            if len(payload) != count * layout.row_bytes:
                raise EOFError(f"{path}: truncated array payload at row {lo}")
            yield lo, np.frombuffer(payload, dtype=layout.dtype).reshape(
                count, layout.shape[1]
            )


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: dict[str, int]
    delta_path: Path
    shape: tuple[int, int]

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        layout = read_npy_layout(delta_path)
        if layout.shape != (len(perts), len(genes)):
            raise ValueError(f"{name}: inconsistent donor table")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            shape=layout.shape,
        )


@dataclass(frozen=True)
class AxisMap:
    all_source: np.ndarray
    unique_source: np.ndarray
    unique_target: np.ndarray


def map_bank_to_axis(bank: Bank, official: Sequence[str]) -> AxisMap:
    official_pos = {gene: index for index, gene in enumerate(official)}
    all_source: list[int] = []
    last_source: dict[int, int] = {}
    for source, gene in enumerate(bank.genes):
        destination = official_pos.get(gene)
        if destination is None:
            continue
        all_source.append(source)
        last_source[destination] = source
    unique_target = np.asarray(sorted(last_source), dtype=np.int64)
    unique_source = np.asarray(
        [last_source[position] for position in unique_target], dtype=np.int64
    )
    return AxisMap(
        all_source=np.asarray(all_source, dtype=np.int64),
        unique_source=unique_source,
        unique_target=unique_target,
    )


def stream_column_mean(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    total = np.zeros(len(columns), dtype=np.float64)
    for lo, rows in iter_npy_row_blocks(bank.delta_path, chunk_rows):
        block = np.asarray(rows[:, columns], dtype=np.float32)
        total += block.sum(axis=0, dtype=np.float64)
        del block
        enforce_memory_limit(max_rss_gib, f"{bank.name} centering")
    return total / bank.shape[0]


def median_centered_norm(
    bank: Bank,
    columns: np.ndarray,
    mean: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> float:
    norms = np.empty(bank.shape[0], dtype=np.float32)
    mean32 = mean.astype(np.float32)
    for lo, rows in iter_npy_row_blocks(bank.delta_path, chunk_rows):
        centered = np.asarray(rows[:, columns], dtype=np.float32)
        centered -= mean32
        norms[lo : lo + len(rows)] = np.linalg.norm(centered, axis=1)
        del centered
        enforce_memory_limit(max_rss_gib, f"{bank.name} norm")
    return float(np.median(norms)) or 1.0


@dataclass
class PreparedSource:
    raw_truth: np.ndarray
    available: np.ndarray
    mapped_genes: int
    raw_median_norm: float

    @property
    def effect_scale(self) -> float:
        return EFFECT_SCALE / self.raw_median_norm


def prepare_source(
    bank: Bank,
    targets: Sequence[str],
    official: Sequence[str],
    panel_positions: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> PreparedSource:
    mapping = map_bank_to_axis(bank, official)
    mean_all = stream_column_mean(
        bank,
        mapping.all_source,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median = median_centered_norm(
        bank,
        mapping.all_source,
        mean_all,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    all_position = {
        int(source): index
        for index, source in enumerate(mapping.all_source.tolist())
    }
    unique_mean = np.asarray(
        [mean_all[all_position[int(source)]] for source in mapping.unique_source],
        dtype=np.float32,
    )
    raw = np.zeros((len(targets), len(official)), dtype=np.float32)
    available = np.zeros(len(targets), dtype=bool)
    layout = read_npy_layout(bank.delta_path)
    with bank.delta_path.open("rb") as handle:
        for index, target in enumerate(targets):
            row = bank.pert_pos.get(target)
            if row is None:
                continue
            handle.seek(layout.data_offset + row * layout.row_bytes)
            payload = handle.read(layout.row_bytes)
            if len(payload) != layout.row_bytes:
                raise EOFError(
                    f"{bank.delta_path}: truncated payload at row {row}"
                )
            full_row = np.frombuffer(payload, dtype=layout.dtype)
            values = np.asarray(
                full_row[mapping.unique_source], dtype=np.float32
            )
            raw[index, mapping.unique_target] = values - unique_mean
            available[index] = True
    raw[:, panel_positions] = 0.0
    return PreparedSource(
        raw_truth=raw,
        available=available,
        mapped_genes=len(mapping.unique_target),
        raw_median_norm=raw_median,
    )


def union_columns(
    source_names: Sequence[str],
    banks: dict[str, Bank],
    official: Sequence[str],
    panel: set[str],
) -> np.ndarray:
    union = set().union(*(set(banks[name].genes) for name in source_names))
    return np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in union and gene not in panel
        ],
        dtype=np.int64,
    )


@dataclass(frozen=True)
class StoredNpyLayout:
    path: Path
    shape: tuple[int, ...]
    dtype: np.dtype
    data_offset: int
    order: str


def stored_npy_layout(
    path: Path,
    member: str,
) -> StoredNpyLayout:
    """Locate an uncompressed NPY member stored inside an NPZ archive."""
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(member)
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"{path}:{member} is compressed and cannot be mmaped")
        header_offset = info.header_offset
    with path.open("rb") as handle:
        handle.seek(header_offset)
        local_header = handle.read(30)
        fields = struct.unpack("<IHHHHHIIIHH", local_header)
        if fields[0] != 0x04034B50:
            raise ValueError(f"{path}:{member} has an invalid ZIP header")
        name_length, extra_length = fields[-2:]
        handle.seek(header_offset + 30 + name_length + extra_length)
        version = np.lib.format.read_magic(handle)
        if version == (1, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_1_0(handle)
        elif version == (2, 0):
            shape, fortran, dtype = np.lib.format.read_array_header_2_0(handle)
        else:
            raise ValueError(f"{path}:{member}: unsupported NPY version {version}")
        data_offset = handle.tell()
    return StoredNpyLayout(
        path=path,
        shape=tuple(map(int, shape)),
        dtype=np.dtype(dtype),
        data_offset=data_offset,
        order="F" if fortran else "C",
    )


def load_stored_npy_member(path: Path, member: str) -> np.ndarray:
    layout = stored_npy_layout(path, member)
    count = math.prod(layout.shape)
    with path.open("rb") as handle:
        handle.seek(layout.data_offset)
        payload = handle.read(count * layout.dtype.itemsize)
    if len(payload) != count * layout.dtype.itemsize:
        raise EOFError(f"{path}:{member}: truncated array payload")
    return np.frombuffer(payload, dtype=layout.dtype).reshape(
        layout.shape, order=layout.order
    ).copy()


def iter_operator_row_blocks(
    operator: np.ndarray | StoredNpyLayout,
    block_rows: int,
) -> Iterator[tuple[int, np.ndarray]]:
    if isinstance(operator, StoredNpyLayout):
        if operator.order != "C" or len(operator.shape) != 2:
            raise ValueError("stored operator must be a C-order matrix")
        row_bytes = operator.shape[1] * operator.dtype.itemsize
        with operator.path.open("rb") as handle:
            handle.seek(operator.data_offset)
            for lo in range(0, operator.shape[0], block_rows):
                count = min(block_rows, operator.shape[0] - lo)
                payload = handle.read(count * row_bytes)
                if len(payload) != count * row_bytes:
                    raise EOFError(
                        f"{operator.path}: truncated operator at row {lo}"
                    )
                yield lo, np.frombuffer(
                    payload, dtype=operator.dtype
                ).reshape(count, operator.shape[1])
        return
    for lo in range(0, operator.shape[0], block_rows):
        hi = min(operator.shape[0], lo + block_rows)
        yield lo, np.asarray(operator[lo:hi])


def apply_operator_rows(
    matrix: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray | StoredNpyLayout,
    *,
    block_rows: int,
    max_rss_gib: float,
    label: str,
) -> np.ndarray:
    """Compute matrix @ V @ diag(lambda^-gamma/2) @ V.T in bounded memory."""
    x = np.asarray(matrix, dtype=np.float32)
    shape = evecs.shape
    if x.ndim != 2 or x.shape[1] != shape[0]:
        raise ValueError(f"{label}: matrix/operator shape mismatch")
    if shape[0] != shape[1] or len(evals) != shape[0]:
        raise ValueError(f"{label}: invalid eigensystem")
    latent = np.zeros((x.shape[0], x.shape[1]), dtype=np.float32)
    temporary = np.empty_like(latent)
    for lo, rows in iter_operator_row_blocks(evecs, block_rows):
        hi = lo + len(rows)
        block = np.array(rows, dtype=np.float32, copy=True)
        np.matmul(x[:, lo:hi], block, out=temporary)
        latent += temporary
        del block
        enforce_memory_limit(max_rss_gib, f"{label} forward {hi}")
    del temporary
    weights = np.power(
        np.asarray(evals, dtype=np.float64), -WHITEN_GAMMA / 2.0
    ).astype(np.float32)
    latent *= weights
    output = np.empty_like(x)
    for lo, rows in iter_operator_row_blocks(evecs, block_rows):
        hi = lo + len(rows)
        block = np.array(rows, dtype=np.float32, copy=True)
        destination = np.empty((x.shape[0], hi - lo), dtype=np.float32)
        np.matmul(latent, block.T, out=destination)
        output[:, lo:hi] = destination
        del block, destination
        enforce_memory_limit(max_rss_gib, f"{label} reverse {hi}")
    return output


def row_normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, np.float32(1e-12))


def pds_from_scores(scores: np.ndarray) -> tuple[float, np.ndarray]:
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    ranks = (scores > own).sum(axis=1)
    ranks = ranks + 0.5 * ((scores == own).sum(axis=1) - 1)
    per_target = 1.0 - ranks / (n - 1)
    return float(np.mean(per_target)), per_target


@dataclass(frozen=True)
class Candidate:
    method: str
    placement: str
    eta: float

    @property
    def key(self) -> str:
        if self.eta == 0:
            return "identity:eta=0"
        return f"{self.method}:{self.placement}:eta={self.eta:g}"

    @property
    def complexity(self) -> tuple[float, int, int]:
        return (
            self.eta,
            0 if self.method == "sd" else 1,
            0 if self.placement == "pre" else 1,
        )


def candidate_grid() -> tuple[Candidate, ...]:
    return (
        Candidate("identity", "none", 0.0),
        *(
            Candidate(method, placement, eta)
            for method in METHODS
            for placement in PLACEMENTS
            for eta in ETAS
        ),
    )


def variance_ratio(
    controls: dict[str, np.ndarray],
    target: str,
    sources: Sequence[str],
    method: str,
    replicate: str,
    columns: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    target_scale = controls[f"{target}_{replicate}_{method}"][columns]
    source_values = np.vstack(
        [controls[f"{name}_{replicate}_{method}"][columns] for name in sources]
    )
    valid_source = np.isfinite(source_values) & (source_values > 0)
    log_values = np.zeros_like(source_values, dtype=np.float64)
    np.log(source_values, out=log_values, where=valid_source)
    counts = valid_source.sum(axis=0)
    source_scale = np.full(len(columns), np.nan, dtype=np.float64)
    valid_any = counts > 0
    source_scale[valid_any] = (
        np.exp(log_values[:, valid_any].sum(axis=0) / counts[valid_any])
    )
    valid = (
        np.isfinite(target_scale)
        & (target_scale > 0)
        & np.isfinite(source_scale)
        & (source_scale > 0)
    )
    ratio = np.ones(len(columns), dtype=np.float32)
    ratio[valid] = np.clip(
        target_scale[valid] / source_scale[valid],
        RATIO_CLIP[0],
        RATIO_CLIP[1],
    ).astype(np.float32)
    return ratio, {
        "valid_genes": int(valid.sum()),
        "identity_genes": int((~valid).sum()),
        "clipped_low": int(np.sum(valid & (target_scale / source_scale < RATIO_CLIP[0]))),
        "clipped_high": int(np.sum(valid & (target_scale / source_scale > RATIO_CLIP[1]))),
        "ratio_median": float(np.median(ratio[valid])) if valid.any() else 1.0,
        "ratio_p10": float(np.percentile(ratio[valid], 10)) if valid.any() else 1.0,
        "ratio_p90": float(np.percentile(ratio[valid], 90)) if valid.any() else 1.0,
    }


@dataclass
class Role:
    name: str
    target: str
    sources: tuple[str, ...]
    genes: np.ndarray
    prediction: np.ndarray
    truth_unit: np.ndarray
    prediction_whitened: np.ndarray | None = None
    truth_unit_whitened: np.ndarray | None = None


def evaluate_role_grid(
    role: Role,
    controls: dict[str, np.ndarray],
    candidates: Sequence[Candidate],
    replicate: str,
) -> tuple[dict[str, dict[str, object]], dict[str, np.ndarray], dict[str, dict]]:
    if role.prediction_whitened is None or role.truth_unit_whitened is None:
        raise RuntimeError(f"{role.name}: whitening outputs are absent")
    baseline_scores = role.prediction_whitened @ role.truth_unit.T
    baseline, baseline_per_target = pds_from_scores(baseline_scores)
    rows: dict[str, dict[str, object]] = {}
    per_target: dict[str, np.ndarray] = {}
    ratio_meta: dict[str, dict] = {}
    for candidate in candidates:
        if candidate.eta == 0:
            score = baseline
            values = baseline_per_target
        else:
            ratio, metadata = variance_ratio(
                controls,
                role.target,
                role.sources,
                candidate.method,
                replicate,
                role.genes,
            )
            ratio_meta[candidate.method] = metadata
            factor = np.power(ratio, candidate.eta).astype(np.float32)
            if candidate.placement == "pre":
                similarities = (
                    role.prediction * factor
                ) @ role.truth_unit_whitened.T
            elif candidate.placement == "post":
                similarities = (
                    role.prediction_whitened * factor
                ) @ role.truth_unit.T
            else:
                raise ValueError(candidate.placement)
            score, values = pds_from_scores(similarities)
        rows[candidate.key] = {
            "pds": score,
            "baseline_pds": baseline,
            "delta": score - baseline,
        }
        per_target[candidate.key] = values
    return rows, per_target, ratio_meta


def choose_candidate(
    inner_results: Sequence[dict[str, dict[str, object]]],
    candidates: Sequence[Candidate],
    tie_tolerance: float,
) -> dict[str, object]:
    by_key = {candidate.key: candidate for candidate in candidates}
    summaries = []
    for candidate in candidates:
        deltas = [
            float(result[candidate.key]["delta"]) for result in inner_results
        ]
        summaries.append(
            {
                "candidate": candidate.key,
                "worst_delta": float(min(deltas)),
                "mean_delta": float(np.mean(deltas)),
                "direction_deltas": deltas,
            }
        )
    best_worst = max(row["worst_delta"] for row in summaries)
    near_worst = [
        row for row in summaries
        if row["worst_delta"] >= best_worst - tie_tolerance
    ]
    best_mean = max(row["mean_delta"] for row in near_worst)
    near = [
        row for row in near_worst
        if row["mean_delta"] >= best_mean - tie_tolerance
    ]
    selected = min(
        near,
        key=lambda row: by_key[row["candidate"]].complexity,
    )
    return {
        **selected,
        "selection_rule": (
            "maximize worst source-to-source delta, then mean delta; "
            f"within {tie_tolerance:g}, prefer lower eta, SD, then pre"
        ),
        "grid": summaries,
    }


def permuted_ratio_null(
    role: Role,
    controls: dict[str, np.ndarray],
    candidate: Candidate,
    replicate: str,
    baseline: float,
    *,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    if candidate.eta == 0:
        return {
            "permutations": permutations,
            "observed_delta": 0.0,
            "null_mean": 0.0,
            "null_p95": 0.0,
            "empirical_p": 1.0,
        }
    ratio, _ = variance_ratio(
        controls,
        role.target,
        role.sources,
        candidate.method,
        replicate,
        role.genes,
    )
    movable = np.flatnonzero(np.abs(ratio - 1.0) > 1e-8)
    rng = np.random.default_rng(seed)

    def score_factor(factor: np.ndarray) -> float:
        if candidate.placement == "pre":
            scores = (
                role.prediction * factor
            ) @ role.truth_unit_whitened.T
        else:
            scores = (
                role.prediction_whitened * factor
            ) @ role.truth_unit.T
        return pds_from_scores(scores)[0]

    observed = score_factor(np.power(ratio, candidate.eta).astype(np.float32))
    null = []
    for _ in range(permutations):
        shuffled = ratio.copy()
        shuffled[movable] = ratio[rng.permutation(movable)]
        null.append(
            score_factor(np.power(shuffled, candidate.eta).astype(np.float32))
            - baseline
        )
    null_array = np.asarray(null)
    observed_delta = observed - baseline
    return {
        "permutations": permutations,
        "movable_genes": int(len(movable)),
        "observed_delta": observed_delta,
        "null_mean": float(np.mean(null_array)),
        "null_p05": float(np.percentile(null_array, 5)),
        "null_p95": float(np.percentile(null_array, 95)),
        "empirical_p": float(
            (1 + np.sum(null_array >= observed_delta)) / (permutations + 1)
        ),
    }


def write_per_target(rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    with PER_TARGET_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def evaluate(
    *,
    sample_per_rep: int,
    n_strata: int,
    selection_targets: int,
    chunk_rows: int,
    operator_block: int,
    permutations: int,
    tie_tolerance: float,
    max_rss_gib: float,
) -> dict[str, object]:
    controls, control_metadata = load_control_cache(
        sample_per_rep, n_strata
    )
    official = read_csv_column(
        ROOT / "data/challenge/gene_names.csv", "gene_name"
    )
    panel = read_csv_column(
        ROOT / "data/challenge/pert_counts.csv", "target_gene"
    )
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[gene] for gene in panel], dtype=np.int64
    )
    banks = {
        name: Bank.load(name, path) for name, path in BANK_PATHS.items()
    }
    targets = tuple(
        gene
        for gene in panel
        if all(gene in banks[name].pert_pos for name in MAIN_NAMES)
    )
    if len(targets) != 259:
        raise ValueError(f"expected 259 shared VCC targets, got {len(targets)}")
    prepared = {}
    for name in MAIN_NAMES:
        print(f"preparing effects: {name}", flush=True)
        prepared[name] = prepare_source(
            banks[name],
            targets,
            official,
            panel_positions,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        enforce_memory_limit(max_rss_gib, f"{name} prepared")
    candidates = candidate_grid()
    by_key = {candidate.key: candidate for candidate in candidates}
    selection_rng = np.random.default_rng(670067)
    selection_indices = np.sort(
        selection_rng.choice(
            len(targets), min(selection_targets, len(targets)), replace=False
        )
    )
    folds: dict[str, dict[str, object]] = {}
    per_target_rows: list[dict[str, object]] = []

    for fold_index, held in enumerate(MAIN_NAMES):
        sources = tuple(name for name in MAIN_NAMES if name != held)
        print(f"\nfold held={held} sources={sources}", flush=True)
        wcols = union_columns(sources, banks, official, panel_set)
        operator_path = WHITENER_DIR / f"held_{held}.npz"
        operator_metadata = json.loads(
            operator_path.with_suffix(".json").read_text()
        )
        if tuple(operator_metadata["donors"]) != sources:
            raise ValueError(
                f"{held}: whitener donors {operator_metadata['donors']} != {sources}"
            )
        if held in operator_metadata["donors"]:
            raise AssertionError(f"{held}: held donor leaked into whitener")
        if int(operator_metadata["n_cols"]) != len(wcols):
            raise ValueError(f"{held}: whitener axis length mismatch")
        if operator_metadata["wcols_sha256"] != sha256_bytes(wcols.tobytes()):
            raise ValueError(f"{held}: whitener axis hash mismatch")
        evals = load_stored_npy_member(operator_path, "evals.npy")
        evecs = stored_npy_layout(operator_path, "evecs.npy")
        if evecs.shape != (len(wcols), len(wcols)):
            raise ValueError(f"{held}: whitener shape {evecs.shape}")
        axis_genes = np.asarray(official, dtype=object)[wcols]

        inner_roles = []
        for target, predictor in ((sources[0], sources[1]), (sources[1], sources[0])):
            inner_roles.append(
                Role(
                    name=f"{predictor}_to_{target}",
                    target=target,
                    sources=(predictor,),
                    genes=wcols,
                    prediction=(
                        prepared[predictor].raw_truth[
                            np.ix_(selection_indices, wcols)
                        ]
                        * np.float32(prepared[predictor].effect_scale)
                    ),
                    truth_unit=row_normalize(
                        prepared[target].raw_truth[
                            np.ix_(selection_indices, wcols)
                        ]
                    ),
                )
            )
        inner_input = np.vstack(
            [role.prediction for role in inner_roles]
            + [role.truth_unit for role in inner_roles]
        ).astype(np.float32, copy=False)
        print(
            f"  whitening inner selection rows={len(inner_input)} "
            f"axis={len(wcols)}",
            flush=True,
        )
        inner_output = apply_operator_rows(
            inner_input,
            evals,
            evecs,
            block_rows=operator_block,
            max_rss_gib=max_rss_gib,
            label=f"{held} inner",
        )
        cursor = 0
        for role in inner_roles:
            count = len(role.prediction)
            role.prediction_whitened = inner_output[cursor : cursor + count]
            cursor += count
        for role in inner_roles:
            count = len(role.truth_unit)
            role.truth_unit_whitened = inner_output[cursor : cursor + count]
            cursor += count
        inner_by_rep = {}
        for replicate in CONTROL_REPLICATES:
            inner_results = []
            inner_details = {}
            for role in inner_roles:
                grid, _, ratio_meta = evaluate_role_grid(
                    role, controls, candidates, replicate
                )
                inner_results.append(grid)
                inner_details[role.name] = {
                    "grid": grid,
                    "ratio_metadata": ratio_meta,
                }
            selection = choose_candidate(
                inner_results, candidates, tie_tolerance
            )
            inner_by_rep[replicate] = {
                "selection": selection,
                "directions": inner_details,
            }
        selected_key = inner_by_rep["a"]["selection"]["candidate"]
        selected_candidate = by_key[selected_key]
        del inner_input, inner_output, inner_roles
        gc.collect()
        enforce_memory_limit(max_rss_gib, f"{held} selection complete")
        print(
            f"  source-only selection A={selected_key}; "
            f"B={inner_by_rep['b']['selection']['candidate']}",
            flush=True,
        )

        source_prediction = np.asarray(
            prepared[sources[0]].raw_truth[:, wcols],
            dtype=np.float32,
        )
        source_prediction *= np.float32(prepared[sources[0]].effect_scale)
        source_prediction += (
            np.asarray(
                prepared[sources[1]].raw_truth[:, wcols],
                dtype=np.float32,
            )
            * np.float32(prepared[sources[1]].effect_scale)
        )
        source_prediction *= np.float32(0.5)
        outer_role = Role(
            name=f"outer_{held}",
            target=held,
            sources=sources,
            genes=wcols,
            prediction=source_prediction,
            truth_unit=row_normalize(
                prepared[held].raw_truth[:, wcols]
            ),
        )
        outer_input = np.vstack(
            [outer_role.prediction, outer_role.truth_unit]
        ).astype(np.float32, copy=False)
        print(
            f"  whitening outer rows={len(outer_input)} axis={len(wcols)}",
            flush=True,
        )
        outer_output = apply_operator_rows(
            outer_input,
            evals,
            evecs,
            block_rows=operator_block,
            max_rss_gib=max_rss_gib,
            label=f"{held} outer",
        )
        n_outer = len(targets)
        outer_role.prediction_whitened = outer_output[:n_outer]
        outer_role.truth_unit_whitened = outer_output[n_outer:]
        outer_by_rep = {}
        outer_per_target = {}
        for replicate in CONTROL_REPLICATES:
            grid, per_target, ratio_meta = evaluate_role_grid(
                outer_role, controls, candidates, replicate
            )
            outer_by_rep[replicate] = {
                "grid": grid,
                "ratio_metadata": ratio_meta,
                "selected": {
                    "candidate": selected_key,
                    **grid[selected_key],
                },
            }
            outer_per_target[replicate] = per_target
        baseline = float(
            outer_by_rep["a"]["grid"]["identity:eta=0"]["pds"]
        )
        symmetry_error = float(
            np.max(
                np.abs(
                    outer_role.prediction_whitened
                    @ outer_role.truth_unit.T
                    - outer_role.prediction
                    @ outer_role.truth_unit_whitened.T
                )
            )
        )
        null = permuted_ratio_null(
            outer_role,
            controls,
            selected_candidate,
            "a",
            baseline,
            permutations=permutations,
            seed=6700 + fold_index,
        )
        for target_index, target in enumerate(targets):
            for replicate in CONTROL_REPLICATES:
                reported_candidates = (
                    ("identity:eta=0",)
                    if selected_key == "identity:eta=0"
                    else ("identity:eta=0", selected_key)
                )
                for key in reported_candidates:
                    values = outer_per_target[replicate][key]
                    per_target_rows.append(
                        {
                            "held": held,
                            "target_gene": target,
                            "control_replicate": replicate,
                            "candidate": key,
                            "pds": float(values[target_index]),
                            "delta_vs_baseline": float(
                                values[target_index]
                                - outer_per_target[replicate][
                                    "identity:eta=0"
                                ][target_index]
                            ),
                        }
                    )
        folds[held] = {
            "sources": list(sources),
            "selection_target_count": int(len(selection_indices)),
            "selection_target_sha256": sha256_bytes(
                "\n".join(targets[index] for index in selection_indices).encode()
            ),
            "held_effects_used_during_selection": False,
            "operator": {
                "path": str(operator_path.resolve()),
                "sha256": sha256_file(operator_path),
                "metadata_sha256": sha256_file(
                    operator_path.with_suffix(".json")
                ),
                "donors": operator_metadata["donors"],
                "genes": len(wcols),
                "retrieval_genes": int(
                    sum(gene in set(banks[held].genes) for gene in axis_genes)
                ),
            },
            "inner": inner_by_rep,
            "selected_candidate": selected_key,
            "outer": outer_by_rep,
            "permuted_gene_ratio_null": null,
            "symmetry_max_abs_error": symmetry_error,
        }
        print(
            f"  outer baseline={baseline:.9f} "
            f"selected A={outer_by_rep['a']['selected']['pds']:.9f} "
            f"delta={outer_by_rep['a']['selected']['delta']:+.9f} "
            f"stability B delta={outer_by_rep['b']['selected']['delta']:+.9f} "
            f"null_p={null['empirical_p']:.4f}",
            flush=True,
        )
        del (
            source_prediction,
            outer_input,
            outer_output,
            outer_role,
            evals,
            evecs,
        )
        gc.collect()
        enforce_memory_limit(max_rss_gib, f"{held} fold complete")

    selected = [folds[name]["selected_candidate"] for name in MAIN_NAMES]
    primary_a = [
        float(folds[name]["outer"]["a"]["selected"]["delta"])
        for name in MAIN_NAMES
    ]
    stability_b = [
        float(folds[name]["outer"]["b"]["selected"]["delta"])
        for name in MAIN_NAMES
    ]
    null_passes = [
        float(folds[name]["permuted_gene_ratio_null"]["observed_delta"])
        > float(folds[name]["permuted_gene_ratio_null"]["null_p95"])
        for name in MAIN_NAMES
    ]
    unanimous = len(set(selected)) == 1 and selected[0] != "identity:eta=0"
    passed = (
        unanimous
        and all(value > 0 for value in primary_a)
        and all(value > 0 for value in stability_b)
        and float(np.mean(primary_a)) >= 0.002
        and sum(null_passes) >= 2
    )

    fixed_grid = {}
    for candidate in candidates:
        key = candidate.key
        deltas = [
            float(folds[name]["outer"]["a"]["grid"][key]["delta"])
            for name in MAIN_NAMES
        ]
        fixed_grid[key] = {
            **{name: delta for name, delta in zip(MAIN_NAMES, deltas, strict=True)},
            "macro_delta": float(np.mean(deltas)),
            "positive_folds": int(np.sum(np.asarray(deltas) > 0)),
            "warning": "exploratory; evaluated on the outer folds",
        }
    payload = {
        "experiment": 67,
        "method": "context-control variance transport",
        "design": {
            "outer_donors": list(MAIN_NAMES),
            "outer_targets": len(targets),
            "selection_targets_per_fold": int(len(selection_indices)),
            "control_estimators": list(METHODS),
            "etas": [0.0, *ETAS],
            "placements": list(PLACEMENTS),
            "ratio_clip": list(RATIO_CLIP),
            "whitening_gamma": WHITEN_GAMMA,
            "held_controls_allowed": True,
            "held_perturbation_effects_used_for_selection": False,
            "cell_eval2": False,
        },
        "control_cache": control_metadata,
        "prepared_sources": {
            name: {
                "mapped_genes": prepared[name].mapped_genes,
                "raw_median_norm": prepared[name].raw_median_norm,
            }
            for name in MAIN_NAMES
        },
        "folds": folds,
        "primary": {
            "selected_candidates": selected,
            "unanimous_nonidentity_selection": unanimous,
            "replicate_a_fold_deltas": primary_a,
            "replicate_a_macro_delta": float(np.mean(primary_a)),
            "replicate_b_fold_deltas": stability_b,
            "replicate_b_macro_delta": float(np.mean(stability_b)),
            "null_exceedances": int(sum(null_passes)),
            "ship_gate": (
                "same nonidentity source-selected candidate in all folds; "
                "positive on every outer fold in both control replicates; "
                "replicate-A macro delta >= 0.002; observed delta above "
                "permuted-ratio p95 in at least two folds"
            ),
            "passed": passed,
            "decision": "ship" if passed else "reject",
        },
        "fixed_grid_outer_optimistic": fixed_grid,
        "leakage_audit": {
            "selection_uses_only_two_source_donor_effects": True,
            "outer_held_effects_accessed_by_evaluation_only_after_selection": True,
            "outer_held_controls_used": True,
            "all_300_vcc_target_coordinates_excluded": True,
            "fold_whiteners_exclude_held_donor": True,
            "operator_axis_hashes_verified": True,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    write_per_target(per_target_rows)
    RESULTS_JSON.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"\nwrote {RESULTS_JSON.relative_to(ROOT)} and "
        f"{PER_TARGET_CSV.relative_to(ROOT)}; "
        f"decision={payload['primary']['decision']} "
        f"peak RSS={peak_rss_gib():.3f} GiB",
        flush=True,
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=(
            "orion-workers",
            "controls",
            "evaluate",
            "all",
            "orion-worker",
        ),
        default="all",
    )
    parser.add_argument("--sample-per-replicate", type=int, default=512)
    parser.add_argument("--strata", type=int, default=32)
    parser.add_argument("--selection-targets", type=int, default=96)
    parser.add_argument("--chunk-rows", type=int, default=64)
    parser.add_argument("--operator-block", type=int, default=64)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--tie-tolerance", type=float, default=0.00025)
    parser.add_argument(
        "--max-rss-gib",
        type=float,
        default=0.70,
        help="internal stop below the external 0.75 GiB process-tree cap",
    )
    parser.add_argument("--orion-file", type=Path)
    parser.add_argument("--orion-output", type=Path)
    parser.add_argument("--orion-wanted", type=int)
    args = parser.parse_args()

    if args.stage == "orion-worker":
        if (
            args.orion_file is None
            or args.orion_output is None
            or args.orion_wanted is None
        ):
            parser.error(
                "--orion-file, --orion-output, and --orion-wanted are "
                "required for --stage orion-worker"
            )
        extract_orion_control_file(
            args.orion_file,
            args.orion_output,
            wanted=args.orion_wanted,
            max_rss_gib=args.max_rss_gib,
        )
        return
    if args.stage in ("orion-workers", "all"):
        build_orion_worker_cache(
            "HEK293T",
            sample_per_rep=args.sample_per_replicate,
            n_strata=args.strata,
            max_rss_gib=args.max_rss_gib,
        )
        if args.stage == "orion-workers":
            return
    if args.stage in ("controls", "all"):
        build_control_cache(
            sample_per_rep=args.sample_per_replicate,
            n_strata=args.strata,
            max_rss_gib=args.max_rss_gib,
        )
    if args.stage in ("evaluate", "all"):
        evaluate(
            sample_per_rep=args.sample_per_replicate,
            n_strata=args.strata,
            selection_targets=args.selection_targets,
            chunk_rows=args.chunk_rows,
            operator_block=args.operator_block,
            permutations=args.permutations,
            tie_tolerance=args.tie_tolerance,
            max_rss_gib=args.max_rss_gib,
        )


if __name__ == "__main__":
    main()
