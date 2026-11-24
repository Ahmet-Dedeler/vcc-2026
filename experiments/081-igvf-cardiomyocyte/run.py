#!/usr/bin/env python3
"""Leakage-safe IGVF cardiomyocyte sparse-direction transfer experiment.

The script is intentionally self-contained and low-memory. It reads large NPY
donor matrices and uncompressed NPY members inside the frozen NPZ whiteners
with ``os.pread`` so file-backed pages cannot accumulate in process RSS.
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import io
import json
import math
import os
import resource
import struct
import sys
import zipfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

import numpy as np
from numpy.lib import format as npformat
from scipy import sparse as sp
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DATA = ROOT / "data/external/igvf-cardiomyocyte"

MAIN_NAMES = ("k562_gwps", "orion_hct116", "orion_hek293t")
EFFECT_SCALE = 6.817
WHITEN_GAMMA = 0.65
REPRESENTATIONS = ("log2fc", "signed_sqrt", "sign")
SCALE_MODES = ("global_l2", "per_target_l2")
FRACTIONS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35)
SPLIT_SEED = "exp081:"
BOOTSTRAP_REPS = 20_000
NULL_REPS = 500
RNG_SEED = 8_102_026
MAX_RSS_GIB = 0.90

EXPECTED_IGVF = {
    "IGVFFI0830FXFI.csv.gz": {
        "size": 1_887_277,
        "md5": "39097b9eb95bb1e43e245480d242fea9",
        "sha256": "e4a68a15db8609b4f57b0e1387da0848f7447346fec2fbf19cc7c9ffdc41477f",
    },
    "IGVFFI8270UPKB.csv.gz": {
        "size": 348_103,
        "md5": "766b1fd06fe7d641104d2569a735a63f",
        "sha256": "dcd5eb865e1caa61557cf7c39ca39b42745f1ada4f8465bcfb53b91fc1c4925b",
    },
}

EXPECTED_WHITENERS = {
    "k562_gwps": "13d988fcc3b129d25d254524db73e8f8253143fa493ee64d9654f0f61bc3ac55",
    "orion_hct116": "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc",
    "orion_hek293t": "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1",
}


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_rss(stage: str, limit_gib: float = MAX_RSS_GIB) -> None:
    peak = peak_rss_gib()
    if peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.4f} GiB >= limit {limit_gib:.4f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing required column {column!r}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def assert_output_scope(path: Path) -> None:
    resolved = path.expanduser().resolve()
    if resolved != HERE and HERE not in resolved.parents:
        raise ValueError(f"refusing write outside experiment 081: {resolved}")


@dataclass(frozen=True)
class NpyHeader:
    shape: tuple[int, ...]
    fortran_order: bool
    dtype: np.dtype
    data_offset: int


def read_npy_header(path: Path) -> NpyHeader:
    with path.open("rb") as handle:
        version = npformat.read_magic(handle)
        readers = {
            (1, 0): npformat.read_array_header_1_0,
            (2, 0): npformat.read_array_header_2_0,
            (3, 0): npformat.read_array_header_2_0,
        }
        if version not in readers:
            raise ValueError(f"{path}: unsupported NPY version {version}")
        shape, fortran_order, dtype = readers[version](handle)
        return NpyHeader(
            shape=tuple(map(int, shape)),
            fortran_order=bool(fortran_order),
            dtype=np.dtype(dtype),
            data_offset=handle.tell(),
        )


class COrderNpyReader:
    """Bounded row reader for a C-order NPY file."""

    def __init__(self, path: Path):
        self.path = path.resolve()
        self.header = read_npy_header(self.path)
        if self.header.fortran_order or len(self.header.shape) != 2:
            raise ValueError(f"{path}: expected a 2-D C-order NPY")
        if self.header.dtype.hasobject:
            raise TypeError(f"{path}: object arrays are unsupported")
        self.row_bytes = self.header.shape[1] * self.header.dtype.itemsize

    def iter_blocks(self, rows_per_block: int) -> Iterable[tuple[int, np.ndarray]]:
        fd = os.open(self.path, os.O_RDONLY)
        try:
            for lo in range(0, self.header.shape[0], rows_per_block):
                count = min(rows_per_block, self.header.shape[0] - lo)
                raw = os.pread(
                    fd,
                    count * self.row_bytes,
                    self.header.data_offset + lo * self.row_bytes,
                )
                if len(raw) != count * self.row_bytes:
                    raise EOFError(f"{self.path}: truncated at row {lo}")
                yield lo, np.frombuffer(raw, dtype=self.header.dtype).reshape(
                    count, self.header.shape[1]
                )
        finally:
            os.close(fd)

    def read_selected(
        self, rows: Sequence[int], columns: np.ndarray
    ) -> np.ndarray:
        columns = np.asarray(columns, dtype=np.int64)
        result = np.empty((len(rows), len(columns)), dtype=np.float32)
        fd = os.open(self.path, os.O_RDONLY)
        try:
            for destination, row in enumerate(rows):
                if row < 0 or row >= self.header.shape[0]:
                    raise IndexError(f"{self.path}: row {row} out of range")
                raw = os.pread(
                    fd,
                    self.row_bytes,
                    self.header.data_offset + row * self.row_bytes,
                )
                if len(raw) != self.row_bytes:
                    raise EOFError(f"{self.path}: truncated at row {row}")
                whole = np.frombuffer(raw, dtype=self.header.dtype)
                result[destination] = whole[columns]
        finally:
            os.close(fd)
        return result


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: Mapping[str, int]
    reader: COrderNpyReader

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        path = path.resolve()
        genes = read_lines(path / "genes.txt", unique=False)
        perts = read_lines(path / "perts.txt")
        reader = COrderNpyReader(path / "delta.npy")
        if reader.header.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: NPY shape {reader.header.shape} != labels "
                f"{(len(perts), len(genes))}"
            )
        if reader.header.dtype != np.dtype("float32"):
            raise TypeError(f"{name}: expected float32, got {reader.header.dtype}")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={target: index for index, target in enumerate(perts)},
            reader=reader,
        )


@dataclass(frozen=True)
class BankStats:
    all_source_columns: np.ndarray
    unique_source_columns: np.ndarray
    unique_official_positions: np.ndarray
    unique_mean: np.ndarray
    raw_median_norm: float
    duplicate_symbols: tuple[str, ...]


def map_bank(bank: Bank, official: Sequence[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    official_pos = {gene: index for index, gene in enumerate(official)}
    all_source: list[int] = []
    last_source: dict[int, int] = {}
    duplicates: set[str] = set()
    seen: set[int] = set()
    for source_column, gene in enumerate(bank.genes):
        destination = official_pos.get(gene)
        if destination is None:
            continue
        all_source.append(source_column)
        if destination in seen:
            duplicates.add(gene)
        seen.add(destination)
        last_source[destination] = source_column
    if not all_source:
        raise ValueError(f"{bank.name}: no genes map to official axis")
    unique_official = np.asarray(sorted(last_source), dtype=np.int64)
    unique_source = np.asarray(
        [last_source[position] for position in unique_official], dtype=np.int64
    )
    return (
        np.asarray(all_source, dtype=np.int64),
        unique_source,
        unique_official,
        tuple(sorted(duplicates)),
    )


def compute_bank_stats(
    bank: Bank,
    official: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> BankStats:
    all_source, unique_source, unique_official, duplicates = map_bank(bank, official)
    total = np.zeros(len(all_source), dtype=np.float64)
    for _, block in bank.reader.iter_blocks(chunk_rows):
        selected = np.asarray(block[:, all_source], dtype=np.float64)
        total += selected.sum(axis=0, dtype=np.float64)
        enforce_rss(f"{bank.name} mean", max_rss_gib)
    mean_all = total / bank.reader.header.shape[0]
    source_mean = {
        int(source): float(mean_all[index])
        for index, source in enumerate(all_source.tolist())
    }
    unique_mean = np.asarray(
        [source_mean[int(source)] for source in unique_source], dtype=np.float32
    )
    norms = np.empty(bank.reader.header.shape[0], dtype=np.float32)
    for lo, block in bank.reader.iter_blocks(chunk_rows):
        selected = np.asarray(block[:, all_source], dtype=np.float64)
        selected -= mean_all
        norms[lo : lo + len(block)] = np.linalg.norm(selected, axis=1)
        enforce_rss(f"{bank.name} scale", max_rss_gib)
    raw_median = float(np.median(norms)) or 1.0
    return BankStats(
        all_source_columns=all_source,
        unique_source_columns=unique_source,
        unique_official_positions=unique_official,
        unique_mean=unique_mean,
        raw_median_norm=raw_median,
        duplicate_symbols=duplicates,
    )


def prepare_bank_rows(
    bank: Bank,
    stats: BankStats,
    targets: Sequence[str],
    selected_official_positions: np.ndarray,
    excluded_official_positions: set[int],
    *,
    scaled: bool,
) -> np.ndarray:
    official_to_selected = {
        int(position): index
        for index, position in enumerate(selected_official_positions.tolist())
    }
    source_columns: list[int] = []
    destination_columns: list[int] = []
    means: list[float] = []
    for source, position, mean in zip(
        stats.unique_source_columns,
        stats.unique_official_positions,
        stats.unique_mean,
        strict=True,
    ):
        destination = official_to_selected.get(int(position))
        if destination is None or int(position) in excluded_official_positions:
            continue
        source_columns.append(int(source))
        destination_columns.append(destination)
        means.append(float(mean))
    rows = [bank.pert_pos[target] for target in targets]
    values = bank.reader.read_selected(rows, np.asarray(source_columns, dtype=np.int64))
    values -= np.asarray(means, dtype=np.float32)
    if scaled:
        values *= np.float32(EFFECT_SCALE / stats.raw_median_norm)
    output = np.zeros(
        (len(targets), len(selected_official_positions)), dtype=np.float32
    )
    output[:, np.asarray(destination_columns, dtype=np.int64)] = values
    return output


def npz_member_header(path: Path, member: str) -> tuple[NpyHeader, int]:
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(member)
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"{path}:{member}: member must be uncompressed")
        local_header_offset = info.header_offset
    with path.open("rb") as handle:
        handle.seek(local_header_offset)
        header = handle.read(30)
        if len(header) != 30:
            raise EOFError(f"{path}: truncated ZIP local header")
        fields = struct.unpack("<IHHHHHIIIHH", header)
        if fields[0] != 0x04034B50:
            raise ValueError(f"{path}: invalid ZIP local header signature")
        name_len, extra_len = fields[-2], fields[-1]
        npy_start = local_header_offset + 30 + name_len + extra_len
        handle.seek(npy_start)
        version = npformat.read_magic(handle)
        readers = {
            (1, 0): npformat.read_array_header_1_0,
            (2, 0): npformat.read_array_header_2_0,
            (3, 0): npformat.read_array_header_2_0,
        }
        shape, fortran_order, dtype = readers[version](handle)
        data_offset = handle.tell()
    return (
        NpyHeader(
            shape=tuple(map(int, shape)),
            fortran_order=bool(fortran_order),
            dtype=np.dtype(dtype),
            data_offset=data_offset,
        ),
        npy_start,
    )


class FrozenOperator:
    def __init__(
        self,
        path: Path,
        metadata_path: Path,
        expected_sha256: str,
        expected_donors: Sequence[str],
        expected_wcols: np.ndarray,
    ):
        self.path = path.resolve()
        self.metadata_path = metadata_path.resolve()
        actual_sha = sha256_file(self.path)
        if actual_sha != expected_sha256:
            raise ValueError(
                f"{self.path}: SHA-256 {actual_sha} != {expected_sha256}"
            )
        self.metadata = json.loads(self.metadata_path.read_text())
        if self.metadata.get("donors") != list(expected_donors):
            raise ValueError(
                f"{self.path}: donors {self.metadata.get('donors')} != "
                f"{list(expected_donors)}"
            )
        if float(self.metadata.get("gamma", -1)) != WHITEN_GAMMA:
            raise ValueError(f"{self.path}: unexpected gamma")
        expected_hash = hashlib.sha256(expected_wcols.tobytes()).hexdigest()
        if self.metadata.get("wcols_sha256") != expected_hash:
            raise ValueError(f"{self.path}: operator-axis hash mismatch")
        if int(self.metadata.get("n_cols", -1)) != len(expected_wcols):
            raise ValueError(f"{self.path}: operator-axis length mismatch")
        with zipfile.ZipFile(self.path) as archive:
            evals_bytes = archive.read("evals.npy")
        self.evals = np.asarray(np.load(io.BytesIO(evals_bytes)), dtype=np.float32)
        header, _ = npz_member_header(self.path, "evecs.npy")
        if header.fortran_order or header.dtype != np.dtype("float32"):
            raise ValueError(f"{self.path}: unexpected evec layout")
        if header.shape != (len(expected_wcols), len(expected_wcols)):
            raise ValueError(f"{self.path}: evec shape {header.shape} is invalid")
        self.header = header
        self.n = len(expected_wcols)
        self.row_bytes = self.n * self.header.dtype.itemsize

    def read_rows(self, lo: int, hi: int) -> np.ndarray:
        fd = os.open(self.path, os.O_RDONLY)
        try:
            raw = os.pread(
                fd,
                (hi - lo) * self.row_bytes,
                self.header.data_offset + lo * self.row_bytes,
            )
        finally:
            os.close(fd)
        if len(raw) != (hi - lo) * self.row_bytes:
            raise EOFError(f"{self.path}: truncated eigenvector member")
        return np.frombuffer(raw, dtype=np.float32).reshape(hi - lo, self.n)

    def apply(
        self,
        matrix: np.ndarray,
        *,
        block_rows: int,
        max_rss_gib: float,
    ) -> tuple[np.ndarray, dict[str, float]]:
        matrix = np.ascontiguousarray(matrix, dtype=np.float32)
        if matrix.ndim != 2 or matrix.shape[1] != self.n:
            raise ValueError(
                f"operator input shape {matrix.shape} != (*, {self.n})"
            )
        before = float(
            np.median(np.linalg.norm(matrix.astype(np.float64), axis=1))
        ) or 1.0
        stage = np.zeros_like(matrix)
        temporary = np.empty_like(matrix)
        for lo in range(0, self.n, block_rows):
            hi = min(lo + block_rows, self.n)
            evec_rows = self.read_rows(lo, hi)
            if not np.isfinite(evec_rows).all():
                raise ValueError(f"{self.path}: non-finite eigenvector block")
            with np.errstate(all="ignore"):
                np.matmul(matrix[:, lo:hi], evec_rows, out=temporary)
            if not np.isfinite(temporary).all():
                raise FloatingPointError(
                    f"{self.path}: non-finite first-pass product {lo}:{hi}"
                )
            stage += temporary
            enforce_rss(f"operator first pass {lo}:{hi}", max_rss_gib)
        weights = np.power(
            np.maximum(self.evals, float(self.evals.max()) * 1e-12),
            -WHITEN_GAMMA / 2.0,
            dtype=np.float32,
        )
        stage *= weights
        output = np.zeros_like(matrix)
        for lo in range(0, self.n, block_rows):
            hi = min(lo + block_rows, self.n)
            evec_rows = self.read_rows(lo, hi)
            with np.errstate(all="ignore"):
                output[:, lo:hi] = stage @ evec_rows.T
            if not np.isfinite(output[:, lo:hi]).all():
                raise FloatingPointError(
                    f"{self.path}: non-finite second-pass product {lo}:{hi}"
                )
            enforce_rss(f"operator second pass {lo}:{hi}", max_rss_gib)
        after = float(
            np.median(np.linalg.norm(output.astype(np.float64), axis=1))
        ) or 1.0
        scale = before / after
        output *= np.float32(scale)
        return output, {
            "before_operator_median_norm": before,
            "after_operator_median_norm": after,
            "applied_scale": scale,
        }


@dataclass(frozen=True)
class IGVFData:
    effects: Mapping[str, Mapping[int, float]]
    num_cells: Mapping[str, int]
    row_counts: Mapping[str, int]
    de_targets: frozenset[str]
    guide_targets: frozenset[str]
    de_rows: int
    guide_rows: int
    response_genes: frozenset[str]
    duplicate_pairs: int


def verify_igvf_files(data_dir: Path) -> dict[str, dict[str, object]]:
    observed: dict[str, dict[str, object]] = {}
    for name, expected in EXPECTED_IGVF.items():
        path = data_dir / name
        if not path.is_file():
            raise FileNotFoundError(path)
        size = path.stat().st_size
        md5 = md5_file(path)
        sha256 = sha256_file(path)
        if size != expected["size"]:
            raise ValueError(f"{path}: size {size} != {expected['size']}")
        if md5 != expected["md5"]:
            raise ValueError(f"{path}: MD5 {md5} != {expected['md5']}")
        if sha256 != expected["sha256"]:
            raise ValueError(f"{path}: SHA-256 {sha256} != {expected['sha256']}")
        observed[name] = {"size": size, "md5": md5, "sha256": sha256}
    return observed


def load_igvf(
    de_path: Path,
    guide_path: Path,
    official: Sequence[str],
) -> IGVFData:
    official_pos = {gene: index for index, gene in enumerate(official)}
    weighted_sum: dict[tuple[str, int], float] = defaultdict(float)
    weight_sum: dict[tuple[str, int], float] = defaultdict(float)
    target_cells: dict[str, int] = {}
    target_rows: dict[str, int] = defaultdict(int)
    de_targets: set[str] = set()
    response_genes: set[str] = set()
    de_rows = 0
    with gzip.open(de_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "gene_names",
            "intended_target_name",
            "num_cell",
            "log2fc",
            "pval-empirical",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"{de_path}: unexpected schema {reader.fieldnames}")
        for row in reader:
            de_rows += 1
            target = row["intended_target_name"].strip()
            gene = row["gene_names"].strip()
            if not target or not gene:
                continue
            value = float(row["log2fc"])
            cells = int(float(row["num_cell"]))
            if not math.isfinite(value) or cells <= 0:
                raise ValueError(f"{de_path}: invalid row {de_rows}")
            de_targets.add(target)
            response_genes.add(gene)
            target_cells[target] = max(cells, target_cells.get(target, 0))
            target_rows[target] += 1
            destination = official_pos.get(gene)
            if destination is None:
                continue
            key = (target, destination)
            weighted_sum[key] += value * cells
            weight_sum[key] += cells
    effects: dict[str, dict[int, float]] = defaultdict(dict)
    duplicate_pairs = 0
    for (target, destination), total in weighted_sum.items():
        effects[target][destination] = total / weight_sum[(target, destination)]
    pair_counts: dict[tuple[str, str], int] = defaultdict(int)
    with gzip.open(de_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            target = row["intended_target_name"].strip()
            gene = row["gene_names"].strip()
            if target and gene:
                pair_counts[(target, gene)] += 1
    duplicate_pairs = sum(count > 1 for count in pair_counts.values())

    guide_targets: set[str] = set()
    guide_rows = 0
    with gzip.open(guide_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"guide_id", "targeting", "type", "gene_name"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"{guide_path}: unexpected schema {reader.fieldnames}")
        for row in reader:
            guide_rows += 1
            gene = row["gene_name"].strip()
            if gene:
                guide_targets.add(gene)

    return IGVFData(
        effects={target: dict(values) for target, values in effects.items()},
        num_cells=target_cells,
        row_counts=dict(target_rows),
        de_targets=frozenset(de_targets),
        guide_targets=frozenset(guide_targets),
        de_rows=de_rows,
        guide_rows=guide_rows,
        response_genes=frozenset(response_genes),
        duplicate_pairs=duplicate_pairs,
    )


def split_targets(targets: Sequence[str]) -> dict[str, tuple[str, ...]]:
    parts: dict[str, list[str]] = {"fit": [], "tune": [], "test": []}
    for target in sorted(targets):
        digest = hashlib.sha256((SPLIT_SEED + target).encode()).digest()
        bucket = int.from_bytes(digest[:8], "big") % 10
        part = "fit" if bucket < 6 else "tune" if bucket < 8 else "test"
        parts[part].append(target)
    return {name: tuple(values) for name, values in parts.items()}


def operator_columns(
    source_banks: Sequence[Bank],
    official: Sequence[str],
    panel: set[str],
) -> np.ndarray:
    source_union = set().union(*(set(bank.genes) for bank in source_banks))
    return np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in source_union and gene not in panel
        ],
        dtype=np.int64,
    )


def igvf_matrix(
    igvf: IGVFData,
    targets: Sequence[str],
    selected_official_positions: np.ndarray,
    excluded_official_positions: set[int],
    representation: str,
    clip_abs: float,
) -> sp.csr_matrix:
    if representation not in REPRESENTATIONS:
        raise ValueError(f"unknown representation {representation!r}")
    selected = {
        int(position): index
        for index, position in enumerate(selected_official_positions.tolist())
    }
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for row_index, target in enumerate(targets):
        for official_position, raw_value in igvf.effects.get(target, {}).items():
            column = selected.get(int(official_position))
            if column is None or official_position in excluded_official_positions:
                continue
            value = float(np.clip(raw_value, -clip_abs, clip_abs))
            if representation == "signed_sqrt":
                value = math.copysign(math.sqrt(abs(value)), value)
            elif representation == "sign":
                value = math.copysign(1.0, value)
            if value == 0.0:
                continue
            rows.append(row_index)
            cols.append(column)
            values.append(value)
    matrix = sp.csr_matrix(
        (
            np.asarray(values, dtype=np.float32),
            (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32)),
        ),
        shape=(len(targets), len(selected_official_positions)),
        dtype=np.float32,
    )
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix


@dataclass(frozen=True)
class ScaleParams:
    mode: str
    global_factor: float
    row_factor_low: float
    row_factor_high: float


def fit_scale(
    raw: sp.csr_matrix,
    baseline_norm: np.ndarray,
    mode: str,
) -> ScaleParams:
    residual_norm = np.sqrt(np.asarray(raw.multiply(raw).sum(axis=1)).ravel())
    baseline_norm = np.asarray(baseline_norm, dtype=np.float64)
    if baseline_norm.shape != (raw.shape[0],):
        raise ValueError(
            f"baseline norm shape {baseline_norm.shape} != {(raw.shape[0],)}"
        )
    valid = residual_norm > 0
    if int(valid.sum()) < 100:
        raise ValueError("too few nonzero IGVF fit directions for scale fitting")
    ratios = baseline_norm[valid] / residual_norm[valid]
    global_factor = float(np.median(baseline_norm[valid])) / (
        float(np.median(residual_norm[valid])) or 1.0
    )
    low, high = np.percentile(ratios, [5.0, 95.0])
    if mode not in SCALE_MODES:
        raise ValueError(f"unknown scale mode {mode!r}")
    return ScaleParams(
        mode=mode,
        global_factor=global_factor,
        row_factor_low=float(low),
        row_factor_high=float(high),
    )


def scale_sparse(
    raw: sp.csr_matrix,
    baseline_pre_norm: np.ndarray,
    params: ScaleParams,
) -> sp.csr_matrix:
    residual_norm = np.sqrt(np.asarray(raw.multiply(raw).sum(axis=1)).ravel())
    if params.mode == "global_l2":
        factors = np.full(raw.shape[0], params.global_factor, dtype=np.float32)
        factors[residual_norm == 0] = 0.0
    elif params.mode == "per_target_l2":
        factors = np.zeros(raw.shape[0], dtype=np.float32)
        valid = residual_norm > 0
        factors[valid] = np.clip(
            baseline_pre_norm[valid] / residual_norm[valid],
            params.row_factor_low,
            params.row_factor_high,
        )
    else:
        raise ValueError(params.mode)
    return raw.multiply(factors[:, None]).tocsr()


def pds_from_similarities(similarities: np.ndarray) -> np.ndarray:
    own = np.diag(similarities)[:, None]
    n = similarities.shape[0]
    if n < 2:
        raise ValueError("PDS requires at least two targets")
    ranks = (similarities > own).sum(axis=1)
    ranks = ranks + 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


@dataclass
class FoldPack:
    held: str
    targets: tuple[str, ...]
    wcols: np.ndarray
    retrieval_positions: np.ndarray
    retrieval_official_positions: np.ndarray
    baseline_pre_norm: np.ndarray
    baseline_retrieval: np.ndarray
    truth_unit: np.ndarray
    base_dot: np.ndarray
    base_norm_sq: np.ndarray
    baseline_pds: np.ndarray
    source_similarity: np.ndarray
    raw_by_representation: Mapping[str, sp.csr_matrix]
    whitening: Mapping[str, float]
    source_names: tuple[str, ...]
    whitener_sha256: str


def evaluate_residual(
    pack: FoldPack,
    residual: sp.csr_matrix,
    fraction: float,
) -> np.ndarray:
    residual = residual[:, pack.retrieval_positions].tocsr()
    residual_dot = np.asarray(residual @ pack.truth_unit.T, dtype=np.float32)
    cross = np.asarray(
        residual.multiply(pack.baseline_retrieval).sum(axis=1)
    ).ravel()
    residual_norm_sq = np.asarray(residual.multiply(residual).sum(axis=1)).ravel()
    denominator = np.sqrt(
        np.maximum(
            pack.base_norm_sq
            + 2.0 * fraction * cross
            + (fraction * fraction) * residual_norm_sq,
            1e-20,
        )
    )
    similarities = (
        pack.base_dot + np.float32(fraction) * residual_dot
    ) / denominator[:, None]
    return pds_from_similarities(similarities)


def build_fold_pack(
    held: str,
    targets: Sequence[str],
    banks: Mapping[str, Bank],
    stats: Mapping[str, BankStats],
    official: Sequence[str],
    panel: set[str],
    excluded_positions: set[int],
    igvf: IGVFData,
    clip_abs: float,
    whitener_dir: Path,
    *,
    operator_block_rows: int,
    max_rss_gib: float,
) -> FoldPack:
    source_names = tuple(name for name in MAIN_NAMES if name != held)
    source_banks = tuple(banks[name] for name in source_names)
    wcols = operator_columns(source_banks, official, panel)
    metadata_path = whitener_dir / f"held_{held}.json"
    whitener_path = whitener_dir / f"held_{held}.npz"
    operator = FrozenOperator(
        whitener_path,
        metadata_path,
        EXPECTED_WHITENERS[held],
        source_names,
        wcols,
    )
    baseline_pre = np.zeros((len(targets), len(wcols)), dtype=np.float32)
    for name in source_names:
        source_matrix = prepare_bank_rows(
            banks[name],
            stats[name],
            targets,
            wcols,
            excluded_positions,
            scaled=True,
        )
        baseline_pre += source_matrix
        del source_matrix
    baseline_pre /= np.float32(len(source_names))
    baseline_pre_norm = np.linalg.norm(
        baseline_pre.astype(np.float64), axis=1
    ).astype(np.float32)
    baseline_post, whitening = operator.apply(
        baseline_pre,
        block_rows=operator_block_rows,
        max_rss_gib=max_rss_gib,
    )
    del baseline_pre

    held_gene_set = set(banks[held].genes)
    retrieval_positions = np.asarray(
        [
            index
            for index, official_position in enumerate(wcols.tolist())
            if official[official_position] in held_gene_set
            and int(official_position) not in excluded_positions
        ],
        dtype=np.int64,
    )
    retrieval_official = wcols[retrieval_positions]
    truth = prepare_bank_rows(
        banks[held],
        stats[held],
        targets,
        retrieval_official,
        excluded_positions,
        scaled=False,
    )
    baseline_retrieval = np.ascontiguousarray(
        baseline_post[:, retrieval_positions], dtype=np.float32
    )
    del baseline_post
    truth_norm = np.maximum(
        np.linalg.norm(truth.astype(np.float64), axis=1), 1e-12
    )
    truth_unit = np.asarray(truth / truth_norm[:, None], dtype=np.float32)
    del truth
    base_norm_sq = np.asarray(
        np.sum(
            baseline_retrieval.astype(np.float64) ** 2,
            axis=1,
        ),
        dtype=np.float64,
    )
    with np.errstate(all="ignore"):
        base_dot = np.asarray(baseline_retrieval @ truth_unit.T, dtype=np.float32)
    if not np.isfinite(base_dot).all():
        raise FloatingPointError(f"{held}: non-finite baseline/truth similarities")
    similarities = base_dot / np.sqrt(np.maximum(base_norm_sq, 1e-20))[:, None]
    baseline_pds = pds_from_similarities(similarities)
    base_unit = baseline_retrieval / np.sqrt(
        np.maximum(base_norm_sq, 1e-20)
    )[:, None]
    with np.errstate(all="ignore"):
        source_similarity = np.asarray(base_unit @ base_unit.T, dtype=np.float32)
    if not np.isfinite(source_similarity).all():
        raise FloatingPointError(f"{held}: non-finite source similarities")
    raw_by_representation = {
        representation: igvf_matrix(
            igvf,
            targets,
            wcols,
            excluded_positions,
            representation,
            clip_abs,
        )
        for representation in REPRESENTATIONS
    }
    enforce_rss(f"{held} fold pack", max_rss_gib)
    return FoldPack(
        held=held,
        targets=tuple(targets),
        wcols=wcols,
        retrieval_positions=retrieval_positions,
        retrieval_official_positions=retrieval_official,
        baseline_pre_norm=baseline_pre_norm,
        baseline_retrieval=baseline_retrieval,
        truth_unit=truth_unit,
        base_dot=base_dot,
        base_norm_sq=base_norm_sq,
        baseline_pds=baseline_pds,
        source_similarity=source_similarity,
        raw_by_representation=raw_by_representation,
        whitening=whitening,
        source_names=source_names,
        whitener_sha256=EXPECTED_WHITENERS[held],
    )


def fit_scale_parameters(
    targets: Sequence[str],
    banks: Mapping[str, Bank],
    stats: Mapping[str, BankStats],
    official: Sequence[str],
    panel: set[str],
    proxy_excluded_positions: set[int],
    igvf: IGVFData,
    clip_abs: float,
    *,
    target_chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, dict[str, ScaleParams]]:
    result: dict[str, dict[str, ScaleParams]] = {}
    for held in MAIN_NAMES:
        source_names = tuple(name for name in MAIN_NAMES if name != held)
        wcols = operator_columns(tuple(banks[name] for name in source_names), official, panel)
        baseline_norm = np.empty(len(targets), dtype=np.float32)
        for lo in range(0, len(targets), target_chunk_rows):
            hi = min(lo + target_chunk_rows, len(targets))
            target_chunk = targets[lo:hi]
            baseline_pre = np.zeros(
                (len(target_chunk), len(wcols)), dtype=np.float32
            )
            for source_name in source_names:
                source_matrix = prepare_bank_rows(
                    banks[source_name],
                    stats[source_name],
                    target_chunk,
                    wcols,
                    proxy_excluded_positions,
                    scaled=True,
                )
                baseline_pre += source_matrix
                del source_matrix
            baseline_pre /= np.float32(len(source_names))
            baseline_norm[lo:hi] = np.linalg.norm(
                baseline_pre.astype(np.float64), axis=1
            )
            del baseline_pre
            enforce_rss(f"{held} scale-fit chunk", max_rss_gib)
        by_representation: dict[str, ScaleParams] = {}
        for representation in REPRESENTATIONS:
            raw = igvf_matrix(
                igvf,
                targets,
                wcols,
                proxy_excluded_positions,
                representation,
                clip_abs,
            )
            for mode in SCALE_MODES:
                by_representation[f"{representation}:{mode}"] = fit_scale(
                    raw, baseline_norm, mode
                )
        result[held] = by_representation
        del baseline_norm
        gc.collect()
        enforce_rss(f"{held} scale fit", max_rss_gib)
    return result


def score_configuration(
    packs: Mapping[str, FoldPack],
    scale_params: Mapping[str, Mapping[str, ScaleParams]],
    representation: str,
    scale_mode: str,
    fraction: float,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    per_fold: dict[str, np.ndarray] = {}
    scores: dict[str, float] = {}
    key = f"{representation}:{scale_mode}"
    for held, pack in packs.items():
        raw = pack.raw_by_representation[representation]
        residual = scale_sparse(
            raw, pack.baseline_pre_norm, scale_params[held][key]
        )
        pds = (
            pack.baseline_pds.copy()
            if fraction == 0.0
            else evaluate_residual(pack, residual, fraction)
        )
        per_fold[held] = pds
        scores[held] = float(np.mean(pds))
    return per_fold, scores


def tune_grid(
    packs: Mapping[str, FoldPack],
    scale_params: Mapping[str, Mapping[str, ScaleParams]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    baseline = {
        held: float(np.mean(pack.baseline_pds)) for held, pack in packs.items()
    }
    configs = [("log2fc", "global_l2", 0.0)]
    configs.extend(
        (representation, scale_mode, fraction)
        for representation in REPRESENTATIONS
        for scale_mode in SCALE_MODES
        for fraction in FRACTIONS
        if fraction > 0
    )
    for representation, scale_mode, fraction in configs:
        _, scores = score_configuration(
            packs, scale_params, representation, scale_mode, fraction
        )
        deltas = {held: scores[held] - baseline[held] for held in MAIN_NAMES}
        rows.append(
            {
                "representation": representation,
                "scale_mode": scale_mode,
                "fraction": fraction,
                **{f"pds_{held}": scores[held] for held in MAIN_NAMES},
                **{f"delta_{held}": deltas[held] for held in MAIN_NAMES},
                "worst_fold_delta": min(deltas.values()),
                "macro_delta": float(np.mean(list(deltas.values()))),
            }
        )
    selected = max(
        rows,
        key=lambda row: (
            float(row["worst_fold_delta"]),
            float(row["macro_delta"]),
            -float(row["fraction"]),
            row["representation"] == "log2fc",
            row["scale_mode"] == "global_l2",
        ),
    )
    return {
        "representation": selected["representation"],
        "scale_mode": selected["scale_mode"],
        "fraction": selected["fraction"],
        "selection_rule": (
            "maximize worst held-fold delta, then macro delta, then smaller "
            "fraction, then log2fc, then global_l2; identity is mandatory"
        ),
        "selected_tune_metrics": selected,
    }, rows


def candidate_grid() -> list[tuple[str, str, float]]:
    configs = [("log2fc", "global_l2", 0.0)]
    configs.extend(
        (representation, scale_mode, fraction)
        for representation in REPRESENTATIONS
        for scale_mode in SCALE_MODES
        for fraction in FRACTIONS
        if fraction > 0
    )
    return configs


def tune_grid_streaming(
    targets: Sequence[str],
    banks: Mapping[str, Bank],
    stats: Mapping[str, BankStats],
    official: Sequence[str],
    panel: set[str],
    excluded_positions: set[int],
    igvf: IGVFData,
    clip_abs: float,
    whitener_dir: Path,
    scale_params: Mapping[str, Mapping[str, ScaleParams]],
    *,
    operator_block_rows: int,
    max_rss_gib: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    configs = candidate_grid()
    values: dict[tuple[str, str, float], dict[str, float]] = {
        config: {} for config in configs
    }
    baselines: dict[str, float] = {}
    for held in MAIN_NAMES:
        pack = build_fold_pack(
            held,
            targets,
            banks,
            stats,
            official,
            panel,
            excluded_positions,
            igvf,
            clip_abs,
            whitener_dir,
            operator_block_rows=operator_block_rows,
            max_rss_gib=max_rss_gib,
        )
        baselines[held] = float(np.mean(pack.baseline_pds))
        for representation, scale_mode, fraction in configs:
            key = f"{representation}:{scale_mode}"
            residual = scale_sparse(
                pack.raw_by_representation[representation],
                pack.baseline_pre_norm,
                scale_params[held][key],
            )
            pds = (
                pack.baseline_pds
                if fraction == 0.0
                else evaluate_residual(pack, residual, fraction)
            )
            values[(representation, scale_mode, fraction)][held] = float(
                np.mean(pds)
            )
        del pack
        gc.collect()
        enforce_rss(f"{held} tune fold released", max_rss_gib)

    rows: list[dict[str, object]] = []
    for representation, scale_mode, fraction in configs:
        scores = values[(representation, scale_mode, fraction)]
        deltas = {held: scores[held] - baselines[held] for held in MAIN_NAMES}
        rows.append(
            {
                "representation": representation,
                "scale_mode": scale_mode,
                "fraction": fraction,
                **{f"pds_{held}": scores[held] for held in MAIN_NAMES},
                **{f"delta_{held}": deltas[held] for held in MAIN_NAMES},
                "worst_fold_delta": min(deltas.values()),
                "macro_delta": float(np.mean(list(deltas.values()))),
            }
        )
    selected = max(
        rows,
        key=lambda row: (
            float(row["worst_fold_delta"]),
            float(row["macro_delta"]),
            -float(row["fraction"]),
            row["representation"] == "log2fc",
            row["scale_mode"] == "global_l2",
        ),
    )
    return {
        "representation": selected["representation"],
        "scale_mode": selected["scale_mode"],
        "fraction": selected["fraction"],
        "selection_rule": (
            "maximize worst held-fold delta, then macro delta, then smaller "
            "fraction, then log2fc, then global_l2; identity is mandatory"
        ),
        "selected_tune_metrics": selected,
    }, rows


def target_strata(
    igvf: IGVFData,
    targets: Sequence[str],
    excluded_positions: set[int],
    clip_abs: float,
) -> np.ndarray:
    nnz = np.asarray(
        [
            sum(
                position not in excluded_positions
                for position in igvf.effects.get(target, {})
            )
            for target in targets
        ],
        dtype=np.float64,
    )
    norms = np.asarray(
        [
            math.sqrt(
                sum(
                    float(np.clip(value, -clip_abs, clip_abs)) ** 2
                    for position, value in igvf.effects.get(target, {}).items()
                    if position not in excluded_positions
                )
            )
            for target in targets
        ],
        dtype=np.float64,
    )
    cells = np.asarray(
        [igvf.num_cells.get(target, 0) for target in targets], dtype=np.float64
    )

    def bins(values: np.ndarray) -> np.ndarray:
        positive = values > 0
        output = np.zeros(len(values), dtype=np.int16)
        if positive.sum() < 4:
            output[positive] = 1
            return output
        edges = np.unique(np.quantile(values[positive], [0.25, 0.5, 0.75]))
        output[positive] = np.digitize(values[positive], edges, right=True) + 1
        return output

    nnz_bin, norm_bin, cell_bin = bins(nnz), bins(norms), bins(cells)
    return (
        nnz_bin.astype(np.int64) * 100
        + norm_bin.astype(np.int64) * 10
        + cell_bin.astype(np.int64)
    )


def matched_permutation(strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    permutation = np.arange(len(strata), dtype=np.int64)
    for value in np.unique(strata):
        positions = np.flatnonzero(strata == value)
        if len(positions) > 1:
            permutation[positions] = rng.permutation(positions)
    return permutation


def sign_randomize(raw: sp.csr_matrix, rng: np.random.Generator) -> sp.csr_matrix:
    output = raw.copy()
    signs = rng.integers(0, 2, size=len(output.data), dtype=np.int8)
    output.data *= np.where(signs == 0, -1.0, 1.0).astype(np.float32)
    return output


def derive_target_clusters(
    source_similarities: Sequence[np.ndarray],
    *,
    mean_cluster_size: int = 8,
) -> np.ndarray:
    similarities = np.mean(
        np.stack(source_similarities, axis=0),
        axis=0,
    )
    similarities = np.clip((similarities + similarities.T) * 0.5, -1.0, 1.0)
    np.fill_diagonal(similarities, 1.0)
    distance = np.clip(1.0 - similarities, 0.0, 2.0)
    condensed = squareform(distance, checks=False)
    linkage = hierarchy.linkage(condensed, method="average")
    n_clusters = max(10, math.ceil(len(similarities) / mean_cluster_size))
    return hierarchy.fcluster(linkage, t=n_clusters, criterion="maxclust") - 1


def cluster_bootstrap(
    per_fold_delta: Mapping[str, np.ndarray],
    clusters: np.ndarray,
    *,
    reps: int,
    seed: int,
) -> dict[str, object]:
    matrix = np.stack([per_fold_delta[held] for held in MAIN_NAMES], axis=1)
    target_delta = matrix.mean(axis=1)
    unique = np.unique(clusters)
    members = [np.flatnonzero(clusters == cluster) for cluster in unique]
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for rep in range(reps):
        selected = rng.integers(0, len(members), size=len(members))
        indices = np.concatenate([members[index] for index in selected])
        draws[rep] = float(np.mean(target_delta[indices]))
    return {
        "method": (
            "paired target-cluster bootstrap; clusters from candidate-independent "
            "incumbent prediction cosine geometry"
        ),
        "reps": reps,
        "seed": seed,
        "clusters": int(len(unique)),
        "max_cluster_size": int(max(map(len, members))),
        "estimate": float(target_delta.mean()),
        "ci95": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "positive_probability": float(np.mean(draws > 0)),
    }


def fold_null_draws(
    pack: FoldPack,
    params: ScaleParams,
    representation: str,
    fraction: float,
    label_permutations: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    sign_draws = np.empty(len(label_permutations), dtype=np.float64)
    label_draws = np.empty(len(label_permutations), dtype=np.float64)
    raw = pack.raw_by_representation[representation]
    for rep, permutation in enumerate(label_permutations):
        signed = sign_randomize(raw, rng)
        signed_scaled = scale_sparse(signed, pack.baseline_pre_norm, params)
        signed_pds = evaluate_residual(pack, signed_scaled, fraction)
        sign_draws[rep] = float(np.mean(signed_pds - pack.baseline_pds))

        permuted = raw[permutation]
        label_scaled = scale_sparse(permuted, pack.baseline_pre_norm, params)
        label_pds = evaluate_residual(pack, label_scaled, fraction)
        label_draws[rep] = float(np.mean(label_pds - pack.baseline_pds))
    return sign_draws, label_draws


def evaluate_stage(
    name: str,
    targets: Sequence[str],
    banks: Mapping[str, Bank],
    stats: Mapping[str, BankStats],
    official: Sequence[str],
    panel: set[str],
    excluded_positions: set[int],
    igvf: IGVFData,
    clip_abs: float,
    whitener_dir: Path,
    scale_params: Mapping[str, Mapping[str, ScaleParams]],
    selection: Mapping[str, object],
    *,
    operator_block_rows: int,
    max_rss_gib: float,
    bootstrap_reps: int,
    null_reps: int,
    seed: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    representation = str(selection["representation"])
    scale_mode = str(selection["scale_mode"])
    fraction = float(selection["fraction"])
    key = f"{representation}:{scale_mode}"
    strata = target_strata(igvf, targets, excluded_positions, clip_abs)
    permutation_rng = np.random.default_rng(seed + 1)
    label_permutations = np.stack(
        [matched_permutation(strata, permutation_rng) for _ in range(null_reps)],
        axis=0,
    ).astype(np.int32)

    baseline_by_fold: dict[str, np.ndarray] = {}
    candidate_by_fold: dict[str, np.ndarray] = {}
    per_fold_delta: dict[str, np.ndarray] = {}
    baseline_scores: dict[str, float] = {}
    scores: dict[str, float] = {}
    deltas: dict[str, float] = {}
    source_similarities: list[np.ndarray] = []
    sign_null_by_fold: list[np.ndarray] = []
    label_null_by_fold: list[np.ndarray] = []
    fold_details: dict[str, object] = {}

    for fold_index, held in enumerate(MAIN_NAMES):
        pack = build_fold_pack(
            held,
            targets,
            banks,
            stats,
            official,
            panel,
            excluded_positions,
            igvf,
            clip_abs,
            whitener_dir,
            operator_block_rows=operator_block_rows,
            max_rss_gib=max_rss_gib,
        )
        residual = scale_sparse(
            pack.raw_by_representation[representation],
            pack.baseline_pre_norm,
            scale_params[held][key],
        )
        candidate_pds = (
            pack.baseline_pds.copy()
            if fraction == 0.0
            else evaluate_residual(pack, residual, fraction)
        )
        baseline_by_fold[held] = pack.baseline_pds.copy()
        candidate_by_fold[held] = candidate_pds
        per_fold_delta[held] = candidate_pds - pack.baseline_pds
        baseline_scores[held] = float(np.mean(pack.baseline_pds))
        scores[held] = float(np.mean(candidate_pds))
        deltas[held] = scores[held] - baseline_scores[held]
        source_similarities.append(pack.source_similarity.copy())

        sign_draws, label_draws = fold_null_draws(
            pack,
            scale_params[held][key],
            representation,
            fraction,
            label_permutations,
            seed=seed + 100 + fold_index,
        )
        sign_null_by_fold.append(sign_draws)
        label_null_by_fold.append(label_draws)
        fold_details[held] = {
            "source_names": list(pack.source_names),
            "operator_genes": int(len(pack.wcols)),
            "retrieval_genes": int(len(pack.retrieval_positions)),
            "whitener_sha256": pack.whitener_sha256,
            "whitening": dict(pack.whitening),
        }
        del pack
        gc.collect()
        enforce_rss(f"{name} {held} fold released", max_rss_gib)

    clusters = derive_target_clusters(source_similarities)
    bootstrap = cluster_bootstrap(
        per_fold_delta,
        clusters,
        reps=bootstrap_reps,
        seed=seed,
    )
    sign_draws = np.mean(np.stack(sign_null_by_fold, axis=0), axis=0)
    label_draws = np.mean(np.stack(label_null_by_fold, axis=0), axis=0)
    nulls = {
        "reps": null_reps,
        "seed": seed + 1,
        "exact_sparsity_sign_null": {
            "definition": (
                "independent random sign for every released nonzero edge; exact "
                "target/gene support, magnitudes, and nnz preserved"
            ),
            "mean": float(np.mean(sign_draws)),
            "p95": float(np.percentile(sign_draws, 95)),
            "p99": float(np.percentile(sign_draws, 99)),
        },
        "matched_label_null": {
            "definition": (
                "whole IGVF target rows permuted within nnz, L2-norm, and "
                "cell-count strata; the same permutation is used in every fold"
            ),
            "mean": float(np.mean(label_draws)),
            "p95": float(np.percentile(label_draws, 95)),
            "p99": float(np.percentile(label_draws, 99)),
            "strata": int(len(np.unique(strata))),
        },
    }
    macro_delta = float(np.mean(list(deltas.values())))
    covered = np.asarray(
        [
            any(
                position not in excluded_positions
                for position in igvf.effects.get(target, {})
            )
            for target in targets
        ],
        dtype=bool,
    )
    if np.any(
        np.concatenate(
            [
                per_fold_delta[held][~covered]
                for held in MAIN_NAMES
            ]
        )
        != 0
    ):
        raise AssertionError(f"{name}: missing IGVF targets changed")
    sign_p95 = float(nulls["exact_sparsity_sign_null"]["p95"])
    label_p95 = float(nulls["matched_label_null"]["p95"])
    gate = {
        "fraction_nonzero": float(selection["fraction"]) > 0,
        "every_fold_nonnegative": min(deltas.values()) >= 0,
        "macro_delta_at_least_0.003": macro_delta >= 0.003,
        "cluster_bootstrap_ci_low_positive": float(bootstrap["ci95"][0]) > 0,
        "beats_sign_null_p95": macro_delta > sign_p95,
        "beats_matched_label_null_p95": macro_delta > label_p95,
        "complete_gallery": True,
        "missing_targets_exactly_preserved": True,
    }
    gate["passed"] = all(gate.values())

    rows: list[dict[str, object]] = []
    for index, target in enumerate(targets):
        row: dict[str, object] = {
            "stage": name,
            "target_gene": target,
            "cluster": int(clusters[index]),
            "igvf_covered": bool(covered[index]),
            "igvf_filtered_rows": int(igvf.row_counts.get(target, 0)),
            "igvf_num_cells": int(igvf.num_cells.get(target, 0)),
        }
        for held in MAIN_NAMES:
            row[f"baseline_pds_{held}"] = float(baseline_by_fold[held][index])
            row[f"candidate_pds_{held}"] = float(candidate_by_fold[held][index])
            row[f"delta_{held}"] = float(per_fold_delta[held][index])
        row["macro_delta"] = float(
            np.mean([per_fold_delta[held][index] for held in MAIN_NAMES])
        )
        rows.append(row)
    result = {
        "stage": name,
        "targets": len(targets),
        "target_order_sha256": sha256_lines(targets),
        "igvf_covered_targets": int(covered.sum()),
        "igvf_missing_targets_preserved": int((~covered).sum()),
        "baseline_pds": baseline_scores,
        "candidate_pds": scores,
        "delta": deltas,
        "macro_baseline_pds": float(np.mean(list(baseline_scores.values()))),
        "macro_candidate_pds": float(np.mean(list(scores.values()))),
        "macro_delta": macro_delta,
        "bootstrap": bootstrap,
        "nulls": nulls,
        "gate": gate,
        "folds": fold_details,
    }
    enforce_rss(f"{name} complete", max_rss_gib)
    return result, rows


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    assert_output_scope(path)
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    assert_output_scope(args.output)
    official = read_csv_column(args.genes_csv, "gene_name")
    panel = read_csv_column(args.panel_csv, "target_gene")
    if len(official) != 18_533 or len(panel) != 300:
        raise ValueError("unexpected VCC gene or target count")
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_set = set(panel)
    panel_positions = {official_pos[gene] for gene in panel}

    downloads = verify_igvf_files(args.igvf_dir)
    igvf = load_igvf(
        args.igvf_dir / "IGVFFI0830FXFI.csv.gz",
        args.igvf_dir / "IGVFFI8270UPKB.csv.gz",
        official,
    )
    if igvf.de_rows != 27_578 or igvf.guide_rows != 14_151:
        raise ValueError("IGVF row counts differ from the audited release")
    if len(igvf.de_targets) != 1_889 or len(igvf.guide_targets) != 2_062:
        raise ValueError("IGVF target counts differ from the audited release")
    if len(igvf.de_targets & panel_set) != 87:
        raise ValueError("IGVF DE panel overlap differs from 87")
    if len(igvf.guide_targets & panel_set) != 91:
        raise ValueError("IGVF guide panel overlap differs from 91")

    banks = {
        name: Bank.load(name, args.donor_root / name) for name in MAIN_NAMES
    }
    all_three_targets = set.intersection(*(set(bank.perts) for bank in banks.values()))
    nonpanel_raw_shared = sorted((set(igvf.de_targets) & all_three_targets) - panel_set)
    nonpanel_shared = sorted((set(igvf.effects) & all_three_targets) - panel_set)
    panel_gallery = tuple(target for target in panel if target in all_three_targets)
    panel_raw_igvf = set(panel_gallery) & set(igvf.de_targets)
    panel_igvf = set(panel_gallery) & set(igvf.effects)
    if len(nonpanel_shared) < args.min_nonpanel_targets:
        raise ValueError(
            f"no-go: only {len(nonpanel_shared)} non-panel IGVF/all-donor targets"
        )
    if len(panel_gallery) < args.min_panel_gallery:
        raise ValueError(f"no-go: only {len(panel_gallery)} complete panel targets")
    if len(panel_igvf) < args.min_panel_igvf_targets:
        raise ValueError(f"no-go: only {len(panel_igvf)} IGVF-covered panel targets")
    if not set(nonpanel_shared) <= igvf.guide_targets:
        missing = sorted(set(nonpanel_shared) - igvf.guide_targets)
        raise ValueError(f"DE targets absent from guide manifest: {missing[:20]}")

    splits = split_targets(nonpanel_shared)
    if min(map(len, splits.values())) < args.min_split_targets:
        raise ValueError(
            f"no-go: split sizes { {key: len(value) for key, value in splits.items()} }"
        )
    proxy_target_positions = {
        official_pos[target]
        for target in nonpanel_shared
        if target in official_pos
    }
    proxy_excluded = panel_positions | proxy_target_positions
    fit_values = [
        abs(value)
        for target in splits["fit"]
        for position, value in igvf.effects.get(target, {}).items()
        if position not in proxy_excluded
    ]
    if len(fit_values) < 5_000:
        raise ValueError(f"no-go: only {len(fit_values)} fit target-response rows")
    clip_abs = float(np.percentile(np.asarray(fit_values), 99.0))

    stats: dict[str, BankStats] = {}
    for name, bank in banks.items():
        stats[name] = compute_bank_stats(
            bank,
            official,
            chunk_rows=args.donor_chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        print(
            f"[source] {name} median_norm={stats[name].raw_median_norm:.6f}",
            flush=True,
        )

    scale_params = fit_scale_parameters(
        splits["fit"],
        banks,
        stats,
        official,
        panel_set,
        proxy_excluded,
        igvf,
        clip_abs,
        target_chunk_rows=args.scale_target_chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    enforce_rss("scale parameters frozen", args.max_rss_gib)

    selection, tune_rows = tune_grid_streaming(
        splits["tune"],
        banks,
        stats,
        official,
        panel_set,
        proxy_excluded,
        igvf,
        clip_abs,
        args.whitener_dir,
        scale_params,
        operator_block_rows=args.operator_block_rows,
        max_rss_gib=args.max_rss_gib,
    )
    tune_summary = dict(selection["selected_tune_metrics"])
    gc.collect()
    enforce_rss("tune selection frozen", args.max_rss_gib)
    print(f"[selection] {json.dumps(selection, sort_keys=True)}", flush=True)

    test_result, test_rows = evaluate_stage(
        "nonpanel_test",
        splits["test"],
        banks,
        stats,
        official,
        panel_set,
        proxy_excluded,
        igvf,
        clip_abs,
        args.whitener_dir,
        scale_params,
        selection,
        operator_block_rows=args.operator_block_rows,
        max_rss_gib=args.max_rss_gib,
        bootstrap_reps=args.bootstrap_reps,
        null_reps=args.null_reps,
        seed=args.seed,
    )

    # This is the one and only point at which panel perturbation effects are read.
    panel_result, panel_rows = evaluate_stage(
        "vcc_panel_final",
        panel_gallery,
        banks,
        stats,
        official,
        panel_set,
        panel_positions,
        igvf,
        clip_abs,
        args.whitener_dir,
        scale_params,
        selection,
        operator_block_rows=args.operator_block_rows,
        max_rss_gib=args.max_rss_gib,
        bootstrap_reps=args.bootstrap_reps,
        null_reps=args.null_reps,
        seed=args.seed + 10_000,
    )

    promotion = {
        "nonpanel_test_gate": bool(test_result["gate"]["passed"]),
        "vcc_panel_final_gate": bool(panel_result["gate"]["passed"]),
    }
    promotion["promote"] = all(promotion.values())
    decision = "PROMOTE" if promotion["promote"] else "NO-GO"

    result: dict[str, object] = {
        "experiment": "081-igvf-cardiomyocyte",
        "status": "completed",
        "decision": decision,
        "promotion": promotion,
        "selection": selection,
        "tune_selected_summary": tune_summary,
        "nonpanel_test": test_result,
        "vcc_panel_final": panel_result,
        "coverage": {
            "de_rows": igvf.de_rows,
            "de_targets": len(igvf.de_targets),
            "guide_rows": igvf.guide_rows,
            "guide_targets": len(igvf.guide_targets),
            "response_genes": len(igvf.response_genes),
            "duplicate_target_response_pairs_aggregated": igvf.duplicate_pairs,
            "nonpanel_raw_de_targets_shared": len(nonpanel_raw_shared),
            "nonpanel_usable_targets_shared": len(nonpanel_shared),
            "split_sizes": {key: len(value) for key, value in splits.items()},
            "panel_complete_gallery": len(panel_gallery),
            "panel_raw_de_targets": len(panel_raw_igvf),
            "panel_official_response_targets": len(panel_igvf),
            "panel_usable_after_coordinate_exclusion": int(
                panel_result["igvf_covered_targets"]
            ),
            "panel_incumbent_fallback": int(
                panel_result["igvf_missing_targets_preserved"]
            ),
            "fit_rows_for_clip": len(fit_values),
            "fit_clip_abs_log2fc_p99": clip_abs,
        },
        "protocol": {
            "incumbent": (
                "equal-source held-donor LODO baseline with frozen source-only "
                "union whitening at gamma=0.65"
            ),
            "adapter": (
                "post-whitening sparse residual; unreported targets and genes "
                "remain exactly incumbent"
            ),
            "representations": list(REPRESENTATIONS),
            "scale_modes": list(SCALE_MODES),
            "fractions": list(FRACTIONS),
            "split": "SHA256 target split 60/20/20, seed prefix exp081:",
            "selection_data": "non-panel fit/tune targets only",
            "nonpanel_test_opened_after_selection": True,
            "panel_truth_opened_once_after_selection_and_nonpanel_test": True,
            "complete_galleries": True,
            "all_300_panel_coordinates_excluded": True,
            "all_proxy_target_coordinates_excluded_in_nonpanel_evaluation": True,
            "count_emission": False,
            "cell_eval2": False,
            "cloud": False,
            "submission": False,
        },
        "scale_parameters": {
            held: {
                key: asdict(value) for key, value in values.items()
            }
            for held, values in scale_params.items()
        },
        "input_hashes": {
            "downloads": downloads,
            "run_py": sha256_file(Path(__file__)),
            "test_run_py": sha256_file(HERE / "test_run.py"),
            "provenance_json": sha256_file(HERE / "PROVENANCE.json"),
            "panel_csv": sha256_file(args.panel_csv),
            "genes_csv": sha256_file(args.genes_csv),
            "whiteners": EXPECTED_WHITENERS,
            "donor_metadata": {
                name: {
                    "genes_txt": sha256_file(bank.path / "genes.txt"),
                    "perts_txt": sha256_file(bank.path / "perts.txt"),
                    "meta_json": sha256_file(bank.path / "meta.json"),
                }
                for name, bank in banks.items()
            },
        },
        "resource": {
            "max_allowed_rss_gib": args.max_rss_gib,
            "peak_rss_gib": peak_rss_gib(),
            "operator_block_rows": args.operator_block_rows,
            "donor_chunk_rows": args.donor_chunk_rows,
            "scale_target_chunk_rows": args.scale_target_chunk_rows,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(HERE / "tune_grid.csv", tune_rows)
    write_csv(HERE / "per_target_pds.csv", [*test_rows, *panel_rows])
    print(
        json.dumps(
            {
                "decision": decision,
                "selection": selection,
                "test_macro_delta": test_result["macro_delta"],
                "panel_macro_delta": panel_result["macro_delta"],
                "peak_rss_gib": peak_rss_gib(),
            },
            indent=2,
        ),
        flush=True,
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-root", type=Path, default=ROOT / "data/donors")
    parser.add_argument("--igvf-dir", type=Path, default=DATA)
    parser.add_argument(
        "--whitener-dir",
        type=Path,
        default=(
            ROOT
            / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
        ),
    )
    parser.add_argument(
        "--genes-csv", type=Path, default=ROOT / "data/challenge/gene_names.csv"
    )
    parser.add_argument(
        "--panel-csv", type=Path, default=ROOT / "data/challenge/pert_counts.csv"
    )
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument("--donor-chunk-rows", type=int, default=32)
    parser.add_argument("--operator-block-rows", type=int, default=256)
    parser.add_argument("--scale-target-chunk-rows", type=int, default=64)
    parser.add_argument("--max-rss-gib", type=float, default=MAX_RSS_GIB)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--null-reps", type=int, default=NULL_REPS)
    parser.add_argument("--seed", type=int, default=RNG_SEED)
    parser.add_argument("--min-nonpanel-targets", type=int, default=1_000)
    parser.add_argument("--min-split-targets", type=int, default=250)
    parser.add_argument("--min-panel-gallery", type=int, default=250)
    parser.add_argument("--min-panel-igvf-targets", type=int, default=60)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing unguarded real-data execution")
    if not 0 < args.max_rss_gib < 1.0:
        raise SystemExit("--max-rss-gib must be positive and below 1.0")
    if (
        args.donor_chunk_rows < 1
        or args.operator_block_rows < 1
        or args.scale_target_chunk_rows < 1
    ):
        raise SystemExit("chunk sizes must be positive")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
